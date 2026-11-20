"""Step 4: is SLC30A1-10 expression linked to survival within each cancer type?

For each gene x cancer type x endpoint (OS = overall survival, PFI = progression-free
interval) we fit a Cox model: expression (z-scored within that cancer type, so the
hazard ratio is per 1 SD) adjusted for age at diagnosis. BH-FDR is applied separately
per endpoint. Primary tumours only (barcode type 01), one sample per patient.

Rules (same spirit as Steps 2-3):
- gene-cancer pairs with >20% missing expression are excluded, never imputed
- at least 20 events needed (about 10 events per variable for 2 variables)
- genes with almost no variation in a cancer type are skipped
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
GENES = [f"SLC30A{i}" for i in range(1, 11)]
ENDPOINTS = ["OS", "PFI"]
MIN_EVENTS = 20
MAX_MISSING = 0.20

Path("results/figures").mkdir(parents=True, exist_ok=True)
df = pd.read_csv("data/processed/znt_expression_survival.csv", index_col=0)
df = df[df.index.str[13:15] == "01"]                      # primary tumours only
df = df.rename(columns={"age_at_initial_pathologic_diagnosis": "age"})

rows = []
for ct, sub in df.groupby("cancer_type"):
    for ep in ENDPOINTS:
        base = sub[[ep, f"{ep}.time", "age"]].dropna()
        base = base[base[f"{ep}.time"] > 0]
        for g in GENES:
            row = dict(cancer_type=ct, endpoint=ep, gene=g, tested=False,
                       frac_missing=round(sub[g].isna().mean(), 3))
            d = base.join(sub[g]).dropna()
            row.update(n=len(d), events=int(d[ep].sum()) if len(d) else 0)
            if row["frac_missing"] > MAX_MISSING:
                row["reason"] = ">20% missing"
            elif row["events"] < MIN_EVENTS:
                row["reason"] = f"<{MIN_EVENTS} events"
            elif d[g].std() < 0.1:
                row["reason"] = "no variation"
            else:
                d = d.assign(expr_z=(d[g] - d[g].mean()) / d[g].std())
                try:
                    cph = CoxPHFitter().fit(d[[f"{ep}.time", ep, "expr_z", "age"]],
                                            duration_col=f"{ep}.time", event_col=ep)
                    s = cph.summary.loc["expr_z"]
                    row.update(tested=True, HR=s["exp(coef)"], HR_low=s["exp(coef) lower 95%"],
                               HR_high=s["exp(coef) upper 95%"], p=s["p"], reason="")
                except Exception as e:                    # convergence failure etc.
                    row["reason"] = f"fit failed: {type(e).__name__}"
            rows.append(row)

res = pd.DataFrame(rows)
for ep in ENDPOINTS:
    m = res["tested"] & (res["endpoint"] == ep)
    res.loc[m, "fdr"] = multipletests(res.loc[m, "p"], method="fdr_bh")[1]
res.to_csv("results/survival_cox.csv", index=False)

# Heatmaps: log2(HR) per 1 SD, * = FDR < 0.05, grey = not tested
fig, axes = plt.subplots(2, 1, figsize=(14, 9))
for ax, ep in zip(axes, ENDPOINTS):
    r = res[res["endpoint"] == ep]
    hr = np.log2(r.pivot(index="gene", columns="cancer_type", values="HR")).reindex(GENES)
    fdr = r.pivot(index="gene", columns="cancer_type", values="fdr").reindex(GENES)
    stars = fdr.map(lambda x: "*" if pd.notna(x) and x < 0.05 else "")
    ax.set_facecolor("lightgrey")
    sns.heatmap(hr, cmap="RdBu_r", center=0, vmin=-1, vmax=1, annot=stars.values, fmt="",
                linewidths=0.5, cbar_kws={"label": "log2 hazard ratio per 1 SD"}, ax=ax)
    ax.set_title(f"{ep}: red = higher expression, worse outcome; * FDR < 0.05; grey = not tested")
    ax.set_xlabel(""); ax.set_ylabel("")
plt.tight_layout()
plt.savefig("results/figures/survival_cox_heatmap.png", dpi=200)

# Summary
for ep in ENDPOINTS:
    r = res[res["endpoint"] == ep]
    sig = r[r["fdr"] < 0.05]
    print(f"\n=== {ep} ===  tested: {int(r['tested'].sum())} of {len(r)} | FDR<0.05: {len(sig)}"
          f" (worse with high expression: {(sig['HR'] > 1).sum()}, better: {(sig['HR'] < 1).sum()})")
    print(sig.sort_values("p")[["cancer_type", "gene", "n", "events", "HR", "HR_low", "HR_high", "fdr"]]
          .round(3).to_string(index=False))
print("\nNot tested, by reason:\n", res.loc[~res["tested"], "reason"].value_counts().to_string())
