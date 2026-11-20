"""Step 2: tumour vs normal expression of SLC30A1-10 across TCGA cancer types."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

GENES = [f"SLC30A{i}" for i in range(1, 11)]
MIN_N = 10            # min samples per group (after dropping missing values)
MAX_MISSING = 0.20    # exclude gene-cancer pairs with >20% missing expression

Path("results/figures").mkdir(parents=True, exist_ok=True)
df = pd.read_csv("data/processed/znt_expression_all_samples.csv", index_col=0)
df = df[df["group"].isin(["tumour", "normal"])]

rows = []
for ct, sub in df.groupby("cancer_type"):
    if (sub["group"] == "normal").sum() < MIN_N:
        continue
    for g in GENES:
        t = sub.loc[sub.group == "tumour", g].dropna()
        n = sub.loc[sub.group == "normal", g].dropna()
        frac_missing = sub[g].isna().mean()
        row = dict(cancer_type=ct, gene=g, n_tumour=len(t), n_normal=len(n),
                   frac_missing=round(frac_missing, 3), tested=False)
        if frac_missing <= MAX_MISSING and len(t) >= MIN_N and len(n) >= MIN_N:
            _, p = mannwhitneyu(t, n, alternative="two-sided")
            row.update(tested=True, median_tumour=t.median(), median_normal=n.median(),
                       delta_log2=t.median() - n.median(), p=p)   # data are log2 scale
        rows.append(row)

res = pd.DataFrame(rows)
mask = res["tested"]
res.loc[mask, "fdr"] = multipletests(res.loc[mask, "p"], method="fdr_bh")[1]
res.to_csv("results/tumour_vs_normal.csv", index=False)

# Heatmap: genes x cancer types, colour = median log2 difference (tumour - normal), * = FDR < 0.05
delta = res.pivot(index="gene", columns="cancer_type", values="delta_log2").loc[GENES]
fdr = res.pivot(index="gene", columns="cancer_type", values="fdr").loc[GENES]
stars = fdr.map(lambda x: "*" if pd.notna(x) and x < 0.05 else "")
lim = float(np.nanmax(np.abs(delta.values)))
fig, ax = plt.subplots(figsize=(11, 5))
ax.set_facecolor("lightgrey")   # grey = not tested
sns.heatmap(delta, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=stars.values, fmt="",
            linewidths=0.5, cbar_kws={"label": "median log2 difference (tumour - normal)"}, ax=ax)
ax.set_title("SLC30A (ZnT) expression: tumour vs normal, TCGA\n* FDR < 0.05; grey = not tested")
ax.set_xlabel(""); ax.set_ylabel("")
plt.tight_layout()
plt.savefig("results/figures/tumour_vs_normal_heatmap.png", dpi=200)

# Summary
print("Gene-cancer pairs total:", len(res), "| tested:", int(mask.sum()))
print("Significant (FDR<0.05):", int((res["fdr"] < 0.05).sum()))
print("\nNot tested:\n", res.loc[~mask, ["cancer_type", "gene", "frac_missing"]].to_string(index=False))
sig = res[res["fdr"] < 0.05].copy()
sig["abs_delta"] = sig["delta_log2"].abs()
print("\nTop 10 by effect size:\n",
      sig.sort_values("abs_delta", ascending=False).head(10)[["cancer_type", "gene", "n_tumour", "n_normal", "delta_log2", "fdr"]].round(3).to_string(index=False))
