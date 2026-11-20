"""Step 3: paired robustness check for the tumour vs normal results.

Many TCGA patients donated both a tumour (barcode type 01) and an adjacent normal
sample (type 11). Comparing each patient's tumour with their own normal removes
between-person differences. If a Step 2 hit is real, the paired test should agree
in direction. Same rules as Step 2: >20% missing = excluded, no imputation.
"""
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

GENES = [f"SLC30A{i}" for i in range(1, 11)]
MIN_PAIRS = 10
MAX_MISSING = 0.20

df = pd.read_csv("data/processed/znt_expression_all_samples.csv", index_col=0)
df = df[df["group"].isin(["tumour", "normal"])].copy()
df["patient"] = df.index.str[:12]          # TCGA-XX-XXXX = one patient

unpaired = pd.read_csv("results/tumour_vs_normal.csv")
tested = unpaired[unpaired["tested"]]

rows = []
for _, r in tested.iterrows():
    sub = df[df["cancer_type"] == r["cancer_type"]]
    if sub[r["gene"]].isna().mean() > MAX_MISSING:
        continue
    # one value per patient per group (mean if a patient has several aliquots)
    wide = sub.groupby(["patient", "group"])[r["gene"]].mean().unstack()
    if not {"tumour", "normal"} <= set(wide.columns):
        continue
    wide = wide.dropna(subset=["tumour", "normal"])
    row = dict(cancer_type=r["cancer_type"], gene=r["gene"], n_pairs=len(wide),
               unpaired_delta=r["delta_log2"], unpaired_fdr=r["fdr"])
    if len(wide) >= MIN_PAIRS:
        diff = wide["tumour"] - wide["normal"]
        p = wilcoxon(diff).pvalue if (diff != 0).any() else 1.0
        row.update(paired_median_diff=diff.median(), frac_pairs_down=(diff < 0).mean(), p_paired=p)
    rows.append(row)

res = pd.DataFrame(rows)
m = res["p_paired"].notna()
res.loc[m, "fdr_paired"] = multipletests(res.loc[m, "p_paired"], method="fdr_bh")[1]
same_dir = (res["paired_median_diff"] * res["unpaired_delta"]) > 0
res["confirmed"] = m & (res["unpaired_fdr"] < 0.05) & (res["fdr_paired"] < 0.05) & same_dir
res.to_csv("results/tumour_vs_normal_paired.csv", index=False)

big = res[(res["unpaired_fdr"] < 0.05) & (res["unpaired_delta"].abs() >= 1)]
print("Pairs with enough matched patients:", int(m.sum()), "of", len(res))
print("Unpaired hits (FDR<0.05):", int((res["unpaired_fdr"] < 0.05).sum()),
      "| confirmed by paired test:", int(res["confirmed"].sum()))
print("\nLarge-effect hits (|delta| >= 1 log2) and whether the paired test agrees:\n")
print(big.sort_values("unpaired_delta")[["cancer_type", "gene", "n_pairs", "unpaired_delta",
      "paired_median_diff", "frac_pairs_down", "fdr_paired", "confirmed"]]
      .round(3).to_string(index=False))
