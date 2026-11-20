"""Step 1: extract SLC30A1-10 from the TCGA PanCanAtlas expression matrix,
label samples (tumour/normal, cancer type) and merge with survival data."""
import numpy as np
import pandas as pd
from pathlib import Path

RAW = Path("data/raw")
OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)

GENES = [f"SLC30A{i}" for i in range(1, 11)]

# 1. Stream the big genes x samples matrix, keeping only the ZnT rows
chunks = []
for chunk in pd.read_csv(RAW / "expression.xena.gz", sep="\t", index_col=0, chunksize=2000):
    keep = chunk.loc[chunk.index.isin(GENES)]
    if not keep.empty:
        chunks.append(keep)
expr = pd.concat(chunks).T[GENES]      # samples x genes, log2(norm_value+1)
expr.index.name = "sample"

# 2. Survival / clinical table (one row per tumour sample)
surv = pd.read_csv(RAW / "survival.tsv", sep="\t")
surv = surv.rename(columns={"cancer type abbreviation": "cancer_type"})

# 3. Label samples from the TCGA barcode: TCGA-XX-XXXX-01 -> 01 = primary tumour, 11 = solid normal
expr["patient"] = expr.index.str[:12]
expr["sample_type_code"] = expr.index.str[13:15]
expr["group"] = expr["sample_type_code"].map({"01": "tumour", "11": "normal"}).fillna("other")
patient_to_cancer = surv.drop_duplicates("_PATIENT").set_index("_PATIENT")["cancer_type"]
expr["cancer_type"] = expr["patient"].map(patient_to_cancer)
expr.to_csv(OUT / "znt_expression_all_samples.csv")

# 4. Merge expression with survival (tumour samples that have outcome data)
cols = ["sample", "_PATIENT", "cancer_type", "age_at_initial_pathologic_diagnosis", "gender",
        "ajcc_pathologic_tumor_stage", "OS", "OS.time", "DSS", "DSS.time",
        "DFI", "DFI.time", "PFI", "PFI.time"]
merged = surv[cols].merge(expr[GENES], left_on="sample", right_index=True, how="inner")
merged.to_csv(OUT / "znt_expression_survival.csv", index=False)

# 5. Quick sanity report
print("Expression table:", expr.shape)
print("\nSamples per group:\n", expr["group"].value_counts())
print("\nNormal samples per cancer type:\n", expr.loc[expr.group == "normal", "cancer_type"].value_counts().to_string())
print("\nMerged expression + survival:", merged.shape)
print("Cancer types with survival data:", merged["cancer_type"].nunique())
print("\nMissing values in merged table:\n", merged[GENES + ["OS", "OS.time"]].isna().sum().to_string())
