"""Step 5 (AI layer): do ZnT genes add prognostic information beyond clinical factors?

Cancer types: KIRC (kidney clear cell) and LGG (lower-grade glioma), the two types
with the strongest Step 4 signals and enough deaths to model.

Part A - penalised (elastic-net) Cox model, overall survival.
  Compares three feature sets with repeated 5-fold cross-validation:
  clinical only | ZnT genes only | clinical + ZnT. Metric: C-index on held-out
  patients (0.5 = coin flip, 1.0 = perfect ranking of who dies first).
Part B - random forest predicting death within 3 years, explained with SHAP.
  Patients censored before 3 years are excluded (their 3-year status is unknown).
  Cross-validated ROC AUC is reported; SHAP shows which features drive predictions.

Clinical features: age, sex, and AJCC stage (KIRC only; LGG has no AJCC stage).
Missing data: genes with >20% missing in a cancer type are dropped; remaining rows
with any missing value are dropped (complete-case). Nothing is imputed.
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
from pathlib import Path
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict, StratifiedKFold

warnings.filterwarnings("ignore")
GENES = [f"SLC30A{i}" for i in range(1, 11)]
CANCERS = ["KIRC", "LGG"]
MAX_MISSING = 0.20
SEED = 42
YEARS = 3

Path("results/figures").mkdir(parents=True, exist_ok=True)
df = pd.read_csv("data/processed/znt_expression_survival.csv", index_col=0)
df = df[df.index.str[13:15] == "01"]                      # primary tumours only
stage_map = {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4}


def prepare(ct):
    sub = df[df["cancer_type"] == ct].copy()
    genes = [g for g in GENES if sub[g].isna().mean() <= MAX_MISSING]
    sub["age"] = sub["age_at_initial_pathologic_diagnosis"]
    sub["male"] = (sub["gender"] == "MALE").astype(int)
    clinical = ["age", "male"]
    if ct == "KIRC":
        sub["stage"] = sub["ajcc_pathologic_tumor_stage"].map(stage_map)
        clinical.append("stage")
    sub = sub[sub["OS.time"] > 0]
    sub = sub[["OS", "OS.time"] + clinical + genes].dropna()
    return sub, clinical, genes


def cv_cindex(d, features, n_repeats=5):
    """Repeated stratified 5-fold CV C-index for an elastic-net Cox model."""
    scores = []
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=SEED)
    for tr, te in cv.split(d, d["OS"]):
        train, test = d.iloc[tr], d.iloc[te]
        mu, sd = train[features].mean(), train[features].std().replace(0, 1)
        Xtr = ((train[features] - mu) / sd).assign(T=train["OS.time"], E=train["OS"])
        Xte = (test[features] - mu) / sd
        cph = CoxPHFitter(penalizer=0.1, l1_ratio=0.5).fit(Xtr, "T", "E")
        risk = cph.predict_partial_hazard(Xte)
        scores.append(concordance_index(test["OS.time"], -risk, test["OS"]))
    return np.mean(scores), np.std(scores)


rows, coef_tables = [], []
for ct in CANCERS:
    d, clinical, genes = prepare(ct)
    print(f"\n===== {ct}: {len(d)} patients, {int(d['OS'].sum())} deaths, genes used: {len(genes)} =====")

    # ---- Part A: elastic-net Cox, three feature sets ----
    for name, feats in [("clinical only", clinical), ("ZnT genes only", genes),
                        ("clinical + ZnT", clinical + genes)]:
        m, s = cv_cindex(d, feats)
        rows.append(dict(cancer_type=ct, model="elastic-net Cox", features=name,
                         metric="C-index", mean=m, sd=s, n=len(d)))
        print(f"  Cox C-index  {name:<16} {m:.3f} ± {s:.3f}")

    # coefficients from a fit on all patients (for the README table)
    feats = clinical + genes
    X = (d[feats] - d[feats].mean()) / d[feats].std()
    cph = CoxPHFitter(penalizer=0.1, l1_ratio=0.5).fit(X.assign(T=d["OS.time"], E=d["OS"]), "T", "E")
    coef_tables.append(cph.params_.rename("coef_per_SD").to_frame().assign(cancer_type=ct))

    # ---- Part B: random forest, death within 3 years, + SHAP ----
    horizon = 365.25 * YEARS
    b = d[(d["OS.time"] >= horizon) | (d["OS"] == 1)].copy()
    b["dead_3y"] = ((b["OS"] == 1) & (b["OS.time"] < horizon)).astype(int)
    y = b["dead_3y"]
    print(f"  3-year task: {len(b)} patients ({y.sum()} died within {YEARS} years)")
    for name, feats in [("clinical only", clinical), ("ZnT genes only", genes),
                        ("clinical + ZnT", clinical + genes)]:
        aucs = []
        for rep in range(5):
            rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=5,
                                        class_weight="balanced", random_state=SEED + rep, n_jobs=-1)
            cv = StratifiedKFold(5, shuffle=True, random_state=SEED + rep)
            p = cross_val_predict(rf, b[feats], y, cv=cv, method="predict_proba")[:, 1]
            aucs.append(roc_auc_score(y, p))
        rows.append(dict(cancer_type=ct, model="random forest (3-yr death)", features=name,
                         metric="ROC AUC", mean=np.mean(aucs), sd=np.std(aucs), n=len(b)))
        print(f"  RF AUC       {name:<16} {np.mean(aucs):.3f} ± {np.std(aucs):.3f}")

    feats = clinical + genes
    rf = RandomForestClassifier(n_estimators=500, min_samples_leaf=5, class_weight="balanced",
                                random_state=SEED, n_jobs=-1).fit(b[feats], y)
    sv = shap.TreeExplainer(rf).shap_values(b[feats])
    sv = sv[1] if isinstance(sv, list) else (sv[..., 1] if np.ndim(sv) == 3 else sv)
    plt.figure()
    shap.summary_plot(sv, b[feats], show=False, max_display=len(feats))
    plt.title(f"{ct}: what drives predicted 3-year death (random forest, SHAP)")
    plt.tight_layout()
    plt.savefig(f"results/figures/shap_{ct}.png", dpi=200, bbox_inches="tight")
    plt.close()
    imp = pd.Series(np.abs(sv).mean(0), index=feats).sort_values(ascending=False)
    print("  Mean |SHAP| ranking:", ", ".join(f"{k} {v:.3f}" for k, v in imp.items()))

res = pd.DataFrame(rows)
res.to_csv("results/ml_model_comparison.csv", index=False)
pd.concat(coef_tables).to_csv("results/elastic_net_cox_coefficients.csv")

# Comparison figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, (model, metric) in zip(axes, [("elastic-net Cox", "C-index"), ("random forest (3-yr death)", "ROC AUC")]):
    r = res[res["model"] == model]
    piv = r.pivot(index="features", columns="cancer_type", values="mean").loc[
        ["clinical only", "ZnT genes only", "clinical + ZnT"]]
    err = r.pivot(index="features", columns="cancer_type", values="sd").loc[piv.index]
    piv.plot.bar(yerr=err, ax=ax, rot=0, capsize=3)
    ax.axhline(0.5, color="grey", ls="--", lw=1)
    ax.set_ylim(0.4, 1.0); ax.set_ylabel(f"cross-validated {metric}"); ax.set_xlabel("")
    ax.set_title(model)
plt.tight_layout()
plt.savefig("results/figures/ml_model_comparison.png", dpi=200)
print("\nSaved: results/ml_model_comparison.csv, results/elastic_net_cox_coefficients.csv, figures.")
