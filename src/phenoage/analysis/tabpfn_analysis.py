"""TabPFN foundation-model analysis for mortality prediction.

TabPFN is a tabular foundation model (Prior-Data Fitted Network) that
performs in-context learning — it sees the entire training set as "context"
and predicts on new samples in a single forward pass, with no gradient-based
training.  This makes it especially strong on small datasets where tree
ensembles and neural networks tend to overfit.

Approach
--------
TabPFN is a classifier, not a native survival model.  We frame the
prediction task the same way Levine's Gompertz model does: predict
10-year mortality (binary: died within 120 months vs alive/censored
at ≥120 months).  This lets us compare AUC and calibration directly.

We also compute the concordance index (C-index) on the full survival
data using the predicted probability as the risk score, which gives a
fair apples-to-apples comparison with the Gompertz linear model and
RSF.

Models compared
---------------
1. Linear PhenoAge (Levine Gompertz, published coefficients)
2. TabPFN (PA9 features — same 10 features as PhenoAge)
3. TabPFN (PA9 + sex)
4. TabPFN (extended features — PA9 + sex + BMI, BP, lipids, etc.)

Usage:
    uv run python -m src.phenoage.tabpfn_analysis
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score

from ..constants import (
    FEATURE_COLS,
    LEVINE_BIOAGE,
    LEVINE_GAMMA,
    PHENOAGE_INTERCEPT,
    PHENOAGE_LOG_SCALE,
    PHENOAGE_RATE,
)
from ..model import compute_mortality_score

OUT_DIR = Path("output/phenoage_analysis")

# Feature sets (same as rsf_expanded.py)
PA9_FEATURES = (
    FEATURE_COLS  # albumin_gL, creat_umol, glucose_mmol, lncrp, lymph, mcv, rdw, alp, wbc, age
)
PA9_SEX_FEATURES = PA9_FEATURES + ["sex"]
EXTENDED_FEATURES = PA9_FEATURES + [
    "sex",
    "bmi",
    "sbp",
    "dbp",
    "total_chol",
    "triglycerides",
    "hdl",
    "hba1c",
    "bun",
    "uric_acid",
    "ggt",
    "total_bilirubin",
]

FEATURE_LABELS = {
    "albumin_gL": "Albumin (g/L)",
    "creat_umol": "Creatinine (µmol/L)",
    "glucose_mmol": "Glucose (mmol/L)",
    "lncrp": "ln(CRP) (mg/dL)",
    "lymph": "Lymphocyte %",
    "mcv": "MCV (fL)",
    "rdw": "RDW (%)",
    "alp": "ALP (U/L)",
    "wbc": "WBC (10⁹/L)",
    "age": "Age (years)",
    "sex": "Sex (0=M, 1=F)",
    "bmi": "BMI (kg/m²)",
    "sbp": "Systolic BP (mmHg)",
    "dbp": "Diastolic BP (mmHg)",
    "total_chol": "Total Cholesterol (mg/dL)",
    "triglycerides": "Triglycerides (mg/dL)",
    "hdl": "HDL (mg/dL)",
    "hba1c": "HbA1c (%)",
    "bun": "BUN (mg/dL)",
    "uric_acid": "Uric Acid (mg/dL)",
    "ggt": "GGT (U/L)",
    "total_bilirubin": "Total Bilirubin (mg/dL)",
}


# ═══════════════════════════════════════════════════════════════════
# DATA LOADING (reuse from existing modules)
# ═══════════════════════════════════════════════════════════════════


def load_train_data():
    """Load NHANES III training data (n ≈ 9,926) with extended features."""
    from .validate_models import _load_train

    df = _load_train()
    # Mortality definition: exclude accidental/external deaths, cap at 20 yr
    df.loc[df.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    df.loc[(df.mortstat == 1) & (df.permth_exm > 240), "mortstat"] = 0
    df.loc[df.permth_exm > 240, "permth_exm"] = 240
    # Sex feature
    df["sex"] = (df["HSSEX"] == 2).astype(float)
    # Map NHANES III variable names to canonical extended feature names
    iii_rename = {
        "BMPBMI": "bmi",
        "PEPMNK1R": "sbp",
        "PEPMNK5R": "dbp",
        "TCP": "total_chol",
        "TGP": "triglycerides",
        "HDP": "hdl",
        "GHP": "hba1c",
        "BUP": "bun",
        "UAP": "uric_acid",
        "GGPSI": "ggt",
        "TBP": "total_bilirubin",
    }
    for old, new in iii_rename.items():
        if old in df.columns:
            df[new] = df[old]
    return df


def load_test_data():
    """Load NHANES IV 1999-2010 validation data with expanded features."""
    from .rsf_expanded import build_validation_set, load_nhanes_iv_expanded

    raw = load_nhanes_iv_expanded()
    df = build_validation_set(raw, fasting_hours=8)
    return df


# ═══════════════════════════════════════════════════════════════════
# BINARY CLASSIFICATION TARGETS
# ═══════════════════════════════════════════════════════════════════


def make_10yr_target(df: pd.DataFrame, t_months: int = 120):
    """Create binary 10-year mortality target.

    Include only subjects who are:
      - dead within t_months (label=1), OR
      - alive with follow-up ≥ t_months (label=0)

    Subjects censored before t_months are EXCLUDED (we don't know
    their true outcome).  This is the standard approach for turning
    survival data into a classification problem at a fixed horizon.

    Returns (mask, y):
        mask: boolean array indexing into df
        y: int array of {0, 1}
    """
    died_before = (df.mortstat == 1) & (df.permth_exm <= t_months)
    alive_after = df.permth_exm >= t_months  # includes censored ≥ t_months
    mask = died_before | alive_after
    y = died_before.astype(int)
    return mask.values, y.values


# ═══════════════════════════════════════════════════════════════════
# LINEAR PHENOAGE BASELINE
# ═══════════════════════════════════════════════════════════════════


def linear_risk(df: pd.DataFrame) -> np.ndarray:
    """Compute linear predictor xb from Levine published coefficients."""
    xb = np.full(len(df), LEVINE_BIOAGE["intercept"])
    for feat in FEATURE_COLS:
        xb += LEVINE_BIOAGE[feat] * df[feat].values
    return xb


def linear_mortality_prob(df: pd.DataFrame, t_months: int = 120) -> np.ndarray:
    """Compute Levine Gompertz 10-year mortality probability."""
    xb = linear_risk(df)
    return compute_mortality_score(xb, gamma=LEVINE_GAMMA, t_months=t_months)


# ═══════════════════════════════════════════════════════════════════
# TABPFN WRAPPER
# ═══════════════════════════════════════════════════════════════════


def fit_tabpfn(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    n_estimators: int = 4,
    model_version: str = "v2",
    max_train: int = 2000,
    n_bags: int = 5,
):
    """Fit a TabPFN classifier via bagged ensemble.

    TabPFN inference cost scales super-linearly with training set size.
    To keep wall-clock time reasonable on CPU we subsample to `max_train`
    and average over `n_bags` independent subsamples (stratified so each
    bag has the same positive rate as the full dataset).

    Returns a list of fitted classifiers.
    """
    from tabpfn import TabPFNClassifier
    from tabpfn.settings import ModelVersion, settings

    # Use v2 (freely available, no gated access needed)
    version_map = {"v2": ModelVersion.V2, "v2.5": ModelVersion.V2_5, "v2.6": ModelVersion.V2_6}
    settings.tabpfn.model_version = version_map.get(model_version, ModelVersion.V2)

    X_arr = X_train.values if hasattr(X_train, "values") else np.asarray(X_train)
    eff_train = min(len(X_arr), max_train)

    print(
        f"\n  Fitting TabPFN ({model_version}): {n_bags} bags × "
        f"{eff_train} samples, p={X_arr.shape[1]}, "
        f"n_estimators={n_estimators}"
    )
    print(f"  Class balance: {y_train.mean():.3f} positive rate")

    rng = np.random.RandomState(42)
    clfs = []
    for i in range(n_bags):
        # Stratified subsample
        pos_idx = np.where(y_train == 1)[0]
        neg_idx = np.where(y_train == 0)[0]
        n_pos = min(len(pos_idx), max(1, int(eff_train * y_train.mean())))
        n_neg = min(len(neg_idx), eff_train - n_pos)
        sub_pos = rng.choice(pos_idx, n_pos, replace=len(pos_idx) < n_pos)
        sub_neg = rng.choice(neg_idx, n_neg, replace=len(neg_idx) < n_neg)
        sub_idx = np.concatenate([sub_pos, sub_neg])
        rng.shuffle(sub_idx)

        clf = TabPFNClassifier(
            n_estimators=n_estimators,
            ignore_pretraining_limits=True,
            device="cpu",
        )
        clf.fit(X_arr[sub_idx], y_train[sub_idx])
        clfs.append(clf)
        print(f"    Bag {i+1}/{n_bags}: n={len(sub_idx)}, " f"pos={y_train[sub_idx].sum()}")

    print(f"  TabPFN ensemble fit complete ({n_bags} bags).")
    return clfs


def predict_tabpfn(clfs: list, X: pd.DataFrame) -> np.ndarray:
    """Predict mortality probability with TabPFN (ensembled).

    Predictions are made in one call per classifier (no batching needed
    since train set is subsampled to a manageable size).
    """
    X_arr = X.values if hasattr(X, "values") else np.asarray(X)
    all_probs = np.zeros(len(X_arr))

    for clf in clfs:
        p = clf.predict_proba(X_arr)
        all_probs += p[:, 1]

    return all_probs / len(clfs)  # average across bags


# ═══════════════════════════════════════════════════════════════════
# EVALUATION
# ═══════════════════════════════════════════════════════════════════


def evaluate_classification(y_true: np.ndarray, y_prob: np.ndarray, label: str) -> dict:
    """Evaluate binary classification metrics."""
    auc = roc_auc_score(y_true, y_prob)
    brier = brier_score_loss(y_true, y_prob)
    # Clip for log_loss
    y_prob_clipped = np.clip(y_prob, 1e-10, 1 - 1e-10)
    ll = log_loss(y_true, y_prob_clipped)

    print(f"  {label:45s}  AUC={auc:.4f}  Brier={brier:.4f}  LogLoss={ll:.4f}")
    return {"label": label, "auc": auc, "brier": brier, "logloss": ll}


def evaluate_cindex(df: pd.DataFrame, risk: np.ndarray, label: str) -> dict:
    """Evaluate C-index on full survival data (including censored)."""
    from sksurv.metrics import concordance_index_censored

    event = df.mortstat.values.astype(bool)
    time = df.permth_exm.values.astype(float)
    c, concordant, discordant, tied, _ = concordance_index_censored(event, time, risk)
    print(f"  {label:45s}  C-index={c:.4f}")
    return {"label": label, "cindex": c, "concordant": concordant, "discordant": discordant}


# ═══════════════════════════════════════════════════════════════════
# CALIBRATION PLOT
# ═══════════════════════════════════════════════════════════════════


def plot_calibration(
    results: list[dict], test_df: pd.DataFrame, mask_10yr: np.ndarray, y_10yr: np.ndarray
):
    """Create calibration plots comparing models."""
    import altair as alt

    n_bins = 10
    rows = []
    for r in results:
        if "y_prob" not in r or "y_true" not in r:
            continue
        probs = r["y_prob"]
        y_true = r["y_true"]
        # Bin predictions into deciles
        bins = pd.qcut(probs, n_bins, duplicates="drop")
        for bin_label in bins.unique():
            mask = bins == bin_label
            rows.append(
                {
                    "model": r["label"],
                    "predicted": probs[mask].mean(),
                    "observed": y_true[mask].mean(),
                    "n": mask.sum(),
                }
            )

    cal_df = pd.DataFrame(rows)

    # Perfect calibration line
    line_df = pd.DataFrame({"x": [0, 0.5], "y": [0, 0.5]})
    ref_line = (
        alt.Chart(line_df).mark_line(color="gray", strokeDash=[4, 4]).encode(x="x:Q", y="y:Q")
    )

    points = (
        alt.Chart(cal_df)
        .mark_point(filled=True, size=80)
        .encode(
            x=alt.X(
                "predicted:Q",
                title="Predicted 10-yr Mortality Probability",
                scale=alt.Scale(domain=[0, 0.5]),
            ),
            y=alt.Y(
                "observed:Q",
                title="Observed 10-yr Mortality Rate",
                scale=alt.Scale(domain=[0, 0.5]),
            ),
            color=alt.Color("model:N", title="Model"),
            tooltip=["model:N", "predicted:Q", "observed:Q", "n:Q"],
        )
    )

    lines = (
        alt.Chart(cal_df)
        .mark_line()
        .encode(
            x="predicted:Q",
            y="observed:Q",
            color="model:N",
        )
    )

    chart = (ref_line + lines + points).properties(
        width=450, height=400, title="Calibration: Predicted vs Observed 10-Year Mortality"
    )
    return chart


# ═══════════════════════════════════════════════════════════════════
# ROC CURVE PLOT
# ═══════════════════════════════════════════════════════════════════


def plot_roc_curves(results: list[dict], y_true: np.ndarray):
    """Create ROC curves comparing models."""
    import altair as alt
    from sklearn.metrics import roc_curve

    rows = []
    for r in results:
        if "y_prob" not in r or "y_true" not in r:
            # Use default y_true only if lengths match
            if "y_prob" not in r:
                continue
            y = r.get("y_true", y_true)
            if len(y) != len(r["y_prob"]):
                continue
        else:
            y = r["y_true"]
        fpr, tpr, _ = roc_curve(y, r["y_prob"])
        # Subsample for reasonable chart size
        idx = np.linspace(0, len(fpr) - 1, 200).astype(int)
        for i in idx:
            rows.append(
                {
                    "model": f"{r['label']} (AUC={r['auc']:.3f})",
                    "fpr": fpr[i],
                    "tpr": tpr[i],
                }
            )

    roc_df = pd.DataFrame(rows)

    # Diagonal reference
    diag = pd.DataFrame({"x": [0, 1], "y": [0, 1]})
    ref = alt.Chart(diag).mark_line(color="gray", strokeDash=[4, 4]).encode(x="x:Q", y="y:Q")

    curves = (
        alt.Chart(roc_df)
        .mark_line(strokeWidth=2)
        .encode(
            x=alt.X("fpr:Q", title="False Positive Rate"),
            y=alt.Y("tpr:Q", title="True Positive Rate"),
            color=alt.Color("model:N", title="Model"),
        )
    )

    chart = (ref + curves).properties(
        width=450, height=400, title="ROC Curves: 10-Year Mortality Prediction"
    )
    return chart


# ═══════════════════════════════════════════════════════════════════
# FEATURE IMPORTANCE VIA PERMUTATION
# ═══════════════════════════════════════════════════════════════════


def tabpfn_permutation_importance(
    clf, X: pd.DataFrame, y: np.ndarray, n_repeats: int = 5
) -> pd.DataFrame:
    """Compute permutation importance for TabPFN."""
    print(f"\n  Computing permutation importance ({n_repeats} repeats)...")

    base_auc = roc_auc_score(y, predict_tabpfn(clf, X))
    results = []

    for j, feat in enumerate(X.columns):
        aucs = []
        for _ in range(n_repeats):
            X_perm = X.copy()
            X_perm[feat] = np.random.permutation(X_perm[feat].values)
            perm_probs = predict_tabpfn(clf, X_perm)
            aucs.append(roc_auc_score(y, perm_probs))
        drop = base_auc - np.mean(aucs)
        results.append(
            {
                "feature": feat,
                "feature_label": FEATURE_LABELS.get(feat, feat),
                "importance": drop,
                "importance_std": np.std([base_auc - a for a in aucs]),
            }
        )
        print(f"    {feat:20s}: ΔAUC = {drop:+.4f}")

    return pd.DataFrame(results).sort_values("importance", ascending=False)


def plot_feature_importance(imp_df: pd.DataFrame, title: str = "TabPFN Feature Importance"):
    """Bar chart of permutation importance."""
    import altair as alt

    chart = (
        alt.Chart(imp_df)
        .mark_bar()
        .encode(
            x=alt.X("importance:Q", title="ΔAUC (permutation drop)"),
            y=alt.Y("feature_label:N", sort="-x", title=""),
            color=alt.condition(
                alt.datum.importance > 0,
                alt.value("steelblue"),
                alt.value("salmon"),
            ),
            tooltip=["feature_label:N", "importance:Q", "importance_std:Q"],
        )
        .properties(width=450, height=max(200, len(imp_df) * 25), title=title)
    )

    return chart


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("  TABPFN FOUNDATION MODEL ANALYSIS")
    print("  Mortality prediction: NHANES III → NHANES IV")
    print("=" * 72)

    # ── Load data ──
    train_df = load_train_data()
    test_df = load_test_data()

    print(
        f"\n  Train (NHANES III): n={len(train_df):,d}, "
        f"deaths={int(train_df.mortstat.sum()):,d}"
    )
    print(f"  Test  (NHANES IV):  n={len(test_df):,d}, " f"deaths={int(test_df.mortstat.sum()):,d}")

    # ── 10-year binary classification target ──
    mask_train, y_train = make_10yr_target(train_df)
    mask_test, y_test = make_10yr_target(test_df)

    train_cls = train_df[mask_train].copy()
    test_cls = test_df[mask_test].copy()

    print(f"\n  10-year classification subsets (censored < 120 mo excluded):")
    print(
        f"    Train: n={len(train_cls):,d}, "
        f"deaths={y_train[mask_train].sum():,d} "
        f"({y_train[mask_train].mean():.1%} positive)"
    )
    print(
        f"    Test:  n={len(test_cls):,d}, "
        f"deaths={y_test[mask_test].sum():,d} "
        f"({y_test[mask_test].mean():.1%} positive)"
    )

    y_train_cls = y_train[mask_train]
    y_test_cls = y_test[mask_test]

    all_results = []

    # ══════════════════════════════════════════════════════════════
    # MODEL 1: Linear PhenoAge (Levine Gompertz)
    # ══════════════════════════════════════════════════════════════
    print("\n" + "─" * 72)
    print("  MODEL 1: Linear PhenoAge (Levine Gompertz)")
    print("─" * 72)

    levine_prob_test = linear_mortality_prob(test_cls)
    r_levine = evaluate_classification(y_test_cls, levine_prob_test, "Levine Gompertz")
    r_levine["y_prob"] = levine_prob_test
    r_levine["y_true"] = y_test_cls

    # C-index on full test set (including censored)
    levine_risk_full = linear_risk(test_df)
    ci_levine = evaluate_cindex(test_df, levine_risk_full, "Levine Gompertz")
    r_levine.update(ci_levine)
    all_results.append(r_levine)

    # ══════════════════════════════════════════════════════════════
    # MODEL 2: TabPFN — PA9 features (same 10 as PhenoAge)
    # ══════════════════════════════════════════════════════════════
    print("\n" + "─" * 72)
    print("  MODEL 2: TabPFN (PA9 features)")
    print("─" * 72)

    X_train_pa9 = train_cls[PA9_FEATURES]
    X_test_pa9 = test_cls[PA9_FEATURES]

    clf_pa9 = fit_tabpfn(X_train_pa9, y_train_cls)
    prob_pa9 = predict_tabpfn(clf_pa9, X_test_pa9)

    r_pa9 = evaluate_classification(y_test_cls, prob_pa9, "TabPFN (PA9)")
    r_pa9["y_prob"] = prob_pa9
    r_pa9["y_true"] = y_test_cls

    # C-index: predict on full test (includes censored subjects)
    risk_pa9_full = predict_tabpfn(clf_pa9, test_df[PA9_FEATURES])
    ci_pa9 = evaluate_cindex(test_df, risk_pa9_full, "TabPFN (PA9)")
    r_pa9.update(ci_pa9)
    all_results.append(r_pa9)

    # ══════════════════════════════════════════════════════════════
    # MODEL 3: TabPFN — PA9 + sex
    # ══════════════════════════════════════════════════════════════
    print("\n" + "─" * 72)
    print("  MODEL 3: TabPFN (PA9 + sex)")
    print("─" * 72)

    X_train_sex = train_cls[PA9_SEX_FEATURES]
    X_test_sex = test_cls[PA9_SEX_FEATURES]

    clf_sex = fit_tabpfn(X_train_sex, y_train_cls)
    prob_sex = predict_tabpfn(clf_sex, X_test_sex)

    r_sex = evaluate_classification(y_test_cls, prob_sex, "TabPFN (PA9 + sex)")
    r_sex["y_prob"] = prob_sex
    r_sex["y_true"] = y_test_cls

    risk_sex_full = predict_tabpfn(clf_sex, test_df[PA9_SEX_FEATURES])
    ci_sex = evaluate_cindex(test_df, risk_sex_full, "TabPFN (PA9 + sex)")
    r_sex.update(ci_sex)
    all_results.append(r_sex)

    # ══════════════════════════════════════════════════════════════
    # MODEL 4: TabPFN — Extended features
    # ══════════════════════════════════════════════════════════════
    print("\n" + "─" * 72)
    print("  MODEL 4: TabPFN (extended features)")
    print("─" * 72)

    available_ext = []
    for feat in EXTENDED_FEATURES:
        in_train = feat in train_cls.columns and train_cls[feat].notna().mean() > 0.5
        in_test = feat in test_cls.columns and test_cls[feat].notna().mean() > 0.5
        if in_train and in_test:
            available_ext.append(feat)
        elif in_test and not in_train:
            print(f"  ⚠ {feat}: in test but not train — skipping")

    print(f"  Available extended features: {len(available_ext)}")
    for f in available_ext:
        print(f"    {f}")

    clf_ext = None
    if len(available_ext) > len(PA9_SEX_FEATURES):
        train_ext = train_cls.dropna(subset=available_ext)
        test_ext = test_cls.dropna(subset=available_ext)
        y_train_ext = ((train_ext.mortstat == 1) & (train_ext.permth_exm <= 120)).astype(int).values
        y_test_ext = ((test_ext.mortstat == 1) & (test_ext.permth_exm <= 120)).astype(int).values

        print(f"  Extended train: n={len(train_ext):,d}, " f"deaths={y_train_ext.sum():,d}")
        print(f"  Extended test:  n={len(test_ext):,d}, " f"deaths={y_test_ext.sum():,d}")

        X_train_ext = train_ext[available_ext]
        X_test_ext = test_ext[available_ext]

        clf_ext = fit_tabpfn(X_train_ext, y_train_ext)
        prob_ext = predict_tabpfn(clf_ext, X_test_ext)

        r_ext = evaluate_classification(
            y_test_ext, prob_ext, f"TabPFN ({len(available_ext)} feats)"
        )
        r_ext["y_prob"] = prob_ext
        r_ext["y_true"] = y_test_ext

        # C-index on complete-extended subset
        test_ext_full = test_df.dropna(subset=available_ext)
        risk_ext_full = predict_tabpfn(clf_ext, test_ext_full[available_ext])
        ci_ext = evaluate_cindex(
            test_ext_full, risk_ext_full, f"TabPFN ({len(available_ext)} feats)"
        )
        r_ext.update(ci_ext)
        all_results.append(r_ext)
    else:
        print("  No additional features available — skipping extended model.")

    # ══════════════════════════════════════════════════════════════
    # MODEL 5: RSF baseline (C-index only — skip slow survival func)
    # ══════════════════════════════════════════════════════════════
    print("\n" + "─" * 72)
    print("  MODEL 5: RSF baseline (PA9 + sex) — C-index only")
    print("─" * 72)

    try:
        from sksurv.ensemble import RandomSurvivalForest
        from sksurv.metrics import concordance_index_censored

        def _struct_y(df):
            return np.array(
                [(bool(e), t) for e, t in zip(df.mortstat.values, df.permth_exm.values)],
                dtype=[("event", bool), ("time", float)],
            )

        rsf = RandomSurvivalForest(
            n_estimators=200,
            min_samples_leaf=15,
            n_jobs=-1,
            random_state=42,
        )
        rsf.fit(train_df[PA9_SEX_FEATURES], _struct_y(train_df))
        rsf_risk = rsf.predict(test_df[PA9_SEX_FEATURES])

        ci_rsf = evaluate_cindex(test_df, rsf_risk, "RSF (PA9 + sex)")
        r_rsf = {"label": "RSF (PA9 + sex)"}
        r_rsf.update(ci_rsf)
        all_results.append(r_rsf)
    except Exception as e:
        print(f"  RSF failed: {e}")

    # ══════════════════════════════════════════════════════════════
    # SUMMARY TABLE
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 72)
    print("  SUMMARY: MODEL COMPARISON")
    print("=" * 72)
    print(f"  {'Model':<45s} {'AUC':>7s} {'Brier':>7s} {'C-idx':>7s}")
    print(f"  {'-'*45} {'-'*7} {'-'*7} {'-'*7}")
    for r in all_results:
        auc_s = f"{r['auc']:.4f}" if "auc" in r else "   —  "
        brier_s = f"{r['brier']:.4f}" if "brier" in r else "   —  "
        ci_s = f"{r['cindex']:.4f}" if "cindex" in r else "   —  "
        print(f"  {r['label']:<45s} {auc_s:>7s} {brier_s:>7s} {ci_s:>7s}")

    # ══════════════════════════════════════════════════════════════
    # PLOTS
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 72)
    print("  GENERATING PLOTS")
    print("=" * 72)

    # ROC curves
    plot_results = [r for r in all_results if "y_prob" in r]
    roc_chart = plot_roc_curves(plot_results, y_test_cls)
    roc_chart.save(str(OUT_DIR / "tabpfn_roc_curves.html"))
    print(f"  Saved: {OUT_DIR}/tabpfn_roc_curves.html")

    # Calibration
    cal_chart = plot_calibration(plot_results, test_cls, mask_test, y_test_cls)
    cal_chart.save(str(OUT_DIR / "tabpfn_calibration.html"))
    print(f"  Saved: {OUT_DIR}/tabpfn_calibration.html")

    # ── Save summary JSON ──
    summary = []
    for r in all_results:
        summary.append(
            {k: v for k, v in r.items() if k != "y_prob" and not isinstance(v, np.ndarray)}
        )
    with open(OUT_DIR / "tabpfn_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  Saved: {OUT_DIR}/tabpfn_summary.json")

    print("\n" + "=" * 72)
    print("  DONE")
    print("=" * 72)


if __name__ == "__main__":
    main()
