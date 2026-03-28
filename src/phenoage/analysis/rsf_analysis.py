"""Random Survival Forest analysis + single-variable impact charts.

Builds on the PhenoAge pipeline:
  1. Single-variable impact charts for the LINEAR PhenoAge model
  2. Fits a Random Survival Forest (RSF) via scikit-survival
  3. SHAP analysis for the RSF
  4. "Optimal profile" analysis: best achievable mortality score by age/sex

Usage:
    uv run python -m src.phenoage.rsf_analysis
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

# ── Local imports ──
from ..constants import (
    FEATURE_COLS,
    GOMPERTZ_GAMMA,
    LEVINE_BIOAGE,
    LEVINE_COEFFICIENTS,
    LEVINE_GAMMA,
    PHENOAGE_INTERCEPT,
    PHENOAGE_LOG_SCALE,
    PHENOAGE_RATE,
)
from ..constants import gompertz_nll as _gompertz_nll
from ..model import compute_mortality_score, compute_phenoage

OUT_DIR = Path("output/phenoage_analysis")

# Pretty names for features
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
}

BIOMARKER_COLS = [f for f in FEATURE_COLS if f != "age"]

# RSF gets sex too (PhenoAge linear model does not use sex)
RSF_FEATURES = FEATURE_COLS + ["sex"]

# Healthy reference ranges (approximate, for display)
HEALTHY_RANGES = {
    "albumin_gL": (35, 50),
    "creat_umol": (53, 97),
    "glucose_mmol": (3.9, 5.6),
    "lncrp": (-4.6, -0.1),  # CRP 0.01–0.9 mg/dL
    "lymph": (20, 40),
    "mcv": (80, 100),
    "rdw": (11.5, 14.5),
    "alp": (44, 147),
    "wbc": (4.5, 11.0),
}


# ═══════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════


def load_train_data():
    """Load the n=9,926 NHANES III training sample."""
    from .validate_models import _load_train

    df = _load_train()
    return df


def load_test_data():
    """Load NHANES IV 1999-2010 test data with mortality."""
    from .validate_models import load_nhanes_iv_test

    df = load_nhanes_iv_test()
    return df


# ═══════════════════════════════════════════════════════════════════
# PART 1: LINEAR MODEL — SINGLE-VARIABLE IMPACT
# ═══════════════════════════════════════════════════════════════════


def linear_single_variable_impact(train_df: pd.DataFrame) -> pd.DataFrame:
    """For each biomarker, sweep its value while holding others at median.

    Returns a long-form DataFrame for charting.
    """
    # Median profile
    medians = train_df[FEATURE_COLS].median()
    intercept = LEVINE_BIOAGE["intercept"]
    gamma = LEVINE_GAMMA

    rows = []
    for feat in FEATURE_COLS:
        lo, hi = train_df[feat].quantile(0.01), train_df[feat].quantile(0.99)
        vals = np.linspace(lo, hi, 200)
        for v in vals:
            profile = medians.copy()
            profile[feat] = v
            xb = intercept + sum(LEVINE_BIOAGE[f] * profile[f] for f in FEATURE_COLS)
            M = float(compute_mortality_score(np.array([xb]), gamma=gamma)[0])
            pa = float(compute_phenoage(np.array([M]))[0])
            rows.append(
                {
                    "feature": FEATURE_LABELS.get(feat, feat),
                    "feature_key": feat,
                    "value": v,
                    "phenoage": pa,
                    "mortality_score": M,
                    "xb": xb,
                }
            )

    return pd.DataFrame(rows)


def plot_linear_impact(impact_df: pd.DataFrame, train_df: pd.DataFrame):
    """Create Altair charts for single-variable impact."""
    import altair as alt

    median_age = float(train_df["age"].median())

    charts = []
    for feat in FEATURE_COLS:
        label = FEATURE_LABELS.get(feat, feat)
        sub = impact_df[impact_df.feature_key == feat].copy()

        # Add a reference line at median age (PhenoAge = chronological age)
        base = (
            alt.Chart(sub)
            .mark_line(strokeWidth=2)
            .encode(
                x=alt.X("value:Q", title=label),
                y=alt.Y("phenoage:Q", title="PhenoAge (years)"),
                tooltip=["value:Q", "phenoage:Q", "mortality_score:Q"],
            )
            .properties(width=280, height=200, title=label)
        )

        # Horizontal line at median age
        rule = (
            alt.Chart(pd.DataFrame({"y": [median_age]}))
            .mark_rule(color="gray", strokeDash=[4, 4])
            .encode(y="y:Q")
        )

        # Healthy range shading if available
        if feat in HEALTHY_RANGES:
            lo, hi = HEALTHY_RANGES[feat]
            band = (
                alt.Chart(pd.DataFrame({"x": [lo], "x2": [hi]}))
                .mark_rect(color="green", opacity=0.08)
                .encode(x="x:Q", x2="x2:Q")
            )
            charts.append(band + base + rule)
        else:
            charts.append(base + rule)

    # Arrange in a 2x5 grid with independent y-scales per chart
    rows = []
    for i in range(0, len(charts), 3):
        row_charts = charts[i : i + 3]
        row = row_charts[0]
        for c in row_charts[1:]:
            row = row | c
        rows.append(row)
    combined = rows[0]
    for r in rows[1:]:
        combined = combined & r

    return combined


# ═══════════════════════════════════════════════════════════════════
# PART 2: RANDOM SURVIVAL FOREST
# ═══════════════════════════════════════════════════════════════════


def prepare_survival_data(df: pd.DataFrame, cap_months: int = 240, use_sex: bool = True):
    """Prepare data for scikit-survival RSF."""
    df = df.copy()

    # Exclude age-related deaths
    df.loc[df.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    df.loc[(df.mortstat == 1) & (df.permth_exm > cap_months), "mortstat"] = 0
    df.loc[df.permth_exm > cap_months, "permth_exm"] = cap_months

    # Add sex feature (0=Male, 1=Female for the RSF)
    if use_sex:
        sex_col = "HSSEX" if "HSSEX" in df.columns else "sex_code"
        if sex_col in df.columns:
            df["sex"] = (df[sex_col] == 2).astype(float)
        else:
            df["sex"] = 0.0  # fallback
        feature_cols = RSF_FEATURES
    else:
        feature_cols = FEATURE_COLS

    X = df[feature_cols].copy()
    # structured array for sksurv
    y = np.array(
        [(bool(e), t) for e, t in zip(df.mortstat.values, df.permth_exm.values)],
        dtype=[("event", bool), ("time", float)],
    )
    return X, y, df


def fit_rsf(
    X_train,
    y_train,
    n_estimators=500,
    max_depth=None,
    min_samples_leaf=15,
    n_jobs=-1,
    random_state=42,
):
    """Fit a Random Survival Forest using scikit-survival."""
    from sksurv.ensemble import RandomSurvivalForest

    print(
        f"\n  Fitting RSF: n={len(X_train)}, trees={n_estimators}, " f"min_leaf={min_samples_leaf}"
    )

    rsf = RandomSurvivalForest(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_leaf=min_samples_leaf,
        n_jobs=n_jobs,
        random_state=random_state,
    )
    rsf.fit(X_train, y_train)
    print("  RSF fit complete.")
    return rsf


def evaluate_rsf(rsf, X_test, y_test, label="RSF"):
    """Evaluate RSF with C-index."""
    from sksurv.metrics import concordance_index_censored

    # Risk score: negative of predicted survival at 120 months
    # Or use the built-in predict method which returns risk scores
    risk = rsf.predict(X_test)
    c, concordant, discordant, tied, _ = concordance_index_censored(
        y_test["event"], y_test["time"], risk
    )
    print(
        f"  {label} C-index: {c:.4f} " f"(concordant={concordant:,d}, discordant={discordant:,d})"
    )
    return {"cindex": c, "risk": risk}


# ═══════════════════════════════════════════════════════════════════
# PART 3: SHAP ANALYSIS
# ═══════════════════════════════════════════════════════════════════


def compute_shap_values(rsf, X_sample: pd.DataFrame, n_background: int = 100):
    """Compute SHAP-like values via permutation importance per sample.

    scikit-survival RSF is not supported by TreeExplainer, and KernelExplainer
    is too slow. We use a fast marginal contribution approach: for each sample,
    permute each feature and measure the change in prediction.
    """
    print(f"\n  Computing permutation-based feature attributions " f"on {len(X_sample)} samples...")

    np.random.seed(42)
    n = len(X_sample)
    n_features = X_sample.shape[1]

    # Baseline predictions
    base_preds = rsf.predict(X_sample)

    # For each feature, shuffle it and measure the change per sample
    shap_values = np.zeros((n, n_features))
    n_repeats = 5

    for j in range(n_features):
        feat = X_sample.columns[j]
        diffs = np.zeros(n)
        for _ in range(n_repeats):
            X_perm = X_sample.copy()
            X_perm.iloc[:, j] = np.random.permutation(X_perm.iloc[:, j].values)
            perm_preds = rsf.predict(X_perm)
            diffs += base_preds - perm_preds
        shap_values[:, j] = diffs / n_repeats
        print(f"    {feat}: mean attribution = {np.mean(np.abs(shap_values[:, j])):.4f}")

    print(f"  Attributions computed: shape={shap_values.shape}")
    return shap_values, X_sample


def plot_shap_importance(shap_values, X_sample):
    """Create SHAP summary using Altair (beeswarm-style)."""
    import altair as alt

    # Mean absolute SHAP per feature
    mean_abs = np.abs(shap_values).mean(axis=0)
    importance = pd.DataFrame(
        {
            "feature": [FEATURE_LABELS.get(f, f) for f in X_sample.columns],
            "mean_abs_shap": mean_abs,
        }
    ).sort_values("mean_abs_shap", ascending=True)

    chart = (
        alt.Chart(importance)
        .mark_bar()
        .encode(
            x=alt.X("mean_abs_shap:Q", title="Mean |SHAP value|"),
            y=alt.Y("feature:N", sort="-x", title=""),
        )
        .properties(width=400, height=300, title="RSF Feature Importance (SHAP)")
    )

    return chart


def plot_shap_dependence(shap_values, X_sample, train_df):
    """SHAP dependence plots for each feature — equivalent of linear impact."""
    import altair as alt

    charts = []
    for i, feat in enumerate(X_sample.columns):
        label = FEATURE_LABELS.get(feat, feat)
        dep_df = pd.DataFrame(
            {
                "value": X_sample[feat].values,
                "shap": shap_values[:, i],
            }
        )

        base = (
            alt.Chart(dep_df)
            .mark_circle(size=4, opacity=0.3)
            .encode(
                x=alt.X("value:Q", title=label),
                y=alt.Y("shap:Q", title="SHAP value"),
                tooltip=["value:Q", "shap:Q"],
            )
            .properties(width=280, height=200, title=label)
        )

        # Zero line
        rule = (
            alt.Chart(pd.DataFrame({"y": [0]}))
            .mark_rule(color="gray", strokeDash=[4, 4])
            .encode(y="y:Q")
        )

        # Healthy range
        if feat in HEALTHY_RANGES:
            lo, hi = HEALTHY_RANGES[feat]
            band = (
                alt.Chart(pd.DataFrame({"x": [lo], "x2": [hi]}))
                .mark_rect(color="green", opacity=0.08)
                .encode(x="x:Q", x2="x2:Q")
            )
            charts.append(band + base + rule)
        else:
            charts.append(base + rule)

    rows = []
    for i in range(0, len(charts), 2):
        row = charts[i]
        if i + 1 < len(charts):
            row = row | charts[i + 1]
        rows.append(row)
    combined = rows[0]
    for r in rows[1:]:
        combined = combined & r

    combined = combined.properties(title="RSF: SHAP Dependence Plots (Single-Variable Impact)")
    return combined


# ═══════════════════════════════════════════════════════════════════
# PART 4: OPTIMAL PROFILE ANALYSIS
# ═══════════════════════════════════════════════════════════════════


def optimal_linear_profile(train_df: pd.DataFrame):
    """Find the best/worst PhenoAge from the linear model.

    Since the model is linear, optimal = push each variable to its
    extreme in the direction of lower xb.
    """
    intercept = LEVINE_BIOAGE["intercept"]
    gamma = LEVINE_GAMMA

    results = []
    for age in range(20, 85, 5):
        for sex_label in ["Both"]:
            # Best: minimize xb by pushing each biomarker to its healthiest extreme
            # Worst: maximize xb
            profile_best = {}
            profile_worst = {}
            for feat in BIOMARKER_COLS:
                coef = LEVINE_BIOAGE[feat]
                lo = train_df[feat].quantile(0.01)
                hi = train_df[feat].quantile(0.99)
                if coef > 0:
                    profile_best[feat] = lo  # lower value → lower xb
                    profile_worst[feat] = hi
                else:
                    profile_best[feat] = hi  # higher value → lower xb (neg coef)
                    profile_worst[feat] = lo

            for label, profile in [
                ("Best", profile_best),
                ("Worst", profile_worst),
                ("Median", train_df[BIOMARKER_COLS].median().to_dict()),
            ]:
                profile["age"] = age
                xb = intercept + sum(LEVINE_BIOAGE[f] * profile[f] for f in FEATURE_COLS)
                M = float(compute_mortality_score(np.array([xb]), gamma=gamma)[0])
                pa = float(compute_phenoage(np.array([M]))[0])
                results.append(
                    {
                        "age": age,
                        "profile": label,
                        "phenoage": pa,
                        "phenoage_accel": pa - age,
                        "mortality_score": M,
                    }
                )

    return pd.DataFrame(results)


def optimal_rsf_profile(rsf, train_df: pd.DataFrame):
    """Find optimal biomarker profiles per age/sex using the RSF.

    For each age×sex, try the 5th and 95th percentile of each biomarker
    in a grid, then pick the combination with lowest predicted risk.
    Since full grid is 2^9=512, this is tractable.
    """
    from itertools import product

    results = []
    # Compute per-feature lo/hi
    extremes = {}
    for feat in BIOMARKER_COLS:
        extremes[feat] = (
            train_df[feat].quantile(0.05),
            train_df[feat].quantile(0.95),
        )

    # Generate all 2^9 = 512 biomarker combinations
    combos = list(product(*[extremes[f] for f in BIOMARKER_COLS]))

    # RSF uses sex — detect from feature names
    rsf_features = (
        list(rsf.feature_names_in_) if hasattr(rsf, "feature_names_in_") else FEATURE_COLS
    )
    has_sex = "sex" in rsf_features

    sex_vals = [("Male", 0.0), ("Female", 1.0)] if has_sex else [("Both", 0.0)]
    print(f"\n  Searching {len(combos)} biomarker combos × 13 ages × {len(sex_vals)} sex...")

    for age in range(20, 85, 5):
        for sex_label, sex_val in sex_vals:
            # Build X matrix for all combos at this age/sex
            rows = []
            for combo in combos:
                row = {}
                for j, feat in enumerate(BIOMARKER_COLS):
                    row[feat] = combo[j]
                row["age"] = age
                if has_sex:
                    row["sex"] = sex_val
                rows.append(row)
            X_df = pd.DataFrame(rows)[rsf_features]
            risk = rsf.predict(X_df)

            # Best (lowest risk) and worst (highest risk)
            best_idx = np.argmin(risk)
            worst_idx = np.argmax(risk)

            for label, idx in [("Best", best_idx), ("Worst", worst_idx)]:
                results.append(
                    {
                        "age": age,
                        "sex": sex_label,
                        "profile": label,
                        "risk_score": risk[idx],
                        **{f"opt_{f}": X_df.iloc[idx][f] for f in BIOMARKER_COLS},
                    }
                )

            # Median
            med_profile = train_df[BIOMARKER_COLS].median().to_dict()
            med_profile["age"] = age
            if has_sex:
                med_profile["sex"] = sex_val
            X_med = pd.DataFrame([med_profile])[rsf_features]
            med_risk = rsf.predict(X_med)[0]
            results.append(
                {
                    "age": age,
                    "sex": sex_label,
                    "profile": "Median",
                    "risk_score": med_risk,
                    **{f"opt_{f}": med_profile[f] for f in BIOMARKER_COLS},
                }
            )

    return pd.DataFrame(results)


def plot_optimal_profiles(linear_df, rsf_df):
    """Chart optimal profiles for linear and RSF models."""
    import altair as alt

    # Linear model: PhenoAge vs age
    linear_chart = (
        alt.Chart(linear_df)
        .mark_line(point=True)
        .encode(
            x=alt.X("age:Q", title="Chronological Age"),
            y=alt.Y("phenoage:Q", title="PhenoAge"),
            color=alt.Color(
                "profile:N",
                title="Profile",
                scale=alt.Scale(
                    domain=["Best", "Median", "Worst"], range=["green", "steelblue", "red"]
                ),
            ),
            strokeDash=alt.StrokeDash("profile:N"),
            tooltip=["age:Q", "profile:N", "phenoage:Q", "phenoage_accel:Q"],
        )
        .properties(width=400, height=300, title="Linear PhenoAge: Optimal vs Worst Profiles")
    )

    # Reference line: PhenoAge = chronological age
    ref = (
        alt.Chart(pd.DataFrame({"x": [20, 85], "y": [20, 85]}))
        .mark_line(color="gray", strokeDash=[4, 4])
        .encode(x="x:Q", y="y:Q")
    )

    linear_combined = ref + linear_chart

    # RSF model: risk score vs age, faceted by sex
    rsf_chart = (
        alt.Chart(rsf_df)
        .mark_line(point=True)
        .encode(
            x=alt.X("age:Q", title="Chronological Age"),
            y=alt.Y("risk_score:Q", title="RSF Risk Score"),
            color=alt.Color(
                "profile:N",
                title="Profile",
                scale=alt.Scale(
                    domain=["Best", "Median", "Worst"], range=["green", "steelblue", "red"]
                ),
            ),
            strokeDash=alt.StrokeDash("profile:N"),
            tooltip=["age:Q", "sex:N", "profile:N", "risk_score:Q"],
        )
        .properties(width=300, height=250)
        .facet(column=alt.Column("sex:N", title="Sex"))
        .properties(title="RSF: Optimal vs Worst Profiles by Sex")
    )

    return linear_combined & rsf_chart


def plot_optimal_biomarker_table(rsf_df):
    """Show what biomarker values the RSF picks as 'optimal' per age/sex."""
    import altair as alt

    best = rsf_df[rsf_df.profile == "Best"].copy()
    # Melt to long form
    opt_cols = [c for c in best.columns if c.startswith("opt_")]
    id_vars = ["age"]
    if "sex" in best.columns:
        id_vars.append("sex")
    melted = best.melt(
        id_vars=id_vars,
        value_vars=opt_cols,
        var_name="feature",
        value_name="optimal_value",
    )
    melted["feature"] = melted["feature"].str.replace("opt_", "")
    melted["feature_label"] = melted["feature"].map(lambda f: FEATURE_LABELS.get(f, f))

    base = (
        alt.Chart(melted)
        .mark_rect()
        .encode(
            x=alt.X("age:O", title="Age"),
            y=alt.Y("feature_label:N", title="", sort=[FEATURE_LABELS[f] for f in BIOMARKER_COLS]),
            color=alt.Color(
                "optimal_value:Q", scale=alt.Scale(scheme="viridis"), title="Optimal Value"
            ),
            tooltip=["age:O", "feature_label:N", "optimal_value:Q"]
            + (["sex:N"] if "sex" in melted.columns else []),
        )
        .properties(width=350, height=250)
    )

    if "sex" in melted.columns:
        chart = base.facet(column=alt.Column("sex:N", title="Sex")).properties(
            title="RSF: Optimal Biomarker Values by Age & Sex"
        )
    else:
        chart = base.properties(title="RSF: Optimal Biomarker Values by Age (Best Profile)")

    return chart


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("  PHENOAGE & RSF ANALYSIS")
    print("=" * 72)

    # ── Load data ──
    train_df = load_train_data()
    test_df = load_test_data()

    # Test data needs mortality-definition filtering
    test_df.loc[test_df.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    test_df.loc[(test_df.mortstat == 1) & (test_df.permth_exm > 240), "mortstat"] = 0
    test_df.loc[test_df.permth_exm > 240, "permth_exm"] = 240

    # ─────────────────────────────────────
    # PART 1: Linear single-variable impact
    # ─────────────────────────────────────
    print("\n" + "=" * 72)
    print("  PART 1: LINEAR MODEL — SINGLE-VARIABLE IMPACT")
    print("=" * 72)

    impact_df = linear_single_variable_impact(train_df)
    impact_df.to_parquet(OUT_DIR / "linear_impact.parquet", index=False)

    linear_chart = plot_linear_impact(impact_df, train_df)
    linear_chart.save(str(OUT_DIR / "linear_single_variable_impact.html"))
    print(f"  Saved: {OUT_DIR}/linear_single_variable_impact.html")

    # Optimal profiles (linear)
    linear_opt = optimal_linear_profile(train_df)
    linear_opt.to_csv(OUT_DIR / "linear_optimal_profiles.csv", index=False)
    print(f"\n  Linear optimal profiles:")
    for _, row in linear_opt.iterrows():
        print(
            f"    Age {row.age:3.0f} {row.profile:6s}: "
            f"PhenoAge={row.phenoage:.1f} (accel={row.phenoage_accel:+.1f})"
        )

    # ─────────────────────────────────────
    # PART 2: Random Survival Forest
    # ─────────────────────────────────────
    print("\n" + "=" * 72)
    print("  PART 2: RANDOM SURVIVAL FOREST")
    print("=" * 72)

    X_train, y_train, train_surv = prepare_survival_data(train_df)
    X_test, y_test, test_surv = prepare_survival_data(test_df)

    rsf = fit_rsf(X_train, y_train, n_estimators=500, min_samples_leaf=15)

    # Evaluate
    train_result = evaluate_rsf(rsf, X_train, y_train, "RSF (train)")
    test_result = evaluate_rsf(rsf, X_test, y_test, "RSF (test)")

    # Compare with linear model C-index on same data
    from sksurv.metrics import concordance_index_censored

    # Linear xb as risk score
    xb_test = np.full(len(X_test), LEVINE_BIOAGE["intercept"])
    for feat in FEATURE_COLS:
        xb_test += LEVINE_BIOAGE[feat] * X_test[feat].values
    c_linear, _, _, _, _ = concordance_index_censored(y_test["event"], y_test["time"], xb_test)
    print(f"  Linear Levine (test) C-index: {c_linear:.4f}")
    print(f"  RSF improvement: {test_result['cindex'] - c_linear:+.4f}")

    # ─────────────────────────────────────
    # PART 3: SHAP
    # ─────────────────────────────────────
    print("\n" + "=" * 72)
    print("  PART 3: SHAP ANALYSIS")
    print("=" * 72)

    # Use a subsample for SHAP (speed)
    np.random.seed(42)
    shap_n = min(2000, len(X_test))
    shap_idx = np.random.choice(len(X_test), shap_n, replace=False)
    X_shap = X_test.iloc[shap_idx].copy()

    shap_values, X_shap_out = compute_shap_values(rsf, X_shap, n_background=200)

    # SHAP importance chart
    importance_chart = plot_shap_importance(shap_values, X_shap_out)
    importance_chart.save(str(OUT_DIR / "rsf_shap_importance.html"))
    print(f"  Saved: {OUT_DIR}/rsf_shap_importance.html")

    # SHAP dependence charts
    dep_chart = plot_shap_dependence(shap_values, X_shap_out, train_df)
    dep_chart.save(str(OUT_DIR / "rsf_shap_dependence.html"))
    print(f"  Saved: {OUT_DIR}/rsf_shap_dependence.html")

    # ─────────────────────────────────────
    # PART 4: Optimal profiles
    # ─────────────────────────────────────
    print("\n" + "=" * 72)
    print("  PART 4: OPTIMAL PROFILES")
    print("=" * 72)

    rsf_opt = optimal_rsf_profile(rsf, train_df)
    rsf_opt.to_csv(OUT_DIR / "rsf_optimal_profiles.csv", index=False)

    # Chart
    opt_chart = plot_optimal_profiles(linear_opt, rsf_opt)
    opt_chart.save(str(OUT_DIR / "optimal_profiles.html"))
    print(f"  Saved: {OUT_DIR}/optimal_profiles.html")

    # Heatmap of RSF optimal values
    heat_chart = plot_optimal_biomarker_table(rsf_opt)
    heat_chart.save(str(OUT_DIR / "rsf_optimal_heatmap.html"))
    print(f"  Saved: {OUT_DIR}/rsf_optimal_heatmap.html")

    # Print RSF optimal profiles
    print(f"\n  RSF optimal biomarker profiles (Best):")
    best = rsf_opt[rsf_opt.profile == "Best"]
    for _, row in best.iterrows():
        sex_str = f" {row.sex:6s}" if "sex" in row.index else ""
        vals = ", ".join(
            f"{FEATURE_LABELS[f].split('(')[0].strip()}={row[f'opt_{f}']:.1f}"
            for f in BIOMARKER_COLS
        )
        print(f"    Age {row.age:3.0f}{sex_str}: risk={row.risk_score:.3f}  {vals}")

    # ─────────────────────────────────────
    # Summary
    # ─────────────────────────────────────
    print("\n" + "=" * 72)
    print("  SUMMARY")
    print("=" * 72)
    print(f"  Training: n={len(X_train)}, deaths={y_train['event'].sum()}")
    print(f"  Test: n={len(X_test)}, deaths={y_test['event'].sum()}")
    print(f"  Linear C-index (test): {c_linear:.4f}")
    print(f"  RSF C-index (test): {test_result['cindex']:.4f}")
    print(f"  RSF C-index (train): {train_result['cindex']:.4f}")
    print(f"\n  Output files in: {OUT_DIR}/")
    print(f"    linear_single_variable_impact.html")
    print(f"    rsf_shap_importance.html")
    print(f"    rsf_shap_dependence.html")
    print(f"    optimal_profiles.html")
    print(f"    rsf_optimal_heatmap.html")


if __name__ == "__main__":
    main()
