"""Compare age-prediction RSF vs survival RSF on mortality C-index.

Question: If we train RSF to predict chronological age (much more data,
no mortality follow-up needed), does the resulting risk score predict
mortality as well as a survival-trained model?

Approach:
  1. Train RSF-age on ALL NHANES III+IV PA9-complete data (~33k, predict age)
  2. Train RSF-surv on NHANES III mortality-linked data (~10k, predict survival)
  3. Evaluate BOTH on the same held-out NHANES IV validation set (mortality C-index)
  4. Also compare linear PhenoAge

Usage:
    uv run python -m src.phenoage.age_vs_survival
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored

from ..constants import (
    FEATURE_COLS,
    LEVINE_BIOAGE,
    LEVINE_GAMMA,
    PHENOAGE_INTERCEPT,
    PHENOAGE_LOG_SCALE,
    PHENOAGE_RATE,
)
from ..model import compute_mortality_score, compute_phenoage

OUT_DIR = Path("output/phenoage_analysis")

BIOMARKER_COLS = [f for f in FEATURE_COLS if f != "age"]


# ═══════════════════════════════════════════════════════════════════
# DATA
# ═══════════════════════════════════════════════════════════════════


def load_nhanes_iii_full():
    """Load ALL NHANES III with PA9 biomarkers (no mortality requirement)."""
    from ..reproduction.find_42 import (
        apply_5sd_outlier,
        apply_creat_calibration,
        load_all_biomarkers,
    )

    df = load_all_biomarkers()
    df = apply_creat_calibration(df)
    df = df[(df.HSAGEIR >= 20) & (df.HSAGEIR <= 84)].copy()

    # Derived features
    df["albumin_gL"] = df["AMP"] * 10
    df["creat_umol"] = df["CEP"] * 88.4017
    df["glucose_mmol"] = df["SGP"] * 0.0555
    df["lncrp"] = np.where(df["CRP"] > 0, np.log(df["CRP"]), np.nan)
    df = df.rename(
        columns={
            "LMPPCNT": "lymph",
            "MVPSI": "mcv",
            "RWP": "rdw",
            "APPSI": "alp",
            "WCPSI": "wbc",
            "HSAGEIR": "age",
        }
    )
    df["sex"] = (df["HSSEX"] == 2).astype(float)
    df["source"] = "NHANES_III"

    # Keep only PA9 complete
    df = df.dropna(subset=FEATURE_COLS)
    return df


def load_nhanes_iv_full():
    """Load ALL NHANES IV with PA9 biomarkers (no mortality requirement)."""
    from .rsf_expanded import load_nhanes_iv_expanded

    df = load_nhanes_iv_expanded()
    df = df[(df.age >= 20) & (df.age <= 84)].copy()
    df = df.dropna(subset=FEATURE_COLS + ["cycle"])
    df["source"] = "NHANES_IV"
    return df


def load_validation_set():
    """Load Levine-style validation: NHANES IV, fasting≥8h, with mortality."""
    from .rsf_expanded import build_validation_set, load_nhanes_iv_expanded

    test_raw = load_nhanes_iv_expanded()
    return build_validation_set(test_raw, fasting_hours=8)


def make_structured_y(df):
    return np.array(
        [(bool(e), t) for e, t in zip(df.mortstat.values, df.permth_exm.values)],
        dtype=[("event", bool), ("time", float)],
    )


def compute_linear_risk(df):
    xb = np.full(len(df), LEVINE_BIOAGE["intercept"])
    for feat in FEATURE_COLS:
        xb += LEVINE_BIOAGE[feat] * df[feat].values
    return xb


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("  AGE-PREDICTION vs SURVIVAL-PREDICTION RSF")
    print("=" * 72)

    # ── 1. Load all data ──
    print("\n  Loading data...")
    iii_full = load_nhanes_iii_full()
    print(f"  NHANES III (PA9 complete): n={len(iii_full):,d}")

    iv_full = load_nhanes_iv_full()
    print(f"  NHANES IV (PA9 complete): n={len(iv_full):,d}")

    print(f"\n  Loading validation set...")
    val_df = load_validation_set()
    val_y = make_structured_y(val_df)

    # Remove validation set SEQNs from IV training pool
    val_seqns = set(val_df.SEQN.values)
    iv_train = iv_full[~iv_full.SEQN.isin(val_seqns)].copy()
    print(f"  IV after removing validation: n={len(iv_train):,d}")

    # Combined age-prediction training set
    age_cols = FEATURE_COLS + ["sex"]  # age is already in FEATURE_COLS
    age_train = pd.concat(
        [
            iii_full[age_cols],
            iv_train[age_cols],
        ],
        ignore_index=True,
    ).dropna()
    print(f"  Combined age-prediction training: n={len(age_train):,d}")

    # Survival training (NHANES III with mortality)
    from .validate_models import _load_train

    surv_train = _load_train()
    surv_train.loc[surv_train.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    surv_train.loc[(surv_train.mortstat == 1) & (surv_train.permth_exm > 240), "mortstat"] = 0
    surv_train.loc[surv_train.permth_exm > 240, "permth_exm"] = 240
    surv_train["sex"] = (surv_train["HSSEX"] == 2).astype(float)
    print(
        f"  Survival training (NHANES III): n={len(surv_train):,d}, "
        f"deaths={int(surv_train.mortstat.sum()):,d}"
    )

    # ── 2. Train models ──
    features = FEATURE_COLS  # 9 biomarkers + age (no sex for fair comparison with linear)
    features_sex = FEATURE_COLS + ["sex"]
    biomarker_sex = BIOMARKER_COLS + ["sex"]  # NO age — for bio-age prediction

    print("\n" + "=" * 72)
    print("  FITTING MODELS")
    print("=" * 72)

    # Model A: Linear PhenoAge
    print("\n  Model A: Linear PhenoAge (Levine published coefficients)")
    xb_val = compute_linear_risk(val_df)

    # Model B: RSF-survival (standard)
    print("\n  Model B: RSF-survival (NHANES III, n=9,926)")
    surv_X = surv_train[features].copy()
    surv_y = make_structured_y(surv_train)
    rsf_surv = RandomSurvivalForest(
        n_estimators=500, min_samples_leaf=15, n_jobs=-1, random_state=42
    )
    rsf_surv.fit(surv_X, surv_y)
    risk_surv = rsf_surv.predict(val_df[features])

    # Model C: RF-age (predict chronological age, ALL data)
    print(f"\n  Model C: RF-age (predict age, n={len(age_train):,d})")
    rf_age = RandomForestRegressor(
        n_estimators=500, min_samples_leaf=15, n_jobs=-1, random_state=42
    )
    rf_age.fit(age_train[features], age_train["age"])
    predicted_age = rf_age.predict(val_df[features])
    # Use predicted age as risk score (higher predicted age = higher risk)
    # But also compute "age acceleration" = predicted_age - actual_age
    age_accel = predicted_age - val_df["age"].values

    # Model D: RF-age with sex
    print(f"\n  Model D: RF-age+sex (predict age, n={len(age_train):,d})")
    rf_age_sex = RandomForestRegressor(
        n_estimators=500, min_samples_leaf=15, n_jobs=-1, random_state=42
    )
    rf_age_sex.fit(age_train[features_sex], age_train["age"])
    predicted_age_sex = rf_age_sex.predict(val_df[features_sex])
    age_accel_sex = predicted_age_sex - val_df["age"].values

    # Model E: RSF-survival with sex
    print(f"\n  Model E: RSF-surv+sex (NHANES III, n=9,926)")
    surv_X_sex = surv_train[features_sex].copy()
    rsf_surv_sex = RandomSurvivalForest(
        n_estimators=500, min_samples_leaf=15, n_jobs=-1, random_state=42
    )
    rsf_surv_sex.fit(surv_X_sex, surv_y)
    risk_surv_sex = rsf_surv_sex.predict(val_df[features_sex])

    # Model F: RF-age trained only on biomarkers (predicts age FROM biomarkers)
    # Then the residual (predicted_age - actual_age) = "biological age acceleration"
    # This is the classic "biological age" approach
    print(
        f"\n  Model F: RF-biomarker-age (predict age from 9 biomarkers only, n={len(age_train):,d})"
    )
    rf_bio = RandomForestRegressor(
        n_estimators=500, min_samples_leaf=15, n_jobs=-1, random_state=42
    )
    rf_bio.fit(age_train[BIOMARKER_COLS + ["sex"]], age_train["age"])
    bio_age = rf_bio.predict(val_df[BIOMARKER_COLS + ["sex"]])
    bio_accel = bio_age - val_df["age"].values

    # Model G: Linear combination of bio_age + chronological age (like PhenoAge concept)
    # Use a simple Cox-like approach: risk = alpha * bio_age + beta * real_age
    # Fit alpha, beta on validation to see the best linear combo
    print(f"\n  Model G: Bio-age + chrono-age linear combination")
    # Actually, just use RSF-survival with bio_age as an additional feature
    # More interestingly: RSF-survival trained on 19k samples using bio_age features

    # Model H: RSF-survival trained on bio-age features from 19k
    # Train a survival model using bio_accel + age + sex on the 10k III data
    print(f"\n  Model H: RSF-surv with bio-age-accel + age + sex")
    # Predict bio_age on III training data
    bio_age_train = rf_bio.predict(surv_train[biomarker_sex])
    surv_train_enhanced = surv_train[["age"]].copy()
    surv_train_enhanced["bio_age"] = bio_age_train
    surv_train_enhanced["bio_accel"] = bio_age_train - surv_train["age"].values
    surv_train_enhanced["sex"] = surv_train["sex"].values

    val_enhanced = val_df[["age"]].copy()
    val_enhanced["bio_age"] = bio_age
    val_enhanced["bio_accel"] = bio_accel
    val_enhanced["sex"] = val_df["sex"].values

    rsf_enhanced = RandomSurvivalForest(
        n_estimators=500, min_samples_leaf=15, n_jobs=-1, random_state=42
    )
    rsf_enhanced.fit(surv_train_enhanced, surv_y)
    risk_enhanced = rsf_enhanced.predict(val_enhanced)

    # ── 3. Evaluate on mortality ──
    print("\n" + "=" * 72)
    print("  MORTALITY C-INDEX ON VALIDATION SET")
    print("=" * 72)
    print(f"  Validation: n={len(val_df):,d}, deaths={int(val_y['event'].sum())}")

    results = {}

    def eval_cindex(risk, label):
        c, conc, disc, tied, _ = concordance_index_censored(val_y["event"], val_y["time"], risk)
        results[label] = c
        return c

    c_linear = eval_cindex(xb_val, "Linear PhenoAge")
    c_rsf_surv = eval_cindex(risk_surv, "RSF-survival (10 feats)")
    c_rsf_surv_sex = eval_cindex(risk_surv_sex, "RSF-survival+sex (11 feats)")
    c_rf_age = eval_cindex(predicted_age, "RF-age (predicted age as risk)")
    c_rf_age_accel = eval_cindex(age_accel, "RF-age (age acceleration as risk)")
    c_rf_age_sex = eval_cindex(predicted_age_sex, "RF-age+sex (predicted age)")
    c_rf_age_accel_sex = eval_cindex(age_accel_sex, "RF-age+sex (age accel)")
    c_bio_age = eval_cindex(bio_age, "RF-biomarker-age (bio age)")
    c_bio_accel = eval_cindex(bio_accel, "RF-biomarker-age (bio accel)")
    c_chrono = eval_cindex(val_df["age"].values, "Chronological age alone")
    c_enhanced = eval_cindex(risk_enhanced, "RSF-surv(bio_accel+age+sex)")

    # ── Summary ──
    print(f"\n  {'Model':<50s} {'C-index':>8s} {'Train n':>8s}")
    print(f"  {'-'*50} {'-'*8} {'-'*8}")
    print(f"  {'Chronological age alone':<50s} {c_chrono:>8.4f} {'—':>8s}")
    print(f"  {'─'*50} {'─'*8} {'─'*8}")
    print(f"  {'Linear PhenoAge (Levine)':<50s} {c_linear:>8.4f} {'9,926':>8s}")
    print(f"  {'─'*50} {'─'*8} {'─'*8}")
    print(f"  {'RSF-survival (10 feats, III only)':<50s} {c_rsf_surv:>8.4f} {'9,926':>8s}")
    print(f"  {'RSF-survival+sex (11 feats, III only)':<50s} {c_rsf_surv_sex:>8.4f} {'9,926':>8s}")
    print(f"  {'─'*50} {'─'*8} {'─'*8}")
    print(
        f"  {'RF-age: predicted age as risk (10 feats)':<50s} {c_rf_age:>8.4f} {len(age_train):>8,d}"
    )
    print(
        f"  {'RF-age: age acceleration as risk':<50s} {c_rf_age_accel:>8.4f} {len(age_train):>8,d}"
    )
    print(
        f"  {'RF-age+sex: predicted age as risk':<50s} {c_rf_age_sex:>8.4f} {len(age_train):>8,d}"
    )
    print(
        f"  {'RF-age+sex: age acceleration as risk':<50s} {c_rf_age_accel_sex:>8.4f} {len(age_train):>8,d}"
    )
    print(f"  {'─'*50} {'─'*8} {'─'*8}")
    print(
        f"  {'RF-biomarker-age: bio age as risk (9+sex)':<50s} {c_bio_age:>8.4f} {len(age_train):>8,d}"
    )
    print(
        f"  {'RF-biomarker-age: bio accel as risk':<50s} {c_bio_accel:>8.4f} {len(age_train):>8,d}"
    )
    print(f"  {'─'*50} {'─'*8} {'─'*8}")
    print(f"  {'RSF-surv(bio_accel+age+sex, III only)':<50s} {c_enhanced:>8.4f} {'9,926':>8s}")

    # ── Age prediction quality ──
    print(f"\n  Age prediction quality:")
    from sklearn.metrics import mean_absolute_error, r2_score

    print(f"  {'Model':<45s} {'R²':>6s} {'MAE':>6s}")
    print(f"  {'-'*45} {'-'*6} {'-'*6}")

    for label, pred in [
        ("RF-age (10 feats)", predicted_age),
        ("RF-age+sex (11 feats)", predicted_age_sex),
        ("RF-biomarker-age (9+sex)", bio_age),
    ]:
        r2 = r2_score(val_df["age"].values, pred)
        mae = mean_absolute_error(val_df["age"].values, pred)
        print(f"  {label:<45s} {r2:>6.3f} {mae:>6.2f}")

    # PhenoAge for comparison
    M = compute_mortality_score(xb_val, gamma=LEVINE_GAMMA)
    phenoage = compute_phenoage(M)
    r2_pa = r2_score(val_df["age"].values, phenoage)
    mae_pa = mean_absolute_error(val_df["age"].values, phenoage)
    print(f"  {'Linear PhenoAge':<45s} {r2_pa:>6.3f} {mae_pa:>6.2f}")


if __name__ == "__main__":
    main()
