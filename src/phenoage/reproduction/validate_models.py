"""Validate PhenoAge models on NHANES IV test data with linked mortality.

Compares three models:
  1. Levine's published coefficients (from the paper)
  2. Our "naive" MLE (no winsorizing, just 5SD + completeness filter)
  3. Our "best" MLE (winsorize 2/98 + cap glucose 6.0 + lymph 20/80)

Performance metrics:
  - C-index (concordance) for mortality prediction
  - PhenoAge-Age correlation
  - PhenoAge acceleration vs mortality (HR per year of acceleration)
  - AUC for 10-year mortality

Usage:
    uv run python -m src.phenoage.validate_models
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

# ── Local imports ──
from ..constants import (
    FEATURE_COLS,
    LEVINE_BIOAGE,
    LEVINE_GAMMA,
    PHENOAGE_INTERCEPT,
    PHENOAGE_LOG_SCALE,
    PHENOAGE_RATE,
)
from ..constants import gompertz_nll as _gompertz_nll
from ..model import compute_mortality_score, compute_phenoage
from .find_42 import (
    apply_5sd_outlier,
    apply_creat_calibration,
    build_candidate_list,
    load_all_biomarkers,
)

# ═══════════════════════════════════════════════════════════════════
# STEP 1 — Train three models on NHANES III
# ═══════════════════════════════════════════════════════════════════


def _load_train():
    """Load and prepare the n=9,926 training sample with 2015 mortality."""
    df = load_all_biomarkers()
    df = df.drop(columns=["mortstat", "ucod_leading", "permth_exm"], errors="ignore")
    mort = pd.read_parquet("data/processed/LMF_Files/LMF_all_MORT_2015.parquet")
    mort = mort[mort.cycle == "NHANES_III"][
        ["SEQN", "eligstat", "mortstat", "ucod_leading", "permth_exm"]
    ].copy()
    df = df.merge(mort.rename(columns={"eligstat": "elig15"}), on="SEQN", how="left")
    df = df[(df.HSAGEIR >= 20) & (df.HSAGEIR <= 84) & (df.elig15 == 1)].copy()
    df = apply_creat_calibration(df)

    # Mortality filters
    df.loc[df.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    df.loc[(df.mortstat == 1) & (df.permth_exm > 240), "mortstat"] = 0
    df.loc[df.permth_exm > 240, "permth_exm"] = 240

    # 5SD
    candidates = build_candidate_list()
    valid = {n: c for n, c in candidates.items() if c in df.columns}
    df = apply_5sd_outlier(df, list(valid.values()))

    # Completeness: PA9 + core6
    pa9 = [
        valid[n]
        for n in ["albumin", "alp", "creat", "glucose", "crp", "lymph", "mcv", "rdw", "wbc"]
    ]
    core6 = [valid[n] for n in ["ggt", "fev", "waist", "vitaminC", "cadmium", "trig"]]
    df = df.dropna(subset=pa9 + core6 + ["permth_exm", "mortstat"])

    # Derived features
    df["albumin_gL"] = df["AMP"] * 10
    df["creat_umol"] = df["CEP"] * 88.4017
    df["glucose_mmol"] = df["SGP"] * 0.0555
    df["lncrp"] = np.log(df["CRP"])
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
    df = df.dropna(subset=["lncrp"])
    return df


def _winsorize(x, lo=0.02, hi=0.98):
    q = np.quantile(x[~np.isnan(x)], [lo, hi])
    return np.clip(x, q[0], q[1])


def _fit(X, time, event, label):
    """Fit Gompertz PH and return (gamma, intercept, betas_dict)."""
    x0 = np.zeros(2 + X.shape[1])
    x0[0] = np.log(LEVINE_GAMMA)
    x0[1] = LEVINE_BIOAGE["intercept"]
    for i, col in enumerate(FEATURE_COLS):
        x0[2 + i] = LEVINE_BIOAGE[col]
    res = minimize(
        _gompertz_nll,
        x0,
        args=(X, time, event),
        method="L-BFGS-B",
        options={"maxiter": 10000, "ftol": 1e-15},
    )
    gamma = np.exp(res.x[0])
    intercept = res.x[1]
    betas = dict(zip(FEATURE_COLS, res.x[2:]))
    print(f"  {label}: gamma={gamma:.7f}, n={len(X)}, deaths={int(event.sum())}")
    return gamma, intercept, betas


def train_models():
    """Train and return three coefficient sets."""
    print("=" * 72)
    print("  TRAINING MODELS ON NHANES III (n=9,926, 2015 mortality)")
    print("=" * 72)

    df = _load_train()
    print(f"  Training sample: n={len(df)}, deaths={int(df.mortstat.sum())}")

    # ── Model 1: Levine published ──
    levine_gamma = LEVINE_GAMMA
    levine_intercept = LEVINE_BIOAGE["intercept"]
    levine_betas = {k: v for k, v in LEVINE_BIOAGE.items() if k != "intercept"}
    print(f"  Levine published: gamma={levine_gamma:.7f}")

    # ── Model 2: Naive MLE (no winsorizing) ──
    X = df[FEATURE_COLS].values.astype(float)
    time = df.permth_exm.values.astype(float)
    event = df.mortstat.values.astype(float)
    naive_gamma, naive_intercept, naive_betas = _fit(X, time, event, "Naive MLE")

    # ── Model 3: Best MLE (winsorize + caps) ──
    df_w = df.copy()
    for v in [
        "albumin_gL",
        "creat_umol",
        "glucose_mmol",
        "lncrp",
        "lymph",
        "mcv",
        "rdw",
        "alp",
        "wbc",
    ]:
        df_w[v] = _winsorize(df_w[v].values, 0.02, 0.98)
    df_w["glucose_mmol"] = np.minimum(df_w["glucose_mmol"], 6.0)
    df_w["lymph"] = _winsorize(df_w["lymph"].values, 0.20, 0.80)

    X_w = df_w[FEATURE_COLS].values.astype(float)
    best_gamma, best_intercept, best_betas = _fit(X_w, time, event, "Best MLE (winsorized)")

    models = {
        "Levine published": (levine_gamma, levine_intercept, levine_betas),
        "Naive MLE": (naive_gamma, naive_intercept, naive_betas),
        "Best MLE (winsorized)": (best_gamma, best_intercept, best_betas),
    }
    return models


# ═══════════════════════════════════════════════════════════════════
# STEP 2 — Load NHANES IV test data with mortality
# ═══════════════════════════════════════════════════════════════════


def load_nhanes_iv_test():
    """Load NHANES IV 1999–2010 with biomarkers + mortality for validation."""
    from ..load_data import load_nhanes_iv
    from .constants import LEVINE_VALIDATION_CYCLES

    print("\n" + "=" * 72)
    print("  LOADING NHANES IV TEST DATA (1999-2010 with mortality)")
    print("=" * 72)

    df = load_nhanes_iv(cycles=LEVINE_VALIDATION_CYCLES)

    # Merge mortality (2015 vintage where available, else 2019)
    mort_2015 = pd.read_parquet("data/processed/LMF_Files/LMF_all_MORT_2015.parquet")
    mort_2015 = mort_2015[mort_2015.cycle != "NHANES_III"]

    # Use cycles that have 2015 mortality
    df = df.merge(
        mort_2015[["SEQN", "mortstat", "ucod_leading", "permth_exm"]],
        on="SEQN",
        how="inner",
    )
    df = df[df.permth_exm.notna() & df.mortstat.notna()].copy()

    # Age filter
    df = df[(df.age >= 20) & (df.age <= 84)].copy()

    # Create PhenoAge features
    if "glucose_mmol" not in df.columns and "glucose_mgdl" in df.columns:
        df["glucose_mmol"] = df["glucose_mgdl"] / 18.016
    if "log_crp" not in df.columns and "crp_mgdl" in df.columns:
        df["log_crp"] = np.where(df["crp_mgdl"] > 0, np.log(df["crp_mgdl"]), np.nan)
    if "albumin_gL" not in df.columns and "albumin" in df.columns:
        df["albumin_gL"] = df["albumin"] * 10
    if "creat_umol" not in df.columns and "creatinine" in df.columns:
        df["creat_umol"] = df["creatinine"] * 88.4017

    # Rename to match training feature names
    rename_map = {
        "lymphocyte_pct": "lymph",
        "log_crp": "lncrp",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    # Drop incomplete
    feats = [
        "albumin_gL",
        "creat_umol",
        "glucose_mmol",
        "lncrp",
        "lymph",
        "mcv",
        "rdw",
        "alp",
        "wbc",
        "age",
    ]
    df = df.dropna(subset=feats + ["permth_exm", "mortstat"])

    print(f"  Test sample: n={len(df):,d}, deaths={int(df.mortstat.sum()):,d}")
    print(f"  Cycles: {df.cycle.value_counts().to_dict()}" if "cycle" in df.columns else "")
    return df


# ═══════════════════════════════════════════════════════════════════
# STEP 3 — Score and evaluate
# ═══════════════════════════════════════════════════════════════════


def concordance_index(time, event, risk_score):
    """Compute Harrell's C-index."""
    n = len(time)
    concordant = 0
    discordant = 0
    tied = 0
    # Only compare pairs where at least one had an event
    for i in range(n):
        if event[i] == 0:
            continue
        for j in range(n):
            if i == j:
                continue
            if time[j] > time[i]:  # j survived longer than i who died
                if risk_score[j] < risk_score[i]:
                    concordant += 1
                elif risk_score[j] > risk_score[i]:
                    discordant += 1
                else:
                    tied += 1
    total = concordant + discordant + tied
    if total == 0:
        return 0.5
    return (concordant + 0.5 * tied) / total


def fast_cindex(time, event, risk):
    """Fast approximate C-index using sorted comparison."""
    # For large datasets, use lifelines or sksurv if available
    try:
        from lifelines.utils import concordance_index as li_ci

        return li_ci(time, -risk, event)
    except ImportError:
        pass
    try:
        from sksurv.metrics import concordance_index_censored

        c, _, _, _, _ = concordance_index_censored(event.astype(bool), time, risk)
        return c
    except ImportError:
        pass
    # Fallback: sample-based approximation for speed
    idx = np.where(event == 1)[0]
    if len(idx) == 0:
        return 0.5
    np.random.seed(42)
    if len(idx) > 500:
        idx = np.random.choice(idx, 500, replace=False)
    conc = disc = tied = 0
    for i in idx:
        mask = time > time[i]
        r_i = risk[i]
        r_j = risk[mask]
        conc += (r_j < r_i).sum()
        disc += (r_j > r_i).sum()
        tied += (r_j == r_i).sum()
    total = conc + disc + tied
    return (conc + 0.5 * tied) / total if total > 0 else 0.5


def evaluate_model(df, gamma, intercept, betas, label, t_months=120):
    """Score PhenoAge and compute performance metrics."""
    # Compute xb
    xb = np.full(len(df), intercept, dtype=float)
    for feat, coef in betas.items():
        if feat in df.columns:
            xb += coef * df[feat].values

    # Mortality score
    M = compute_mortality_score(xb, gamma=gamma, t_months=t_months)

    # PhenoAge
    M_clipped = np.clip(M, 1e-10, 1 - 1e-10)
    phenoage = (
        PHENOAGE_INTERCEPT + np.log(-PHENOAGE_LOG_SCALE * np.log(1 - M_clipped)) / PHENOAGE_RATE
    )
    phenoage_accel = phenoage - df["age"].values

    # Metrics
    age_corr = np.corrcoef(phenoage, df["age"].values)[0, 1]
    mean_accel = np.nanmean(phenoage_accel)

    # C-index (using xb as risk score — higher xb = higher risk)
    time = df["permth_exm"].values
    event = df["mortstat"].values
    cindex = fast_cindex(time, event, xb)

    # 10-year mortality AUC (if enough follow-up)
    has_10yr = (time >= 120) | (event == 1)
    died_10yr = (event == 1) & (time <= 120)
    if has_10yr.sum() > 100:
        from sklearn.metrics import roc_auc_score

        try:
            auc_10yr = roc_auc_score(died_10yr[has_10yr], M[has_10yr])
        except Exception:
            auc_10yr = np.nan
    else:
        auc_10yr = np.nan

    return {
        "label": label,
        "n": len(df),
        "deaths": int(event.sum()),
        "age_corr": age_corr,
        "mean_accel": mean_accel,
        "cindex": cindex,
        "auc_10yr": auc_10yr,
        "phenoage": phenoage,
        "xb": xb,
    }


def main():
    # Train
    models = train_models()

    # Load test
    test_df = load_nhanes_iv_test()

    # Evaluate each model
    print("\n" + "=" * 72)
    print("  VALIDATION ON NHANES IV")
    print("=" * 72)

    results = []
    for label, (gamma, intercept, betas) in models.items():
        r = evaluate_model(test_df, gamma, intercept, betas, label)
        results.append(r)
        print(f"\n  {label}:")
        print(f"    PhenoAge-Age correlation: {r['age_corr']:.4f}")
        print(f"    Mean PhenoAge acceleration: {r['mean_accel']:.2f} years")
        print(f"    C-index (mortality): {r['cindex']:.4f}")
        print(
            f"    10-year mortality AUC: {r['auc_10yr']:.4f}"
            if not np.isnan(r["auc_10yr"])
            else "    10-year mortality AUC: N/A"
        )

    # Correlation between models' PhenoAge
    print(f"\n  PhenoAge correlation between models:")
    for i, r1 in enumerate(results):
        for r2 in results[i + 1 :]:
            corr = np.corrcoef(r1["phenoage"], r2["phenoage"])[0, 1]
            print(f"    {r1['label']} vs {r2['label']}: r={corr:.4f}")

    # Summary table
    print(f"\n  {'='*72}")
    print(f"  {'Model':<30s} {'Age r':>7s} {'C-idx':>7s} {'AUC':>7s} {'Accel':>7s}")
    print(f"  {'-'*30} {'-'*7} {'-'*7} {'-'*7} {'-'*7}")
    for r in results:
        auc_str = f"{r['auc_10yr']:.4f}" if not np.isnan(r["auc_10yr"]) else "N/A"
        print(
            f"  {r['label']:<30s} {r['age_corr']:>7.4f} {r['cindex']:>7.4f} {auc_str:>7s} {r['mean_accel']:>+7.2f}"
        )


if __name__ == "__main__":
    main()
