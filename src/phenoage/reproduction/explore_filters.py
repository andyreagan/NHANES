"""Exploration script to match Levine 2018 NHANES III training sample.

Levine 2018 reports n=9,926 with 1,846 deaths for the NHANES III training set.
This script lets you iterate on filter configurations and compare the resulting
Gompertz PH coefficients to Levine's published values.

Usage:
    uv run python -m src.phenoage.explore_filters

Adjust the CONFIGS list at the bottom to try different filter combinations.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from ..constants import (
    FEATURE_COLS,
    LEVINE_BIOAGE,
    LEVINE_COEFFICIENTS,
    LEVINE_GAMMA,
)
from ..constants import gompertz_nll as _gompertz_nll

# Column positions in NHANES III lab.dat (1-indexed inclusive)
LAB_COLSPECS = [
    ("SEQN", 1, 5),
    ("HSSEX", 15, 15),
    ("HSAGEIR", 16, 17),
    ("SDPPHASE", 40, 40),
    ("WTPFEX6", 59, 67),
    ("WTPFSD6", 95, 103),
    ("PHPFAST", 1263, 1267),
    ("WCPSI", 1278, 1282),
    ("LMPPCNT", 1283, 1287),
    ("MVPSI", 1340, 1344),
    ("RWP", 1360, 1364),
    ("CRP", 1667, 1671),
    ("SGP", 1758, 1760),
    ("CEP", 1784, 1787),
    ("APPSI", 1835, 1838),
    ("AMP", 1846, 1848),
    ("G1P", 1866, 1870),
]


# ============================================================================
# Data loading helpers
# ============================================================================
def load_nhanes_iii_raw() -> pd.DataFrame:
    """Parse NHANES III lab.dat with sentinel-value cleanup."""
    colspecs = [(s - 1, e) for _, s, e in LAB_COLSPECS]
    names = [n for n, _, _ in LAB_COLSPECS]
    df = pd.read_fwf(
        "data/raw/NHANES_III/lab.dat",
        colspecs=colspecs,
        header=None,
        names=names,
        na_values=[".", ""],
    )
    # Sentinel values  (see NHANES III lab.sas documentation)
    for col in ["CRP"]:
        df.loc[df[col] >= 8888.80, col] = np.nan
    for col in ["CEP"]:
        df.loc[df[col] >= 8888, col] = np.nan
    for col in ["SGP", "AMP"]:
        df.loc[df[col] >= 888, col] = np.nan
    for col in ["APPSI"]:
        df.loc[df[col] >= 8888, col] = np.nan
    for col in ["LMPPCNT", "WCPSI", "MVPSI", "RWP"]:
        df.loc[df[col] >= 88888, col] = np.nan
    for col in ["G1P"]:
        df.loc[df[col] >= 888, col] = np.nan
    df.loc[df.PHPFAST >= 88888, "PHPFAST"] = np.nan
    return df


def load_mortality() -> pd.DataFrame:
    """Load NHANES III linked mortality file."""
    from ..load_data import load_all_mortality

    return load_all_mortality(["NHANES_III"])


# ============================================================================
# Core: build sample from a config dict and (optionally) fit Gompertz
# ============================================================================
def build_sample(
    raw: pd.DataFrame,
    mort: pd.DataFrame,
    *,
    fasting_min: float | None = 8.0,
    creat_calibrate: bool = True,
    outlier_sd: float | None = 5.0,
    mort_window_months: int = 240,
    mort_exclude_age_related: bool = True,
    crp_transform: str = "log",  # "log" or "log1p"
    glucose_source: str = "SGP",  # "SGP" (serum) or "G1P" (fasting plasma)
    phase: int | str = "both",  # 1, 2, or "both"
    age_min: int = 20,
    age_max: int = 84,
) -> pd.DataFrame:
    """Apply a filter configuration and return the analysis-ready sample.

    Returns a DataFrame with the 10 model features renamed to FEATURE_COLS,
    plus ``permth_exm`` and ``mortstat`` for survival modelling.
    """
    df = raw.merge(
        mort[["SEQN", "eligstat", "mortstat", "permth_exm", "ucod_leading"]],
        on="SEQN",
        how="left",
    )

    # Outlier removal is done on the broadest population (age ≥ 20, no upper cap)
    # BEFORE subsetting to age_max, matching the BioAge R code ordering.
    df = df[df.HSAGEIR >= age_min].copy()

    # Creatinine assay calibration (Selvin et al.)
    if creat_calibrate:
        df["CEP"] = df["CEP"] * 0.960 - 0.184

    # 5-SD outlier removal, per gender
    if outlier_sd is not None:
        bio = ["AMP", "APPSI", "CRP", "CEP", glucose_source, "LMPPCNT", "WCPSI", "MVPSI", "RWP"]
        for col in bio:
            for g in [1, 2]:
                mask = df.HSSEX == g
                m, s = df.loc[mask, col].mean(), df.loc[mask, col].std()
                if pd.notna(m) and pd.notna(s) and s > 0:
                    df.loc[
                        mask & ((df[col] > m + outlier_sd * s) | (df[col] < m - outlier_sd * s)),
                        col,
                    ] = np.nan

    # Age upper bound
    df = df[df.HSAGEIR <= age_max].copy()

    # Phase filter
    if phase != "both":
        df = df[df.SDPPHASE == int(phase)]

    # Fasting filter
    if fasting_min is not None:
        df = df[df.PHPFAST >= fasting_min]

    # Mortality processing
    if mort_exclude_age_related:
        # ucod_leading 4=diabetes, 8=Alzheimer's, 10=external causes
        df.loc[df.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    if mort_window_months < 9999:
        df.loc[(df.mortstat == 1) & (df.permth_exm > mort_window_months), "mortstat"] = 0
        df["permth_exm"] = np.minimum(df["permth_exm"], mort_window_months)

    # Derive model features (BioAge units)
    df["albumin_gL"] = df["AMP"] * 10
    df["creat_umol"] = df["CEP"] * 88.4017
    df["glucose_mmol"] = df[glucose_source] * 0.0555
    if crp_transform == "log":
        df["lncrp"] = np.where(df["CRP"] > 0, np.log(df["CRP"]), np.nan)
    else:
        df["lncrp"] = np.log(1 + df["CRP"])
    df.loc[df["CRP"].isna(), "lncrp"] = np.nan

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

    # Complete cases
    needed = FEATURE_COLS + ["permth_exm", "mortstat"]
    return df.dropna(subset=needed).copy()


def fit_and_compare(sample: pd.DataFrame, label: str = "") -> dict:
    """Fit Gompertz on *sample* and compare coefficients to Levine's.

    Returns a summary dict with n, deaths, gamma, per-feature fitted/levine/pct,
    and an overall ``avg_pct_diff``.
    """
    X = sample[FEATURE_COLS].values.astype(float)
    time = sample["permth_exm"].values.astype(float)
    event = sample["mortstat"].values.astype(float)

    # Warm-start from Levine's values
    x0 = np.zeros(2 + len(FEATURE_COLS))
    x0[0] = np.log(LEVINE_GAMMA)
    x0[1] = LEVINE_BIOAGE["intercept"]
    for i, col in enumerate(FEATURE_COLS):
        x0[2 + i] = LEVINE_BIOAGE[col]

    res = minimize(
        _gompertz_nll,
        x0,
        args=(X, time, event),
        method="L-BFGS-B",
        options={"maxiter": 10_000, "ftol": 1e-15},
    )

    gamma = np.exp(res.x[0])
    intercept = res.x[1]
    betas = res.x[2:]
    n = len(sample)
    deaths = int(event.sum())

    rows = []
    for i, col in enumerate(FEATURE_COLS):
        lev = LEVINE_BIOAGE[col]
        pct = abs(betas[i] - lev) / abs(lev) * 100 if lev != 0 else 0
        rows.append({"feature": col, "fitted": betas[i], "levine": lev, "pct_diff": pct})
    avg_pct = np.mean([r["pct_diff"] for r in rows])

    return {
        "label": label,
        "n": n,
        "deaths": deaths,
        "gamma": gamma,
        "gamma_pct": abs(gamma - LEVINE_GAMMA) / LEVINE_GAMMA * 100,
        "intercept": intercept,
        "avg_pct_diff": avg_pct,
        "coefficients": rows,
        "converged": res.success,
    }


def print_result(r: dict) -> None:
    """Pretty-print the output of ``fit_and_compare``."""
    print(f"\n{'=' * 72}")
    print(f"  {r['label']}")
    print(f"  n = {r['n']:,d} (target 9,926, diff {r['n']-9926:+,d})")
    print(f"  deaths = {r['deaths']:,d} (target 1,846, diff {r['deaths']-1846:+,d})")
    print(f"  gamma = {r['gamma']:.7f}  (Levine {LEVINE_GAMMA:.7f}, " f"{r['gamma_pct']:.1f}% off)")
    print(f"  avg |%diff| across betas = {r['avg_pct_diff']:.1f}%")
    print(f"  {'feature':<14s} {'fitted':>12s} {'levine':>12s} {'%diff':>8s}")
    print(f"  {'-'*14} {'-'*12} {'-'*12} {'-'*8}")
    for c in r["coefficients"]:
        flag = " ✓" if c["pct_diff"] < 5 else " <<<" if c["pct_diff"] > 30 else ""
        print(
            f"  {c['feature']:<14s} {c['fitted']:>12.7f} {c['levine']:>12.7f} "
            f"{c['pct_diff']:>7.1f}%{flag}"
        )
    print(f"{'=' * 72}")


# ============================================================================
# Main: define configs and run
# ============================================================================
def main():
    print("Loading NHANES III lab data...")
    raw = load_nhanes_iii_raw()
    print(f"  {raw.shape[0]:,d} rows")

    print("Loading mortality data...")
    mort = load_mortality()

    # ------------------------------------------------------------------
    # CONFIGS — edit this list to explore different filter combinations.
    #
    # Each entry is (label, kwargs_dict).  The kwargs are passed to
    # ``build_sample``.  Defaults (if a key is omitted):
    #   fasting_min=8.0, creat_calibrate=True, outlier_sd=5.0,
    #   mort_window_months=240, mort_exclude_age_related=True,
    #   crp_transform="log", glucose_source="SGP", phase="both",
    #   age_min=20, age_max=84
    # ------------------------------------------------------------------
    CONFIGS = [
        # -- Best coefficient match so far --
        (
            "A: fast>=8, 5SD, 20yr excl, log(crp), cal",
            dict(
                fasting_min=8,
                crp_transform="log",
            ),
        ),
        # -- Closest to n=9,926 (fasting ~7.2h) --
        (
            "B: fast>=7.2, 5SD, 20yr excl, log(crp), cal",
            dict(
                fasting_min=7.2,
                crp_transform="log",
            ),
        ),
        # -- No fasting baseline --
        (
            "C: no fast, 5SD, 20yr excl, log(crp), cal",
            dict(
                fasting_min=None,
                crp_transform="log",
            ),
        ),
        # -- log1p comparison --
        (
            "D: fast>=8, 5SD, 20yr excl, log1p(crp), cal",
            dict(
                fasting_min=8,
                crp_transform="log1p",
            ),
        ),
        # -- Simplest: 10yr, no processing --
        (
            "E: fast>=8, 10yr, no outlier, no cal, log(crp)",
            dict(
                fasting_min=8,
                crp_transform="log",
                mort_window_months=120,
                mort_exclude_age_related=False,
                creat_calibrate=False,
                outlier_sd=None,
            ),
        ),
    ]

    results = []
    for label, kwargs in CONFIGS:
        sample = build_sample(raw, mort, **kwargs)
        r = fit_and_compare(sample, label)
        print_result(r)
        results.append(r)

    # Summary table
    print(f"\n{'=' * 72}")
    print(f"  SUMMARY  (sorted by avg coefficient %diff)")
    print(f"{'=' * 72}")
    results.sort(key=lambda x: x["avg_pct_diff"])
    print(f"  {'avg%':>6s} {'γ%':>5s} {'n':>7s} {'D':>6s} | label")
    print(f"  {'-'*6} {'-'*5} {'-'*7} {'-'*6}-+{'-'*40}")
    for r in results:
        print(
            f"  {r['avg_pct_diff']:>5.1f}% {r['gamma_pct']:>4.1f}% "
            f"{r['n']:>7,d} {r['deaths']:>6,d} | {r['label']}"
        )


if __name__ == "__main__":
    main()
