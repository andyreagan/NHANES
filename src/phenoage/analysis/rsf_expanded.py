"""Expanded RSF analysis: proper Levine validation set + more features.

1. Matches Levine/Liu validation: NHANES IV 1999-2010, ages 20-84,
   fasting ≥8h, complete PA9 biomarkers, with mortality follow-up.
2. Compares Linear PhenoAge vs RSF on identical held-out data.
3. Expands RSF with additional NHANES features (lipids, BMI, BP, HbA1c, etc.)

Usage:
    uv run python -m src.phenoage.rsf_expanded
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored

from ..constants import (
    FEATURE_COLS,
    LEVINE_BIOAGE,
    LEVINE_GAMMA,
    LEVINE_VALIDATION_CYCLES,
    NHANES_IV_CYCLES,
    PHENOAGE_INTERCEPT,
    PHENOAGE_LOG_SCALE,
    PHENOAGE_RATE,
)
from ..model import compute_mortality_score, compute_phenoage

OUT_DIR = Path("output/phenoage_analysis")

# ── Feature sets ──

# PA9: the 9 PhenoAge biomarkers (in BioAge units) + age
PA9_FEATURES = (
    FEATURE_COLS  # albumin_gL, creat_umol, glucose_mmol, lncrp, lymph, mcv, rdw, alp, wbc, age
)

# Extended: add sex, BMI, systolic/diastolic BP, total cholesterol,
# triglycerides, HDL, LDL, HbA1c, BUN, uric acid, GGT, total bilirubin
EXTENDED_FEATURES = PA9_FEATURES + [
    "sex",  # 0=M, 1=F
    "bmi",  # kg/m²
    "sbp",  # mmHg
    "dbp",  # mmHg
    "total_chol",  # mg/dL
    "triglycerides",  # mg/dL
    "hdl",  # mg/dL
    "hba1c",  # %
    "bun",  # mg/dL
    "uric_acid",  # mg/dL
    "ggt",  # U/L
    "total_bilirubin",  # mg/dL
]


# ═══════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════


def load_nhanes_iv_expanded():
    """Load NHANES IV 1999-2010 with PhenoAge biomarkers + extended features + fasting + mortality."""
    print("  Loading NHANES IV cycles with expanded features...")

    all_dfs = []
    for cycle_label, suffix, start_year in LEVINE_VALIDATION_CYCLES:
        df = _load_cycle_expanded(cycle_label, suffix, start_year)
        if df is not None:
            all_dfs.append(df)

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"  Combined: {len(combined):,d} rows")
    return combined


def _load_cycle_expanded(cycle: str, suffix: str, start_year: int):
    """Load a single NHANES IV cycle with all available features."""
    from ..load_data import _load_xpt

    raw_dir = Path(f"data/raw/{cycle}")
    if not raw_dir.exists():
        # Try alternate naming
        raw_dir = Path(f"data/raw/{start_year}-{start_year+1}")
        if not raw_dir.exists():
            print(f"    ⚠ No data dir for {cycle}")
            return None

    # File name mapping
    file_map = {
        "DEMO": {
            "1999-2000": "DEMO",
            "2001-2002": "DEMO_B",
            "2003-2004": "DEMO_C",
            "2005-2006": "DEMO_D",
            "2007-2008": "DEMO_E",
            "2009-2010": "DEMO_F",
        },
        "BIOPRO": {
            "1999-2000": "LAB18",
            "2001-2002": "L40_B",
            "2003-2004": "L40_C",
            "2005-2006": "BIOPRO_D",
            "2007-2008": "BIOPRO_E",
            "2009-2010": "BIOPRO_F",
        },
        "CBC": {
            "1999-2000": "LAB25",
            "2001-2002": "L25_B",
            "2003-2004": "L25_C",
            "2005-2006": "CBC_D",
            "2007-2008": "CBC_E",
            "2009-2010": "CBC_F",
        },
        "CRP": {
            "1999-2000": "LAB11",
            "2001-2002": "L11_B",
            "2003-2004": "L11_C",
            "2005-2006": "CRP_D",
            "2007-2008": "CRP_E",
            "2009-2010": "CRP_F",
        },
        "GLU": {
            "1999-2000": "LAB10AM",
            "2001-2002": "L10AM_B",
            "2003-2004": "L10AM_C",
            "2005-2006": "GLU_D",
            "2007-2008": "GLU_E",
            "2009-2010": "GLU_F",
        },
        # Additional files for extended features
        "BMX": {
            "1999-2000": "BMX",
            "2001-2002": "BMX_B",
            "2003-2004": "BMX_C",
            "2005-2006": "BMX_D",
            "2007-2008": "BMX_E",
            "2009-2010": "BMX_F",
        },
        "BPX": {
            "1999-2000": "BPX",
            "2001-2002": "BPX_B",
            "2003-2004": "BPX_C",
            "2005-2006": "BPX_D",
            "2007-2008": "BPX_E",
            "2009-2010": "BPX_F",
        },
        "TCHOL": {
            "1999-2000": "LAB13",
            "2001-2002": "L13_B",
            "2003-2004": "L13_C",
            "2005-2006": "TCHOL_D",
            "2007-2008": "TCHOL_E",
            "2009-2010": "TCHOL_F",
        },
        "TRIGLY": {
            "1999-2000": "LAB13AM",
            "2001-2002": "L13AM_B",
            "2003-2004": "L13AM_C",
            "2005-2006": "TRIGLY_D",
            "2007-2008": "TRIGLY_E",
            "2009-2010": "TRIGLY_F",
        },
        "HDL": {
            "1999-2000": "LAB13",
            "2001-2002": "L13_B",
            "2003-2004": "L13_C",
            "2005-2006": "HDL_D",
            "2007-2008": "HDL_E",
            "2009-2010": "HDL_F",
        },
        "GHB": {
            "1999-2000": "LAB10",
            "2001-2002": "L10_B",
            "2003-2004": "L10_C",
            "2005-2006": "GHB_D",
            "2007-2008": "GHB_E",
            "2009-2010": "GHB_F",
        },
    }

    # 1. DEMO
    demo_file = file_map["DEMO"].get(cycle, f"DEMO_{suffix}")
    demo = _load_xpt(raw_dir / f"{demo_file}.XPT")
    if demo is None:
        return None
    df = demo[["SEQN"]].copy()
    for col in ["RIDAGEYR", "RIAGENDR", "WTMEC2YR", "RIDRETH1"]:
        if col in demo.columns:
            df[col] = demo[col]

    # 2. BIOPRO (albumin, creatinine, ALP, BUN, uric acid, GGT, total bilirubin)
    biopro_file = file_map["BIOPRO"].get(cycle)
    biopro = _load_xpt(raw_dir / f"{biopro_file}.XPT") if biopro_file else None
    if biopro is not None:
        # Standard columns
        for col in ["LBXSAL", "LBXSAPSI", "LBXSBU", "LBXSUA", "LBXSGTSI", "LBXSTB"]:
            if col in biopro.columns:
                df = df.merge(biopro[["SEQN", col]], on="SEQN", how="left")
        # Creatinine - handle name variation
        creat_col = "LBXSCR" if "LBXSCR" in biopro.columns else "LBDSCR"
        if creat_col in biopro.columns:
            biopro_cr = biopro[["SEQN", creat_col]].rename(columns={creat_col: "LBXSCR"})
            df = df.merge(biopro_cr, on="SEQN", how="left")

    # 3. CBC
    cbc_file = file_map["CBC"].get(cycle)
    cbc = _load_xpt(raw_dir / f"{cbc_file}.XPT") if cbc_file else None
    if cbc is not None:
        for col in ["LBXWBCSI", "LBXLYPCT", "LBXMCVSI", "LBXRDW"]:
            if col in cbc.columns:
                df = df.merge(cbc[["SEQN", col]], on="SEQN", how="left")

    # 4. CRP
    crp_file = file_map["CRP"].get(cycle)
    crp = _load_xpt(raw_dir / f"{crp_file}.XPT") if crp_file else None
    if crp is not None:
        if "LBXCRP" in crp.columns:
            df = df.merge(crp[["SEQN", "LBXCRP"]], on="SEQN", how="left")

    # 5. GLU (includes fasting hours)
    glu_file = file_map["GLU"].get(cycle)
    glu = _load_xpt(raw_dir / f"{glu_file}.XPT") if glu_file else None
    if glu is not None:
        glu_cols = ["SEQN", "LBXGLU"]
        if "PHAFSTHR" in glu.columns:
            glu_cols.append("PHAFSTHR")
        df = df.merge(glu[glu_cols], on="SEQN", how="left")

    # 6. BMX (BMI)
    bmx_file = file_map["BMX"].get(cycle)
    bmx = _load_xpt(raw_dir / f"{bmx_file}.XPT") if bmx_file else None
    if bmx is not None and "BMXBMI" in bmx.columns:
        df = df.merge(bmx[["SEQN", "BMXBMI"]], on="SEQN", how="left")

    # 7. BPX (blood pressure)
    bpx_file = file_map["BPX"].get(cycle)
    bpx = _load_xpt(raw_dir / f"{bpx_file}.XPT") if bpx_file else None
    if bpx is not None:
        # Average of up to 4 readings
        sbp_cols = [c for c in bpx.columns if c.startswith("BPXSY") and c[-1].isdigit()]
        dbp_cols = [c for c in bpx.columns if c.startswith("BPXDI") and c[-1].isdigit()]
        if sbp_cols:
            bpx["mean_sbp"] = bpx[sbp_cols].mean(axis=1)
        if dbp_cols:
            bpx["mean_dbp"] = bpx[dbp_cols].mean(axis=1)
        bp_merge = ["SEQN"]
        if "mean_sbp" in bpx.columns:
            bp_merge.append("mean_sbp")
        if "mean_dbp" in bpx.columns:
            bp_merge.append("mean_dbp")
        if len(bp_merge) > 1:
            df = df.merge(bpx[bp_merge], on="SEQN", how="left")

    # 8. Total cholesterol
    tchol_file = file_map["TCHOL"].get(cycle)
    tchol = _load_xpt(raw_dir / f"{tchol_file}.XPT") if tchol_file else None
    if tchol is not None and "LBXTC" in tchol.columns:
        df = df.merge(tchol[["SEQN", "LBXTC"]], on="SEQN", how="left")

    # 9. Triglycerides
    trigly_file = file_map["TRIGLY"].get(cycle)
    trigly = _load_xpt(raw_dir / f"{trigly_file}.XPT") if trigly_file else None
    if trigly is not None and "LBXTR" in trigly.columns:
        df = df.merge(trigly[["SEQN", "LBXTR"]], on="SEQN", how="left")

    # 10. HDL
    hdl_file = file_map["HDL"].get(cycle)
    # For 1999-2004, HDL is in the same file as total cholesterol (LAB13/L13)
    if hdl_file == tchol_file and tchol is not None:
        hdl = tchol
    else:
        hdl = _load_xpt(raw_dir / f"{hdl_file}.XPT") if hdl_file else None
    if hdl is not None and "LBDHDL" in hdl.columns:
        if "LBDHDL" not in df.columns:
            df = df.merge(hdl[["SEQN", "LBDHDL"]], on="SEQN", how="left")
    elif hdl is not None and "LBDHDD" in hdl.columns:
        df = df.merge(
            hdl[["SEQN", "LBDHDD"]].rename(columns={"LBDHDD": "LBDHDL"}), on="SEQN", how="left"
        )

    # 11. HbA1c
    ghb_file = file_map["GHB"].get(cycle)
    ghb = _load_xpt(raw_dir / f"{ghb_file}.XPT") if ghb_file else None
    if ghb is not None and "LBXGH" in ghb.columns:
        df = df.merge(ghb[["SEQN", "LBXGH"]], on="SEQN", how="left")

    # ── Harmonize names ──
    df["cycle"] = cycle
    rename = {
        "RIDAGEYR": "age",
        "RIAGENDR": "sex_code",
        "WTMEC2YR": "exam_weight",
        "LBXSAL": "albumin",  # g/dL
        "LBXSCR": "creatinine",  # mg/dL
        "LBXSAPSI": "alp",  # U/L
        "LBXWBCSI": "wbc",  # 10^9/L
        "LBXLYPCT": "lymph",  # %
        "LBXMCVSI": "mcv",  # fL
        "LBXRDW": "rdw",  # %
        "LBXCRP": "crp_mgdl",  # mg/dL
        "LBXGLU": "glucose_mgdl",  # mg/dL
        "BMXBMI": "bmi",
        "mean_sbp": "sbp",
        "mean_dbp": "dbp",
        "LBXTC": "total_chol",
        "LBXTR": "triglycerides",
        "LBDHDL": "hdl",
        "LBXGH": "hba1c",
        "LBXSBU": "bun",
        "LBXSUA": "uric_acid",
        "LBXSGTSI": "ggt",
        "LBXSTB": "total_bilirubin",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    # Derived features
    df["glucose_mmol"] = df.get("glucose_mgdl", pd.Series(dtype=float)) / 18.016
    if "crp_mgdl" in df.columns:
        df["lncrp"] = np.where(df.crp_mgdl > 0, np.log(df.crp_mgdl), np.nan)
    if "albumin" in df.columns:
        df["albumin_gL"] = df.albumin * 10
    if "creatinine" in df.columns:
        df["creat_umol"] = df.creatinine * 88.4017
    if "sex_code" in df.columns:
        df["sex"] = (df.sex_code == 2).astype(float)

    return df


# ═══════════════════════════════════════════════════════════════════
# VALIDATION SET CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════


def build_validation_set(df: pd.DataFrame, fasting_hours: int = 8):
    """Build Levine/Liu-style validation set."""
    # Merge mortality
    mort = pd.read_parquet("data/processed/LMF_Files/LMF_all_MORT_2015.parquet")
    mort = mort[mort.cycle != "NHANES_III"]
    df = df.merge(
        mort[["SEQN", "mortstat", "ucod_leading", "permth_exm"]],
        on="SEQN",
        how="inner",
    )

    # Age 20-84
    df = df[(df.age >= 20) & (df.age <= 84)].copy()

    # Fasting ≥ 8 hours
    if "PHAFSTHR" in df.columns:
        n_before = len(df)
        df = df[df.PHAFSTHR >= fasting_hours].copy()
        print(f"  Fasting ≥{fasting_hours}h filter: {n_before:,d} → {len(df):,d}")
    else:
        print("  ⚠ No fasting hours variable — skipping fasting filter")

    # Complete PA9 biomarkers
    pa9_raw = [
        "albumin_gL",
        "creat_umol",
        "glucose_mmol",
        "lncrp",
        "lymph",
        "mcv",
        "rdw",
        "alp",
        "wbc",
    ]
    df = df.dropna(subset=pa9_raw + ["age", "permth_exm", "mortstat"])

    # Mortality definition
    df.loc[df.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    df.loc[(df.mortstat == 1) & (df.permth_exm > 240), "mortstat"] = 0
    df.loc[df.permth_exm > 240, "permth_exm"] = 240

    print(f"  Validation set: n={len(df):,d}, deaths={int(df.mortstat.sum()):,d}")
    return df


# ═══════════════════════════════════════════════════════════════════
# MODEL COMPARISON
# ═══════════════════════════════════════════════════════════════════


def compute_linear_risk(df: pd.DataFrame):
    """Compute linear PhenoAge risk score (xb)."""
    xb = np.full(len(df), LEVINE_BIOAGE["intercept"])
    for feat in FEATURE_COLS:
        xb += LEVINE_BIOAGE[feat] * df[feat].values
    return xb


def fit_and_eval_rsf(X_train, y_train, X_test, y_test, label="RSF", n_trees=500, min_leaf=15):
    """Fit RSF and return C-index on train and test."""
    print(
        f"\n  Fitting {label}: {X_train.shape[1]} features, "
        f"n_train={len(X_train)}, n_test={len(X_test)}"
    )

    rsf = RandomSurvivalForest(
        n_estimators=n_trees,
        min_samples_leaf=min_leaf,
        n_jobs=-1,
        random_state=42,
    )
    rsf.fit(X_train, y_train)

    # Train C-index
    risk_train = rsf.predict(X_train)
    c_train, _, _, _, _ = concordance_index_censored(y_train["event"], y_train["time"], risk_train)

    # Test C-index
    risk_test = rsf.predict(X_test)
    c_test, _, _, _, _ = concordance_index_censored(y_test["event"], y_test["time"], risk_test)

    print(f"  {label} C-index: train={c_train:.4f}, test={c_test:.4f}")
    return rsf, c_train, c_test, risk_test


def make_structured_y(df):
    """Create structured array for sksurv."""
    return np.array(
        [(bool(e), t) for e, t in zip(df.mortstat.values, df.permth_exm.values)],
        dtype=[("event", bool), ("time", float)],
    )


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("  EXPANDED RSF ANALYSIS: PROPER VALIDATION + MORE FEATURES")
    print("=" * 72)

    # ── Load training data (NHANES III) ──
    print("\n  Loading NHANES III training data...")
    from .validate_models import _load_train

    train_df = _load_train()
    train_df.loc[train_df.ucod_leading.isin([4, 8, 10]), "mortstat"] = 0
    train_df.loc[(train_df.mortstat == 1) & (train_df.permth_exm > 240), "mortstat"] = 0
    train_df.loc[train_df.permth_exm > 240, "permth_exm"] = 240
    train_df["sex"] = (train_df["HSSEX"] == 2).astype(float)

    # Map NHANES III extra variables to canonical names
    iii_rename = {
        "BMPBMI": "bmi",
        "PEPMNK1R": "sbp",  # systolic BP (1st reading)
        "PEPMNK5R": "dbp",  # diastolic BP (5th reading)
        "TCP": "total_chol",  # total cholesterol mg/dL
        "TGP": "triglycerides",  # triglycerides mg/dL
        "HDP": "hdl",  # HDL mg/dL
        "GHP": "hba1c",  # glycohemoglobin %
        "BUP": "bun",  # blood urea nitrogen mg/dL
        "UAP": "uric_acid",  # uric acid mg/dL
        "GGPSI": "ggt",  # GGT U/L (SI units)
        "TBP": "total_bilirubin",  # total bilirubin mg/dL
    }
    for old, new in iii_rename.items():
        if old in train_df.columns:
            train_df[new] = train_df[old]

    print(f"  Training: n={len(train_df):,d}, deaths={int(train_df.mortstat.sum()):,d}")

    # ── Load NHANES IV with expanded features ──
    print()
    test_raw = load_nhanes_iv_expanded()

    # Build validation set (Levine-style: fasting ≥8h)
    test_df = build_validation_set(test_raw, fasting_hours=8)

    # Also build the relaxed set (no fasting filter) for comparison
    test_relaxed = build_validation_set(test_raw.copy(), fasting_hours=0)

    # ── Print feature availability ──
    print(f"\n  Extended feature availability (validation set n={len(test_df)}):")
    for feat in EXTENDED_FEATURES:
        n_avail = test_df[feat].notna().sum() if feat in test_df.columns else 0
        pct = n_avail / len(test_df) * 100
        print(f"    {feat:20s}: {n_avail:6,d} ({pct:5.1f}%)")

    # ── Prepare survival data ──
    y_train = make_structured_y(train_df)
    y_test = make_structured_y(test_df)

    # ── Model 1: Linear PhenoAge ──
    print("\n" + "=" * 72)
    print("  MODEL COMPARISONS")
    print("=" * 72)

    xb_test = compute_linear_risk(test_df)
    c_linear, _, _, _, _ = concordance_index_censored(y_test["event"], y_test["time"], xb_test)
    print(f"\n  Linear PhenoAge C-index (test): {c_linear:.4f}")

    # ── Model 2: RSF with PA9 features only ──
    X_train_pa9 = train_df[PA9_FEATURES].copy()
    X_test_pa9 = test_df[PA9_FEATURES].copy()
    _, c_pa9_train, c_pa9_test, _ = fit_and_eval_rsf(
        X_train_pa9, y_train, X_test_pa9, y_test, "RSF (PA9 only)"
    )

    # ── Model 3: RSF with PA9 + sex ──
    pa9_sex = PA9_FEATURES + ["sex"]
    X_train_sex = train_df[pa9_sex].copy()
    X_test_sex = test_df[pa9_sex].copy()
    _, c_sex_train, c_sex_test, _ = fit_and_eval_rsf(
        X_train_sex, y_train, X_test_sex, y_test, "RSF (PA9 + sex)"
    )

    # ── Model 4: RSF with extended features ──
    # Find features available in BOTH train and test
    # Training data doesn't have all extended features — need to check
    available_extended = []
    for feat in EXTENDED_FEATURES:
        in_train = feat in train_df.columns and train_df[feat].notna().sum() > len(train_df) * 0.5
        in_test = feat in test_df.columns and test_df[feat].notna().sum() > len(test_df) * 0.5
        if in_train and in_test:
            available_extended.append(feat)
        elif in_test and not in_train:
            print(f"  ⚠ {feat}: available in test but not train (skipping)")

    # For features only in NHANES IV, train on NHANES IV early cycles, test on later
    # But for now, only use features available in both train and test
    print(f"\n  Available extended features ({len(available_extended)}):")
    for f in available_extended:
        print(f"    {f}")

    if len(available_extended) > len(PA9_FEATURES):
        X_train_ext = train_df[available_extended].dropna()
        y_train_ext = make_structured_y(train_df.loc[X_train_ext.index])
        X_test_ext = test_df[available_extended].dropna()
        y_test_ext = make_structured_y(test_df.loc[X_test_ext.index])
        _, c_ext_train, c_ext_test, _ = fit_and_eval_rsf(
            X_train_ext,
            y_train_ext,
            X_test_ext,
            y_test_ext,
            f"RSF (extended, {len(available_extended)} feats)",
        )
    else:
        c_ext_train = c_ext_test = float("nan")

    # ── Model 5: RSF trained on NHANES IV itself (cross-validated) ──
    # Use 1999-2006 as train, 2007-2010 as test (temporal split)
    print(f"\n  Temporal split: train on 1999-2006, test on 2007-2010...")
    early_cycles = ["1999-2000", "2001-2002", "2003-2004", "2005-2006"]
    late_cycles = ["2007-2008", "2009-2010"]
    iv_train = test_df[test_df.cycle.isin(early_cycles)].copy()
    iv_test = test_df[test_df.cycle.isin(late_cycles)].copy()
    print(f"  IV-train: n={len(iv_train):,d}, IV-test: n={len(iv_test):,d}")

    # Extended features available in NHANES IV
    iv_extended = []
    for feat in EXTENDED_FEATURES:
        in_iv = feat in test_df.columns and test_df[feat].notna().sum() > len(test_df) * 0.5
        if in_iv:
            iv_extended.append(feat)

    iv_train_ext = iv_train[iv_extended].dropna()
    y_iv_train = make_structured_y(iv_train.loc[iv_train_ext.index])
    iv_test_ext = iv_test[iv_extended].dropna()
    y_iv_test = make_structured_y(iv_test.loc[iv_test_ext.index])

    print(f"  IV extended features: {len(iv_extended)}")
    print(f"  IV-train (complete): n={len(iv_train_ext):,d}")
    print(f"  IV-test (complete): n={len(iv_test_ext):,d}")

    # PA9 only on IV temporal split
    iv_train_pa9 = iv_train[PA9_FEATURES].dropna()
    y_iv_train_pa9 = make_structured_y(iv_train.loc[iv_train_pa9.index])
    iv_test_pa9 = iv_test[PA9_FEATURES].dropna()
    y_iv_test_pa9 = make_structured_y(iv_test.loc[iv_test_pa9.index])

    _, c_iv_pa9_train, c_iv_pa9_test, _ = fit_and_eval_rsf(
        iv_train_pa9, y_iv_train_pa9, iv_test_pa9, y_iv_test_pa9, "RSF IV temporal (PA9)"
    )

    if len(iv_extended) > len(PA9_FEATURES):
        _, c_iv_ext_train, c_iv_ext_test, _ = fit_and_eval_rsf(
            iv_train_ext,
            y_iv_train,
            iv_test_ext,
            y_iv_test,
            f"RSF IV temporal ({len(iv_extended)} feats)",
        )
    else:
        c_iv_ext_train = c_iv_ext_test = float("nan")

    # Linear on IV temporal test
    xb_iv_test = compute_linear_risk(iv_test)
    c_iv_linear, _, _, _, _ = concordance_index_censored(
        y_iv_test_pa9["event"], y_iv_test_pa9["time"], xb_iv_test[: len(iv_test_pa9)]
    )

    # ── Summary ──
    print("\n" + "=" * 72)
    print("  SUMMARY: C-INDEX COMPARISON")
    print("=" * 72)
    print(f"  {'Model':<45s} {'Train':>7s} {'Test':>7s}")
    print(f"  {'-'*45} {'-'*7} {'-'*7}")
    print(f"  {'Linear PhenoAge (III→IV)':<45s} {'—':>7s} {c_linear:>7.4f}")
    print(f"  {'RSF PA9 only (III→IV)':<45s} {c_pa9_train:>7.4f} {c_pa9_test:>7.4f}")
    print(f"  {'RSF PA9+sex (III→IV)':<45s} {c_sex_train:>7.4f} {c_sex_test:>7.4f}")
    if not np.isnan(c_ext_test):
        n_ext = len(available_extended)
        print(
            f"  {f'RSF extended {n_ext} feats (III→IV)':<45s} {c_ext_train:>7.4f} {c_ext_test:>7.4f}"
        )
    print(f"  {'—'*45} {'—'*7} {'—'*7}")
    print(f"  {'Linear PhenoAge (IV temporal)':<45s} {'—':>7s} {c_iv_linear:>7.4f}")
    print(
        f"  {'RSF PA9 (IV temporal: 99-06→07-10)':<45s} {c_iv_pa9_train:>7.4f} {c_iv_pa9_test:>7.4f}"
    )
    if not np.isnan(c_iv_ext_test):
        print(
            f"  {f'RSF {len(iv_extended)} feats (IV temporal: 99-06→07-10)':<45s} {c_iv_ext_train:>7.4f} {c_iv_ext_test:>7.4f}"
        )

    print(
        f"\n  Validation set: n={len(test_df):,d} " f"(fasting≥8h, ages 20-84, NHANES IV 1999-2010)"
    )
    print(f"  Levine reported: n=6,209 (different exclusion criteria)")


if __name__ == "__main__":
    main()
