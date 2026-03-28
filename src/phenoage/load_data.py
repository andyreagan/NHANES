"""Data loading utilities for PhenoAge analysis.

Handles:
1. Downloading and parsing NHANES III fixed-width lab data
2. Downloading and loading NHANES IV (continuous) XPT files
3. Downloading and parsing linked mortality files for all cycles
4. Harmonizing variable names across all cycles into a common schema

Downloads use the Makefile rules (same approach as src/m3s_variables/main.py):
  - NHANES IV XPTs:  `make data/raw/{cycle}/{FILE}.XPT`
  - NHANES III .dat: `make data/raw/NHANES_III/lab.dat`
  - Mortality .dat:  `make data/raw/LMF_Files/{FILE}.dat`
"""

from pathlib import Path
from subprocess import run

import numpy as np
import pandas as pd

from .constants import (
    CRP_MGL_TO_MGDL,
    CRP_UNAVAILABLE_CYCLES,
    GLUCOSE_MGDL_TO_MMOL,
    LEVINE_VALIDATION_CYCLES,
    MORTALITY_FILE_COLSPECS,
    MORTALITY_FILE_COLUMNS,
    MORTALITY_FILES,
    NHANES_III_LAB_COLSPECS,
    NHANES_III_TO_PHENOAGE,
    NHANES_IV_CYCLES,
    NHANES_IV_TO_PHENOAGE,
    PHENOAGE_FILE_NAMES,
    PHENOAGE_NHANES_IV_VARS,
    PHENOAGE_VAR_MAPPING,
    get_file_name,
)

DATA_DIR = Path("data")
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"


def make_file(filename: Path) -> bool:
    """Ensure a raw data file exists, downloading via Make if needed.

    This mirrors the pattern in src/m3s_variables/main.py — it delegates
    downloading to the Makefile rules (data/raw/%.XPT, data/raw/NHANES_III/%.dat,
    data/raw/LMF_Files/%.dat).
    """
    if not filename.exists():
        filename.parent.mkdir(parents=True, exist_ok=True)
        run(["make", str(filename)])
    if filename.exists():
        if filename.stat().st_size < 100:
            return False
        return True
    return False


# =============================================================================
# Linked Mortality Files
# =============================================================================


PROCESSED_MORTALITY_PATH = PROCESSED_DIR / "LMF_Files" / "LMF_all_MORT_2019.parquet"


def load_all_mortality(cycles: list[str] | None = None) -> pd.DataFrame:
    """Load linked mortality data for specified cycles.

    First tries the pre-processed parquet at
    ``data/processed/LMF_Files/LMF_all_MORT_2019.parquet``
    (produced by ``make mort_data_processed`` / ``src.parse_lmf.parse_lmf``).
    Falls back to parsing raw .dat files if the parquet is missing.

    Args:
        cycles: List of cycle labels (e.g., ["1999-2000", "NHANES_III"]).
                If None, loads all available cycles.

    Returns:
        DataFrame with SEQN, mortality variables, and 'cycle' column.
    """
    if cycles is None:
        cycles = list(MORTALITY_FILES.keys())

    # Try pre-processed parquet first
    if PROCESSED_MORTALITY_PATH.exists():
        print(f"  Loading pre-processed mortality from {PROCESSED_MORTALITY_PATH}")
        df = pd.read_parquet(PROCESSED_MORTALITY_PATH)
        # Filter to requested cycles
        df = df[df.cycle.isin(cycles)].copy()
        for cycle in cycles:
            cycle_df = df[df.cycle == cycle]
            if cycle_df.shape[0] == 0:
                print(f"  ⚠ No mortality data for {cycle} in processed file")
            else:
                print(
                    f"  Loaded mortality for {cycle}: {cycle_df.shape[0]:,d} rows, "
                    f"{(cycle_df.mortstat == 1).sum():,d} deaths"
                )
        if df.shape[0] > 0:
            return df.reset_index(drop=True)
        print("  ⚠ Processed file empty for requested cycles, falling back to raw .dat files")

    # Fall back to parsing raw .dat files
    all_dfs: list[pd.DataFrame] = []
    for cycle in cycles:
        if cycle not in MORTALITY_FILES:
            print(f"  ⚠ No mortality file for {cycle}")
            continue
        fname = MORTALITY_FILES[cycle]
        dest = RAW_DIR / "LMF_Files" / fname

        if not make_file(dest):
            print(f"  ⚠ Could not download mortality file for {cycle}")
            continue

        df = pd.read_fwf(
            dest,
            colspecs=MORTALITY_FILE_COLSPECS,
            skiprows=0,
            header=None,
            na_values=["."],
        )
        df.columns = MORTALITY_FILE_COLUMNS
        df["cycle"] = cycle

        # For NHANES III, derive the year from the cycle label
        if cycle == "NHANES_III":
            df["year"] = 1991  # midpoint of 1988-1994
        else:
            df["year"] = int(cycle.split("-")[0])

        print(
            f"  Loaded mortality for {cycle}: {df.shape[0]:,d} rows, "
            f"{(df.mortstat == 1).sum():,d} deaths"
        )
        all_dfs.append(df)

    if not all_dfs:
        raise ValueError("No mortality data loaded!")
    return pd.concat(all_dfs, ignore_index=True)


# =============================================================================
# NHANES III data
# =============================================================================


def load_nhanes_iii() -> pd.DataFrame:
    """Download and parse NHANES III lab data for PhenoAge training.

    NHANES III (1988-1994) distributes data as large fixed-width .dat files
    with companion .sas layout programs — a different format from continuous
    NHANES (1999+) which uses per-topic .XPT SAS transport files.

    The lab.dat file (58 MB, ~29k rows, record length 1977) contains all
    laboratory results in a single file.  We parse only the columns needed
    for PhenoAge using the column specs derived from lab.sas.

    Downloads via: ``make data/raw/NHANES_III/lab.dat``

    Returns:
        DataFrame with harmonized PhenoAge variable names.
    """
    lab_dat = RAW_DIR / "NHANES_III" / "lab.dat"

    if not make_file(lab_dat):
        raise FileNotFoundError(
            f"Could not obtain NHANES III lab data at {lab_dat}.  "
            "Try running: make data/raw/NHANES_III/lab.dat"
        )

    print(f"  Parsing NHANES III lab.dat ({lab_dat.stat().st_size / 1e6:.1f} MB)...")

    # Build colspecs for pd.read_fwf (0-indexed, half-open intervals)
    colspecs = [(start - 1, end) for _, start, end in NHANES_III_LAB_COLSPECS]
    names = [name for name, _, _ in NHANES_III_LAB_COLSPECS]

    df = pd.read_fwf(
        lab_dat,
        colspecs=colspecs,
        header=None,
        names=names,
        na_values=[".", ""],
    )
    print(f"  NHANES III lab data: {df.shape[0]:,d} rows, {df.shape[1]} columns")

    # NHANES III uses sentinel values for missing/blank data:
    # 88888, 8888, 888 = blank/missing; replace with NaN
    biomarker_cols = [
        "WCP",
        "WCPSI",
        "LMPPCNT",
        "MVPSI",
        "RWP",
        "CRP",
        "SGP",
        "CEP",
        "APPSI",
        "AMP",
        "G1P",
        "GHP",
    ]
    for col in biomarker_cols:
        if col in df.columns:
            n_sentinel = (df[col] >= 888).sum()
            if n_sentinel > 0:
                df.loc[df[col] >= 888, col] = np.nan
                print(f"    Replaced {n_sentinel:,d} sentinel values (>=888) in {col}")

    # Rename to canonical PhenoAge names
    rename_map = {k: v for k, v in NHANES_III_TO_PHENOAGE.items() if k in df.columns}
    df = df.rename(columns=rename_map)

    # Unit conversions:
    # Glucose: mg/dL → mmol/L (Levine model expects mmol/L)
    if "glucose_mgdl" in df.columns:
        df["glucose_mmol"] = df["glucose_mgdl"] * GLUCOSE_MGDL_TO_MMOL

    # CRP: NHANES III CRP is in mg/dL — compute log directly
    # (Levine model expects log(CRP in mg/dL))
    if "crp_mgdl" in df.columns:
        df["log_crp"] = np.where(df["crp_mgdl"] > 0, np.log(df["crp_mgdl"]), np.nan)

    df["cycle"] = "NHANES_III"
    df["year"] = 1991

    # Keep only relevant columns
    keep_cols = ["SEQN", "cycle", "year"] + [
        v for v in NHANES_III_TO_PHENOAGE.values() if v in df.columns
    ]
    # Add derived columns
    for col in ["glucose_mmol", "log_crp"]:
        if col in df.columns:
            keep_cols.append(col)
    # Also keep SDPPHASE, SDPPSU6, SDPSTRA6 for survey design
    for col in ["SDPPHASE", "SDPPSU6", "SDPSTRA6", "PHPFAST"]:
        if col in df.columns:
            keep_cols.append(col)

    # Deduplicate keep_cols while preserving order
    seen: set[str] = set()
    keep_cols = [c for c in keep_cols if c in df.columns and not (c in seen or seen.add(c))]  # type: ignore[func-returns-value]
    df = df.loc[:, keep_cols].copy()

    print(f"  NHANES III after column selection: {df.shape[0]:,d} rows, {df.shape[1]} columns")
    return df


# =============================================================================
# NHANES IV (continuous) data
# =============================================================================


def _load_xpt(path: Path) -> pd.DataFrame | None:
    """Load a SAS XPT file, returning None if file doesn't exist or is invalid."""
    if not path.exists() or path.stat().st_size < 100:
        return None
    try:
        return pd.read_sas(path)
    except Exception as e:
        print(f"  ⚠ Error reading {path}: {e}")
        return None


def load_nhanes_iv_cycle(
    cycle: str,
    suffix: str,
    start_year: int,
) -> pd.DataFrame | None:
    """Load all PhenoAge-relevant data for a single NHANES IV cycle.

    Downloads XPT files via Make as needed, loads them, and merges on SEQN.

    Args:
        cycle: Cycle label (e.g., "2003-2004")
        suffix: File suffix letter (e.g., "C")
        start_year: Start year of cycle (e.g., 2003)

    Returns:
        Merged DataFrame with SEQN and all PhenoAge biomarkers, or None if
        critical files are missing.
    """
    base_path = RAW_DIR / cycle

    # Load DEMO first as the base
    demo_file = get_file_name("DEMO", cycle, suffix)
    demo_path = base_path / f"{demo_file}.XPT"
    if not make_file(demo_path):
        print(f"  ⚠ Could not obtain DEMO for {cycle}")
        return None

    base_df = _load_xpt(demo_path)
    if base_df is None:
        print(f"  ⚠ Could not load DEMO for {cycle}")
        return None

    demo_vars = [v for v in PHENOAGE_NHANES_IV_VARS["DEMO"] if v in base_df.columns]
    result = base_df.loc[:, ["SEQN"] + demo_vars].copy()
    print(f"  {cycle} DEMO: {result.shape[0]:,d} rows")

    # Load each additional file
    for file_key in ["BIOPRO", "CBC", "CRP", "GLU"]:
        vars_needed = PHENOAGE_NHANES_IV_VARS[file_key]

        # Skip CRP for cycles where it's unavailable
        if file_key == "CRP" and cycle in CRP_UNAVAILABLE_CYCLES:
            print(f"  {cycle} {file_key}: SKIPPED (not available)")
            for v in vars_needed:
                result[v] = np.nan
            continue

        actual_file = get_file_name(file_key, cycle, suffix)
        file_path = base_path / f"{actual_file}.XPT"

        if not make_file(file_path):
            print(f"  ⚠ Could not obtain {file_key} ({actual_file}) for {cycle}")
            for v in vars_needed:
                result[v] = np.nan
            continue

        df = _load_xpt(file_path)
        if df is None:
            print(f"  ⚠ Could not load {file_key} for {cycle}")
            for v in vars_needed:
                result[v] = np.nan
            continue

        # Handle variable name mappings (e.g., LBXHSCRP -> LBXCRP)
        var_mappings = PHENOAGE_VAR_MAPPING.get(file_key, {})
        for canonical_name, alt_name in var_mappings.items():
            if canonical_name not in df.columns and alt_name in df.columns:
                df[canonical_name] = df[alt_name]

        available_vars = [v for v in vars_needed if v in df.columns]
        missing_vars = set(vars_needed) - set(df.columns)
        if missing_vars:
            print(f"  ⚠ {cycle} {file_key}: missing vars {missing_vars}")
            for v in missing_vars:
                df[v] = np.nan
            available_vars = vars_needed

        merge_df = df.loc[:, ["SEQN"] + available_vars].copy()
        # Ensure unique SEQNs
        merge_df = merge_df.drop_duplicates(subset=["SEQN"])
        result = result.merge(merge_df, on="SEQN", how="left")
        n_non_null = {v: (~result[v].isna()).sum() for v in available_vars}
        print(
            f"  {cycle} {file_key} ({actual_file}): "
            f"{merge_df.shape[0]:,d} rows, non-null: {n_non_null}"
        )

    # Handle CRP unit conversion for HSCRP files (2015+)
    # HSCRP files report CRP in mg/L; standard CRP files report in mg/dL
    # Convert HSCRP mg/L → mg/dL by dividing by 10
    if cycle in {"2015-2016", "2017-2018"} and "LBXCRP" in result.columns:
        result["LBXCRP"] = result["LBXCRP"] * CRP_MGL_TO_MGDL

    result["cycle"] = cycle
    result["year"] = start_year

    return result


def load_nhanes_iv(
    cycles: list[tuple[str, str, int]] | None = None,
) -> pd.DataFrame:
    """Load and combine all NHANES IV cycles.

    Args:
        cycles: List of (cycle_label, suffix, start_year) tuples.
                If None, uses all defined cycles.

    Returns:
        Combined DataFrame with harmonized PhenoAge variable names.
    """
    if cycles is None:
        cycles = NHANES_IV_CYCLES

    all_dfs: list[pd.DataFrame] = []
    for cycle, suffix, start_year in cycles:
        print(f"\n{'='*60}")
        print(f"Loading NHANES IV cycle: {cycle}")
        print(f"{'='*60}")
        df = load_nhanes_iv_cycle(cycle, suffix, start_year)
        if df is not None:
            all_dfs.append(df)
        else:
            print(f"  ⚠ Skipping cycle {cycle}")

    if not all_dfs:
        raise ValueError("No NHANES IV data loaded!")

    combined = pd.concat(all_dfs, ignore_index=True)

    # Harmonize to canonical PhenoAge names
    rename_map = {k: v for k, v in NHANES_IV_TO_PHENOAGE.items() if k in combined.columns}
    combined = combined.rename(columns=rename_map)

    print(f"\nCombined NHANES IV: {combined.shape[0]:,d} rows, {combined.shape[1]} columns")
    print(f"Cycles: {combined.cycle.value_counts().to_dict()}")

    return combined


# =============================================================================
# Combined data assembly
# =============================================================================


def prepare_phenoage_data(
    df: pd.DataFrame,
    mortality: pd.DataFrame,
    min_age: int = 20,
    max_age: int = 84,
) -> pd.DataFrame:
    """Merge biomarker data with mortality and prepare for PhenoAge modeling.

    Args:
        df: DataFrame with biomarker data (from load_nhanes_iii or load_nhanes_iv)
        mortality: DataFrame from load_all_mortality
        min_age: Minimum age to include (Levine used 20)
        max_age: Maximum age to include (Levine used 84)

    Returns:
        DataFrame ready for PhenoAge model fitting or scoring, with:
        - log_crp computed
        - Age filtered
        - Mortality merged
        - 10-year mortality indicator created
    """
    df = df.copy()

    # Merge mortality data
    mort_cols = ["SEQN", "eligstat", "mortstat", "permth_exm", "permth_int", "ucod_leading"]
    mort_available = [c for c in mort_cols if c in mortality.columns]
    df = df.merge(mortality.loc[:, mort_available], on="SEQN", how="left")

    # Filter to eligible for mortality follow-up
    n_before = df.shape[0]
    df = df.loc[df.eligstat == 1].copy()
    print(f"  Filtered to mortality-eligible: {n_before:,d} -> {df.shape[0]:,d}")

    # Filter age range
    if "age" in df.columns:
        n_before = df.shape[0]
        df = df.loc[(df.age >= min_age) & (df.age <= max_age)].copy()
        print(f"  Filtered to age {min_age}-{max_age}: {n_before:,d} -> {df.shape[0]:,d}")

    # Unit conversions for NHANES IV:
    # Glucose: mg/dL → mmol/L (Levine model expects mmol/L)
    if "glucose_mgdl" in df.columns and "glucose_mmol" not in df.columns:
        df["glucose_mmol"] = df["glucose_mgdl"] * GLUCOSE_MGDL_TO_MMOL

    # Create log(CRP) — CRP should be in mg/dL at this point
    if "crp_mgdl" in df.columns and "log_crp" not in df.columns:
        df["log_crp"] = np.where(df["crp_mgdl"] > 0, np.log(df["crp_mgdl"]), np.nan)

    # Create 10-year (120-month) mortality indicator
    # Use permth_exm (months from MEC exam to death/censor)
    exposure_col = "permth_exm"
    if exposure_col not in df.columns:
        exposure_col = "permth_int"  # fallback to interview date

    if exposure_col in df.columns:
        df["mort_10yr"] = np.where(
            (df.mortstat == 1) & (df[exposure_col] <= 120),
            1,
            np.where(
                df[exposure_col].isna(),
                np.nan,
                0,
            ),
        )
        # Exposure capped at 120 months
        df["exposure_10yr"] = np.minimum(df[exposure_col], 120)

        n_deaths = (df.mort_10yr == 1).sum()
        n_alive = (df.mort_10yr == 0).sum()
        n_na = df.mort_10yr.isna().sum()
        print(
            f"  10-year mortality: {n_deaths:,d} deaths, {n_alive:,d} alive, "
            f"{n_na:,d} missing exposure"
        )
    else:
        print("  ⚠ No exposure column found; cannot create mortality indicator")

    return df


def load_complete_dataset(
    include_nhanes_iii: bool = True,
    nhanes_iv_cycles: list[tuple[str, str, int]] | None = None,
) -> tuple[pd.DataFrame | None, pd.DataFrame, pd.DataFrame]:
    """Load the complete dataset for PhenoAge analysis.

    Returns:
        Tuple of (nhanes_iii_df, nhanes_iv_df, mortality_df)
    """
    # Determine which mortality files we need
    mort_cycles: list[str] = []
    if include_nhanes_iii:
        mort_cycles.append("NHANES_III")
    iv_cycles = nhanes_iv_cycles or NHANES_IV_CYCLES
    mort_cycles.extend([c[0] for c in iv_cycles])

    print("=" * 60)
    print("Loading mortality data")
    print("=" * 60)
    mortality = load_all_mortality(mort_cycles)

    nhanes_iii_df = None
    if include_nhanes_iii:
        print("\n" + "=" * 60)
        print("Loading NHANES III")
        print("=" * 60)
        nhanes_iii_df = load_nhanes_iii()
        nhanes_iii_df = prepare_phenoage_data(nhanes_iii_df, mortality)

    print("\n" + "=" * 60)
    print("Loading NHANES IV")
    print("=" * 60)
    nhanes_iv_df = load_nhanes_iv(iv_cycles)
    nhanes_iv_df = prepare_phenoage_data(nhanes_iv_df, mortality)

    return nhanes_iii_df, nhanes_iv_df, mortality
