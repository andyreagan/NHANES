"""Parse all NHANES Linked Mortality Files (LMF) into parquet files.

The CDC publishes fixed-width .dat mortality follow-up files for each NHANES
cycle. Multiple **vintages** exist, each extending follow-up further:

  - **2011 vintage**: follow-up through Dec 31, 2011 (only 2003-2004 and
    2005-2006 recovered, from the rnhanesdata R package on GitHub)
  - **2015 vintage**: follow-up through Dec 31, 2015 (9 files, NHANES III–2013-2014,
    recovered from the Internet Archive)
  - **2019 vintage**: follow-up through Dec 31, 2019 (11 files, NHANES III–2017-2018,
    current CDC FTP release)

The 2011 vintage has a **different column layout** from the 2015/2019 vintages —
it includes extra fields (causeavl, mortality source flags) and uses different
column positions. The layout was recovered from the rnhanesdata R package's
``process_mort()`` function.

See ``references/linked_mortality_vintages.md`` for full documentation of all
vintages, layouts, and provenance.

Output files:
    data/processed/LMF_Files/LMF_all_MORT_2019.parquet
    data/processed/LMF_Files/LMF_all_MORT_2015.parquet
    data/processed/LMF_Files/LMF_all_MORT_2011.parquet
    data/processed/LMF_Files/LMF__2003-2004__2005-2006__MORT_2019.parquet  (legacy)

Usage:
    uv run python -m src.parse_lmf.parse_lmf
"""

from pathlib import Path

import pandas as pd

DATA_DIR = Path("data")

# ============================================================================
# Column specifications for the 2015/2019 vintages
# Matches CDC SAS_ReadInProgramAllSurveys.sas (NHANES version).
# 0-indexed, half-open intervals for pd.read_fwf.
# ============================================================================

COLSPECS_2015_2019: list[tuple[int, int]] = [
    (0, 6),  # SEQN          (SAS cols 1-6)
    (14, 15),  # ELIGSTAT      (SAS col 15)
    (15, 16),  # MORTSTAT      (SAS col 16)
    (16, 19),  # UCOD_LEADING  (SAS cols 17-19)
    (19, 20),  # DIABETES      (SAS col 20)
    (20, 21),  # HYPERTEN      (SAS col 21)
    (42, 45),  # PERMTH_INT    (SAS cols 43-45)
    (45, 48),  # PERMTH_EXM    (SAS cols 46-48)
]

COLUMNS_2015_2019: list[str] = [
    "SEQN",
    "eligstat",
    "mortstat",
    "ucod_leading",
    "diabetes",
    "hyperten",
    "permth_int",
    "permth_exm",
]

# ============================================================================
# Column specifications for the 2011 vintage
# Recovered from the rnhanesdata R package's process_mort() function.
# This vintage has additional fields and different column positions.
# Records are up to 54 chars (variable-length lines).
# ============================================================================

COLSPECS_2011: list[tuple[int, int]] = [
    (0, 14),  # SEQN          (R: 1-14, SEQN is in first ~5 chars, rest blank)
    (14, 15),  # ELIGSTAT      (R: 15-15)
    (15, 16),  # MORTSTAT      (R: 16-16)
    (16, 17),  # CAUSEAVL      (R: 17-17) — cause of death data available flag
    (17, 20),  # UCOD_LEADING  (R: 18-20)
    (20, 21),  # DIABETES      (R: 21-21)
    (21, 22),  # HYPERTEN      (R: 22-22)
    (43, 46),  # PERMTH_INT    (R: 44-46)
    (46, 49),  # PERMTH_EXM    (R: 47-49)
    (49, 50),  # MORTSRCE_NDI  (R: 50-50)
    (50, 51),  # MORTSRCE_CMS  (R: 51-51)
    (51, 52),  # MORTSRCE_SSA  (R: 52-52)
    (52, 53),  # MORTSRCE_DC   (R: 53-53)
    (53, 54),  # MORTSRCE_DCL  (R: 54-54)
]

COLUMNS_2011: list[str] = [
    "SEQN",
    "eligstat",
    "mortstat",
    "causeavl",
    "ucod_leading",
    "diabetes",
    "hyperten",
    "permth_int",
    "permth_exm",
    "mortsrce_ndi",
    "mortsrce_cms",
    "mortsrce_ssa",
    "mortsrce_dc",
    "mortsrce_dcl",
]

# ============================================================================
# Value labels (from CDC SAS format statements)
# ============================================================================

ELIGSTAT_LABELS: dict[int, str] = {
    1: "Eligible",
    2: "Under age 18, not available for public release",
    3: "Ineligible",
}

MORTSTAT_LABELS: dict[int, str] = {
    0: "Assumed alive",
    1: "Assumed deceased",
}

UCOD_LEADING_LABELS: dict[str, str] = {
    "001": "Diseases of heart (I00-I09, I11, I13, I20-I51)",
    "002": "Malignant neoplasms (C00-C97)",
    "003": "Chronic lower respiratory diseases (J40-J47)",
    "004": "Accidents (unintentional injuries) (V01-X59, Y85-Y86)",
    "005": "Cerebrovascular diseases (I60-I69)",
    "006": "Alzheimer's disease (G30)",
    "007": "Diabetes mellitus (E10-E14)",
    "008": "Influenza and pneumonia (J09-J18)",
    "009": "Nephritis, nephrotic syndrome and nephrosis (N00-N07, N17-N19, N25-N27)",
    "010": "All other causes (residual)",
}

# ============================================================================
# Vintage definitions
# ============================================================================

VINTAGES: dict[str, dict] = {
    "2019": {
        "raw_dir": DATA_DIR / "raw" / "LMF_Files",
        "colspecs": COLSPECS_2015_2019,
        "columns": COLUMNS_2015_2019,
        "files": {
            "NHANES_III": ("NHANES_III_MORT_2019_PUBLIC.dat", 1991),
            "1999-2000": ("NHANES_1999_2000_MORT_2019_PUBLIC.dat", 1999),
            "2001-2002": ("NHANES_2001_2002_MORT_2019_PUBLIC.dat", 2001),
            "2003-2004": ("NHANES_2003_2004_MORT_2019_PUBLIC.dat", 2003),
            "2005-2006": ("NHANES_2005_2006_MORT_2019_PUBLIC.dat", 2005),
            "2007-2008": ("NHANES_2007_2008_MORT_2019_PUBLIC.dat", 2007),
            "2009-2010": ("NHANES_2009_2010_MORT_2019_PUBLIC.dat", 2009),
            "2011-2012": ("NHANES_2011_2012_MORT_2019_PUBLIC.dat", 2011),
            "2013-2014": ("NHANES_2013_2014_MORT_2019_PUBLIC.dat", 2013),
            "2015-2016": ("NHANES_2015_2016_MORT_2019_PUBLIC.dat", 2015),
            "2017-2018": ("NHANES_2017_2018_MORT_2019_PUBLIC.dat", 2017),
        },
    },
    "2015": {
        "raw_dir": DATA_DIR / "raw" / "LMF_Files_2015",
        "colspecs": COLSPECS_2015_2019,
        "columns": COLUMNS_2015_2019,
        "files": {
            "NHANES_III": ("NHANES_III_MORT_2015_PUBLIC.dat", 1991),
            "1999-2000": ("NHANES_1999_2000_MORT_2015_PUBLIC.dat", 1999),
            "2001-2002": ("NHANES_2001_2002_MORT_2015_PUBLIC.dat", 2001),
            "2003-2004": ("NHANES_2003_2004_MORT_2015_PUBLIC.dat", 2003),
            "2005-2006": ("NHANES_2005_2006_MORT_2015_PUBLIC.dat", 2005),
            "2007-2008": ("NHANES_2007_2008_MORT_2015_PUBLIC.dat", 2007),
            "2009-2010": ("NHANES_2009_2010_MORT_2015_PUBLIC.dat", 2009),
            "2011-2012": ("NHANES_2011_2012_MORT_2015_PUBLIC.dat", 2011),
            "2013-2014": ("NHANES_2013_2014_MORT_2015_PUBLIC.dat", 2013),
        },
    },
    "2011": {
        "raw_dir": DATA_DIR / "raw" / "LMF_Files_2011",
        "colspecs": COLSPECS_2011,
        "columns": COLUMNS_2011,
        "files": {
            # Only 2003-2004 and 2005-2006 recovered (from rnhanesdata R package)
            "2003-2004": ("NHANES_2003_2004_MORT_2011_PUBLIC.dat", 2003),
            "2005-2006": ("NHANES_2005_2006_MORT_2011_PUBLIC.dat", 2005),
        },
    },
}

# Wayback Machine URLs for the 2015 vintage (captured Jan 15, 2021)
WAYBACK_2015_BASE = (
    "https://web.archive.org/web/20210115id_/"
    "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/"
)

OUTPUT_DIR = DATA_DIR / "processed" / "LMF_Files"
LEGACY_OUTPUT_PATH = OUTPUT_DIR / "LMF__2003-2004__2005-2006__MORT_2019.parquet"


def parse_one(
    path: Path,
    cycle: str,
    year: int,
    colspecs: list[tuple[int, int]],
    columns: list[str],
) -> pd.DataFrame:
    """Parse a single linked mortality .dat file."""
    df = pd.read_fwf(
        path,
        colspecs=colspecs,
        skiprows=0,
        header=None,
        na_values=[".", ""],
    )
    df.columns = columns
    df["cycle"] = cycle
    df["year"] = year
    return df


def download_file(url: str, dest: Path) -> bool:
    """Download a file if it doesn't exist. Returns True on success."""
    if dest.exists() and dest.stat().st_size > 100:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    import time
    from urllib.error import HTTPError, URLError
    from urllib.request import Request, urlopen

    for attempt in range(3):
        try:
            req = Request(url, headers={"User-Agent": "NHANES-LMF-download/1.0"})
            with urlopen(req, timeout=60) as resp:
                data = resp.read()
            dest.write_bytes(data)
            return True
        except (HTTPError, URLError, TimeoutError, OSError) as e:
            if attempt < 2:
                time.sleep(2 * (attempt + 1))
            else:
                print(f"    ⚠ Failed to download {url}: {e}")
                return False
    return False


def process_vintage(vintage: str, config: dict) -> pd.DataFrame | None:
    """Parse all files for a single vintage into a combined DataFrame."""
    raw_dir = config["raw_dir"]
    colspecs = config["colspecs"]
    columns = config["columns"]
    files = config["files"]
    output_path = OUTPUT_DIR / f"LMF_all_MORT_{vintage}.parquet"

    print(f"\n{'='*70}")
    print(f"Processing {vintage} vintage")
    print(f"{'='*70}")

    all_dfs: list[pd.DataFrame] = []
    for cycle, (dat_file, year) in files.items():
        dat_path = raw_dir / dat_file
        if not dat_path.exists():
            print(f"  ⚠ Skipping {cycle}: {dat_path} not found")
            continue
        df = parse_one(dat_path, cycle, year, colspecs, columns)
        n_eligible = (df.eligstat == 1).sum()
        n_deaths = (df.mortstat == 1).sum()
        n_with_exm = df.permth_exm.notna().sum()
        n_with_int = df.permth_int.notna().sum()
        print(
            f"  {cycle:15s}  {df.shape[0]:>6,d} rows  "
            f"{n_eligible:>6,d} eligible  "
            f"{n_deaths:>5,d} deaths  "
            f"permth_exm: {n_with_exm:>6,d}  "
            f"permth_int: {n_with_int:>6,d}"
        )
        all_dfs.append(df)

    if not all_dfs:
        print(f"  ⚠ No files found for {vintage} vintage")
        return None

    combined = pd.concat(all_dfs, ignore_index=True)
    extra_cols = [c for c in combined.columns if c not in COLUMNS_2015_2019 + ["cycle", "year"]]
    print(f"\n  Combined: {combined.shape[0]:,d} rows across {len(all_dfs)} cycles")
    print(f"  Total eligible:  {(combined.eligstat == 1).sum():,d}")
    print(f"  Total deaths:    {(combined.mortstat == 1).sum():,d}")
    print(f"  permth_exm max:  {combined.permth_exm.max():.0f} months")
    if extra_cols:
        print(f"  Extra columns:   {extra_cols}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(output_path)
    print(f"\n  Wrote {output_path}")

    return combined


def main() -> None:
    # ------------------------------------------------------------------
    # Download 2015 vintage from Internet Archive (if not already present)
    # ------------------------------------------------------------------
    config_2015 = VINTAGES["2015"]
    missing_2015 = []
    for cycle, (dat_file, _year) in config_2015["files"].items():
        dest = config_2015["raw_dir"] / dat_file
        if not dest.exists() or dest.stat().st_size < 100:
            missing_2015.append((cycle, dat_file))

    if missing_2015:
        print(f"\nDownloading {len(missing_2015)} files for 2015 vintage from Internet Archive...")
        for cycle, dat_file in missing_2015:
            dest = config_2015["raw_dir"] / dat_file
            url = WAYBACK_2015_BASE + dat_file
            print(f"  {cycle}: {dat_file}...", end=" ", flush=True)
            if download_file(url, dest):
                print(f"OK ({dest.stat().st_size:,d} bytes)")
            else:
                print("FAILED")

    # ------------------------------------------------------------------
    # Process all vintages
    # ------------------------------------------------------------------
    combined_2019 = process_vintage("2019", VINTAGES["2019"])
    process_vintage("2015", VINTAGES["2015"])
    process_vintage("2011", VINTAGES["2011"])

    # ------------------------------------------------------------------
    # Legacy output for backwards compatibility
    # ------------------------------------------------------------------
    if combined_2019 is not None:
        legacy = combined_2019[combined_2019.cycle.isin(["2003-2004", "2005-2006"])].copy()
        legacy.to_parquet(LEGACY_OUTPUT_PATH)
        print(f"\n  Wrote {LEGACY_OUTPUT_PATH} (legacy, {legacy.shape[0]:,d} rows)")


if __name__ == "__main__":
    main()
