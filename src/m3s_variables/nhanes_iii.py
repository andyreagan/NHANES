"""Load NHANES III fixed-width data files and map to M3S variable names.

NHANES III (1988-1994) uses a completely different data format from NHANES IV
(continuous, 1999+): three large fixed-width .dat files (adult, exam, lab)
instead of per-topic SAS transport (.XPT) files.

This module:
1. Parses the SAS layout programs (adult.sas, exam.sas, lab.sas) from CDC
2. Extracts only the columns needed for the M3S variable mapping pipeline
3. Renames them to the canonical NHANES IV variable names expected by constants.py

The result is a DataFrame with the same column names as NHANES IV DEMO/BIOPRO/etc.
files, so it can be passed directly into the existing ``process_variables()`` and
``impute_defaults()`` pipeline.

Column provenance (NHANES III → NHANES IV canonical name):
    ─── Demographics (adult.dat) ───
    HSAGEIR   → RIDAGEYR     Age at interview (years)
    HSSEX     → RIAGENDR     Sex (1=Male, 2=Female)
    DMARETHN  → RIDRETH1     Race/ethnicity (recoded; see notes)
    HFA8R     → DMDEDUC2     Highest grade of school (recoded)
    WTPFEX6   → WTMEC2YR     MEC exam weight (both phases)
    DMPPIR    → INDHHINC     Poverty income ratio (≈income proxy)

    ─── Smoking (adult.dat) ───
    HAR1      → SMQ020       Ever smoked 100+ cigarettes
    HAR3      → SMQ040       Do you smoke now (1=every day, 2=some days, 3=not at all)

    ─── Medical conditions (adult.dat) ───
    HAC1A     → MCQ160A      Arthritis
    HAC1C     → MCQ160B      Congestive heart failure (mapped to MCQ160B)
    HAC1D     → MCQ160F      Stroke
    HAC1F     → MCQ160K      Chronic bronchitis
    HAC1G     → MCQ160G      Emphysema
    HAC1O     → MCQ220       Cancer (other than skin)
    HAC1N     → MCQ220_skin  Skin cancer (not separately tracked in NHANES IV)
    HAE2      → DIQ010       Diabetes
    HAC4A     → MCQ300C      Family diabetes
    HAC4B     → MCQ300A      Family heart attack/angina
    HAG5C     → OSQ010c      Spine fracture

    ─── Blood Pressure (exam.dat) ───
    PEP6G1    → BPXSY1       Systolic BP reading 1 (K1, mmHg)
    PEP6H1    → BPXSY2       Systolic BP reading 2
    PEP6I1    → BPXSY3       Systolic BP reading 3
    PEP6G3    → BPXDI1       Diastolic BP reading 1 (K5, mmHg)
    PEP6H3    → BPXDI2       Diastolic BP reading 2
    PEP6I3    → BPXDI3       Diastolic BP reading 3
    PEP6DR    → BPXPLS       Pulse rate (beats/min)

    ─── Body measures (exam.dat) ───
    BMPWT     → BMXWT        Weight (kg)
    BMPHT     → BMXHT        Height (cm)
    BMPBMI    → BMXBMI       BMI

    ─── Self-reported weight (adult.dat) ───
    HAM6S     → WHD020       Self-reported weight (lbs)

    ─── Lab (lab.dat) ───
    AMP       → LBXSAL       Albumin (g/dL)
    AMPSI     → LBDSALSI     Albumin SI (g/L)
    GBP       → LBXSGB       Globulin (g/dL)
    GBPSI     → LBDSGBSI     Globulin SI (g/L)
    APPSI     → LBXSAPSI     Alkaline phosphatase SI (U/L)
    ASPSI     → LBXSASSI     AST (SGOT) SI (U/L)
    ATPSI     → LBXSATSI     ALT (SGPT) SI (U/L)
    BUP       → LBXSBU       Blood urea nitrogen (mg/dL)
    CEP       → LBXSCR       Creatinine (mg/dL)
    TGP       → LBXSTR       Triglycerides (mg/dL)
    TBP       → LBXSTB       Total bilirubin (mg/dL)
    TPP       → LBXSTP       Total protein (g/dL)
    GGPSI     → LBXSGTSI     GGT SI (U/L)
    TCP       → LBXTC        Total cholesterol (mg/dL)
    HDP       → LBDHDD       HDL cholesterol (mg/dL)  [note: NHANES IV pre-2005 = LBXHDD]
    LCP       → LBDLDL       LDL cholesterol (mg/dL)
    GHP       → LBXGH        Glycohemoglobin (%)
    HCP       → LBDHCV       Hepatitis C antibody
    HIVP      → LBDHI        HIV antibody
    URP       → URXUCR       Urinary creatinine (mg/dL)
    UBP       → URXUMA       Urinary albumin (µg/mL → mg/L, see note)

    ─── Depression (exam.dat) ───
    MQPDG43L  → DPQ090       "Feel life not worth living" (≈PHQ-9 Q9)
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path("data")
RAW_DIR = DATA_DIR / "raw" / "NHANES_III"


# ============================================================================
# SAS layout parser
# ============================================================================


def parse_sas_input_block(filepath: str | Path) -> list[tuple[str, int, int]]:
    """Parse a NHANES III SAS layout program and return (name, start, end) tuples.

    Returns 1-indexed inclusive column positions (matching SAS convention).
    """
    with open(filepath) as f:
        text = f.read()
    m = re.search(r"INPUT\s*\n(.*?);\s*\n", text, re.DOTALL)
    if not m:
        raise ValueError(f"Could not find INPUT block in {filepath}")
    specs = []
    for line in m.group(1).strip().split("\n"):
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        name = parts[0]
        pos_str = parts[1]
        if "-" in pos_str:
            start, end = pos_str.split("-")
            specs.append((name, int(start), int(end)))
        else:
            col = int(pos_str)
            specs.append((name, col, col))
    return specs


def download_sas_layout(filename: str) -> Path:
    """Download a NHANES III SAS layout file if not already present."""
    dest = RAW_DIR / filename
    if dest.exists() and dest.stat().st_size > 100:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes3/1a/{filename}"
    from urllib.request import Request, urlopen

    req = Request(url, headers={"User-Agent": "NHANES-download/1.0"})
    with urlopen(req, timeout=60) as resp:
        dest.write_bytes(resp.read())
    print(f"  Downloaded {url} → {dest}")
    return dest


def download_dat_file(filename: str) -> Path:
    """Download a NHANES III .dat file if not already present."""
    dest = RAW_DIR / filename
    if dest.exists() and dest.stat().st_size > 100:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes3/1a/{filename}"
    from urllib.request import Request, urlopen

    print(f"  Downloading {url}...")
    req = Request(url, headers={"User-Agent": "NHANES-download/1.0"})
    with urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())
    print(f"  Downloaded ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def load_fwf(
    dat_path: Path, all_specs: list[tuple[str, int, int]], wanted_vars: list[str]
) -> pd.DataFrame:
    """Load only the requested columns from a NHANES III fixed-width file."""
    spec_map = {name: (start, end) for name, start, end in all_specs}
    keep_specs = []
    keep_names = []
    for var in wanted_vars:
        if var in spec_map:
            start, end = spec_map[var]
            keep_specs.append((start - 1, end))  # convert to 0-indexed half-open
            keep_names.append(var)
        # else: silently skip — caller will handle missing columns

    if not keep_specs:
        return pd.DataFrame()

    df = pd.read_fwf(
        dat_path,
        colspecs=keep_specs,
        header=None,
        names=keep_names,
        na_values=[".", "", " "],
    )
    return df


# ============================================================================
# Variable mapping: NHANES III name → NHANES IV canonical name
# ============================================================================

# Which NHANES III variables to extract from each file, and their NHANES IV name
ADULT_VARS: dict[str, str] = {
    # Demographics
    "SEQN": "SEQN",
    "HSAGEIR": "RIDAGEYR",
    "HSSEX": "RIAGENDR",
    "DMARETHN": "RIDRETH1_raw_iii",  # needs recoding
    "HFA8R": "DMDEDUC2_raw_iii",  # needs recoding
    "WTPFEX6": "WTMEC2YR",
    "DMPPIR": "INDHHINC_raw_iii",  # poverty ratio, not income category
    # Survey design
    "SDPPHASE": "SDPPHASE",
    "SDPPSU6": "SDPPSU6",
    "SDPSTRA6": "SDPSTRA6",
    # Smoking
    "HAR1": "SMQ020",  # ever smoked 100+ cigs (1=yes, 2=no)
    "HAR3": "SMQ040",  # smoke now (1=yes, 2=no → recode to 1=every day,3=not at all)
    # Medical conditions
    "HAC1A": "MCQ160A",  # arthritis
    "HAC1C": "MCQ160B",  # congestive heart failure → MCQ160B
    "HAC1D": "MCQ160F",  # stroke
    "HAC1F": "MCQ160K",  # chronic bronchitis
    "HAC1G": "MCQ160G",  # emphysema
    "HAC1O": "MCQ220",  # cancer (other)
    "HAE2": "DIQ010",  # diabetes
    "HAC4A": "MCQ300C_raw_iii",  # family diabetes (1=yes,2=no; needs recode)
    "HAC4B": "MCQ300A_raw_iii",  # family heart attack (1=yes,2=no)
    "HAG5C": "OSQ010c",  # spine fracture (1=yes, 2=no)
    # Self-reported weight
    "HAM6S": "WHD020",  # weight without clothes (lbs)
    # Depression — not directly comparable; NHANES III used DIS
    # (Diagnostic Interview Schedule) not PHQ-9.
    # MQPDG43L is "felt life not worth living" which is closest to DPQ090
}

EXAM_VARS: dict[str, str] = {
    "SEQN": "SEQN",
    # Blood pressure — use K1 (systolic) and K5 (diastolic) for adults
    "PEP6G1": "BPXSY1",
    "PEP6H1": "BPXSY2",
    "PEP6I1": "BPXSY3",
    "PEP6G3": "BPXDI1",
    "PEP6H3": "BPXDI2",
    "PEP6I3": "BPXDI3",
    "PEP6DR": "BPXPLS",  # pulse rate
    # Body measures
    "BMPWT": "BMXWT",  # weight (kg)
    "BMPHT": "BMXHT",  # height (cm)
    "BMPBMI": "BMXBMI",  # BMI
}

LAB_VARS: dict[str, str] = {
    "SEQN": "SEQN",
    # Biochemistry
    "AMP": "LBXSAL",  # albumin (g/dL)
    "AMPSI": "LBDSALSI",  # albumin SI (g/L)
    "GBP": "LBXSGB",  # globulin (g/dL)
    "GBPSI": "LBDSGBSI",  # globulin SI (g/L)
    "APPSI": "LBXSAPSI",  # alkaline phosphatase SI (U/L)
    "ASPSI": "LBXSASSI",  # AST SI (U/L)
    "ATPSI": "LBXSATSI",  # ALT SI (U/L)
    "BUP": "LBXSBU",  # BUN (mg/dL)
    "CEP": "LBXSCR",  # creatinine (mg/dL)
    "TGP": "LBXSTR",  # triglycerides (mg/dL)
    "TBP": "LBXSTB",  # total bilirubin (mg/dL)
    "TPP": "LBXSTP",  # total protein (g/dL)
    "GGPSI": "LBXSGTSI",  # GGT SI (U/L)
    # Lipids
    "TCP": "LBXTC",  # total cholesterol (mg/dL)
    "HDP": "LBDHDD",  # HDL (mg/dL) — note: canonical is LBDHDD
    "LCP": "LBDLDL",  # LDL (mg/dL)
    # A1C
    "GHP": "LBXGH",  # glycohemoglobin (%)
    # Infection
    "HCP": "LBDHCV",  # hepatitis C antibody (1=pos, 2=neg)
    "HIVP": "LBDHI",  # HIV (1=pos, 2=neg)
    # Urine
    "URP": "URXUCR",  # urinary creatinine (mg/dL)
    "UBP": "URXUMA",  # urinary albumin (µg/mL; NHANES IV uses mg/L = µg/mL)
}


# ============================================================================
# Race/ethnicity recoding
# ============================================================================
# NHANES III DMARETHN: 1=Non-Hispanic white, 2=Non-Hispanic black,
#                      3=Mexican-American, 4=Other
# NHANES IV RIDRETH1:  1=Mexican American, 2=Other Hispanic,
#                      3=Non-Hispanic White, 4=Non-Hispanic Black, 5=Other
RACE_RECODE: dict[int, int] = {
    1: 3,  # NH White → 3
    2: 4,  # NH Black → 4
    3: 1,  # Mexican American → 1
    4: 5,  # Other → 5
}


# Education recoding
# NHANES III HFA8R: grade completed (0-17, where 17=5+ years college)
# NHANES IV DMDEDUC2: 1=<9th grade, 2=9-11th, 3=HS/GED, 4=some college, 5=college+
def recode_education(grade: float) -> float:
    """Convert NHANES III grade (0-17) to NHANES IV DMDEDUC2 categories."""
    if pd.isna(grade) or grade >= 88:
        return np.nan
    grade = int(grade)
    if grade <= 8:
        return 1  # Less than 9th grade
    elif grade <= 11:
        return 2  # 9-11th grade
    elif grade == 12:
        return 3  # HS graduate / GED
    elif grade <= 15:
        return 4  # Some college / AA degree
    else:
        return 5  # College graduate or above


# Smoking recoding
# NHANES III HAR3: 1=yes, 2=no (binary: do you smoke now?)
# NHANES IV SMQ040: 1=every day, 2=some days, 3=not at all
def recode_smoking_now(val: float) -> float:
    """Convert NHANES III HAR3 to NHANES IV SMQ040."""
    if pd.isna(val):
        return np.nan
    if val == 1:
        return 1  # yes → every day (conservative)
    elif val == 2:
        return 3  # no → not at all
    return np.nan


# ============================================================================
# Main loader
# ============================================================================


def load_nhanes_iii_for_m3s() -> pd.DataFrame:
    """Load NHANES III data and return a DataFrame with NHANES IV column names.

    Downloads raw .dat and .sas files from CDC if not present locally.

    Returns:
        DataFrame with columns matching NHANES IV naming conventions,
        ready for the M3S ``process_variables()`` pipeline.
    """
    print("=" * 70)
    print("Loading NHANES III for M3S pipeline")
    print("=" * 70)

    # Download SAS layout files
    adult_sas = download_sas_layout("adult.sas")
    exam_sas = download_sas_layout("exam.sas")
    lab_sas = download_sas_layout("lab.sas")

    # Download data files
    adult_dat = download_dat_file("adult.dat")
    exam_dat = download_dat_file("exam.dat")
    lab_dat = download_dat_file("lab.dat")

    # Parse layouts
    adult_specs = parse_sas_input_block(adult_sas)
    exam_specs = parse_sas_input_block(exam_sas)
    lab_specs = parse_sas_input_block(lab_sas)

    print(
        f"  Parsed layouts: adult={len(adult_specs)}, exam={len(exam_specs)}, lab={len(lab_specs)} variables"
    )

    # Load each file (only requested columns)
    print("  Loading adult.dat...")
    adult_df = load_fwf(adult_dat, adult_specs, list(ADULT_VARS.keys()))
    adult_df = adult_df.rename(columns=ADULT_VARS)
    print(f"    {adult_df.shape[0]:,d} rows, {adult_df.shape[1]} columns")

    print("  Loading exam.dat...")
    exam_df = load_fwf(exam_dat, exam_specs, list(EXAM_VARS.keys()))
    exam_df = exam_df.rename(columns=EXAM_VARS)
    print(f"    {exam_df.shape[0]:,d} rows, {exam_df.shape[1]} columns")

    print("  Loading lab.dat...")
    lab_df = load_fwf(lab_dat, lab_specs, list(LAB_VARS.keys()))
    lab_df = lab_df.rename(columns=LAB_VARS)
    print(f"    {lab_df.shape[0]:,d} rows, {lab_df.shape[1]} columns")

    # Replace sentinel values (≥888 = blank/missing in NHANES III)
    # Only apply to biomarker/questionnaire columns, NOT to weights, ages,
    # design variables, self-reported weight, or other legitimately large values.
    skip_sentinel = {
        "SEQN",
        "WTMEC2YR",
        "RIDAGEYR",
        "SDPPHASE",
        "SDPPSU6",
        "SDPSTRA6",
        "INDHHINC_raw_iii",  # poverty ratio can be > 5
        "BMXWT",
        "BMXHT",
        "BMXBMI",  # weight(kg), height(cm), BMI
        "WHD020",  # self-reported weight (lbs)
        "BPXSY1",
        "BPXSY2",
        "BPXSY3",
        "BPXSY4",  # BP can be > 200 mmHg
        "LBXSTR",  # triglycerides can be > 888
        "LBXTC",  # total cholesterol can be > 888
        "URXUCR",  # urinary creatinine can be > 888
        "LBXSAPSI",  # alkaline phosphatase can be > 888
    }
    for df_to_clean in [adult_df, exam_df, lab_df]:
        for col in df_to_clean.select_dtypes(include=["number"]).columns:
            if col in skip_sentinel:
                continue
            n_sentinel = (df_to_clean[col] >= 888).sum()
            if n_sentinel > 0:
                df_to_clean.loc[df_to_clean[col] >= 888, col] = np.nan

    # Merge on SEQN (adults only: adult.dat has 20,050 adults)
    # exam.dat and lab.dat have all ages (31,311 and 29,314 respectively)
    df = adult_df.merge(exam_df, on="SEQN", how="left")
    df = df.merge(lab_df, on="SEQN", how="left")
    print(f"  Merged: {df.shape[0]:,d} rows, {df.shape[1]} columns")

    # ─── Recoding ───

    # Race/ethnicity
    df["RIDRETH1"] = df["RIDRETH1_raw_iii"].map(RACE_RECODE)

    # Education
    df["DMDEDUC2"] = df["DMDEDUC2_raw_iii"].apply(recode_education)

    # Smoking
    df["SMQ040"] = df["SMQ040"].apply(recode_smoking_now)

    # Family history: NHANES III uses 1=yes, 2=no; NHANES IV uses 1=yes, 2=no, 9=dk
    # These are compatible as-is for the M3S pipeline expressions that check == 1

    # Income: NHANES III DMPPIR is a ratio (0-5+), not a category.
    # NHANES IV INDHHINC is a category (1-15). We keep the raw ratio for now
    # and let the pipeline handle it (income is not required/predictor).
    df["INDHHINC"] = df["INDHHINC_raw_iii"]  # imperfect but not used as predictor

    # Family history recoding: MCQ300C/MCQ300A in pre-2005 are MCQ250A/MCQ250G
    # but our NHANES III values are already 1/2, which matches NHANES IV coding
    df["MCQ300C"] = df.get("MCQ300C_raw_iii")
    df["MCQ300A"] = df.get("MCQ300A_raw_iii")

    # NHANES III doesn't have RIDAGEEX or RIDAGEMN — set to NA
    df["RIDAGEEX"] = np.nan
    df["RIDAGEMN"] = np.nan

    # MCQ160C, MCQ160D, MCQ160E (coronary heart disease, angina, heart attack)
    # not separately available — set to NA (heart_condition will use MCQ160B = CHF)
    for col in ["MCQ160C", "MCQ160D", "MCQ160E"]:
        if col not in df.columns:
            df[col] = np.nan

    # MCQ160L (liver condition) — not directly available
    if "MCQ160L" not in df.columns:
        df["MCQ160L"] = np.nan

    # MCQ053 (blood disorder), BPXPULS (irregular pulse) — not available
    for col in ["MCQ053", "BPXPULS"]:
        if col not in df.columns:
            df[col] = np.nan

    # 4th BP reading not in NHANES III exam protocol
    df["BPXSY4"] = np.nan
    df["BPXDI4"] = np.nan

    # Questionnaire variables not in NHANES III
    for col in [
        "DUQ240",
        "HSQ510",
        "DIQ280",
        "PFD069D",
        "PFD069E",
        "PFD069K",
        "PFQ061B",
        "PFQ061C",
        "CIDDSCOR",
        "MPQ120AB",
        "DPQ090",
        "OSD030ca",
        "SMQ680",
        "SMD070",
        "SMQ620",
        "SMQ660",
        "RHQ360",
        "RDQ134",
        "KIQ022",
        "LBXP1",
    ]:
        if col not in df.columns:
            df[col] = np.nan

    print(f"  After recoding: {df.shape[0]:,d} rows, {df.shape[1]} columns")

    # Drop the raw_iii helper columns
    drop_cols = [c for c in df.columns if c.endswith("_raw_iii")]
    df = df.drop(columns=drop_cols)

    return df
