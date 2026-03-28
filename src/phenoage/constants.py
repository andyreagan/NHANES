"""Constants for PhenoAge computation.

Reproduces Levine et al. 2018:
"An epigenetic biomarker of aging for lifespan and healthspan"
(Aging, 2018, Vol 10, No 4)

Training data: NHANES III (1988-1994), ages 20-84
Validation data: NHANES IV continuous (1999-2010)
Model: Gompertz proportional hazards, 10-year mortality
Biomarkers: 9 clinical chemistry + chronological age
"""

import numpy as np

# =============================================================================
# NHANES cycle definitions
# =============================================================================

# NHANES III: training set for Gompertz model
NHANES_III_YEAR = "NHANES_III"

# NHANES IV (continuous) cycles: validation / application
# (year_label, suffix_letter, start_year)
NHANES_IV_CYCLES: list[tuple[str, str, int]] = [
    ("1999-2000", "", 1999),
    ("2001-2002", "B", 2001),
    ("2003-2004", "C", 2003),
    ("2005-2006", "D", 2005),
    ("2007-2008", "E", 2007),
    ("2009-2010", "F", 2009),
    # CRP not available 2011-2014, so PhenoAge can't be computed
    # Include anyway for completeness; rows will have NA PhenoAge
    ("2011-2012", "G", 2011),
    ("2013-2014", "H", 2013),
    ("2015-2016", "I", 2015),
    ("2017-2018", "J", 2017),
]

# Cycles used by Levine 2018 for validation (CRP available in all)
LEVINE_VALIDATION_CYCLES = [c for c in NHANES_IV_CYCLES if c[2] <= 2010]

# =============================================================================
# NHANES IV file mapping for PhenoAge biomarkers
# =============================================================================

# Map: canonical_file_name -> {cycle_year: actual_file_name}
# For cycles where file names differ from the canonical name
PHENOAGE_FILE_NAMES: dict[str, dict[str, str]] = {
    # Complete Blood Count
    "CBC": {
        "1999-2000": "LAB25",
        "2001-2002": "L25",
        "2003-2004": "L25",
        # 2005+ use CBC
    },
    # C-Reactive Protein
    "CRP": {
        "1999-2000": "LAB11",
        "2001-2002": "L11",
        "2003-2004": "L11",
        # 2005-2010 use CRP
        # 2011-2014: NOT AVAILABLE
        "2015-2016": "HSCRP",
        "2017-2018": "HSCRP",
    },
    # Fasting Glucose (note: L10/LAB10 = GHB/glycohemoglobin, L10AM/LAB10AM = glucose)
    "GLU": {
        "1999-2000": "LAB10AM",
        "2001-2002": "L10AM",
        "2003-2004": "L10AM",
        # 2005+ use GLU
    },
    # Standard Biochemistry Profile (albumin, creatinine, alk phos)
    "BIOPRO": {
        "1999-2000": "LAB18",
        "2001-2002": "L40",
        "2003-2004": "L40",
        # 2005+ use BIOPRO
    },
    # Demographics
    "DEMO": {
        # all cycles use DEMO (with suffix)
    },
}


def get_file_name(canonical: str, cycle: str, suffix: str) -> str:
    """Get the actual XPT file name for a canonical file in a given cycle.

    Args:
        canonical: Canonical file name (e.g., "CBC", "CRP", "GLU", "BIOPRO", "DEMO")
        cycle: Cycle label (e.g., "2003-2004")
        suffix: Suffix letter (e.g., "C", "D", or "" for 1999-2000)

    Returns:
        Full file name like "CBC_D" or "LAB25" (without .XPT extension)
    """
    overrides = PHENOAGE_FILE_NAMES.get(canonical, {})
    base = overrides.get(cycle, canonical)
    if suffix:
        return f"{base}_{suffix}"
    else:
        return base


# =============================================================================
# PhenoAge biomarker variable names in NHANES IV (continuous)
# =============================================================================

# Variables needed from each file for PhenoAge
PHENOAGE_NHANES_IV_VARS: dict[str, list[str]] = {
    "DEMO": ["RIDAGEYR", "RIAGENDR", "RIDRETH1", "WTMEC2YR", "RIDAGEEX", "RIDAGEMN"],
    "BIOPRO": ["LBXSAL", "LBXSCR", "LBXSAPSI"],
    "CBC": ["LBXWBCSI", "LBXLYPCT", "LBXMCVSI", "LBXRDW"],
    "CRP": ["LBXCRP"],  # becomes LBXHSCRP in 2015+
    "GLU": ["LBXGLU"],
}

# Variable name mappings for cycles where names differ
# Maps: canonical_name -> alternate_name_to_look_for
PHENOAGE_VAR_MAPPING: dict[str, dict[str, str]] = {
    # HSCRP files (2015+) use LBXHSCRP instead of LBXCRP
    "CRP": {"LBXCRP": "LBXHSCRP"},
    # 2001-2002 BIOPRO uses LBD prefix instead of LBX for some vars
    "BIOPRO": {
        "LBXSCR": "LBDSCR",
        "LBXSAPSI": "LBDSAPSI",
    },
}

# CRP not available in these cycles
CRP_UNAVAILABLE_CYCLES = {"2011-2012", "2013-2014"}

# =============================================================================
# NHANES III variable names (fixed-width lab.dat file)
# =============================================================================

# Column positions in NHANES III lab.dat (1-indexed, inclusive on both ends)
# Parsed from lab.sas
NHANES_III_LAB_COLSPECS: list[tuple[str, int, int]] = [
    ("SEQN", 1, 5),
    ("HSSEX", 15, 15),
    ("HSAGEIR", 16, 17),
    # Survey design variables
    ("SDPPHASE", 40, 40),
    ("SDPPSU6", 41, 41),
    ("SDPSTRA6", 42, 43),
    # Exam sample weight (MEC examined, both phases combined)
    ("WTPFEX6", 59, 67),
    # Fasting time
    ("PHPFAST", 1263, 1267),
    # CBC
    ("WCP", 1273, 1277),  # White blood cell count (1000 cells/uL -> need to check units)
    ("WCPSI", 1278, 1282),  # WBC: SI (10^9/L) = same as LBXWBCSI
    ("LMPPCNT", 1283, 1287),  # Lymphocyte percent
    ("MVPSI", 1340, 1344),  # Mean cell volume: SI (fL) = same as LBXMCVSI
    ("RWP", 1360, 1364),  # Red cell distribution width (%)
    # Inflammation
    ("CRP", 1667, 1671),  # C-reactive protein (mg/dL)
    # Chemistry
    ("SGP", 1758, 1760),  # Serum glucose (mg/dL)
    ("CEP", 1784, 1787),  # Serum creatinine (mg/dL)
    ("APPSI", 1835, 1838),  # Alkaline phosphatase: SI (U/L)
    ("AMP", 1846, 1848),  # Albumin (g/dL)
    # Glucose tolerance
    ("G1P", 1866, 1870),  # Plasma glucose (mg/dL) — fasting
    ("GHP", 1861, 1864),  # Glycohemoglobin
]

# Map NHANES III variable names → PhenoAge canonical names
# Note: glucose and CRP need unit conversion AFTER renaming
NHANES_III_TO_PHENOAGE: dict[str, str] = {
    "HSAGEIR": "age",
    "HSSEX": "sex_code",  # 1=Male, 2=Female (same as NHANES IV RIAGENDR)
    "AMP": "albumin",  # g/dL (same units as LBXSAL)
    "CEP": "creatinine",  # mg/dL (same units as LBXSCR)
    "SGP": "glucose_mgdl",  # serum glucose mg/dL → will convert to mmol/L
    "CRP": "crp_mgdl",  # mg/dL (same as NHANES IV LBXCRP)
    "LMPPCNT": "lymphocyte_pct",  # % (same units as LBXLYPCT)
    "MVPSI": "mcv",  # fL (same units as LBXMCVSI)
    "RWP": "rdw",  # % (same units as LBXRDW)
    "APPSI": "alp",  # U/L (same units as LBXSAPSI)
    "WCPSI": "wbc",  # 10^9/L (same units as LBXWBCSI)
    "WTPFEX6": "exam_weight",
}

# Map NHANES IV variable names → PhenoAge canonical names
# Note: glucose and CRP need unit conversion AFTER renaming
NHANES_IV_TO_PHENOAGE: dict[str, str] = {
    "RIDAGEYR": "age",
    "RIAGENDR": "sex_code",
    "LBXSAL": "albumin",  # g/dL
    "LBXSCR": "creatinine",  # mg/dL
    "LBXGLU": "glucose_mgdl",  # mg/dL → will convert to mmol/L
    "LBXCRP": "crp_mgdl",  # mg/dL for 1999-2010 cycles
    "LBXLYPCT": "lymphocyte_pct",  # %
    "LBXMCVSI": "mcv",  # fL
    "LBXRDW": "rdw",  # %
    "LBXSAPSI": "alp",  # U/L
    "LBXWBCSI": "wbc",  # 10^9/L
    "WTMEC2YR": "exam_weight",
    "RIDRETH1": "race_ethnicity",
}

# Glucose conversion factor: mg/dL → mmol/L
GLUCOSE_MGDL_TO_MMOL = 1.0 / 18.016

# CRP in HSCRP files (2015+) is in mg/L; divide by 10 to get mg/dL
CRP_MGL_TO_MGDL = 0.1

# The 10 PhenoAge model features (in order used by Levine 2018)
PHENOAGE_FEATURES: list[str] = [
    "albumin",  # g/dL
    "creatinine",  # mg/dL
    "glucose_mmol",  # mmol/L (glucose mg/dL / 18.016)
    "log_crp",  # ln(CRP in mg/dL)
    "lymphocyte_pct",  # %
    "mcv",  # fL
    "rdw",  # %
    "alp",  # U/L
    "wbc",  # 10^9/L (= 1000 cells/µL)
    "age",  # years
]

# =============================================================================
# Published Levine 2018 PhenoAge coefficients
# (From Gompertz proportional hazards fit on NHANES III)
# =============================================================================

# These are the beta coefficients from the Gompertz PH model
# xb = intercept + sum(beta_i * x_i)
# IMPORTANT unit notes (verified by reproducing expected PhenoAge ≈ chronological age
# for healthy individuals):
#   - glucose_mmol is glucose in mmol/L (NHANES glucose mg/dL / 18.016)
#   - log_crp is ln(CRP in mg/dL) (NHANES CRP is in mg/dL for most cycles)
#   - All other biomarkers in conventional US units
LEVINE_COEFFICIENTS: dict[str, float] = {
    "intercept": -19.9067,
    "albumin": -0.0336,  # g/dL
    "creatinine": 0.0095,  # mg/dL
    "glucose_mmol": 0.1953,  # mmol/L
    "log_crp": 0.0954,  # ln(mg/dL)
    "lymphocyte_pct": -0.0120,  # %
    "mcv": 0.0268,  # fL
    "rdw": 0.3306,  # %
    "alp": 0.00188,  # U/L
    "wbc": 0.0554,  # 10^9/L
    "age": 0.0804,  # years
}

# Gompertz model parameters
GOMPERTZ_GAMMA = 0.0076927  # shape parameter (rate of mortality acceleration)
GOMPERTZ_LAMBDA = 0.0553  # scale (not directly used in PhenoAge formula)

# PhenoAge conversion constants
# PhenoAge = PHENOAGE_INTERCEPT + ln(-PHENOAGE_LOG_SCALE * ln(1-M)) / PHENOAGE_RATE
PHENOAGE_INTERCEPT = 141.50225
PHENOAGE_LOG_SCALE = 0.00553
PHENOAGE_RATE = 0.09165

# Months of follow-up for 10-year mortality window
MORTALITY_WINDOW_MONTHS = 120

# =============================================================================
# Linked Mortality File definitions
# =============================================================================

# Column layout from the official CDC SAS/Stata read-in programs
# (SAS_ReadInProgramAllSurveys.sas, NHANES version, May 2022).
# 0-indexed, half-open intervals for pd.read_fwf.
#
# Note: columns 7-14 and 22-42 are blank in NHANES files (those positions
# are used only in the NHIS version for PUBLICID, DODQTR, DODYEAR, WGT_NEW,
# SA_WGT_NEW).
MORTALITY_FILE_COLSPECS: list[tuple[int, int]] = [
    (0, 6),  # SEQN          (SAS cols 1-6)
    (14, 15),  # eligstat      (SAS col 15)
    (15, 16),  # mortstat      (SAS col 16)
    (16, 19),  # ucod_leading  (SAS cols 17-19)
    (19, 20),  # diabetes      (SAS col 20)
    (20, 21),  # hyperten      (SAS col 21)
    (42, 45),  # permth_int    (SAS cols 43-45)
    (45, 48),  # permth_exm    (SAS cols 46-48)
]

MORTALITY_FILE_COLUMNS: list[str] = [
    "SEQN",
    "eligstat",
    "mortstat",
    "ucod_leading",
    "diabetes",
    "hyperten",
    "permth_int",
    "permth_exm",
]

# All NHANES linked mortality files available from CDC
MORTALITY_FILES: dict[str, str] = {
    "NHANES_III": "NHANES_III_MORT_2019_PUBLIC.dat",
    "1999-2000": "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
    "2001-2002": "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
    "2003-2004": "NHANES_2003_2004_MORT_2019_PUBLIC.dat",
    "2005-2006": "NHANES_2005_2006_MORT_2019_PUBLIC.dat",
    "2007-2008": "NHANES_2007_2008_MORT_2019_PUBLIC.dat",
    "2009-2010": "NHANES_2009_2010_MORT_2019_PUBLIC.dat",
    "2011-2012": "NHANES_2011_2012_MORT_2019_PUBLIC.dat",
    "2013-2014": "NHANES_2013_2014_MORT_2019_PUBLIC.dat",
    "2015-2016": "NHANES_2015_2016_MORT_2019_PUBLIC.dat",
    "2017-2018": "NHANES_2017_2018_MORT_2019_PUBLIC.dat",
}


# =============================================================================
# BioAge-unit constants (used by reproduction & analysis modules)
#
# These are the same Levine coefficients as LEVINE_COEFFICIENTS above, but
# keyed by the BioAge R package column names (albumin in g/L, creatinine in
# µmol/L, glucose in mmol/L).  The reproduction and analysis code builds
# data with these column names, so they import these directly.
# =============================================================================

LEVINE_BIOAGE: dict[str, float] = {
    "intercept": -19.90667,
    "albumin_gL": -0.03359355,
    "creat_umol": 0.009506491,
    "glucose_mmol": 0.1953192,
    "lncrp": 0.09536762,
    "lymph": -0.01199984,
    "mcv": 0.02676401,
    "rdw": 0.3306156,
    "alp": 0.001868778,
    "wbc": 0.05542406,
    "age": 0.08035356,
}

LEVINE_GAMMA = 0.0076927

# Feature column names in BioAge units (order matters for coefficient vectors)
FEATURE_COLS: list[str] = [
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


# =============================================================================
# Gompertz negative log-likelihood (shared by fitting & analysis code)
# =============================================================================


def gompertz_nll(params, X, time, event):
    """Gompertz PH negative log-likelihood (flexsurv parameterisation).

    Parameters
    ----------
    params : array-like
        [log(gamma), intercept, beta_1, ..., beta_p]
    X : ndarray, shape (n, p)
        Feature matrix.
    time : ndarray, shape (n,)
        Follow-up time in months.
    event : ndarray, shape (n,)
        Event indicator (1=death, 0=censored).
    """
    gamma = np.exp(params[0])  # keep gamma > 0
    intercept = params[1]
    beta = params[2:]
    xb = intercept + X @ beta
    b = np.exp(np.clip(xb, -500, 500))
    gt = np.clip(gamma * time, -500, 500)
    log_surv = -(b / gamma) * (np.exp(gt) - 1)
    log_haz = xb + gamma * time
    ll = np.sum(event * log_haz) + np.sum(log_surv)
    return -ll if np.isfinite(ll) else 1e20
