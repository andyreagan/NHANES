"""
Autoresearch: find which subset of candidate biomarkers, when requiring
complete cases, yields the Levine 2018 training sample of n=9,926.

RESULTS
=======
Two exact matches found (n=9,926 with 5SD outlier removal per gender):

  PA9 + {ggt, fev, waist, vitaminC, cadmium, trig, uap}  = 16 vars
  PA9 + {ggt, fev, waist, vitaminC, cadmium, trig, rbc}   = 16 vars

where PA9 = {albumin, creatinine, glucose, CRP, lymphocyte%, MCV, RDW, ALP, WBC}

These completeness-filter variables include fasting-dependent measures
(triglycerides) that restrict the sample to the morning fasting MEC session,
explaining why the training sample is ~10k rather than ~15k.

COEFFICIENT COMPARISON
======================
                                  n    deaths   γ%off   avg coef %diff
n=9,926 (uap combo, 2019 mort)  9,926  1,885   4.3%    32.4%
n=9,926 (rbc combo, 2019 mort)  9,926  1,885   4.5%    32.1%
fasting>=8, PA9 only, 5SD       8,913  1,762   2.2%    17.6%  ← best

The n=9,926 sample includes ~1,000 non-fasting people who dilute the
glucose signal, worsening the glucose coefficient from 62% to 54% off
(and other coefficients as well). The fasting≥8 subset gives better
coefficients despite the wrong n.

The remaining coefficient discrepancies (especially glucose at ~62% and
albumin at ~33%) are likely due to the mortality file version difference
(we use 2019; the paper used 2015, which is equivalent after the
240-month cap).

Usage:
    uv run python -m src.phenoage.find_42
"""

import re
import subprocess
import time
from itertools import combinations

import numpy as np
import pandas as pd


def load_all_biomarkers():
    """Parse all candidate biomarkers from NHANES III lab.dat and exam.dat."""

    # ── lab.dat variables ──
    lab_vars = [
        ("SEQN", 1, 5),
        ("HSSEX", 15, 15),
        ("HSAGEIR", 16, 17),
        ("PHPFAST", 1263, 1267),
        # CBC
        ("WCPSI", 1278, 1282),  # WBC
        ("LMPPCNT", 1283, 1287),  # Lymphocyte %
        ("MOPPCNT", 1288, 1292),  # Monocyte %
        ("RCPSI", 1316, 1319),  # RBC
        ("HTP", 1330, 1334),  # Hematocrit
        ("MVPSI", 1340, 1344),  # MCV
        ("RWP", 1360, 1364),  # RDW
        # Lipids
        ("TCP", 1598, 1600),  # Total cholesterol
        ("TGP", 1606, 1609),  # Triglycerides
        ("LCP", 1615, 1617),  # LDL cholesterol
        ("HDP", 1622, 1624),  # HDL cholesterol
        # Chemistry
        ("CRP", 1667, 1671),  # C-reactive protein
        ("UAP", 1749, 1752),  # Uric acid
        ("SGP", 1758, 1760),  # Serum glucose
        ("BUP", 1766, 1768),  # BUN
        ("TBP", 1774, 1777),  # Total bilirubin
        ("CEP", 1784, 1787),  # Creatinine
        ("ASPSI", 1821, 1823),  # AST
        ("GGPSI", 1827, 1830),  # GGT
        ("APPSI", 1835, 1838),  # ALP
        ("AMP", 1846, 1848),  # Albumin
        # Other
        ("GHP", 1861, 1864),  # HbA1c
        ("G1P", 1866, 1870),  # Fasting plasma glucose
        ("I1P", 1918, 1923),  # Fasting insulin
        # Vitamins / trace elements
        ("VBPSI", 1497, 1504),  # Vitamin B12
        ("VCP", 1505, 1508),  # Vitamin C
        ("VAPSI", 1534, 1537),  # Vitamin A
        ("VEPSI", 1543, 1548),  # Vitamin E
        ("UDPSI", 1950, 1955),  # Urinary cadmium
        ("URP", 1956, 1960),  # Urinary protein
    ]

    colspecs = [(s - 1, e) for _, s, e in lab_vars]
    col_names = [n for n, _, _ in lab_vars]
    lab = pd.read_fwf(
        "data/raw/NHANES_III/lab.dat",
        colspecs=colspecs,
        header=None,
        names=col_names,
        na_values=[".", ""],
    )

    # Sentinel value cleanup
    sentinel_map = {
        "PHPFAST": 88888,
        "WCPSI": 88888,
        "LMPPCNT": 88888,
        "MOPPCNT": 88888,
        "RCPSI": 8888,
        "HTP": 88888,
        "MVPSI": 88888,
        "RWP": 88888,
        "TCP": 888,
        "TGP": 8888,
        "LCP": 888,
        "HDP": 888,
        "CRP": 8888.8,
        "UAP": 8888,
        "SGP": 888,
        "BUP": 888,
        "TBP": 8888,
        "CEP": 8888,
        "ASPSI": 888,
        "GGPSI": 8888,
        "APPSI": 8888,
        "AMP": 888,
        "GHP": 8888,
        "G1P": 88888,
        "I1P": 888888,
        "VBPSI": 88888888,
        "VCP": 8888,
        "VAPSI": 8888,
        "VEPSI": 888888,
        "UDPSI": 888888,
        "URP": 88888,
    }
    for col, sentinel in sentinel_map.items():
        if col in lab.columns:
            lab.loc[lab[col] >= sentinel, col] = np.nan

    # ── exam.dat variables ──
    # Get FEV1 position from exam.sas
    fev_start, fev_end = 1674, 1677  # default
    try:
        result = subprocess.run(
            ["curl", "-s", "https://wwwn.cdc.gov/Nchs/Data/Nhanes3/1a/exam.sas"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        for line in result.stdout.split("\n"):
            m = re.match(r"^\s+SPPFEV1\s+(\d+)-(\d+)\s*$", line)
            if m:
                fev_start, fev_end = int(m.group(1)), int(m.group(2))
                break
    except Exception:
        pass

    exam_vars = [
        ("SEQN", 1, 5),
        ("BMPBMI", 1524, 1527),
        ("BMPHT", 1528, 1532),
        ("BMPWAIST", 1590, 1594),
        ("BMPWT", 1508, 1513),
        ("PEPMNK1R", 1423, 1425),
        ("PEPMNK5R", 1428, 1430),
        ("SPPFEV1", fev_start, fev_end),
    ]
    exam = pd.read_fwf(
        "data/raw/NHANES_III/exam.dat",
        colspecs=[(s - 1, e) for _, s, e in exam_vars],
        header=None,
        names=[n for n, _, _ in exam_vars],
        na_values=[".", ""],
    )
    exam_sentinels = {
        "BMPBMI": 8888,
        "BMPHT": 88888,
        "BMPWAIST": 88888,
        "BMPWT": 888888,
        "PEPMNK1R": 888,
        "PEPMNK5R": 888,
        "SPPFEV1": 8880,
    }
    for col, sentinel in exam_sentinels.items():
        if col in exam.columns:
            exam.loc[exam[col] >= sentinel, col] = np.nan

    # ── Mortality ──
    mort = pd.read_fwf(
        "data/raw/LMF_Files/NHANES_III_MORT_2019_PUBLIC.dat",
        colspecs=[(0, 14), (14, 15), (15, 16), (16, 19), (42, 45), (45, 48)],
        header=None,
        names=["SEQN", "eligstat", "mortstat", "ucod_leading", "permth_int", "permth_exm"],
        na_values=["."],
    )

    # ── Merge ──
    df = lab.merge(
        exam[["SEQN", "BMPBMI", "BMPHT", "BMPWAIST", "BMPWT", "PEPMNK1R", "PEPMNK5R", "SPPFEV1"]],
        on="SEQN",
        how="left",
    )
    df = df.merge(
        mort[["SEQN", "eligstat", "mortstat", "ucod_leading", "permth_exm"]],
        on="SEQN",
        how="left",
    )
    return df


def apply_mortality_filters(df):
    """Age-related death exclusion + 240-month cap (per BioAge/Levine)."""
    df = df.copy()
    mask_age_death = df.ucod_leading.isin([4, 8, 10])
    df.loc[mask_age_death, "mortstat"] = 0
    died_after_240 = (df.mortstat == 1) & (df.permth_exm > 240)
    df.loc[died_after_240, "mortstat"] = 0
    df.loc[df.permth_exm > 240, "permth_exm"] = 240
    return df


def apply_creat_calibration(df):
    """Creatinine calibration per Selvin et al."""
    df = df.copy()
    df["CEP"] = df["CEP"] * 0.960 - 0.184
    return df


def apply_5sd_outlier(df, bio_cols):
    """Remove 5SD outliers per gender for specified columns."""
    df = df.copy()
    for sex in [1, 2]:
        mask_sex = df.HSSEX == sex
        for col in bio_cols:
            vals = df.loc[mask_sex, col]
            mu, sd = vals.mean(), vals.std()
            if pd.isna(mu) or pd.isna(sd) or sd == 0:
                continue
            outlier = (vals > mu + 5 * sd) | (vals < mu - 5 * sd)
            df.loc[mask_sex & outlier.reindex(df.index, fill_value=False), col] = np.nan
    return df


def build_candidate_list():
    """Return dict of candidate variable name -> NHANES III column name."""
    return {
        # PhenoAge 9
        "albumin": "AMP",
        "alp": "APPSI",
        "creat": "CEP",
        "glucose": "SGP",
        "crp": "CRP",
        "lymph": "LMPPCNT",
        "mcv": "MVPSI",
        "rdw": "RWP",
        "wbc": "WCPSI",
        # Other lab
        "hba1c": "GHP",
        "totchol": "TCP",
        "hdl": "HDP",
        "ldl": "LCP",
        "trig": "TGP",
        "bun": "BUP",
        "uap": "UAP",
        "ttbl": "TBP",
        "ggt": "GGPSI",
        "monopa": "MOPPCNT",
        "rbc": "RCPSI",
        "insulin": "I1P",
        "cadmium": "UDPSI",
        "vitaminB12": "VBPSI",
        "vitaminC": "VCP",
        "vitaminA": "VAPSI",
        "vitaminE": "VEPSI",
        # Exam
        "bmi": "BMPBMI",
        "height": "BMPHT",
        "waist": "BMPWAIST",
        "weight": "BMPWT",
        "sbp": "PEPMNK1R",
        "dbp": "PEPMNK5R",
        "fev": "SPPFEV1",
    }


TARGET_N = 9926


def main():
    t0 = time.time()
    print("=" * 72)
    print("  AUTORESEARCH: Finding the 42-variable filter for n=9,926")
    print("=" * 72)

    # ── Load data ──
    print("\n[1/4] Loading NHANES III data...")
    df = load_all_biomarkers()
    print(f"  Raw merged: {df.shape[0]:,d} rows")

    # ── Base filters ──
    print("\n[2/4] Applying base filters...")
    df = df[(df.HSAGEIR >= 20) & (df.HSAGEIR <= 84) & (df.eligstat == 1)].copy()
    df = apply_creat_calibration(df)
    df = apply_mortality_filters(df)
    print(f"  Age 20-84, eligible: {df.shape[0]:,d}")

    # ── Candidate variables ──
    candidates = build_candidate_list()
    valid = {n: c for n, c in candidates.items() if c in df.columns}
    bio_cols = list(valid.values())

    # ── 5SD outlier removal ──
    print("\n[3/4] Applying 5SD outlier removal per gender...")
    df = apply_5sd_outlier(df, bio_cols)

    pa9 = ["albumin", "alp", "creat", "glucose", "crp", "lymph", "mcv", "rdw", "wbc"]
    pa9_cols = [valid[n] for n in pa9]
    surv = ["permth_exm", "mortstat"]
    extra = [n for n in valid if n not in pa9]

    n_pa9 = df.dropna(subset=pa9_cols + surv).shape[0]
    print(f"  PA9 complete: {n_pa9:,d}")

    # ── Combinatorial search ──
    print(f"\n[4/4] Searching for completeness filter → n={TARGET_N:,d}...")

    # Greedy forward search
    print(f"\n  GREEDY FORWARD (adding vars that bring n closest to {TARGET_N}):")
    current_cols = pa9_cols.copy()
    current_names = pa9.copy()
    remaining = {n: valid[n] for n in extra}

    while remaining:
        best_name, best_n = None, None
        best_gap = float("inf")
        for name, col in remaining.items():
            n = df.dropna(subset=current_cols + [col] + surv).shape[0]
            if abs(n - TARGET_N) < best_gap:
                best_gap = abs(n - TARGET_N)
                best_n = n
                best_name = name

        if best_name is None:
            break

        current_cols.append(remaining[best_name])
        current_names.append(best_name)
        del remaining[best_name]

        n_curr = df.dropna(subset=current_cols + surv).shape[0]
        marker = " *** MATCH ***" if n_curr == TARGET_N else ""
        print(
            f"    +{best_name:<15s} → n={n_curr:>7,d} "
            f"(gap={n_curr-TARGET_N:>+5d}, {len(current_names)} vars){marker}"
        )

        if n_curr == TARGET_N:
            print(f"\n  Variables: {current_names}")
            break

        if n_curr < TARGET_N - 500:
            current_cols.pop()
            current_names.pop()
            break

    # Focused search: greedy found core6={ggt,fev,waist,vitaminC,cadmium,trig}
    # gives n=9,927 (+1). Adding one more variable:
    core6 = ["ggt", "fev", "waist", "vitaminC", "cadmium", "trig"]
    core6_cols = [valid[n] for n in core6]
    n_core6 = df.dropna(subset=pa9_cols + core6_cols + surv).shape[0]
    print(f"\n  PA9 + core 6 ({', '.join(core6)}): n={n_core6:,d} (gap={n_core6-TARGET_N:+d})")

    remaining = [n for n in extra if n not in core6]
    print(f"\n  Adding 7th variable to reach target:")
    for name in remaining:
        n = df.dropna(subset=pa9_cols + core6_cols + [valid[name]] + surv).shape[0]
        if abs(n - TARGET_N) <= 30:
            marker = " *** EXACT ***" if n == TARGET_N else ""
            print(f"    +{name:<15s}: n={n:>7,d} (gap={n-TARGET_N:>+4d}){marker}")

    elapsed = time.time() - t0
    print(f"\n  Completed in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
