# All NHANES Cohorts: Data Availability for M3S Pipeline

## Overview

The National Health and Nutrition Examination Survey has been conducted in
several waves since 1960. This document catalogs every cohort, what biomarker
data is available, whether mortality follow-up exists, and feasibility for the
M3S variable mapping pipeline.

## Cohort Summary

| Cohort | Years | N (approx) | Biomarkers | Mortality | M3S Feasibility |
|--------|-------|-----------|------------|-----------|-----------------|
| **NHES I** | 1960-1962 | 7,710 | Minimal | ❌ None public | ❌ Not feasible |
| **NHES II** | 1963-1965 | 7,417 | Minimal | ❌ None public | ❌ Not feasible |
| **NHES III** | 1966-1970 | 7,514 | Minimal | ❌ None public | ❌ Not feasible |
| **NHANES I** | 1971-1975 | 20,749 | Partial | ⚠️ NHEFS (through 1992) | ⚠️ ~20 of 45 predictors |
| **NHANES II** | 1976-1980 | 20,322 | Partial | ⚠️ Archived (2004 release) | ⚠️ ~22 of 45 predictors |
| **NHANES III** | 1988-1994 | 33,994 | **Full** | ✅ 2019 vintage | ✅ **In pipeline** |
| **NHANES IV** | 1999-2018 | 101,316 | **Full** | ✅ 2019 vintage | ✅ **In pipeline** |
| **NHANES 2021-2023** | 2021-2023 | 11,933 | **Full** | ❌ Not yet released | ⚠️ Data loads, no mortality |

## Detailed Biomarker Availability

### M3S Pipeline Predictors (45 total)

| Predictor | NHANES I | NHANES II | NHANES III+ |
|-----------|----------|-----------|-------------|
| **age_5** | ✅ | ✅ | ✅ |
| **sex** | ✅ | ✅ | ✅ |
| **albumin** | ✅ N1LB0232 | ✅ | ✅ |
| **alkaline_phosphatase** | ✅ N1LB0459 | ✅ | ✅ |
| **blood_pressure_systolic_avg** | ✅ | ✅ | ✅ |
| **blood_pressure_diastolic_avg** | ✅ | ✅ | ✅ |
| **blood_urea_nitrogen_bun** | ✅ N1LB0472 | ✅ | ✅ |
| **bmi** | ✅ | ✅ | ✅ |
| **cholesterol** | ✅ N1LB0237 | ✅ | ✅ |
| **creatinine** | ✅ N1LB0475 | ✅ | ✅ |
| **pulse_standard_at_rest** | ✅ | ✅ | ✅ |
| **weight** | ✅ | ✅ | ✅ |
| **smoking** | ✅ | ✅ | ✅ |
| **cancer** | ✅ | ✅ | ✅ |
| **endocrine_disorder** (diabetes) | ✅ | ✅ | ✅ |
| **heart_condition** | ✅ | ✅ | ✅ |
| **muscular_disorder** (arthritis) | ✅ | ✅ | ✅ |
| **respiratory_disorder** | ✅ | ✅ | ✅ |
| **egfr** (derived) | ✅ (from creatinine) | ✅ | ✅ |
| **high_density_lipoprotein (HDL)** | ❌ | ❌ | ✅ |
| **cholesterol_hdl_ratio** | ❌ | ❌ | ✅ |
| **low_density_lipoprotein_ldl** | ❌ | ❌ | ✅ |
| **triglycerides** | ❌ | ❌ | ✅ |
| **hemoglobin_a1c** | ❌ | ❌ | ✅ |
| **gamma_glutamyltransferase** (GGT) | ❌ | ❌ | ✅ |
| **sgot_ast** | ❌ | ❌ | ✅ |
| **sgpt_alt** | ❌ | ❌ | ✅ |
| **ast_alt_ratio** | ❌ | ❌ | ✅ |
| **total_bilirubin** | ⚠️ dipstick only | ⚠️ | ✅ |
| **globulin** | ❌ | ❌ | ✅ |
| **albumin_globulin_ratio** | ❌ | ❌ | ✅ |
| **total_protein** | ✅ N1LB0227 | ✅ | ✅ |
| **anti_hcv_hepatitis_c** | ❌ | ❌ | ✅ |
| **hiv_1_eia** | ❌ | ❌ | ✅ |
| **cocaine_metabolites** | ❌ | ❌ | ✅ (some cycles) |
| **nicotine_metabolites_urn** | ❌ | ❌ | ✅ (some cycles) |
| **prostate_specific_antigen** | ❌ | ❌ | ✅ (some cycles) |
| **urn_creatinine_low** | ✅ N1LB0544 | ✅ | ✅ |
| **microalbumin_creatinine** | ❌ | ❌ | ✅ |
| **mental_condition** | ⚠️ CES-D only | ⚠️ | ✅ |
| **fam_diabetes** | ✅ | ✅ | ✅ |
| **fam_vascular** | ✅ | ✅ | ✅ |
| **digestive_condition** | ⚠️ | ⚠️ | ✅ |
| **nervous_system_disorder** | ⚠️ | ⚠️ | ✅ |
| **reproductive_disorder** | ⚠️ | ⚠️ | ✅ |
| **urinary_tract** | ⚠️ | ⚠️ | ✅ |

### Predictor count by era

| Era | Available Predictors | Missing Predictors | Notes |
|-----|---------------------|-------------------|-------|
| NHANES I (1971-1975) | ~20 of 45 | HDL, LDL, triglycerides, GGT, AST, ALT, A1C, globulin, HCV, HIV, PSA, urinary albumin | Pre-lipid-panel era; limited CBC |
| NHANES II (1976-1980) | ~22 of 45 | Same as NHANES I, except may have slightly better CBC | Marginally better lab panel |
| NHANES III (1988-1994) | **43 of 45** | Only cocaine and nicotine metabolites sparse | Full modern clinical chemistry |
| NHANES IV (1999-2018) | **45 of 45** | CRP unavailable 2011-2014 (PhenoAge affected, not M3S) | Complete |

## Data Format by Cohort

| Cohort | Format | Key Files | Download |
|--------|--------|-----------|----------|
| NHANES I | Fixed-width .txt + SAS layouts | DU4800 (biochem), DU4081 (medical hx), DU4233 (exam), DU4111 (anthropometry) | `https://wwwn.cdc.gov/Nchs/Data/nhanes1/DU####.txt` |
| NHANES II | Fixed-width .txt + SAS layouts | DU5411 (hematology/biochem), DU5020 (medical hx), DU5302 (physician exam), DU5301 (anthropometry) | `https://wwwn.cdc.gov/Nchs/Data/nhanes2/DU####.txt` |
| NHANES III | Fixed-width .dat + SAS layouts | adult.dat, exam.dat, lab.dat | `https://wwwn.cdc.gov/Nchs/Data/Nhanes3/1a/` |
| NHANES IV | SAS transport .XPT per topic | DEMO, BIOPRO, CBC, BPX, BMX, SMQ, MCQ, etc. | `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/` |

## Mortality Linkage by Cohort

| Cohort | Public Mortality File | Follow-up Through | Source |
|--------|----------------------|-------------------|--------|
| NHANES I | **NHEFS** N92mort.txt | Dec 31, 1992 | CDC: `data/nhefs/N92mort.txt` |
| NHANES I | NDI linkage | Dec 31, 2019 | ⚠️ Restricted-use only (NCHS RDC) |
| NHANES II | 2004 release | ~Dec 31, 2000 | ⚠️ Was on CDC FTP; archived; may need ICPSR |
| NHANES III | 2019 vintage | Dec 31, 2019 | ✅ CDC FTP (in our pipeline) |
| NHANES IV (1999-2018) | 2019 vintage | Dec 31, 2019 | ✅ CDC FTP (in our pipeline) |
| NHANES 2021-2023 | Not yet released | — | Expected ~2027 |

### NHEFS Mortality File Details

The NHANES I Epidemiologic Follow-up Study (NHEFS) conducted follow-up
interviews in 1982-84, 1986, 1987, and 1992, with vital status ascertained
at each wave. The 1992 mortality file (N92mort.txt) contains:

- Vital status through December 31, 1992
- Cause of death (ICD-9 coded)
- Date of death
- ~21 years of follow-up from baseline (1971-1975)

**Key limitation:** The NHEFS mortality file uses a different format and
different variable names than the 2019 NHANES linked mortality files. It
would need a separate parser (similar to what we built for NHANES III).

**Files on CDC:**
- Data: `https://wwwn.cdc.gov/Nchs/Data/nhefs/N92mort.txt`
- Layout: `https://wwwn.cdc.gov/Nchs/Data/nhefs/mort.inputs.labels.txt`
- Documentation: `https://wwwn.cdc.gov/Nchs/Data/nhefs/1992mort.pdf`

### NHANES II Mortality Situation

The 2004 CDC mortality release included NHANES II linkage
(`matching_methodology_nhanes2_final.pdf` is still on the FTP archive).
However, the actual public-use data files are no longer on the CDC FTP.

**Possible sources:**
- ICPSR (Inter-university Consortium for Political and Social Research)
- Internet Archive Wayback Machine
- Direct request to NCHS

## Recommendations

### Current pipeline (implemented)

- **NHANES III (1988-1994):** 13,825 rows, 5,556 deaths — full biomarker panel
- **NHANES IV (1999-2018):** ~52,000 rows, ~7,300 deaths — full biomarker panel
- **Total: ~66,000 rows, ~12,900 deaths**

### Potential additions (not yet implemented)

1. **NHANES I (1971-1975):** Would add ~14,000 adults with NHEFS mortality,
   but only ~20 of 45 M3S predictors available. Would require:
   - Parser for DU4800.txt (biochemistry) and DU4081/4233 (exam/medical history)
   - Parser for NHEFS N92mort.txt mortality file
   - Lowered predictor threshold (from 40 to ~18)
   - Careful unit conversion (1970s lab methods differ)

2. **NHANES II (1976-1980):** Would add ~15,000 adults, ~22 predictors.
   Would require:
   - Parser for DU5411.txt (hematology/biochemistry) and DU5020/5302 (exam)
   - Recovery of mortality linkage file from archives
   - Same threshold and unit conversion caveats as NHANES I

3. **NHANES 2021-2023:** Already loads (11,933 participants). Will become
   useful once CDC releases the mortality linkage (expected ~2027).

### For a reduced-biomarker model

If we built a model using only the ~20 biomarkers available across ALL cohorts
(NHANES I through present), we could pool ~95,000+ participants with mortality
follow-up spanning 1971-2019 — nearly 50 years of data. The available
"universal" biomarkers are:

- Albumin, total protein, cholesterol, creatinine, BUN, alkaline phosphatase
- Blood pressure (systolic, diastolic), pulse, BMI, weight
- WBC, hematocrit/hemoglobin, lymphocyte %
- Glucose, urinary creatinine
- Smoking, diabetes, heart disease, arthritis, cancer (self-report)
- Age, sex, race/ethnicity, education

This is essentially the "Framingham-era" biomarker panel and would be a
scientifically interesting historical comparison.
