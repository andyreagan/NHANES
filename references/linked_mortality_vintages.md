# NHANES Linked Mortality File Vintages

## Summary

The CDC/NCHS has released multiple vintages of NHANES linked mortality files,
each extending follow-up further. This document catalogs all known vintages,
their availability, and file naming conventions.

| Vintage | Follow-up through | NHANES cycles covered | Files | Status |
|---------|-------------------|----------------------|-------|--------|
| **2010** | ~Dec 31, 2006 | NHANES III, 1999–2004 (3 cycles) | 3 | ❌ Lost (see below) |
| **2011** | Dec 31, 2011 | NHANES III, 1999–2010 (7 cycles) | 7 | ⚠️ 2 of 7 recovered (see below) |
| **2015** | Dec 31, 2015 | NHANES III, 1999–2014 (9 cycles) | 9 | ✅ Recovered from Internet Archive |
| **2019** | Dec 31, 2019 | NHANES III, 1999–2018 (11 cycles) | 11 | ✅ Current CDC FTP release |

## 2019 Vintage (current)

**Source:** https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/

**Read-in programs (current):**
- `SAS_ReadInProgramAllSurveys.sas`
- `Stata_ReadInProgramAllSurveys.do`

**File names:**
```
NHANES_III_MORT_2019_PUBLIC.dat
NHANES_1999_2000_MORT_2019_PUBLIC.dat
NHANES_2001_2002_MORT_2019_PUBLIC.dat
NHANES_2003_2004_MORT_2019_PUBLIC.dat
NHANES_2005_2006_MORT_2019_PUBLIC.dat
NHANES_2007_2008_MORT_2019_PUBLIC.dat
NHANES_2009_2010_MORT_2019_PUBLIC.dat
NHANES_2011_2012_MORT_2019_PUBLIC.dat
NHANES_2013_2014_MORT_2019_PUBLIC.dat
NHANES_2015_2016_MORT_2019_PUBLIC.dat
NHANES_2017_2018_MORT_2019_PUBLIC.dat
```

**NHANES layout (48 chars per record):**
```
SEQN            1-6
ELIGSTAT        15
MORTSTAT        16
UCOD_LEADING    17-19
DIABETES        20
HYPERTEN        21
(blank)         22-42
PERMTH_INT      43-45
PERMTH_EXM      46-48
```

## 2015 Vintage (archived)

**Source:** Internet Archive Wayback Machine (captured Jan 15, 2021)
  - Base URL: `https://web.archive.org/web/20210115id_/https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/`

**Read-in programs (still on CDC FTP under archived_files/):**
- `SAS_ReadInProgramAllSurveys_2015.sas`
- `Stata_ReadInProgramAllSurveys_2015.do`
- `R_ReadInProgramAllSurveys_2015.R`

**File names:**
```
NHANES_III_MORT_2015_PUBLIC.dat
NHANES_1999_2000_MORT_2015_PUBLIC.dat
NHANES_2001_2002_MORT_2015_PUBLIC.dat
NHANES_2003_2004_MORT_2015_PUBLIC.dat
NHANES_2005_2006_MORT_2015_PUBLIC.dat
NHANES_2007_2008_MORT_2015_PUBLIC.dat
NHANES_2009_2010_MORT_2015_PUBLIC.dat
NHANES_2011_2012_MORT_2015_PUBLIC.dat
NHANES_2013_2014_MORT_2015_PUBLIC.dat
```

**Layout:** Same 48-char layout as 2019 vintage. The SAS program uses `SEQN 1-5`
(5 chars) vs `SEQN 1-6` in 2019, but bytes 1-6 parse identically since SEQN
values ≤99999 are right-padded with a space in position 6.

**Relevance:** Levine et al. (2018) PhenoAge paper used this vintage.

## 2011 Vintage (partially recovered)

**Source of recovered files:** The `rnhanesdata` R package on GitHub
(https://github.com/andrew-leroux/rnhanesdata) bundled the 2003-2004 and
2005-2006 files as package data in `inst/extdat/mort/`.

**Recovered files (2 of 7):**
```
NHANES_2003_2004_MORT_2011_PUBLIC.dat   (566,832 bytes, 10,122 rows)
NHANES_2005_2006_MORT_2011_PUBLIC.dat   (579,488 bytes, 10,348 rows)
```

**Still missing (5 of 7):**
```
NHANES_III_MORT_2011_PUBLIC.dat
NHANES_1999_2000_MORT_2011_PUBLIC.dat
NHANES_2001_2002_MORT_2011_PUBLIC.dat
NHANES_2007_2008_MORT_2011_PUBLIC.dat
NHANES_2009_2010_MORT_2011_PUBLIC.dat
```

**Layout (recovered from rnhanesdata `process_mort()` R function):**

Records are variable-length (up to 54 chars). Different from 2015/2019!
```
SEQN            1-14    (SEQN in first ~5 chars, rest blank padding)
ELIGSTAT        15
MORTSTAT        16
CAUSEAVL        17      Cause of Death Data Available flag (extra vs 2015/2019)
UCOD_LEADING    18-20
DIABETES        21
HYPERTEN        22
(blank)         23-43
PERMTH_INT      44-46
PERMTH_EXM      47-49
MORTSRCE_NDI    50      Mortality Source: NDI Match (extra vs 2015/2019)
MORTSRCE_CMS    51      Mortality Source: CMS Information (extra)
MORTSRCE_SSA    52      Mortality Source: SSA Information (extra)
MORTSRCE_DC     53      Mortality Source: Death Certificate Match (extra)
MORTSRCE_DCL    54      Mortality Source: Data Collection (extra)
```

**Key differences from 2015/2019:**
- UCOD_LEADING at positions 18-20 (shifted +1 by CAUSEAVL at 17)
- DIABETES at 21, HYPERTEN at 22 (shifted +1)
- PERMTH_INT at 44-46, PERMTH_EXM at 47-49 (shifted +1)
- Extra columns: CAUSEAVL, MORTSRCE_NDI/CMS/SSA/DC/DCL
- Variable-length lines (15-54 chars depending on eligibility/death)

**Where to search for remaining 5 files:**
- Google: `"NHANES_1999_2000_MORT_2011_PUBLIC"` or `"MORT_2011_PUBLIC" NHANES`
- GitHub: search for R packages or data repos that bundled NHANES 2007-2008 or 2009-2010 mortality data
- Researcher hard drives — epidemiologists who published with NHANES mortality data 2013-2018
- ICPSR, OSF.io, Zenodo, Harvard Dataverse, Dryad

## 2010 Vintage (LOST)

**Last known location:** CDC FTP `linked_mortality/` directory

**Evidence:** Wayback Machine CDX shows these files were crawled on Oct 16, 2011,
but every capture returned HTTP 404. The directory listing from that date also
shows them. The read-in programs and README files were moved to `archived_files/`
and some are STILL on the CDC FTP today.

**File names (different naming convention from later vintages!):**
```
NHANES3_MORT_PUBLIC_USE_2010.DAT          (note: uppercase .DAT)
NHANES99_00_MORT_PUBLIC_USE_2010.DAT
NHANES01_02_MORT_PUBLIC_USE_2010.DAT
NHANES03_04_MORT_PUBLIC_USE_2010.DAT
NHANES05_06_MORT_PUBLIC_USE_2010.DAT      (if it existed — only III through 03-04 confirmed)
```

**Layout (NHANES III, from recovered SAS program — 24 chars per record):**
```
SEQN            1-5
ELIGSTAT        6
MORTSTAT        7
MORTSRCE_DC     8       (Mortality Source: Death Certificate Match)
MORTSRCE_NDI    9       (Mortality Source: NDI Match)
MORTSRCE_SSA    10      (Mortality Source: SSA Information)
MORTSRCE_CMS    11      (Mortality Source: CMS Information)
PERMTH_INT      12-14
PERMTH_EXM      15-17
CAUSEAVL        18      (Cause of Death Data Available)
UCOD_113        19-21   (ICD-10 113-group recode, different from UCOD_LEADING)
DIABETES        22
HYPERTEN        23
HIPFRACT        24      (Hip Fracture flag — dropped in later vintages)
```

**Key differences from 2015/2019:**
- Completely different record layout (24 chars vs 48 chars)
- Additional fields: MORTSRCE_DC/NDI/SSA/CMS, CAUSEAVL, HIPFRACT
- Cause of death uses 113-group ICD-10 recode (UCOD_113) instead of 10-group (UCOD_LEADING)
- Different file naming convention

**Still-live documentation on CDC FTP:**
- `archived_files/nhanes3_readme_mortality_2010.txt`
- `archived_files/nhanes3_sample_ascii_pgm_2010.sas`
- `archived_files/nhanes3_sample_stata_ascii_pgm_2010.do`
- `archived_files/nh99+_readme_mortality_2010.txt` (404 — the `+` causes IIS issues)
- `archived_files/nh99+_sample_ascii_pgm_2010.sas` (404)
- `archived_files/nh99+_sample_stata_ascii_pgm_2010.do` (404)

**Where to search:**
- Google: `"NHANES03_04_MORT_PUBLIC_USE_2010"` or `"MORT_PUBLIC_USE_2010" NHANES`
- Google: `"nhanes3_mort_public_use_2010"` (lowercase variant)
- GitHub code search: `MORT_PUBLIC_USE_2010`
- Researcher hard drives — especially epidemiologists who published with NHANES
  mortality data between 2010 and 2013
- ICPSR, OSF.io, Zenodo, Harvard Dataverse
- The `nh99+` README/SAS files may be recoverable if someone can URL-encode the
  `+` correctly or access the CDC FTP via actual FTP protocol rather than HTTP

## Vintage differences matter

Comparing 2003-2004 cycle across all three recovered vintages:
- **2011 vintage:** 732 deaths, max permth_exm=107 months
- **2015 vintage:** 1,094 deaths, max permth_exm=157 months
- **2019 vintage:** 1,420 deaths, max permth_exm=205 months

Comparing 2003-2004 cycle between 2015 and 2019 vintages:
- 376 people alive in 2015 vintage → deceased in 2019 vintage
- Exposure months grow for surviving subjects (150→201 months for SEQN 21005)
- Even deceased subjects can have slightly adjusted timing (~1 month shifts)
- Total deaths: 1,094 (2015) → 1,420 (2019) for the 2003-2004 cycle alone

## GitHub repos searched (March 2026)

The following repos were cloned and thoroughly searched for mortality .dat files,
.rds/.RData/.csv mortality data, and code references to older vintages. None
contained the 5 missing 2011 vintage files or any 2010 vintage data files.

| Repo | Size | Findings |
|------|------|----------|
| `andrew-leroux/rnhanesdata` | 421M | ✅ **2011 vintage: 2003-2004 and 2005-2006 .dat files** + parsing code in `process_mort()` |
| `ehsanx/Reproducible-NHANES-Analysis` | 709M | 2019 vintage .dat files only (all 10 cycles) |
| `ClaireMargaux/rforphysicians` | 159M | 2019 vintage 2011-2012 file (renamed without vintage year) |
| `lilykoff/nhanes_pa` | 128M | 2019 vintage (processed .rds, no raw .dat) |
| `dhicks/obesity` | 5M | Code references 2010+2011 vintages but no data bundled |
| `kshedden/survival_workshop` | 1.5M | Code references 2011 vintage but no data bundled |
| `niehs-prime/factor_interactions` | 311M | Renamed XPT files, no mortality data |
| `slab-itu/HTSS` | 399M | No mortality data |
| `QuQ-ToT-Orz/HUO_Project2` | 13M | Fork of rnhanesdata; process_mort docs only |
| `mel-hsw/longevity_researchpottal` | 2.6M | No mortality data |
| `BioData-Club/nhanes_explore` | 46M | No mortality data |
| `HMS-AgeVSSurvival/NHANES_preprocessing` | 516K | Code only, no data |
| `sanjaybasu/nhanesml` | 276K | Code only, no data |
| `Marowalker/retrain_w2v` | 2.5G | NLP repo, no NHANES mortality data |
| `muschellij2/rscopus` | 6.1M | Scopus R package, unrelated |
| `LincoleJ/stcs6701_final_project` | 33M | No mortality data |
| `Keniajin/US_BMI_Mortality_Trends` | 400M | One DEMO XPT, no mortality .dat |

## HuggingFace datasets searched (March 2026)

| Dataset | Findings |
|---------|----------|
| `nguyenvy/cleaned_nhanes_1988_2018` | Has `mortality_clean.csv` — **2019 vintage only** (max permth_exm=249) |
| `Healome/nhanes-validation-data` | Standard XPT files, no mortality data |
| `HHS-Official/nchs-*-mortality-*` | Aggregate statistics, not individual-level linked files |

Searches for "NHANES mortality", "linked mortality", "NCHS mortality", and "NHANES 2011" returned no additional datasets with older vintage linked mortality files.
