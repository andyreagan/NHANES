# Process NHANES data

This repository contains Python scripts that processes NHANES data.
They focus on the 2003-2006 data that contains activity files,
with utilities to process the accelerometer data.
In addition to activity processing,
there are utilities to process the demographic
and bloodwork data into a format that works with the public lifescore API,
and to reproduce Levine et al. (2018) PhenoAge biological aging measure.

## Activity processing

In `src/paxraw` we provide a series of scripts to process the PAXRAW
activity data files from the 2003-2006 NHANES cohorts.
You can run the `01_pipeline.ipynb` files by setting the year variable
and running through.
The `03_analysis.ipynb` can also run for each year.
Copies of the run notebooks are in the `output/paxraw` folder.

## PhenoAge (Biological Aging)

The code in `src/phenoage/` reproduces Levine et al. (2018) [2] PhenoAge — a biological
age estimate derived from 9 clinical chemistry biomarkers + chronological age.

**Module structure:**

```
src/phenoage/
├── constants.py              # All shared constants, coefficients, column definitions
├── model.py                  # Gompertz PH model fitting & PhenoAge computation
├── load_data.py              # NHANES III + IV data loading & harmonization
├── run_phenoage.py           # Main pipeline (make phenoage)
├── REPRODUCTION_NOTES.md     # Detailed findings from coefficient reproduction
├── reproduction/             # Investigating Levine 2018 coefficient matching
│   ├── explore_filters.py    #   Filter config iteration & Gompertz comparison
│   ├── find_42.py            #   Search for completeness filter yielding n=9,926
│   └── validate_models.py    #   Compare Levine/naive/winsorized on NHANES IV
└── analysis/                 # Downstream model comparisons & response curves
    ├── rsf_analysis.py       #   Random Survival Forest + SHAP + optimal profiles
    ├── rsf_expanded.py       #   RSF with extended feature sets
    ├── age_vs_survival.py    #   Age-prediction vs survival-prediction comparison
    ├── tabpfn_analysis.py    #   TabPFN foundation model for mortality
    ├── single_variable_curves.py  # Multi-model response curves (Gompertz/GAM/RSF/TabPFN)
    ├── gam_analysis.R        #   GAM survival analysis (R mgcv)
    └── rsf_r_comparison.R    #   RSF comparison (R ranger + randomForestSRC)
```

**Data coverage:**
- **Training**: NHANES III (1988-1994) lab data (~15k complete cases) + linked mortality
- **Validation/Scoring**: NHANES IV continuous (1999-2000 through 2017-2018, 10 cycles)
  - PhenoAge requires CRP, which is unavailable in 2011-2014 cycles
  - ~22,500 scored across 8 available cycles

**Pipeline** (`uv run -m src.phenoage.run_phenoage` or `make phenoage`):
1. Downloads NHANES III lab.dat (fixed-width) + all NHANES IV XPT files + linked mortality
2. Harmonizes variable names and units across cycles (glucose→mmol/L, CRP→mg/dL)
3. Fits a Gompertz proportional hazards model on NHANES III (10-year mortality)
4. Compares fitted coefficients with Levine's published values
5. Scores PhenoAge on all NHANES IV cycles using both fitted and published coefficients
6. Saves results to `data/processed/phenoage/`

**Reproduction** (`uv run -m src.phenoage.reproduction.explore_filters`):
Iterates over configurable filter combinations (fasting threshold, outlier removal,
mortality definition, CRP transform, etc.) and fits a Gompertz model on each,
comparing all coefficients to Levine's published values.  Edit the `CONFIGS` list
in `src/phenoage/reproduction/explore_filters.py` to try new combinations.  See
[`src/phenoage/REPRODUCTION_NOTES.md`](src/phenoage/REPRODUCTION_NOTES.md) for
detailed findings and [`research_agents/README.md`](research_agents/README.md) for
the 16-agent parallel investigation log.

**Biomarkers used** (9 + age):
albumin (g/dL), creatinine (mg/dL), glucose (mmol/L), ln(CRP mg/dL),
lymphocyte %, mean cell volume (fL), red cell distribution width (%),
alkaline phosphatase (U/L), white blood cell count (10⁹/L), chronological age

**Key results:**
- PhenoAge–age correlation: 0.93-0.94 across validation cycles (matching Levine)
- Mean PhenoAge acceleration: +2.7 to +3.6 years (1999-2010 cycles)
- 2015+ cycles show ~+7 year acceleration due to RDW assay method change (known issue)
- All model variants (Levine published, naive MLE, winsorized MLE) predict
  mortality equivalently (C-index within 0.002, PhenoAge r > 0.99 pairwise)

**Note on units:** The Levine 2018 coefficients expect glucose in mmol/L (not mg/dL)
and ln(CRP) where CRP is in mg/dL. This was verified by confirming PhenoAge ≈ chronological
age for individuals with population-median biomarker values.

[2] Levine, M. E., Lu, A. T., Quach, A., et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. Aging, 10(4), 573-591.

## Mortality file processing

In `src/parse_lmf` there is a single script that handles the
fixed-width `.dat` format of the mortality followup files.
It writes out a parquet file.

## Lifescore processing

The code in `src/m3s_variables` focuses on formatting available NHANES data into the format
required for the MassMutual Mortality Score,
as made available publicly through myLifeScore.
For more details on the model itself,
see the writeup [here](https://f.hubspotusercontent40.net/hubfs/5627392/LifeScore%20Labs_Med360.pdf) [1].

[1] Maier, M., Carlotto, H., Saperstein, S., Sanchez, F., Balogun, S., & Merritt, S. (2020). The Accuracy and Transparency of Underwriting with Artificial Intelligence to Transform the Life-Insurance Industry.

### Lifescore API

There is a public API at https://api.lifescoremodels.com/general/mls360
that serves the "myLifeScore" tool.
An example is worth 1000 words:

```
➜  ~ curl 'https://api.lifescoremodels.com/general/mls360' \
  -H 'accept: application/json, text/plain, */*' \
  -H 'accept-language: en-US,en;q=0.9' \
  -H 'content-type: application/json; charset=UTF-8' \
  -H 'priority: u=1, i' \
  -H 'referer: https://www.lifescoremodels.com/' \
  -H 'sec-ch-ua: "Google Chrome";v="129", "Not=A?Brand";v="8", "Chromium";v="129"' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'sec-ch-ua-platform: "macOS"' \
  -H 'sec-fetch-dest: empty' \
  -H 'sec-fetch-mode: cors' \
  -H 'sec-fetch-site: same-site' \
  -H 'user-agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/129.0.0.0 Safari/537.36' \
  --data-raw '{
    "data": {
      "current_smoker": "non-smoker",
      "total_cholesterol": "borderline",
      "cardiovascular_history": "non-cardi",
      "cancer_history": "non-cancer",
      "year_of_birth": 1974,
      "blood_pressure_systolic": 100,
      "blood_pressure_diastolic": 80,
      "blood_pressure": "normal-bp",
      "gender": "male",
      "heart_condition": "non-heart",
      "height": 62,
      "weight": 200,
      "albumin": 5.2
    },
    "pagesOrder": [
      "current_smoker",
      "total_cholesterol",
      "cardiovascular_history",
      "cancer_history",
      "year_of_birth",
      "blood_pressure",
      "gender",
      "heart_condition",
      "height",
      "weight"
    ]
  }'
```

Response:

```
{
  "score": 91,
  "contributions": [
    {
      "group": "build",
      "visual": "-",
      "contribution": -4,
      "variables": [
        {"var": "weight", "value": 200, "median": 197},
        {"var": "bmi", "value": 36.5765, "median": 27.47}
      ]
    },
    {
      "group": "bp_pulse",
      "visual": " ",
      "contribution": -1,
      "variables": [
        {"var": "blood_pressure_systolic", "value": 100, "median": 120},
        {"var": "blood_pressure_diastolic", "value": 80, "median": 78},
        {"var": "pulse_standard_at_rest", "value": 66, "median": 66}
      ]
    },
    {
      "group": "blood_protein",
      "visual": " ",
      "contribution": -2,
      "variables": [
        {"var": "albumin", "value": 5.2, "median": 4.6},
        {"var": "globulin", "value": 2.4, "median": 2.4},
        {"var": "albumin_globulin_ratio", "value": 2.1667, "median": 1.92},
        {"var": "total_protein", "value": 7.6, "median": 7}
      ]
    },
    {
      "group": "lipids",
      "visual": " ",
      "contribution": -1,
      "variables": [
        {"var": "cholesterol", "value": 219.5, "median": 196},
        {"var": "high_density_lipoproteinhdl", "value": 50.8, "median": 50.8},
        {"var": "low_density_lipoprotein_ldl", "value": 144.3, "median": 120.8},
        {"var": "cholesterol_hdl_ratio", "value": 4.3209, "median": 3.81}
      ]
    },
    {"group": "other", "contribution": -1}
  ],
  "lifescore_id": "91f6fe2a-0ca5-4b89-8b47-9eae36d1bbcf",
  "model_input": {
    "age": 50,
    "sex": "M",
    "mv_or": "NA",
    "cancer": "NA",
    "height": 62,
    "weight": 200,
    "albumin": 5.2,
    "glucose": "NA",
    "smoking": "NA",
    "globulin": "NA",
    "sgot_ast": "NA",
    "sgpt_alt": "NA",
    "hemolysis": "NA",
    "hiv_1_eia": "NA",
    "creatinine": "NA",
    "fam_cancer": "f",
    "cholesterol": 219.5,
    "fam_diabetes": "NA",
    "fam_vascular": "f",
    "fructosamine": "NA",
    "triglycerides": "NA",
    "urinary_tract": "NA",
    "blood_disorder": "NA",
    "current_smoker": "f",
    "hemoglobin_a1c": "NA",
    "urn_creatinine": "NA",
    "heart_condition": "f",
    "total_bilirubin": "NA",
    "disability_claim": "NA",
    "mental_condition": "NA",
    "urn_microalbumin": "NA",
    "muscular_disorder": "NA",
    "urn_total_protein": "NA",
    "endocrine_disorder": "NA",
    "cocaine_metabolites": "NA",
    "digestive_condition": "NA",
    "alkaline_phosphatase": "NA",
    "anti_hcv_hepatitis_c": "NA",
    "respiratory_disorder": "NA",
    "reproductive_disorder": "NA",
    "pulse_standard_at_rest": "NA",
    "blood_pressure_systolic": 100,
    "blood_urea_nitrogen_bun": "NA",
    "nervous_system_disorder": "NA",
    "blood_pressure_diastolic": 80,
    "nicotine_metabolites_urn": "NA",
    "gamma_glutamyltransferase": "NA",
    "prostate_specific_antigen": "NA",
    "high_density_lipoproteinhdl": "NA"
  }
}
```
