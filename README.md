# NHANES Analysis

Processing and analysis of CDC [NHANES](https://www.cdc.gov/nchs/nhanes/index.htm)
(National Health and Nutrition Examination Survey) data, covering three workstreams:

1. **PhenoAge** — Reproduce Levine et al. (2018) biological aging measure and compare
   against modern ML models (RSF, TabPFN, GAM)
2. **LifeScore** — Format NHANES biomarker data for the MassMutual Mortality Score API
3. **Accelerometry** — Process raw PAXRAW accelerometer data from 2003-2006 cohorts

All data is downloaded from CDC on demand via `make` targets. Nothing is checked into
`data/` — the `.gitignore` excludes it.

## Quick Start

```bash
# Install dependencies
uv sync

# Run the PhenoAge pipeline (downloads data, fits model, scores all cycles)
make phenoage

# Run the LifeScore variable mapping (all NHANES IV cycles)
make all

# Process accelerometer data for a single cohort
make data/processed/2003-2004/paxraw_c_met_worn_bouts_reliable.parquet
```

## Repository Structure

```
├── src/
│   ├── phenoage/               # PhenoAge biological aging (Levine 2018)
│   │   ├── constants.py        #   Shared constants, coefficients, Gompertz NLL
│   │   ├── model.py            #   Gompertz PH fitting & PhenoAge computation
│   │   ├── load_data.py        #   NHANES III + IV data loading & harmonization
│   │   ├── run_phenoage.py     #   Main pipeline entry point
│   │   ├── REPRODUCTION_NOTES.md  # Detailed reproduction findings
│   │   ├── reproduction/       #   Coefficient matching investigation
│   │   │   ├── explore_filters.py  # Filter config iteration & Gompertz comparison
│   │   │   ├── find_42.py     #   Search for completeness filter yielding n=9,926
│   │   │   └── validate_models.py  # Compare model variants on NHANES IV
│   │   └── analysis/           #   Downstream model comparisons
│   │       ├── rsf_analysis.py #   Random Survival Forest + SHAP
│   │       ├── rsf_expanded.py #   RSF with extended feature sets
│   │       ├── age_vs_survival.py  # Age-prediction vs survival-prediction
│   │       ├── tabpfn_analysis.py  # TabPFN foundation model
│   │       ├── single_variable_curves.py  # Multi-model response curves
│   │       ├── gam_analysis.R  #   GAM survival analysis (R mgcv)
│   │       └── rsf_r_comparison.R  # RSF comparison in R
│   ├── m3s_variables/          # LifeScore (MassMutual Mortality Score) formatting
│   ├── parse_lmf/              # CDC linked mortality file parser
│   ├── download/               # Bulk NHANES file discovery & download
│   └── paxraw/                 # Accelerometer data processing (Jupyter notebooks)
├── research_agents/            # 16 AI agent investigation logs (PhenoAge reproduction)
├── references/                 # Papers (Levine 2018, Liu 2018), LifeScore whitepaper,
│                               #   NHANES cohort catalog, mortality file documentation
├── data/                       # (gitignored) Raw XPT/dat files + processed parquet
├── output/                     # (gitignored) HTML charts, CSVs, PDFs from analyses
├── Makefile                    # All build targets with dependency tracking
└── pyproject.toml              # Python dependencies (managed by uv)
```

## PhenoAge (Biological Aging)

Reproduces Levine et al. (2018) [1] PhenoAge — a biological age estimate derived
from 9 clinical chemistry biomarkers + chronological age, trained on NHANES III
mortality data using a Gompertz proportional hazards model.

**Biomarkers** (9 + chronological age):
albumin, creatinine, glucose, ln(CRP), lymphocyte %, mean cell volume,
red cell distribution width, alkaline phosphatase, white blood cell count

### Pipeline

```bash
make phenoage                  # Full pipeline: download → fit → score
# or equivalently:
uv run -m src.phenoage.run_phenoage
```

1. Downloads NHANES III lab data (fixed-width) + NHANES IV XPT files + linked mortality
2. Harmonizes variable names and units across cycles (glucose→mmol/L, CRP→mg/dL)
3. Fits Gompertz PH model on NHANES III (10-year mortality, n≈10k, ~1,800 deaths)
4. Scores PhenoAge on NHANES IV (1999-2018, 10 cycles, ~22k scored)
5. Saves to `data/processed/phenoage/`

### Reproduction investigation

We conducted an exhaustive investigation into matching Levine's published
coefficients, documented in [`src/phenoage/REPRODUCTION_NOTES.md`](src/phenoage/REPRODUCTION_NOTES.md)
and the 16 parallel agent logs in [`research_agents/`](research_agents/README.md).

**Bottom line:** All coefficient variants produce PhenoAge scores that correlate
r > 0.99 with each other and predict mortality equivalently (C-index within 0.002).
The remaining differences are due to irreducible data composition differences in
the NHANES III files, not implementation errors.

### Analysis commands

```bash
make validate-phenoage         # Compare Levine/naive/winsorized models on NHANES IV
make rsf-analysis              # Random Survival Forest + SHAP importance
make tabpfn-analysis           # TabPFN foundation model comparison
make gam-analysis              # GAM survival analysis (R)
make single-variable-curves    # Multi-model biomarker response curves
```

### Key results

- PhenoAge–age correlation: 0.93-0.94 across validation cycles (matching Levine)
- All models (linear Gompertz, RSF, TabPFN) give similar C-index (~0.86)
- Linear Gompertz provides the right inductive bias for NHANES-scale data (~10k samples)
- 2015+ cycles show ~+7 year acceleration artifact due to RDW assay method change

**Unit note:** Glucose in mmol/L (not mg/dL) and ln(CRP) where CRP is in mg/dL.

[1] Levine, M. E., Lu, A. T., Quach, A., et al. (2018). An epigenetic biomarker
of aging for lifespan and healthspan. *Aging*, 10(4), 573-591.

## LifeScore (MassMutual Mortality Score)

Formats NHANES biomarker data for the [myLifeScore](https://www.lifescoremodels.com/)
API, a publicly available mortality prediction model from MassMutual/LifeScore Labs [2].

```bash
make all                       # Process all NHANES IV cycles (1999-2023)
make download-m3s-deps         # Download only the XPT files needed
```

Processes demographics, blood pressure, body measurements, smoking status,
medical conditions, and lab results across all NHANES IV cycles into the format
expected by the LifeScore API. See [`references/LifeScore Labs_Med360.pdf`](references/LifeScore%20Labs_Med360.pdf)
for the model whitepaper.

<details>
<summary>LifeScore API example</summary>

```bash
curl 'https://api.lifescoremodels.com/general/mls360' \
  -H 'content-type: application/json' \
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
    }
  }'
```

Returns a score (0-100) with per-group contributions (build, bp_pulse,
blood_protein, lipids, etc.).

</details>

[2] Maier, M., et al. (2020). The Accuracy and Transparency of Underwriting
with Artificial Intelligence to Transform the Life-Insurance Industry.

## Accelerometry

Processes raw PAXRAW accelerometer data from NHANES 2003-2006 cohorts
(minute-level activity counts for ~15k participants over 7 days).

```bash
make data/processed/2003-2004/paxraw_c_met_worn_bouts_reliable.parquet
make data/processed/2005-2006/paxraw_d_met_worn_bouts_reliable.parquet
```

Pipeline notebooks in `src/paxraw/`:
- `01_pipeline.ipynb` — MET conversion, wear detection, bout identification, reliability filtering
- `02_algorithm_development.ipynb` — Speed benchmarks for the processing algorithms
- `03_analysis.ipynb` — Visualization of activity patterns
- `04_analyze_steps.ipynb` — Step count analysis

Rendered HTML copies of run notebooks are in `output/paxraw/`.

## Data Infrastructure

### NHANES file download

```bash
make manifest                  # Scrape CDC to discover all available XPT files
make download-all              # Download everything (~1,600 files)
make download-lab              # Download only Laboratory component
make download-demo             # Download only Demographics component
```

### Linked mortality files

```bash
make mort_data_all             # Download all mortality .dat files (2019 vintage)
make mort_data_processed       # Parse into parquet (2019 + 2015 vintages)
```

The parser (`src/parse_lmf/`) handles multiple vintages (2011, 2015, 2019) of the
CDC fixed-width mortality follow-up files. The 2015 vintage is recovered from the
Internet Archive for PhenoAge reproduction purposes. See
[`references/linked_mortality_vintages.md`](references/linked_mortality_vintages.md)
for full documentation of all known vintages.
