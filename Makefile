.PHONY: mort_data mort_data_all mort_data_processed all phenoage manifest download-all download-lab download-demo

RAW=data/raw
PROCESSED=data/processed

################################################################################
# Linked mortality files
################################################################################

# Add dependency on the py file, so they stay more freshly downloaded
data/raw/LMF_Files/%.dat : src/parse_lmf/parse_lmf.py
	mkdir -p data/raw/LMF_Files
	wget -O $@ https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/$(notdir $@) --no-use-server-timestamps

mort_data : data/raw/LMF_Files/NHANES_2003_2004_MORT_2019_PUBLIC.dat data/raw/LMF_Files/NHANES_2005_2006_MORT_2019_PUBLIC.dat

# All mortality files (for PhenoAge analysis)
MORT_ALL_FILES := $(foreach y,III 1999_2000 2001_2002 2003_2004 2005_2006 2007_2008 2009_2010 2011_2012 2013_2014 2015_2016 2017_2018,data/raw/LMF_Files/NHANES_$(y)_MORT_2019_PUBLIC.dat)
mort_data_all : $(MORT_ALL_FILES)

# Process ALL mortality files into parquet — both 2019 and 2015 vintages.
# The 2015 vintage is downloaded from the Internet Archive by the script itself.
# Also writes the legacy 2003-2006 subset for backwards compatibility.
data/processed/LMF_Files/LMF_all_MORT_2019.parquet data/processed/LMF_Files/LMF_all_MORT_2015.parquet &: src/parse_lmf/parse_lmf.py $(MORT_ALL_FILES)
	mkdir -p data/processed/LMF_Files
	uv run python -m src.parse_lmf.parse_lmf

# Legacy output (subset of all) is produced as a side-effect of the above
data/processed/LMF_Files/LMF__2003-2004__2005-2006__MORT_2019.parquet : data/processed/LMF_Files/LMF_all_MORT_2019.parquet

mort_data_processed : data/processed/LMF_Files/LMF_all_MORT_2019.parquet

################################################################################
# Raw NHANES files
################################################################################

data/raw/%.XPT :
	# convert the %.XPT which will come like "2003-2004/DEMO_C.XPT"
	# into the new format "2003/DataFiles/DEMO_C.xpt"
	mkdir -p $(dir $@)
	wget -O $@ https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/$(firstword $(subst -, ,$(patsubst %/,%,$(dir $*))))/DataFiles/$(notdir $*).xpt --no-use-server-timestamps

################################################################################
# NHANES III raw files (fixed-width .dat, different distribution from CDC)
################################################################################

data/raw/NHANES_III/%.dat :
	mkdir -p $(dir $@)
	wget -O $@ https://wwwn.cdc.gov/Nchs/Data/Nhanes3/1a/$(notdir $@) --no-use-server-timestamps

################################################################################
# Processed activity
################################################################################

# Wildcard on the year, with the fixed suffixes
# Uses the split notebooks: 01_pipeline for data processing (production)
# 02_algorithm_development for speed comparisons (dev only)
# 03_analysis for visualizations (optional)
data/processed/2003-2004/paxraw_c_met_worn_bouts_reliable.parquet: src/paxraw/01_pipeline.ipynb src/paxraw/utils.py
	cd src/paxraw && uv run --extra paxraw --with papermill papermill 01_pipeline.ipynb /dev/stdout -p year 2003-2004 -k python3 > /dev/null
data/processed/2005-2006/paxraw_d_met_worn_bouts_reliable.parquet: src/paxraw/01_pipeline.ipynb src/paxraw/utils.py
	cd src/paxraw && uv run --extra paxraw --with papermill papermill 01_pipeline.ipynb /dev/stdout -p year 2005-2006 -k python3 > /dev/null

################################################################################
# M3S covariates
################################################################################

# The pipeline now processes ALL NHANES IV cycles (1999-2000 through 2021-2023).
# Individual cycle outputs go to data/processed/{cycle}/nhanes.parquet.
# Combined output goes to data/processed/nhanes.parquet.
data/processed/nhanes.parquet : src/m3s_variables/generate_M3S_NHANES_dataset.py src/m3s_variables/main.py src/m3s_variables/constants.py src/m3s_variables/util.py data/processed/LMF_Files/LMF_all_MORT_2019.parquet
	uv run -m src.m3s_variables.generate_M3S_NHANES_dataset

################################################################################
# PhenoAge (Levine 2018 reproduction)
################################################################################

data/processed/phenoage/nhanes_iv_phenoage_scored.parquet : src/phenoage/run_phenoage.py src/phenoage/model.py src/phenoage/load_data.py src/phenoage/constants.py data/processed/LMF_Files/LMF_all_MORT_2019.parquet
	uv run -m src.phenoage.run_phenoage

phenoage : data/processed/phenoage/nhanes_iv_phenoage_scored.parquet

all: data/processed/nhanes.parquet

################################################################################
# Bulk NHANES download (all cohorts, all components)
################################################################################

# Step 1: Scrape the CDC data pages to discover all available XPT files.
# Produces a CSV manifest listing every downloadable file with its URL.
data/nhanes_file_manifest.csv : src/download/scrape_nhanes_manifest.py
	uv run python -m src.download.scrape_nhanes_manifest

manifest : data/nhanes_file_manifest.csv

# Step 2: Download all files listed in the manifest.
# Skips files that are already present locally. Parallelized (4 workers).
download-all : data/nhanes_file_manifest.csv
	uv run python -m src.download.download_nhanes_files

# Convenience targets for downloading subsets:
download-lab : data/nhanes_file_manifest.csv
	uv run python -m src.download.download_nhanes_files --components Laboratory

download-demo : data/nhanes_file_manifest.csv
	uv run python -m src.download.download_nhanes_files --components Demographics

# Download only the files needed by the M3S variable mapping pipeline.
# This is much smaller than downloading everything (~200 files vs ~1600).
M3S_FILE_BASES := DEMO SMQ MCQ BPX BMX DUQ HSQ DIQ PFQ CIQDEP DPQ OSQ \
	SMQRTU RHQ RDQ KIQ_U WHQ MPQ BIOPRO HEPC TCHOL HDL HIV TRIGLY \
	ALB_CR GHB PSA L40 L02 L13 L13AM L10 L03 L11PSA L16 L24PP \
	LAB18 LAB02 LAB13 LAB13AM LAB10 LAB03 LAB16 LAB11PSA SMQMEC \
	CBC CRP GLU HSCRP BPXO

download-m3s-deps : data/nhanes_file_manifest.csv
	uv run python -m src.download.download_nhanes_files --files $(M3S_FILE_BASES)

.PHONY: validate-phenoage
validate-phenoage : research_agents/phenoage_train_2015.csv
	uv run -m src.phenoage.reproduction.validate_models

.PHONY: rsf-analysis
rsf-analysis :
	uv run -m src.phenoage.analysis.rsf_analysis
	Rscript src/phenoage/analysis/rsf_r_comparison.R

.PHONY: tabpfn-analysis
tabpfn-analysis :
	uv run -m src.phenoage.analysis.tabpfn_analysis

.PHONY: gam-analysis
gam-analysis : research_agents/phenoage_train_2015.csv
	Rscript src/phenoage/analysis/gam_analysis.R

.PHONY: single-variable-curves
single-variable-curves :
	uv run -m src.phenoage.analysis.single_variable_curves
