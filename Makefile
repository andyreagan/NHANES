.PHONY: mort_data all

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

data/processed/LMF_Files/LMF__2003-2004__2005-2006__MORT_2019.parquet : src/parse_lmf/parse_lmf.py mort_data
	mkdir -p data/processed/LMF_Files
	uv run $<

################################################################################
# Raw NHANES files
################################################################################

data/raw/%.XPT :
	# convert the %.XPT which will come like "2003-2004/DEMO_C.XPT"
	# into the new format "2003/DataFiles/DEMO_C.xpt"
	mkdir -p $(dir $@)
	wget -O $@ https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/$(firstword $(subst -, ,$(patsubst %/,%,$(dir $*))))/DataFiles/$(notdir $*).xpt --no-use-server-timestamps

################################################################################
# Processed activity
################################################################################

# Wildcard on the year, with the fixed suffixes
data/raw/%/paxraw_c_met_worn_bouts_reliable.parquet: src/paxraw/main.ipynb
	uv run --with papermill papermill $< $*-run.ipynb -p year $*
	uv run --with jupyter jupyter nbconvert $< --to pdf
data/raw/%/paxraw_d_met_worn_bouts_reliable.parquet: src/paxraw/main.ipynb
	uv run --with papermill papermill $< $*-run.ipynb -p year $*
	uv run --with jupyter jupyter nbconvert $< --to pdf

# Simple versions,
data/processed/%/PAXRAW_C_simple.parquet : src/paxraw/01_simple.py data/raw/%/paxraw_c.xpt
	uv run $< $*
data/processed/%/PAXRAW_D_simple.parquet : src/paxraw/01_simple.py data/raw/%/paxraw_d.xpt
	uv run $< $*

################################################################################
# M3S covariates
################################################################################

data/processed/nhanes.parquet : src/m3s_variables/generate_M3S_NHANES_dataset.py src/m3s_variables/main.py src/m3s_variables/constants.py src/m3s_variables/util.py data/processed/LMF_Files/LMF__2003-2004__2005-2006__MORT_2019.parquet
	uv run -m src.m3s_variables.generate_M3S_NHANES_dataset

all: data/processed/nhanes.parquet
