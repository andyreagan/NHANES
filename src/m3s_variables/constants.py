from itertools import chain

# metadata for which NHANES files have which variables
# this could be built automatically, of course, by crawling all of the NHANES files
# we'll just include the subset that we need
# we need this lookup in the first place so that we only load variables that we need
# NHANES would be too wide otherwise
NHANES_VAR_FILES: dict = {
    "BIOPRO": [
        "LBXSAL",
        "LBDSALSI",
        "LBXSGB",
        "LBDSGBSI",
        "LBXSAPSI",
        "LBXSASSI",
        "LBXSATSI",
        "LBXSBU",
        "LBXSCR",
        "LBXSTR",
        "LBXSTB",
        "LBXSTP",
        "LBXSGTSI",
    ],
    "DEMO": [
        "RIDAGEYR",
        "RIDAGEMN",
        "RIAGENDR",
        "RIDAGEEX",
        "WTMEC2YR",
        "INDHHINC",
        "RIDRETH1",
        "DMDEDUC2",
    ],
    "SMQ": ["SMQ040", "SMD070", "SMQ020", "SMQ620", "SMQ660"],
    "MCQ": [
        "MCQ300A",
        "MCQ300C",
        "MCQ053",
        "MCQ220",
        "MCQ160F",
        "MCQ160L",
        "MCQ160A",
        "MCQ160B",
        "MCQ160C",
        "MCQ160D",
        "MCQ160E",
        "MCQ160G",
        "MCQ160K",
    ],
    "HEPC": ["LBDHCV"],
    "BPX": [f"BPX{t}{i+1}" for t in {"DI", "SY"} for i in range(4)] + ["BPXPLS", "BPXPULS"],
    "BMX": ["BMXBMI", "BMXWT", "BMXHT"],
    "TCHOL": ["LBXTC"],
    "HDL": ["LBDHDD"],
    "DUQ": ["DUQ240"],
    "HSQ": ["HSQ510"],
    "DIQ": [
        "DIQ010",
        "DIQ280",
    ],  # the first is diabetes, both years, latter is A1C, only 2005-2006
    "HIV": ["LBDHI"],
    "TRIGLY": ["LBDLDL"],
    "PFQ": ["PFD069D", "PFD069E", "PFD069K", "PFQ061B", "PFQ061C"],
    # face-to-face depression interview only in 2003:
    "CIQDEP": ["CIDDSCOR"],
    # no pain screen only in 2003:
    "MPQ": ["MPQ120AB"],
    # depression screen only in 2005, use just last question:
    "DPQ": ["DPQ090"],
    "OSQ": ["OSD030ca", "OSQ010c"],
    "SMQRTU": ["SMQ680"],
    "RHQ": ["RHQ360"],
    # keep only the medication question
    # https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/RDQ_D.htm#RDQ134
    "RDQ": ["RDQ134"],
    "KIQ_U": ["KIQ022"],
    "ALB_CR": ["URXUCR", "URXUMA"],  # these are from urine
    "WHQ": ["WHD020"],  # note: this is pounds
    "GHB": ["LBXGH"],
    "PSA": ["LBXP1"],
    # just note the vars that we'll load from the linked mortality file:
    "mortality": ["eligstat", "mortstat", "permth_exm", "ucod_leading"],
    # this gets created by the loading script
    "other": ["year"],
}
FILE_MAPPING_PRE2005: dict = {
    "BIOPRO": "L40",
    "HEPC": "L02",
    "TCHOL": "L13",
    "HDL": "L13",
    "GHB": "L10",
    "HIV": "L03",
    "TRIGLY": "L13AM",
    "SMQMEC": "SMQRTU",
    "KIQ_U": "L24PP",
    "ALB_CR": "L16",
    "L11PSA": "PSA",
}
VAR_MAPPING_PRE2005: dict = {
    "MCQ": {"MCQ300C": "MCQ250A", "MCQ300A": "MCQ250G"},
    "HDL": {"LBDHDD": "LBXHDD"},
    "DUQ": {"DUQ240": "DUQ100"},
    # "DIQ": {"DIQ280": "LBXGH"},
}

# for imputing,
# - we impute all predictors that are not "required" (if required and NA, we'd drop that row)
# - we also go ahead and impute things that have impute=True
# - after these imputations, we could recompute all the combos

TF = {"expr": "ifelse(_ == 1, 't', 'f')", "impute_value": "f", "predictor": True}
RN = {"expr": "ifelse(_ == 1, 'REA', 'NON')", "impute_value": "NON", "predictor": True}
PN = {"expr": "ifelse(_ == 1, 'POS', 'NEG')", "impute_value": "NEG", "predictor": True}

# store the lookup for m3s covariates as data
# so that we can loop through and analyze
LOOKUP: dict = {
    "age_5_hh": {"depends_on": ["RIDAGEYR"], "expr": "5 * round(_ / 5)"},
    "age_5_mn": {"depends_on": ["RIDAGEMN"], "expr": "5 * round(_ / 12 / 5)"},
    "age_5_ex": {"depends_on": ["RIDAGEEX"], "expr": "5 * round(_ / 12 / 5)"},
    "age": {
        "depends_on": ["RIDAGEEX", "RIDAGEMN", "RIDAGEYR"],
        "expr": "coalesce(_[0] / 12, _[1] / 12, _[2])",
        "na_if_dependency_na": "all",
    },
    "age_5": {
        "depends_on": ["age"],
        "expr": "5 * round(_ / 5)",
        "na_if_dependency_na": "all",
        "required": True,
        "predictor": True,
    },
    "income": {"depends_on": ["INDHHINC"]},
    "exam_offset": {"depends_on": ["RIDAGEEX", "RIDAGEYR"], "expr": "_[0] - _[1]*12"},
    "survey_offset": {
        "depends_on": ["RIDAGEMN", "RIDAGEYR"],
        "expr": "_[0] - _[1]*12",
    },
    "race": {
        "depends_on": ["RIDRETH1"],
        "expr": '["Mexican American", "Other Hispanic", "White", "Black", "Other"][int(_) - 1]',
    },
    # need to create a year column for these ones:
    "exposure_2011": {"depends_on": ["year"], "expr": "(2011 - _ + 1)*12"},
    "exposure_2011_offset": {
        "depends_on": ["exposure_2011", "exam_offset", "survey_offset"],
        "na_if_dependency_na_subset": ["exposure_2011"],
        "expr": "_[0] - coalesce(_[1], _[2], 0)",
    },
    "exposure_2015": {"depends_on": ["year"], "expr": "(2015 - _ + 1)*12"},
    "exposure_2015_offset": {
        "depends_on": ["exposure_2015", "exam_offset", "survey_offset"],
        "na_if_dependency_na_subset": ["exposure_2015"],
        "expr": "_[0] - coalesce(_[1], _[2], 0)",
    },
    # need the LMF for these ones, don't impute because we always have it:
    "in_death_file": {
        "depends_on": ["eligstat"],
        "expr": "(_ == 1)*1",
        "required": True,
    },
    "accidental_death": {
        "depends_on": ["ucod_leading"],
        "expr": "(_ == 4)*1",
        "impute": True,
        "impute_value": 0,
    },
    # let this coalesce in the full exposure for the 546 missing values
    # but careful: 172 of those without expsosure have died (we just don't know when)
    "exposure": {
        "depends_on": ["permth_exm", "exposure_2015_offset"],
    },
    # just rename:
    "exp_integral": {
        "depends_on": ["exposure"],
        "required": True,
    },
    # always have this one, but filter to only use it where we have exposure
    # by setting to NA where exposure is NA and there is a death
    # (keep the people that are assumed alive - we can give them full exposure)
    "is_dead": {
        "depends_on": ["permth_exm", "mortstat"],
        "expr": "ifelse(isna(_[0]) and (_[1] == 1), None, _[1])",
        "required": True,
    },
    "is_dead_2011": {
        "depends_on": ["exposure", "exposure_2011_offset", "is_dead"],
        "expr": "ifelse(_[0] > _[1], 0, _[2])",
    },
    "censored_exposure_2011": {
        "depends_on": ["exposure_2011_offset", "exposure"],
        "expr": "min(_[0], _[1])",
    },
    "albumin_siunit": {"depends_on": ["LBDSALSI"]},
    "albumin": {
        "depends_on": ["LBXSAL"],
        "valid_range": [None, 6],
        "predictor": True,
    },
    "sex": {
        "depends_on": ["RIAGENDR"],
        "expr": "ifelse(_ == 1, 'M', 'F')",
        "required": True,
    },
    "smoking_self_report_current": {"depends_on": ["SMQ040"], "expr": "_ <= 2"},
    "smoking_self_report_quit": {
        "depends_on": ["SMQ020", "SMQ040"],
        "expr": "ifelse(_[1] == 1, (_[0] == 3)*1, 0)",
    },
    # since they had to have answered that they smoked at least 100 cigarettes
    # lifetime to get asked the current smoker question,
    # we can just check on whether they aren't currently smoking
    # smoking should just be the current self report
    "smoking": {
        "depends_on": ["SMQ040"],
        "expr": "ifelse(_ <= 2, 't', 'f')",
        "impute_value": "f",
        "predictor": True,
    },
    "fam_diabetes": {"depends_on": ["MCQ300C"], **TF},
    "fam_vascular": {"depends_on": ["MCQ300A"], **TF},
    "alkaline_phosphatase": {"depends_on": ["LBXSAPSI"], "predictor": True, "caps": [0, 1000]},
    "anti_hcv_hepatitis_c": {"depends_on": ["LBDHCV"], **RN, "sample": False},
    "blood_disorder": {"depends_on": ["MCQ053"], **TF},
    **{
        f"blood_pressure_diastolic_{i+1}": {
            "depends_on": [f"BPXDI{i+1}"],
            "valid_range": [30, 130],
            "impute": i < 1,
        }
        for i in range(4)
    },
    **{
        f"blood_pressure_systolic_{i+1}": {
            "depends_on": [f"BPXSY{i+1}"],
            "valid_range": [50, 220],
            "impute": i < 1,
        }
        for i in range(4)
    },
    "blood_pressure_diastolic_avg": {
        "depends_on": [f"blood_pressure_diastolic_{i+1}" for i in range(4)],
        "expr": "nanmean(_)",
        "na_if_dependency_na": "all",
        "valid_range": [30, 130],
        "predictor": True,
        "post_impute": True,
    },
    "blood_pressure_systolic_avg": {
        "depends_on": [f"blood_pressure_systolic_{i+1}" for i in range(4)],
        "expr": "nanmean(_)",
        "na_if_dependency_na": "all",
        "valid_range": [50, 220],
        "predictor": True,
        "post_impute": True,
    },
    "blood_urea_nitrogen_bun": {"depends_on": ["LBXSBU"], "predictor": True},
    "height": {"depends_on": ["BMXHT"], "impute": True, "expr": "_ / 2.54"},
    "weight_self_report": {"depends_on": ["WHD020"], "valid_range": [70, 450]},
    "weight": {
        "depends_on": ["BMXWT"],
        "expr": "_ * 2.2",
        "valid_range": [70, 450],
        "predictor": True,
    },
    "bmi_nhanes": {"depends_on": ["BMXBMI"]},
    "bmi": {
        "depends_on": ["height", "weight"],
        "expr": "_[1] * 703 / (_[0] ** 2)",
        "predictor": True,
        "post_impute": True,
    },
    "cancer": {"depends_on": ["MCQ220"], **TF},
    "stroke": {
        "depends_on": ["MCQ160F"],
        "expr": "ifelse(_ == 1, 't', 'f')",
        "impute_value": "f",
        "predictor": False,
        "impute": True,
    },
    "mobility_limitation": {
        "depends_on": ["PFQ061B", "PFQ061C"],
        "expr": "ifelse(any(c(_[0] == 2, _[0] == 3, _[0] == 4, _[1] == 2, _[1] == 3, _[1] == 4)), 't', 'f')",
        "na_if_dependency_na": "all",
        "impute": True,
        "impute_value": "f",
        "predictor": False,
    },
    "less_than_highschool": {
        "depends_on": ["DMDEDUC2"],
        "expr": "ifelse(any(c(_ == 1, _ == 2)), 't', 'f')",
        "impute_value": "f",
        "predictor": False,
        "impute": True,
    },
    "highschool_inclu_GED": {
        "depends_on": ["DMDEDUC2"],
        "expr": "ifelse(_ == 3, 't', 'f')",
        "impute_value": "f",
        "predictor": False,
        "impute": True,
    },
    "greater_than_highschool": {
        "depends_on": ["DMDEDUC2"],
        "expr": "ifelse(any(c(_ == 4, _ == 5)), 't', 'f')",
        "impute_value": "f",
        "predictor": False,
        "impute": True,
    },
    "non_hispanic_black": {
        "depends_on": ["RIDRETH1"],
        "expr": "ifelse(_ == 4, 't', 'f')",
        "impute_value": "f",
        "predictor": False,
        "impute": True,
    },
    "hispanic": {
        "depends_on": ["RIDRETH1"],
        "expr": "ifelse(any(c(_ == 1, _ == 2)), 't', 'f')",
        "impute_value": "f",
        "predictor": False,
        "impute": True,
    },
    "cholesterol": {
        "depends_on": ["LBXTC"],
        "valid_range": [70, 600],
        "predictor": True,
    },
    "high_density_lipoproteinhdl": {
        "depends_on": ["LBDHDD"],
        "valid_range": [10, None],
        "predictor": True,
    },
    "cholesterol_hdl_ratio": {
        "depends_on": ["cholesterol", "high_density_lipoproteinhdl"],
        "expr": "_[0] / _[1]",
        "predictor": True,
        "post_impute": True,
    },
    "cocaine_metabolites": {"depends_on": ["DUQ240"], **PN, "sample": False},
    "creatinine": {"depends_on": ["LBXSCR"], "predictor": True, "caps": [0, 10]},
    # "HSQ510" is whether stomach sick in past 30 days, way more common (1641 positive)
    # "MCQ160L" is diagnosed liver condition:
    "digestive_condition": {"depends_on": ["MCQ160L"], **TF},
    "disability_claim": None,
    "egfr": {
        "depends_on": ["creatinine", "age_5", "sex"],
        "expr": "175 * _[0]**-1.154 * _[1]**-0.203 * ifelse(_[2] == 'F', 0.742, 1)",
        "post_impute": True,
        "predictor": True,
        "caps": [0, 250],
    },
    # the question being used here is whether they have diabetes...
    "endocrine_disorder": {"depends_on": ["DIQ010"], **TF},
    # could use fam prostate cancer: MCQ265
    "fam_cancer": None,
    "gamma_glutamyltransferase": {"depends_on": ["LBXSGTSI"], "predictor": True, "caps": [0, 1000]},
    "globulin_siunit": {"depends_on": ["LBDSGBSI"]},
    "globulin": {"depends_on": ["LBXSGB"], "predictor": True, "caps": [0, 7]},
    "albumin_globulin_ratio": {
        "depends_on": ["albumin", "globulin"],
        "expr": "_[0] / _[1]",
        "post_impute": True,
        "predictor": True,
        "caps": [0, 6],
    },
    # BPXPULS is irregular pulse...
    # "heart_condition": {"depends_on": ["BPXPULS"]},
    "heart_condition": {
        "depends_on": ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E"],
        "expr": "ifelse(nanmin(_)==1, 't', 'f')",
        "na_if_dependency_na": "all",
        "predictor": True,
        "impute_value": "f",
    },
    #
    "hemoglobin_a1c": {"depends_on": ["LBXGH"], "predictor": True},
    # DIQ280 is self report, has values of 999 for don't know, only a few values
    # and only a single year of data!
    # "hemoglobin_a1c": {"depends_on": ["DIQ280"], "predictor": True},
    "hiv_1_eia": {"depends_on": ["LBDHI"], **RN, "sample": False},
    "triglycerides": {
        "depends_on": ["LBXSTR"],
        "expr": "_",
        "valid_range": [0, None],
        "impute": True,
        "caps": [0, 2000],
    },
    # we have the raw,
    # and it exists in 561 cases where we don't have all the formula inputs
    # and the correlation checks out
    "low_density_lipoprotein_ldl_nhanes": {"depends_on": ["LBDLDL"]},
    # so fill in:
    "low_density_lipoprotein_ldl": {
        "depends_on": [
            "cholesterol",
            "triglycerides",
            "high_density_lipoproteinhdl",
            "low_density_lipoprotein_ldl_nhanes",
        ],
        "expr": "coalesce(_[0] - pmin(_[1]/5, 400) - _[2], _[3])",
        "predictor": True,
        "post_impute": True,
    },
    # depression: PFD069D <= 66666
    # also check year-specific depression screens
    "mental_condition": {
        "depends_on": ["PFD069D", "CIDDSCOR", "DPQ090"],
        "expr": "ifelse(any(c(_[0] <= 66666, _[1] == 1, _[2] == 2, _[2] == 3)), 't', 'f')",
        "na_if_dependency_na": "all",
        "impute_value": "f",
        "predictor": True,
    },
    # this is arthritis...but ok
    "muscular_disorder": {"depends_on": ["MCQ160A"], **TF},
    "mv_or": None,
    # the first is for spine pain, only in first year
    # the second is age at spine fracture
    # "nervous_system_disorder": {"depends_on": ["MPQ120AB", "OSD030ca"]},
    # for spine pain, the values are 37 when true (139 people)
    # "spine_pain": {"depends_on": ["MPQ120AB"], "expr": "(_ == 37)*1"}
    # "spine_pain": {"depends_on": ["MPQ120AB"], "expr": "(_ > 0)*1"}
    # instead let's use the yes/no on spine fracture,
    # so we don't have a cohort bias
    # but, turns out, both OSD030ca and OSQ010c are null on all participants
    "nervous_system_disorder": {"depends_on": ["OSQ010c"], **TF},
    "nicotine_metabolites_urn": {"depends_on": ["SMQ680"], **PN, "sample": False},
    "prostate_specific_antigen": {
        "depends_on": ["LBXP1", "RIAGENDR"],
        # set to 0 for females
        "expr": "ifelse((_[1] == 'F') or (_[0] < 0.5), 'NEG', 'POS')",
        "predictor": True,
    },
    "pulse_standard_at_rest": {
        "depends_on": ["BPXPLS"],
        "predictor": True,
        "valid_range": [40, 150],
    },
    # could impute medians ^ with: "impute_medians": {'age_5': [20, 25, 30], 'sex': 'M', 'pulse_standard_at_rest': 100}},
    # just the one question on endometriosis
    "reproductive_disorder": {"depends_on": NHANES_VAR_FILES["RHQ"], **TF},
    # been told emphysema, chronic bronchitis
    # ["MCQ160G", "MCQ160K"]
    # days have lung/breathing problem
    # PFD069K <= 66666
    # report lung/breating problem:
    # PFQ063A,B,C,D,E == 21
    "respiratory_disorder": {
        "depends_on": ["MCQ160G", "MCQ160K"] + ["PFD069K"] + NHANES_VAR_FILES["RDQ"],
        "expr": "ifelse(any(c(_[0]==1,_[1]==1,_[2]<=66666,_[3]==1)), 't', 'f')",
        "na_if_dependency_na": "all",
        "impute_value": "f",
        "predictor": True,
    },
    "sgot_ast": {"depends_on": ["LBXSASSI"], "predictor": True, "caps": [0, 800]},
    "sgpt_alt": {"depends_on": ["LBXSATSI"], "predictor": True, "caps": [0, 800]},
    "ast_alt_ratio": {
        "depends_on": ["sgot_ast", "sgpt_alt"],
        "expr": "_[0] / _[1]",
        "predictor": True,
        "post_impute": True,
    },
    "total_bilirubin": {"depends_on": ["LBXSTB"], "predictor": True, "caps": [0, 6]},
    # the raw value here exactly matches the summed values
    "total_protein_nhanes": {"depends_on": ["LBXSTP"]},
    "total_protein": {
        "depends_on": ["albumin", "globulin"],
        "expr": "_[0] + _[1]",
        "predictor": True,
    },
    "urinary_tract": {"depends_on": ["KIQ022"], **TF},
    "urn_creatinine": {"depends_on": ["URXUCR"], "impute": True},
    "urn_creatinine_low": {
        "depends_on": ["urn_creatinine"],
        "expr": "ifelse(_ <= 10, 'POS', 'NEG')",
        "impute_value": "NEG",
        "predictor": True,
    },
    "urn_microalbumin": {
        "depends_on": ["URXUMA"],
        "impute": True,
        "expr": "_ / 10",
        "caps": [0, 700],
    },
    # "urn_total_protein": {"depends_on": ["urn_microalbumin"], "predictor": True},
    # "urn_pc_ratio": {
    #     # note, we're just using urine albumin as urine total protein:
    #     "depends_on": ["urn_microalbumin", "urn_creatinine"],
    #     "expr": "_[0] / _[1]",
    #     "predictor": True,
    #     "post_impute": True
    # },
    "urn_total_protein": None,
    "urn_pc_ratio": None,
    "microalbumin_creatinine": {
        "depends_on": ["urn_microalbumin", "urn_creatinine"],
        "expr": "_[0] / _[1]",
        "predictor": True,
        "post_impute": True,
    },
    "population_weight": {"depends_on": ["WTMEC2YR"], "required": True},
}

required_non_predictors = set(
    [
        k
        for k, v in LOOKUP.items()
        if v is not None and (not v.get("predictor", False) and v.get("required", False))
    ]
)
required_predictors = set(
    [
        k
        for k, v in LOOKUP.items()
        if v is not None and (v.get("predictor", False) and v.get("required", False))
    ]
)
required = set([k for k, v in LOOKUP.items() if v is not None and v.get("required", False)])

imputed_non_predictors = set(
    [
        k
        for k, v in LOOKUP.items()
        if v is not None and (not v.get("predictor", False) and v.get("impute", False))
    ]
)

predictors = set([k for k, v in LOOKUP.items() if v is not None and v.get("predictor", False)])

post_imputed_predictors = set(
    [
        k
        for k, v in LOOKUP.items()
        if v is not None and (v.get("predictor", False) and v.get("post_impute", False))
    ]
)

post_imputed_predictor_deps = set(
    chain(*[LOOKUP[predictor]["depends_on"] for predictor in post_imputed_predictors])
)
deps_not_imputed = post_imputed_predictor_deps - predictors - imputed_non_predictors - required

median_but_do_not_sample = [
    var
    for var in (predictors | imputed_non_predictors) - post_imputed_predictors
    if not LOOKUP[var].get("sample", True)
]

keep = predictors | required | imputed_non_predictors

nhanes_vars_mapped = set(
    chain(
        *NHANES_VAR_FILES.values(),
        *map(lambda x: list(x.values()), VAR_MAPPING_PRE2005.values()),
    )
)
lookup_deps_all = set(chain(*[v.get("depends_on", []) for v in LOOKUP.values() if v is not None]))

missing_deps = lookup_deps_all - nhanes_vars_mapped - set(LOOKUP.keys())
assert len(missing_deps) == 0
