import pandas as pd

pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)

lemurs = pd.read_csv("data/lemurs/mm-survey-labs-sample.csv")

# need to map the columns from the lemurs.columns
# (see public https://api.lifescoremodels.com/general/mls360 for variable names)
# - age_5
# - albumin
# - albumin_globulin_ratio
# - alkaline_phosphatase
# - anti_hcv_hepatitis_c
# - ast_alt_ratio
# - blood_disorder
# - blood_pressure_diastolic_avg
# - blood_pressure_systolic_avg
# - blood_urea_nitrogen_bun
# - bmi
# - cancer
# - cholesterol
# - cholesterol_hdl_ratio
# - cocaine_metabolites
# - creatinine
# - digestive_condition
# - egfr
# - endocrine_disorder
# - fam_diabetes
# - fam_vascular
# - gamma_glutamyltransferase
# - globulin
# - heart_condition
# - hemoglobin_a1c
# - high_density_lipoproteinhdl
# - hiv_1_eia
# - low_density_lipoprotein_ldl
# - mental_condition
# - microalbumin_creatinine
# - muscular_disorder
# - nervous_system_disorder
# - nicotine_metabolites_urn
# - prostate_specific_antigen
# - pulse_standard_at_rest
# - reproductive_disorder
# - respiratory_disorder
# - sgot_ast
# - sgpt_alt
# - smoking
# - total_bilirubin
# - total_protein
# - urinary_tract
# - urn_creatinine_low
# - weight
# we were missing these in NHANES but we might have them in LEMURS:
# - urn_pc_ratio
# - urn_total_protein

# from .constants import predictors
# nhanes = pd.read_parquet('../data/processed/nhanes_m3s_vars.parquet')
# nhanes.iloc[:5, :].loc[:, predictors]
# from NHANES, these are example values for the variables that we do have ^

# inputs:

lemurs["weight_kg"] = lemurs["WEIGHT in LBS"] * 0.453592
lemurs["weight"] = lemurs["WEIGHT in LBS"]


def parse_height(x: str) -> float:
    feet, inches = x.rstrip('"').split("'")
    return int(feet) * 12 + float(inches)


def test_parse_height():
    assert parse_height('''5' 7.52"''') == 67.52


test_parse_height()

lemurs["height"] = lemurs["HEIGHT"].apply(parse_height)
lemurs["height_m"] = lemurs["height"] * 0.0254
# conver to M/F from Transgender/Female/Male
lemurs["sex"] = lemurs.apply(
    lambda x: {"Male": "M", "Female": "F", "Transgender": "NB"}[x.Gender], axis=1
)

# create the predictor columns 1 by 1 in the desired format
lemurs["age_5"] = 5 * round(lemurs["AGE"] / 5)
lemurs["albumin"] = lemurs["ALBUMIN"]
lemurs["alkaline_phosphatase"] = lemurs["ALK PHOS"]
lemurs["ast_alt_ratio"] = lemurs["AST"] / lemurs["ALT"]
lemurs["blood_disorder"] = (
    lemurs[
        "In the past 10 years, have you been diagnosed, treated, tested positive for, or been given medical advice by a member of the medical profession for a disease or disorder noted below? Select all that apply (choice=Blood, Spleen, or Immune condition)"
    ]
    != "Unchecked"
)
lemurs["blood_pressure_diastolic_avg"] = lemurs["BP_DIASTOLIC"]
lemurs["blood_pressure_systolic_avg"] = lemurs["BP_SYSTOLIC"]
lemurs["blood_urea_nitrogen_bun"] = lemurs["BUN"]
lemurs["bmi"] = lemurs["weight_kg"] / (lemurs["height_m"] ** 2)
lemurs["cancer"] = (
    lemurs[
        "In the past 10 years, have you been diagnosed, treated, tested positive for, or been given medical advice by a member of the medical profession for a disease or disorder noted below? Select all that apply (choice=Cancer or Tumor)"
    ]
    != "Unchecked"
)
lemurs["cholesterol"] = lemurs["CHOLESTEROL"]
lemurs["cholesterol_hdl_ratio"] = lemurs["CHOL/HDL RATIO"]
lemurs["creatinine"] = lemurs["CREATININE"]
# lemurs['disability_claim'] = lemurs["Have you applied for or received disability payments within the past 12 months for your headaches or migraines?"] == "Checked"  # we have this partially, at least. not in NHANES
# no actual values of the above:
lemurs["disability_claim"] = "f"
lemurs["egfr"] = None  # TODO, this is more complex calculation
# could I map to endocrine from anemia? would need a doctor on that one, maybe
# question: do we only count these if fatal?
# I didn't require this, or have an age cap (cancer > 90 shouldn't count?)
# I counted positive answers for either of the two family questions (the "or" here)
lemurs["fam_cancer"] = (
    lemurs["Condition:"] == "Cancer" or lemurs["Condition:.1"] == "Cancer"
)  # this wasn't in NHANES, but we have it (one of two cases, also disability_claim)
lemurs["fam_diabetes"] = lemurs["Condition:"] == "Diabetes" or lemurs["Condition:.1"] == "Diabetes"
lemurs["fam_vascular"] = (
    lemurs["Condition:"] == "Heart Disease"
    or lemurs["Condition:.1"] == "Heart Disease"
    or lemurs["Condition:"] == "Stroke"
    or lemurs["Condition:.1"] == "Stroke"
)
lemurs["gamma_glutamyltransferase"] = lemurs["GGT"]
lemurs.loc[lemurs["gamma_glutamyltransferase"] == "<13", "gamma_glutamyltransferase"] = 13
# note that we have other things we could include for heart condition, like AFIB, murmer, high BP, etc
# I didn't include them here:
lemurs["heart_condition"] = (
    lemurs[
        "In the past 10 years, have you been diagnosed, treated, tested positive for, or been given medical advice by a member of the medical profession for a disease or disorder noted below? Select all that apply (choice=Cardiovascular condition)"
    ]
    == "Checked"
)
lemurs["hemoglobin_a1c"] = lemurs["HEMOGLOBIN A1C"]
lemurs["high_density_lipoproteinhdl"] = lemurs["HDL"]
lemurs["low_density_lipoprotein_ldl"] = lemurs["LDL"]
lemurs["mental_condition"] = (
    lemurs["Do you have a history of an underlying mental health condition?"] == "Yes"
)
lemurs["mv_or"] = None  # not asked, not in NHANES either
lemurs["nervous_system_disorder"] = (
    lemurs[
        "In the past 10 years, have you been diagnosed, treated, tested positive for, or been given medical advice by a member of the medical profession for a disease or disorder noted below? Select all that apply (choice=Brain, Neurological, or Nervous System condition)"
    ]
    != "Unchecked"
)
lemurs["pulse_standard_at_rest"] = lemurs["PULSE"]
lemurs["sgot_ast"] = lemurs["AST"]
lemurs["sgpt_alt"] = lemurs["ALT"]
lemurs["smoking"] = (
    (lemurs["Cigarettes"] == "Checked")
    or (lemurs["E-cigarettes/Vapes"] == "Checked")
    or (
        lemurs["Other (Pipe, Chewing Tobacco, Cigars, Patch, Gum, Smoking Cessation Aids)"]
        == "Checked"
    )
    or (lemurs["Do you currently use other tobacco products?"] == "Checked")
)
lemurs["total_bilirubin"] = lemurs["BILIRUBIN, TOTAL"]
lemurs.loc[lemurs["total_bilirubin"] == "<0.5", "total_bilirubin"] = 0.5
lemurs["weight"] = lemurs["WEIGHT in LBS"]

lemurs["globulin"] = None  # don't have
lemurs["albumin_globulin_ratio"] = None  # we do not have globulin, only immunoglobulin
lemurs["total_protein"] = None  # no globulin, can't do: lemurs["ALBUMIN"] + lemurs["GLOBULIN"]

lemurs["microalbumin_creatinine"] = None  # no urine values
lemurs["hiv_1_eia"] = None  # not tested (urine?)
lemurs["urn_creatinine_low"] = None  # no urine values
lemurs["urn_pc_ratio"] = None  # no urine values (not in NHANES either)
lemurs["urn_total_protein"] = None  # no urine values (not in NHANES either)
lemurs["cocaine_metabolites"] = None  # no urine values
lemurs["nicotine_metabolites_urn"] = None  # no urine values
lemurs["prostate_specific_antigen"] = None  # don't have
lemurs["anti_hcv_hepatitis_c"] = None  # don't see it (urine?)

lemurs["digestive_condition"] = None  # not asked, not in NHANES either
lemurs["endocrine_disorder"] = None  # not asked
lemurs["muscular_disorder"] = None  # not asked
lemurs["reproductive_disorder"] = None  # not asked
lemurs["respiratory_disorder"] = None  # not asked
lemurs["urinary_tract"] = None  # not asked

lemurs.fillna("NA")

lemurs.to_parquet("data/lemurs/mm-survey-labs-sample-processed.parquet")
