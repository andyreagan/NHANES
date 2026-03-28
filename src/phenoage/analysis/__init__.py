"""PhenoAge downstream analysis — model comparisons and response curves.

Modules:
    rsf_analysis            — Random Survival Forest + single-variable impact charts
    rsf_expanded            — Expanded RSF with additional NHANES features
    age_vs_survival         — Age-prediction RSF vs survival RSF comparison
    tabpfn_analysis         — TabPFN foundation model for mortality prediction
    single_variable_curves  — Multi-model response curves (Gompertz vs GAM vs RSF vs TabPFN)
    gam_analysis.R          — GAM survival analysis (R script)
    rsf_r_comparison.R      — R RSF comparison (ranger + randomForestSRC)
"""
