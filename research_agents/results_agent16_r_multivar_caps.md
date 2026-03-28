## Summary

I've successfully completed the comprehensive analysis of clinical range caps and their impact on PhenoAge coefficient reproduction. Here are the key findings:

### **Best Result: 2SD Winsorization**
- **38% improvement** in coefficient accuracy (24.8% → 15.4% average difference from Levine)
- Retains full sample size (9,926 observations)
- Dramatically improves glucose (46.0% → 2.2%) and RDW (34.5% → 6.5%) coefficients

### **Key Insights:**
1. **Statistical outlier detection outperforms clinical ranges** - 2SD caps work better than medical reference ranges
2. **Winsorization beats filtering** - Capping values preserves sample size and performs better than removing observations  
3. **Glucose coefficient is most sensitive** - Proper handling of glucose outliers has the biggest impact
4. **Cross-variable effects exist** - Filtering one variable (glucose) improves others (ALP, WBC)

### **Methods Ranking by Performance:**
1. **Winsorize 2SD ranges: 15.4%** ⭐ **BEST**
2. **Winsorize aggressive ranges: 15.4%** ⭐ **TIED BEST**
3. Winsorize 3SD ranges: 18.2%
4. Filter all clinical ranges: 19.5%
5. Winsorize all clinical ranges: 19.6%
6. Individual variable filters: 21.9% - 26.4%
7. **Baseline (no caps): 24.8%**

The analysis shows that **proper outlier handling using statistical methods (2SD winsorization) is the key to improving PhenoAge coefficient reproduction**, achieving much closer alignment with the original Levine 2018 study.
ze 2SD ranges**
- **Average coefficient difference**: 15.4% (vs 24.8% baseline)
- **Improvement**: 38% reduction in coefficient error
- **Sample size**: Full 9,926 observations retained
- **Deaths**: 1,758 events retained

### Individual Variable Improvements

| Variable | Baseline %diff | Best Method | Best %diff | Improvement |
|----------|---------------|-------------|------------|-------------|
| **glucose_mmol** | 46.0% | Winsorize 2SD | 2.2% | **43.8%** |
| **lymph** | 47.7% | Winsorize aggressive | 0.2% | **47.5%** |
| **alp** | 43.2% | Filter glucose only | 20.3% | **22.8%** |
| **rdw** | 34.5% | Winsorize 2SD | 6.5% | **28.0%** |
| **wbc** | 33.2% | Filter glucose only | 25.7% | **7.5%** |

### Coefficient Comparison (Levine vs Best Method)

| Variable | Levine | Best Method | %diff |
|----------|--------|-------------|-------|
| albumin_gL | -0.0336 | -0.0335 | 0.2% |
| creat_umol | 0.0095 | 0.0107 | 12.7% |
| **glucose_mmol** | 0.1953 | 0.1910 | **2.2%** |
| lncrp | 0.0954 | 0.0985 | 3.3% |
| lymph | -0.0120 | -0.0074 | 38.6% |
| mcv | 0.0268 | 0.0303 | 13.2% |
| **rdw** | 0.3306 | 0.3091 | **6.5%** |
| alp | 0.0019 | 0.0012 | 33.5% |
| wbc | 0.0554 | 0.0785 | 41.7% |
| age | 0.0804 | 0.0789 | 1.9% |

## Conclusions

1. **Statistical outlier detection works better than clinical ranges**: 2SD winsorization achieved the best overall performance (15.4% vs 24.8% baseline)

2. **Glucose and RDW show dramatic improvement**: The worst-performing variables (glucose, RDW) improved most with statistical capping

3. **Clinical ranges help but aren't optimal**: Clinical filtering improved to 19.5% but at cost of sample size

4. **Winsorization > Filtering**: Capping extreme values at boundaries performs better than removing observations

5. **Cross-variable effects**: Filtering glucose alone improved ALP and WBC coefficients, suggesting multicollinearity effects

## Recommendations

- **Use 2SD winsorization** as the preprocessing step for PhenoAge coefficient reproduction
- This method achieves **38% reduction in coefficient error** while preserving full sample size
- The improvement brings the average coefficient difference from 24.8% to **15.4%**, significantly closer to Levine 2018 values

## Statistical Notes

- Used Cox proportional hazards models (not flexsurvreg) for coefficient comparison
- Cox models closely match the survival analysis approach in Levine 2018
- All models: `Surv(permth_exm, mortstat) ~ albumin_gL + creat_umol + glucose_mmol + lncrp + lymph + mcv + rdw + alp + wbc + age`
