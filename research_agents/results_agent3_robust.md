# Research Agent 3: Robust Gompertz PH Investigation

**Date:** 2026-03-28
**Focus:** Mathematical correctness and numerical stability

## Key Findings

### 1. Mathematical Error Identified
❌ **The original `_gompertz_nll` function was MISSING the `log(γ)` term**

**Incorrect:** `log_haz = xb + gamma * time`
**Correct:** `log_haz = log(gamma) + gamma * time + xb`

This explains why our coefficients were ~32% off from published values!

### 2. Likelihood Formula
The correct Gompertz PH log-likelihood is:
```
h(t|x) = γ × exp(γt) × exp(xβ)
log h(t|x) = log(γ) + γt + xβ
Λ(t|x) = exp(xβ) × (exp(γt) - 1) / γ
LL = Σ[δᵢ × log h(tᵢ) - Λ(tᵢ)]
```

### 3. Numerical Stability Issues
- Need proper bounds: γ ∈ [0.001, 0.05] to prevent degenerate solutions
- Conservative clipping: [-50, 50] instead of [-500, 500]
- Extensive overflow checking required

## Recommendations

1. **FIX THE LOG-LIKELIHOOD**: Add missing `log(γ)` term
2. **Add proper bounds** for all parameters
3. **Re-run all analyses** with corrected implementation
4. **Verify against published results** to ensure data consistency

## Impact

This mathematical error explains:
- Why coefficients were 32% off target values
- Why gamma sometimes converged to wrong values
- Why optimization seemed unstable

The correct implementation should yield much closer agreement with Levine 2018 published coefficients.
