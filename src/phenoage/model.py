"""Gompertz Proportional Hazards model and PhenoAge computation.

Reproduces Levine et al. 2018 PhenoAge methodology:

1. Fit Gompertz PH on NHANES III with 10 biomarkers → get beta coefficients
2. Compute mortality score M from the linear predictor xb
3. Convert M to PhenoAge (biological age estimate)

Gompertz hazard: h(t|x) = b * exp(bt) * exp(xb)
  where b = gamma (shape), xb = X @ beta (linear predictor)

Survival: S(t|x) = exp(-exp(xb) * (exp(bt) - 1) / b)
Mortality probability at time t: M = 1 - S(t|x)

PhenoAge inverts this to find the age at which a reference individual
would have the same mortality risk.
"""

from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm

from .constants import (
    GOMPERTZ_GAMMA,
    LEVINE_COEFFICIENTS,
    MORTALITY_WINDOW_MONTHS,
    PHENOAGE_FEATURES,
    PHENOAGE_INTERCEPT,
    PHENOAGE_LOG_SCALE,
    PHENOAGE_RATE,
)


def gompertz_log_likelihood(
    params: np.ndarray,
    X: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    sample_weight: np.ndarray | None = None,
) -> float:
    """Negative log-likelihood for Gompertz proportional hazards model.

    Parameterization: h(t|x) = gamma * exp(gamma*t + X@beta)
    where params = [gamma, beta_0, beta_1, ..., beta_p]

    For the Gompertz model:
    - Log-hazard: log h(t) = log(gamma) + gamma*t + X@beta
    - Cumulative hazard: H(t) = exp(X@beta) * (exp(gamma*t) - 1) / gamma
    - Log-likelihood: sum_i [delta_i * log(h(t_i)) - H(t_i)]

    Args:
        params: Array of [gamma, beta_0, beta_1, ..., beta_p]
        X: Design matrix (n x p), should include intercept column if desired
        time: Follow-up time in months
        event: Event indicator (1=death, 0=censored)
        sample_weight: Optional survey weights

    Returns:
        Negative log-likelihood (for minimization)
    """
    gamma = params[0]
    beta = params[1:]

    # Ensure gamma > 0 (Gompertz shape parameter)
    if gamma <= 0:
        return 1e15

    xb = X @ beta
    # Clip to prevent overflow
    xb = np.clip(xb, -20, 20)
    gamma_t = gamma * time

    # Log-hazard at event/censoring time
    log_h = np.log(gamma) + gamma_t + xb

    # Cumulative hazard
    exp_gamma_t = np.exp(gamma_t)
    H = np.exp(xb) * (exp_gamma_t - 1.0) / gamma

    # Log-likelihood
    ll = event * log_h - H

    if sample_weight is not None:
        ll = ll * sample_weight

    # Return negative for minimization
    return -np.sum(ll)


def gompertz_log_likelihood_gradient(
    params: np.ndarray,
    X: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    sample_weight: np.ndarray | None = None,
) -> np.ndarray:
    """Gradient of the negative Gompertz log-likelihood.

    Args:
        params: Array of [gamma, beta_0, ..., beta_p]
        X, time, event, sample_weight: Same as gompertz_log_likelihood

    Returns:
        Gradient vector of shape (1 + p,)
    """
    gamma = params[0]
    beta = params[1:]

    if gamma <= 0:
        return np.zeros_like(params)

    n = X.shape[0]
    xb = X @ beta
    xb = np.clip(xb, -20, 20)
    exp_xb = np.exp(xb)
    gamma_t = gamma * time
    exp_gamma_t = np.exp(gamma_t)

    # Gradient w.r.t. gamma
    # d(log h)/d(gamma) = 1/gamma + t
    # d(H)/d(gamma) = exp(xb) * [t*exp(gamma*t)/gamma - (exp(gamma*t)-1)/gamma^2]
    #               = exp(xb)/gamma * [t*exp(gamma*t) - (exp(gamma*t)-1)/gamma]
    dH_dgamma = exp_xb * (time * exp_gamma_t / gamma - (exp_gamma_t - 1.0) / gamma**2)
    dll_dgamma = event * (1.0 / gamma + time) - dH_dgamma

    # Gradient w.r.t. beta
    # d(log h)/d(beta) = X
    # d(H)/d(beta) = exp(xb) * (exp(gamma*t) - 1) / gamma * X = H * X
    H = exp_xb * (exp_gamma_t - 1.0) / gamma
    dll_dbeta = X.T @ (event - H)  # This is sum_i (delta_i - H_i) * x_i

    # Wait, let me redo this more carefully
    # ll_i = delta_i * (log(gamma) + gamma*t_i + xb_i) - exp(xb_i)*(exp(gamma*t_i)-1)/gamma
    # d(ll_i)/d(beta_j) = delta_i * x_ij - x_ij * exp(xb_i) * (exp(gamma*t_i)-1)/gamma
    #                    = x_ij * (delta_i - H_i)

    if sample_weight is not None:
        dll_dgamma = dll_dgamma * sample_weight
        w = sample_weight
        dll_dbeta = X.T @ ((event - H) * w)

    grad = np.zeros_like(params)
    grad[0] = -np.sum(dll_dgamma)
    grad[1:] = -dll_dbeta

    return grad


def fit_gompertz(
    df: pd.DataFrame,
    features: list[str] | None = None,
    time_col: str = "exposure_10yr",
    event_col: str = "mort_10yr",
    weight_col: str | None = "exam_weight",
    include_intercept: bool = True,
    initial_gamma: float = 0.008,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a Gompertz proportional hazards model.

    Args:
        df: DataFrame with biomarker features, time, and event columns
        features: List of feature column names. If None, uses PHENOAGE_FEATURES.
        time_col: Column name for follow-up time (months)
        event_col: Column name for event indicator (1=death, 0=censored)
        weight_col: Column name for survey weights (or None)
        include_intercept: Whether to include an intercept term
        initial_gamma: Starting value for gamma parameter
        verbose: Print progress

    Returns:
        Dictionary with:
        - 'gamma': Estimated Gompertz shape parameter
        - 'beta': Estimated coefficients (including intercept if requested)
        - 'feature_names': Names corresponding to beta coefficients
        - 'coefficients': Dict mapping feature names to coefficients
        - 'n_obs': Number of observations
        - 'n_events': Number of events
        - 'converged': Whether optimization converged
        - 'result': Full scipy OptimizeResult
    """
    if features is None:
        features = PHENOAGE_FEATURES

    # Prepare data: drop rows with any missing values in required columns
    required_cols = features + [time_col, event_col]
    if weight_col and weight_col in df.columns:
        required_cols.append(weight_col)
    else:
        weight_col = None

    df_clean = df.dropna(subset=required_cols).copy()
    if verbose:
        print(f"  Fitting Gompertz PH model:")
        print(f"  {df.shape[0]:,d} rows -> {df_clean.shape[0]:,d} complete cases")
        print(f"  Features: {features}")
        print(f"  Events: {(df_clean[event_col] == 1).sum():,d} / {df_clean.shape[0]:,d}")

    # Build design matrix
    X = df_clean[features].values.astype(np.float64)
    if include_intercept:
        X = np.column_stack([np.ones(X.shape[0]), X])
        feature_names = ["intercept"] + list(features)
    else:
        feature_names = list(features)

    time = df_clean[time_col].values.astype(np.float64)
    event = df_clean[event_col].values.astype(np.float64)

    weights = None
    if weight_col:
        weights = df_clean[weight_col].values.astype(np.float64)
        # Normalize weights to sum to n
        weights = weights / weights.mean()

    # Initial parameters: use Levine's published values as starting point
    n_params = 1 + X.shape[1]  # gamma + betas
    x0 = np.zeros(n_params)
    x0[0] = initial_gamma  # gamma

    # Initialize betas with Levine's published values
    for i, name in enumerate(feature_names):
        if name in LEVINE_COEFFICIENTS:
            x0[i + 1] = LEVINE_COEFFICIENTS[name]

    # Set bounds: gamma > 0
    bounds = [(1e-6, 0.5)] + [(None, None)] * X.shape[1]

    if verbose:
        print(f"  Optimizing with {n_params} parameters...")
        print(f"  Initial gamma: {x0[0]:.6f}")
        ll_init = gompertz_log_likelihood(x0, X, time, event, weights)
        print(f"  Initial -LL: {ll_init:.2f}")

    # Use multiple optimization attempts
    best_result = None
    best_ll = np.inf

    for method in ["L-BFGS-B"]:
        try:
            result = minimize(
                gompertz_log_likelihood,
                x0,
                args=(X, time, event, weights),
                jac=gompertz_log_likelihood_gradient,
                method=method,
                bounds=bounds,
                options={"maxiter": 10000, "ftol": 1e-14, "gtol": 1e-10},
            )
            if result.fun < best_ll:
                best_ll = result.fun
                best_result = result
                if verbose:
                    print(f"  {method}: -LL={result.fun:.2f}, converged={result.success}")
        except Exception as e:
            if verbose:
                print(f"  {method}: failed ({e})")

    # Try again from the best result with tighter tolerances
    if best_result is not None:
        try:
            result2 = minimize(
                gompertz_log_likelihood,
                best_result.x,
                args=(X, time, event, weights),
                jac=gompertz_log_likelihood_gradient,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 10000, "ftol": 1e-15, "gtol": 1e-12},
            )
            if result2.fun < best_ll:
                best_result = result2
                if verbose:
                    print(f"  Refinement: -LL={result2.fun:.2f}, converged={result2.success}")
        except Exception:
            pass

    result = best_result

    gamma_hat = result.x[0]
    beta_hat = result.x[1:]

    coefficients = dict(zip(feature_names, beta_hat))

    if verbose:
        print(f"\n  {'='*50}")
        print(f"  Gompertz PH Model Results")
        print(f"  {'='*50}")
        print(f"  Converged: {result.success}")
        print(f"  Gamma (shape): {gamma_hat:.7f}")
        print(f"  {'Feature':<25s} {'Coefficient':>12s}")
        print(f"  {'-'*25} {'-'*12}")
        for name, coef in zip(feature_names, beta_hat):
            levine_coef = LEVINE_COEFFICIENTS.get(name, None)
            comparison = f" (Levine: {levine_coef:.4f})" if levine_coef is not None else ""
            print(f"  {name:<25s} {coef:>12.6f}{comparison}")
        print(f"\n  Levine gamma: {GOMPERTZ_GAMMA:.7f}, Fitted gamma: {gamma_hat:.7f}")

    return {
        "gamma": gamma_hat,
        "beta": beta_hat,
        "feature_names": feature_names,
        "coefficients": coefficients,
        "n_obs": df_clean.shape[0],
        "n_events": int((event == 1).sum()),
        "converged": result.success,
        "result": result,
    }


def compute_mortality_score(
    xb: np.ndarray | pd.Series,
    gamma: float = GOMPERTZ_GAMMA,
    t_months: int = MORTALITY_WINDOW_MONTHS,
) -> np.ndarray:
    """Compute 10-year mortality probability from linear predictor.

    M = 1 - exp(-exp(xb) * (exp(gamma * t) - 1) / gamma)

    Args:
        xb: Linear predictor values (X @ beta)
        gamma: Gompertz shape parameter
        t_months: Follow-up time in months (120 = 10 years)

    Returns:
        Mortality probability M ∈ [0, 1]
    """
    xb = np.asarray(xb, dtype=np.float64)
    cumhaz = np.exp(xb) * (np.exp(gamma * t_months) - 1.0) / gamma
    M = 1.0 - np.exp(-cumhaz)
    return M


def compute_phenoage(
    M: np.ndarray | pd.Series,
    intercept: float = PHENOAGE_INTERCEPT,
    log_scale: float = PHENOAGE_LOG_SCALE,
    rate: float = PHENOAGE_RATE,
) -> np.ndarray:
    """Convert mortality probability to PhenoAge.

    PhenoAge = intercept + ln(-log_scale * ln(1 - M)) / rate

    Args:
        M: Mortality probability from compute_mortality_score
        intercept, log_scale, rate: PhenoAge conversion constants

    Returns:
        PhenoAge values (biological age estimate)
    """
    M = np.asarray(M, dtype=np.float64)
    # Clip M to avoid log(0) issues
    M = np.clip(M, 1e-10, 1.0 - 1e-10)
    phenoage = intercept + np.log(-log_scale * np.log(1.0 - M)) / rate
    return phenoage


def score_phenoage(
    df: pd.DataFrame,
    coefficients: dict[str, float] | None = None,
    gamma: float | None = None,
    features: list[str] | None = None,
) -> pd.DataFrame:
    """Score a DataFrame with PhenoAge using given or published coefficients.

    Args:
        df: DataFrame with PhenoAge biomarker columns
        coefficients: Dict of feature_name -> coefficient.
                     If None, uses published Levine 2018 coefficients.
        gamma: Gompertz shape parameter. If None, uses Levine's value.
        features: Feature names to use. If None, uses PHENOAGE_FEATURES.

    Returns:
        DataFrame with added columns:
        - 'xb': linear predictor
        - 'mortality_score': 10-year mortality probability
        - 'phenoage': PhenoAge (biological age)
        - 'phenoage_accel': PhenoAge acceleration (PhenoAge - chronological age)
    """
    if coefficients is None:
        coefficients = LEVINE_COEFFICIENTS
    if gamma is None:
        gamma = GOMPERTZ_GAMMA
    if features is None:
        features = PHENOAGE_FEATURES

    df = df.copy()

    # Compute linear predictor
    intercept = coefficients.get("intercept", 0.0)
    xb = np.full(df.shape[0], intercept, dtype=np.float64)

    for feature in features:
        if feature not in df.columns:
            print(f"  ⚠ Missing feature '{feature}'; setting xb to NaN for missing rows")
            xb = np.where(df.get(feature, pd.Series(np.nan)).isna(), np.nan, xb)
            continue
        coef = coefficients.get(feature, 0.0)
        xb = np.where(
            df[feature].isna(),
            np.nan,
            xb + coef * df[feature].values,
        )

    df["xb"] = xb
    df["mortality_score"] = compute_mortality_score(xb, gamma=gamma)
    df["phenoage"] = compute_phenoage(df["mortality_score"])

    if "age" in df.columns:
        df["phenoage_accel"] = df["phenoage"] - df["age"]

    n_scored = (~df["phenoage"].isna()).sum()
    print(
        f"  PhenoAge scored: {n_scored:,d}/{df.shape[0]:,d} " f"({n_scored/df.shape[0]*100:.1f}%)"
    )

    if "age" in df.columns and n_scored > 0:
        valid = df.dropna(subset=["phenoage", "age"])
        corr = valid["phenoage"].corr(valid["age"])
        mean_accel = valid["phenoage_accel"].mean()
        print(f"  PhenoAge-Age correlation: {corr:.3f}")
        print(f"  Mean PhenoAge acceleration: {mean_accel:.2f} years")

    return df
