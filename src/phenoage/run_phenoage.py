"""Main script to reproduce Levine 2018 PhenoAge analysis.

Steps:
1. Load NHANES III lab data + mortality → training set
2. Load NHANES IV (1999-2010) lab data + mortality → validation set
3. Fit Gompertz PH model on NHANES III → estimated coefficients
4. Compare fitted coefficients with Levine's published coefficients
5. Score PhenoAge on NHANES IV using both fitted and published coefficients
6. Save results

Usage:
    uv run -m src.phenoage.run_phenoage
"""

from pathlib import Path

import numpy as np
import pandas as pd

from .constants import (
    LEVINE_COEFFICIENTS,
    LEVINE_VALIDATION_CYCLES,
    NHANES_IV_CYCLES,
    PHENOAGE_FEATURES,
)
from .load_data import load_complete_dataset
from .model import fit_gompertz, score_phenoage

OUTPUT_DIR = Path("data/processed/phenoage")
PLOT_DIR = Path("output/phenoage")


def main(
    fit_model: bool = True,
    score_all_cycles: bool = True,
) -> None:
    """Run the complete PhenoAge reproduction pipeline."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Step 1: Load data
    # =========================================================================
    print("\n" + "#" * 70)
    print("# Step 1: Load data")
    print("#" * 70)

    # For model fitting, use NHANES III (training) and NHANES IV 1999-2010 (validation)
    # For scoring, optionally include all NHANES IV cycles
    iv_cycles = NHANES_IV_CYCLES if score_all_cycles else LEVINE_VALIDATION_CYCLES

    nhanes_iii, nhanes_iv, mortality = load_complete_dataset(
        include_nhanes_iii=fit_model,
        nhanes_iv_cycles=iv_cycles,
    )

    # =========================================================================
    # Step 2: Summarize the data
    # =========================================================================
    print("\n" + "#" * 70)
    print("# Step 2: Data summary")
    print("#" * 70)

    if nhanes_iii is not None:
        print(f"\nNHANES III (training):")
        print(f"  Rows: {nhanes_iii.shape[0]:,d}")
        for feat in PHENOAGE_FEATURES:
            if feat in nhanes_iii.columns:
                n = (~nhanes_iii[feat].isna()).sum()
                print(f"  {feat}: {n:,d} ({n/nhanes_iii.shape[0]*100:.1f}%)")

        # Complete cases for model fitting
        complete = nhanes_iii.dropna(subset=PHENOAGE_FEATURES + ["exposure_10yr", "mort_10yr"])
        print(f"  Complete cases: {complete.shape[0]:,d}")
        print(f"  Deaths (10yr): {(complete.mort_10yr == 1).sum():,d}")
        nhanes_iii.to_parquet(OUTPUT_DIR / "nhanes_iii_phenoage_data.parquet")

    print(f"\nNHANES IV (validation/scoring):")
    print(f"  Rows: {nhanes_iv.shape[0]:,d}")
    for cycle in nhanes_iv.cycle.unique():
        cycle_df = nhanes_iv[nhanes_iv.cycle == cycle]
        n_complete = cycle_df.dropna(subset=PHENOAGE_FEATURES).shape[0]
        print(f"  {cycle}: {cycle_df.shape[0]:,d} rows, " f"{n_complete:,d} complete for PhenoAge")

    nhanes_iv.to_parquet(OUTPUT_DIR / "nhanes_iv_phenoage_data.parquet")

    # =========================================================================
    # Step 3: Fit Gompertz model on NHANES III
    # =========================================================================
    if fit_model and nhanes_iii is not None:
        print("\n" + "#" * 70)
        print("# Step 3: Fit Gompertz PH model on NHANES III")
        print("#" * 70)

        model_result = fit_gompertz(
            nhanes_iii,
            features=PHENOAGE_FEATURES,
            time_col="exposure_10yr",
            event_col="mort_10yr",
            weight_col="exam_weight",
            verbose=True,
        )

        # Save model results
        coef_df = pd.DataFrame(
            {
                "feature": model_result["feature_names"],
                "fitted_coef": model_result["beta"],
                "levine_coef": [
                    LEVINE_COEFFICIENTS.get(f, np.nan) for f in model_result["feature_names"]
                ],
            }
        )
        coef_df["difference"] = coef_df["fitted_coef"] - coef_df["levine_coef"]
        coef_df["pct_difference"] = coef_df["difference"] / coef_df["levine_coef"].abs() * 100

        print("\n  Coefficient comparison:")
        print(coef_df.to_string(index=False))
        coef_df.to_csv(OUTPUT_DIR / "gompertz_coefficients.csv", index=False)

        model_summary = {
            "gamma_fitted": model_result["gamma"],
            "gamma_levine": 0.0076927,
            "n_obs": model_result["n_obs"],
            "n_events": model_result["n_events"],
            "converged": model_result["converged"],
        }
        pd.DataFrame([model_summary]).to_csv(OUTPUT_DIR / "gompertz_model_summary.csv", index=False)

        # =====================================================================
        # Step 4: Score NHANES IV with fitted coefficients
        # =====================================================================
        print("\n" + "#" * 70)
        print("# Step 4: Score NHANES IV with fitted coefficients")
        print("#" * 70)

        nhanes_iv_fitted = score_phenoage(
            nhanes_iv,
            coefficients=model_result["coefficients"],
            gamma=model_result["gamma"],
        )

        # Rename columns to indicate fitted model
        nhanes_iv_fitted = nhanes_iv_fitted.rename(
            columns={
                "xb": "xb_fitted",
                "mortality_score": "mortality_score_fitted",
                "phenoage": "phenoage_fitted",
                "phenoage_accel": "phenoage_accel_fitted",
            }
        )

    # =========================================================================
    # Step 5: Score NHANES IV with published Levine coefficients
    # =========================================================================
    print("\n" + "#" * 70)
    print("# Step 5: Score NHANES IV with published Levine 2018 coefficients")
    print("#" * 70)

    nhanes_iv_levine = score_phenoage(
        nhanes_iv,
        coefficients=LEVINE_COEFFICIENTS,
    )

    # =========================================================================
    # Step 6: Combine and save
    # =========================================================================
    print("\n" + "#" * 70)
    print("# Step 6: Save results")
    print("#" * 70)

    # Combine fitted and published scores
    if fit_model and nhanes_iii is not None:
        for col in [
            "xb_fitted",
            "mortality_score_fitted",
            "phenoage_fitted",
            "phenoage_accel_fitted",
        ]:
            if col in nhanes_iv_fitted.columns:
                nhanes_iv_levine[col] = nhanes_iv_fitted[col]

    # Save final scored dataset
    nhanes_iv_levine.to_parquet(OUTPUT_DIR / "nhanes_iv_phenoage_scored.parquet")

    # Print summary statistics by cycle
    print("\n  Summary by cycle (Levine published coefficients):")
    print(
        f"  {'Cycle':<12s} {'N':>7s} {'N scored':>9s} "
        f"{'Mean Age':>9s} {'Mean PA':>9s} {'Mean Accel':>11s} {'Corr':>6s}"
    )
    print(f"  {'-'*12} {'-'*7} {'-'*9} {'-'*9} {'-'*9} {'-'*11} {'-'*6}")

    for cycle in sorted(nhanes_iv_levine.cycle.unique()):
        cdf = nhanes_iv_levine[nhanes_iv_levine.cycle == cycle]
        valid = cdf.dropna(subset=["phenoage", "age"])
        if valid.shape[0] > 0:
            corr = valid.phenoage.corr(valid.age)
            print(
                f"  {cycle:<12s} {cdf.shape[0]:>7,d} {valid.shape[0]:>9,d} "
                f"{valid.age.mean():>9.1f} {valid.phenoage.mean():>9.1f} "
                f"{valid.phenoage_accel.mean():>11.2f} {corr:>6.3f}"
            )
        else:
            print(
                f"  {cycle:<12s} {cdf.shape[0]:>7,d} {'0':>9s} "
                f"{'N/A':>9s} {'N/A':>9s} {'N/A':>11s} {'N/A':>6s}"
            )

    print(f"\n  Results saved to: {OUTPUT_DIR}")
    print("  Done!")


if __name__ == "__main__":
    main()
