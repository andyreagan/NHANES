from pathlib import Path

import altair as alt
import pandas as pd

from .constants import ALL_NHANES_IV_CYCLES, keep
from .main import main, main_nhanes_iii

# alt.data_transformers.enable('json')
alt.data_transformers.enable("default", max_rows=1000000)


if __name__ == "__main__":
    base_output_path: Path = Path("data", "processed")

    # =========================================================================
    # NHANES III (1988-1994) — fixed-width .dat files
    # =========================================================================

    print(f"\n{'='*80}")
    print("Processing NHANES III (1988-1994)")
    print(f"{'='*80}")
    (base_output_path / "NHANES_III").mkdir(parents=True, exist_ok=True)
    try:
        nhanes_iii_result = main_nhanes_iii(
            output_path=Path("output") / "NHANES_III",
            plot_dists=False,
            plot_missingness_diagnostics=False,
            save_unfiltered_fname=str(base_output_path / "NHANES_III" / "nhanes_all.parquet"),
        )
        nhanes_iii_result.to_parquet(base_output_path / "NHANES_III" / "nhanes.parquet")
    except Exception as e:
        nhanes_iii_result = None
        print(f"⚠ Error processing NHANES III: {e}")
        import traceback

        traceback.print_exc()

    # =========================================================================
    # All NHANES IV cycles (1999-2000 through 2021-2023)
    # =========================================================================

    all_cycles = ALL_NHANES_IV_CYCLES
    nhanes_iv_results: list[pd.DataFrame] = []

    # Process each cycle individually
    year: str
    suffix: str
    start_year: int
    d: pd.DataFrame
    for year, suffix, start_year in all_cycles:
        print(f"\n{'='*80}")
        print(f"Processing cycle: {year} (suffix={suffix!r})")
        print(f"{'='*80}")
        (base_output_path / year).mkdir(parents=True, exist_ok=True)
        try:
            cycle_result = main(
                base_paths=[Path(f"data/raw/{year}")],
                suffixes=[suffix],
                output_path=Path("output") / year,
                plot_dists=False,
                plot_missingness_diagnostics=False,
                save_unfiltered_fname=str(base_output_path / year / "nhanes_all.parquet"),
                cycles=[year],
            )
            cycle_result.to_parquet(base_output_path / year / "nhanes.parquet")
            nhanes_iv_results.append(cycle_result)
        except Exception as e:
            print(f"⚠ Error processing {year}: {e}")
            import traceback

            traceback.print_exc()

    # =========================================================================
    # Combined NHANES IV output (all continuous cycles together)
    # =========================================================================

    print(f"\n{'='*80}")
    print("Processing combined NHANES IV dataset (all continuous cycles)")
    print(f"{'='*80}")
    base_output_path.mkdir(parents=True, exist_ok=True)
    d = main(
        base_paths=[Path(f"data/raw/{year}") for year, _, _ in all_cycles],
        suffixes=[suffix for _, suffix, _ in all_cycles],
        output_path=Path("output") / "combined",
        plot_dists=False,
        plot_missingness_diagnostics=False,
        save_unfiltered_fname=str(base_output_path / "nhanes_all.parquet"),
        cycles=[year for year, _, _ in all_cycles],
    )
    d.to_parquet(base_output_path / "nhanes.parquet")
    d.loc[:, ["SEQN"] + list(keep)].to_parquet(base_output_path / "nhanes_m3s_vars.parquet")

    # =========================================================================
    # Grand combined output (NHANES III + NHANES IV)
    # =========================================================================

    print(f"\n{'='*80}")
    print("Creating grand combined dataset (NHANES III + NHANES IV)")
    print(f"{'='*80}")
    all_dfs = []
    if nhanes_iii_result is not None and nhanes_iii_result.shape[0] > 0:
        all_dfs.append(nhanes_iii_result)
    if d.shape[0] > 0:
        all_dfs.append(d)
    if all_dfs:
        # Use only columns present in all DataFrames (intersection)
        common_cols = set(all_dfs[0].columns)
        for df_part in all_dfs[1:]:
            common_cols &= set(df_part.columns)
        grand = pd.concat([df_part[sorted(common_cols)] for df_part in all_dfs], ignore_index=True)
        grand.to_parquet(base_output_path / "nhanes_all_cohorts.parquet")
        print(f"Grand combined: {grand.shape[0]:,d} rows, {grand.shape[1]} columns")
        print(f"  Deaths: {grand.is_dead.sum():.0f}")
    else:
        print("⚠ No data to combine")
