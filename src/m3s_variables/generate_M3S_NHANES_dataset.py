from pathlib import Path

import altair as alt
import pandas as pd

from .constants import keep
from .main import main

# alt.data_transformers.enable('json')
alt.data_transformers.enable("default", max_rows=1000000)


if __name__ == "__main__":
    base_output_path: Path = Path("data", "processed")
    years: set = {("2003-2004", "C"), ("2005-2006", "D")}
    year: str
    suffix: str
    d: pd.DataFrame
    for year, suffix in years:
        (base_output_path / year).mkdir(parents=True, exist_ok=True)
        main(
            base_paths=[Path(f"data/raw/{year}")],
            suffixes=[suffix],
            output_path=Path("output") / year,
            # plot_dists=True,
            # plot_missingness_diagnostics=True,
            plot_dists=False,
            plot_missingness_diagnostics=False,
            save_unfiltered_fname=str(base_output_path / year / "nhanes_all.parquet"),
        ).to_parquet(base_output_path / year / "nhanes.parquet")
    base_output_path.mkdir(parents=True, exist_ok=True)
    d = main(
        base_paths=[Path(f"data/raw/{year}") for year, _ in years],
        suffixes=[suffix for _, suffix in years],
        output_path=Path("output") / "combined",
        # plot_dists=True,
        # plot_missingness_diagnostics=True,
        plot_dists=False,
        plot_missingness_diagnostics=False,
        save_unfiltered_fname=str(base_output_path / "nhanes_all.parquet"),
    )
    d.to_parquet(base_output_path / "nhanes.parquet")
    d.loc[:, ["SEQN"] + list(keep)].to_parquet(base_output_path / "nhanes_m3s_vars.parquet")
