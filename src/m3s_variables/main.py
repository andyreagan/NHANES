from collections.abc import Callable
from pathlib import Path
from subprocess import run
from typing import Optional, Union

import altair as alt
import numpy as np
import pandas as pd
from altair_saver import save

from .constants import (
    FILE_MAPPING_PRE2005,
    LOOKUP,
    NHANES_VAR_FILES,
    VAR_MAPPING_PRE2005,
    keep,
    predictors,
    required,
)
from .util import bounds_filter, mylambda

WIDTH = 300
HEIGHT = 300


def make_file(filename: Path) -> bool:
    if not filename.exists():
        filename.parent.mkdir(parents=True, exist_ok=True)
        run(["make", filename])
    # return filename.exists()
    # instead, we'll check that the file has a certian size
    if filename.exists():
        if filename.stat().st_size < 100:
            return False
        else:
            return True
    else:
        return False


def load_NHANES_files(
    base_df: pd.DataFrame,
    base_paths: list[Path],
    suffixes: list[str],
    file_lookup: dict,
    file_mapping: dict,
    var_mapping: dict,
) -> pd.DataFrame:
    df_copy = base_df.copy()

    file: str
    vars: list
    for file, vars in file_lookup.items():
        # already loaded this one as the base:
        if file in {"DEMO", "mortality", "other"}:
            continue
        var_files: list = []
        for base_path, suffix in zip(base_paths, suffixes):
            filename: Path = base_path / f"{file}_{suffix}.XPT"
            filename_pre2005: Path
            print(f"Loading {filename}")
            df: pd.DataFrame
            if make_file(filename):
                df = pd.read_sas(filename)
            elif (file in file_mapping) and make_file(
                base_path / f"{file_mapping[file]}_{suffix}.XPT"
            ):
                filename_pre2005 = base_path / f"{file_mapping[file]}_{suffix}.XPT"
                print(
                    f" • Was missing {filename}, found and loading pre-2005 file {filename_pre2005} instead."
                )
                df = pd.read_sas(filename_pre2005)
            elif file in file_mapping:
                filename_pre2005 = base_path / f"{file_mapping[file]}_{suffix}.XPT"
                print(f" • ⚠ Missing both {filename} and {filename_pre2005} ⚠ ")
                # create an empty df:
                df = pd.DataFrame([{"SEQN": 0} | {x: 0 for x in vars}]).iloc[:0, :]
                assert df.shape[0] == 0
            else:
                print(f" • ⚠ Missing {filename}, no alternate ⚠ ")
                # create an empty df:
                df = pd.DataFrame([{"SEQN": 0} | {x: 0 for x in vars}]).iloc[:0, :]
                assert df.shape[0] == 0

            print(f" • Row count from this file is: {df.shape[0]:,d}")
            # make sure it's unique
            assert df.shape[0] == df.SEQN.unique().shape[0]
            # find missing cols from other files
            for col in set(vars) - set(df.columns):
                print(f" • Looking for missing {col=}")
                # see if file has mapped vars
                if file in var_mapping:
                    # see if the col is in those mapped vars
                    if col in var_mapping[file]:
                        print(f" • Found a mapped version for {col=} as {var_mapping[file][col]=}")
                        # see if we have a mapping file
                        # if file in file_mapping:
                        # we should have already loaded that version
                        assert var_mapping[file][col] in df.columns
                        df[col] = df[var_mapping[file][col]]
            missing_cols: set = set(vars) - set(df.columns)
            has_cols: set = set(df.columns) & set(vars)
            if len(missing_cols) > 0:
                print(f" • ⚠ df {has_cols=}, but {missing_cols=}. filling missing with None ⚠")
                for col in missing_cols:
                    df[col] = None
            print(df.columns)
            var_files.append(df)
        dfs: pd.DataFrame = pd.concat([df.loc[:, ["SEQN"] + vars] for df in var_files])
        print(f" • Joining {filename}")
        df_copy = df_copy.merge(dfs, on=["SEQN"], how="left")
        print(f" • Row count is now: {df_copy.shape[0]:,d}")

    return df_copy


def plot_all_dists(df: pd.DataFrame, info: dict, output_path: Path) -> pd.DataFrame:
    nice_names = {"albumin": "Albumin"}
    for var, info in info.items():
        print(" • Creating plot.")
        chart = (
            alt.Chart(df.loc[:, [var]])
            .mark_bar()
            .encode(
                alt.X(var, bin=alt.Bin(maxbins=20), title=nice_names.get(var, var)),
                alt.Y("count()", title="Count of Records"),
            )
        ).properties(width=WIDTH, height=HEIGHT)
        print(f"    Saving plot at {str(output_path / f'{var}.png')}.")
        save(chart, str(output_path / f"{var}.png"), vega_cli_options=["-s 4"])


def test_process_variables():
    from .util import ifelse

    d = pd.DataFrame({"a": [0, 1, np.NaN], "b": [1, 2, 3]})
    all_info = {
        "a_new": {
            "depends_on": ["a", "b"],
            "expr": "ifelse(any(c(_[0] == 0, _[1] == 1)), 'a', 'b')",
        }
    }
    d.loc[:, d.columns].dropna(how="any").apply(lambda x: any([x[0] == 0, x[1] == 1]), axis=1)
    d.loc[:, d.columns].dropna(how="any").apply(
        lambda x: ifelse(any([x[0] == 0, x[1] == 1]), "a", "b"), axis=1
    )
    d.loc[:, d.columns].dropna(how="any").apply(mylambda(all_info["a_new"]["expr"]), axis=1)
    process_variables(d, all_info)


def process_variables(df: pd.DataFrame, all_info: dict) -> pd.DataFrame:
    df_copy = df.copy()

    # TODO: make sure we go in an _ordered_ way here,
    # so we can use columns that we create in subsequent statements
    var: str
    info: Optional[dict]
    for var, info in all_info.items():
        if info is None or "depends_on" not in info or len(info["depends_on"]) == 0:
            print(f"Skipping {var=}, set to NA.")
            df_copy[var] = None
            continue
        print(f"Creating {var=}.")
        missing_cols = set(info["depends_on"]) - set(df_copy.columns)
        if missing_cols:
            print(missing_cols)
        assert len(missing_cols) == 0
        if "expr" in info:
            print(f' • Applying function(_) {{ {info["expr"]} }} to input variable list _.')
            if len(info["depends_on"]) > 1:
                df_var_depends_on_notna = df_copy.loc[:, info["depends_on"]].dropna(
                    how=info.get("na_if_dependency_na", "any"),
                    subset=info.get("na_if_dependency_na_subset", info["depends_on"]),
                )
            else:
                df_var_depends_on_notna = df_copy.loc[
                    :, info["depends_on"][0]
                ].dropna()  # how= is not used, does not accept subset (that wouldn't even make sense!)
            expr = mylambda(info["expr"])
            if df_var_depends_on_notna.shape[0] > 0:
                df_copy[var] = df_var_depends_on_notna.apply(expr, axis=1)
            else:
                print(f" • Attention: no rows with notna depends for {var}")
                df_copy[var] = None  # np.NaN or pd.NA work too
        elif (n := len(info["depends_on"])) > 1:
            print(
                f' • Given no expression, but {n} vars {", ".join(info["depends_on"])}...coalesce them.'
            )
            df_copy[var] = df_copy.loc[:, info["depends_on"]].bfill(axis=1).iloc[:, 0]
        else:
            print(
                f' • Given no expression and only 1 var {info["depends_on"][0]} then we just rename.'
            )
            df_copy[var] = df_copy.loc[:, info["depends_on"][0]]
        print(
            f""" • NA count of {var} is {df_copy[var].isna().sum():,d}, {df_copy[var].isna().sum()/df_copy.shape[0]*100:.1f}%."""
        )
        # stash a raw version
        df_copy[f"{var}_raw"] = df_copy[var]
        if "valid_range" in info:
            print(" • Apply filter")
            var_bounds_filter: Callable = bounds_filter(info["valid_range"])
            # create a valid version
            df_copy[f"{var}_valid"] = df_copy.loc[df_copy[var].apply(var_bounds_filter), var]
            print(
                f""" • Valid count of {var} is {(~df_copy[f'{var}_valid'].isna()).sum():,d}, {(~df_copy[f'{var}_valid'].isna()).sum()/(~df_copy[f'{var}_raw'].isna()).sum()*100:.1f}% of raw values. {(~df_copy[f'{var}_raw'].isna()).sum()-(~df_copy[f'{var}_valid'].isna()).sum():,d} invalids set to NA."""
            )
        else:
            df_copy[f"{var}_valid"] = df_copy[var]
        if (caps := info.get("caps", None)) is not None:
            print(" • Apply caps")
            if caps[0] is not None:
                print(
                    f""" • Capping {(df_copy[f"{var}_valid"] < caps[0]).sum():,d} of {(~df_copy[f'{var}_valid'].isna()).sum():,d} values, or {(df_copy[f"{var}_valid"] < caps[0]).sum()/(~df_copy[f'{var}_valid'].isna()).sum()*100:.1f}%, at cap of {caps[0]}."""
                )
                df_copy.loc[df_copy[f"{var}_valid"] < caps[0], f"{var}_valid"] = caps[0]
            if caps[1] is not None:
                print(
                    f""" • Capping {(df_copy[f"{var}_valid"] > caps[1]).sum():,d} of {(~df_copy[f'{var}_valid'].isna()).sum():,d} values, or {(df_copy[f"{var}_valid"] > caps[1]).sum()/(~df_copy[f'{var}_valid'].isna()).sum()*100:.1f}%, at cap of {caps[1]}."""
                )
                df_copy.loc[df_copy[f"{var}_valid"] > caps[1], f"{var}_valid"] = caps[1]
        # overwrite the base with the valid version
        df_copy[var] = df_copy[f"{var}_valid"]
    return df_copy


def impute_defaults(d: pd.DataFrame, info: dict) -> pd.DataFrame:
    """Impute defaults wherever they are given."""
    d_copy = d.copy()
    for var, var_info in [(k, v) for k, v in info.items() if v is not None]:
        # don't worry about checking post-impute types:
        # they shouldn't have defaults!
        if var_info.get("impute", False) or var_info.get("predictor", False):
            if (val := var_info.get("impute_value", None)) is not None:
                print(f" • Imputing default for {var} to {val}.")
                d_copy[var] = d_copy[var].fillna(val)
            # could, in theory, do the medians here too:
            elif (medians := var_info.get("impute_medians", None)) is not None:
                print(f" • Imputing medians for {var} to {medians}.")
                d_copy[var] = (
                    d_copy.loc[d_copy[var].isna(), ["age_5", "sex"]]
                    .merge(pd.DataFrame(medians), how="left", on=["age_5", "sex"])
                    .loc[:, var]
                )
    return d_copy


def missingness_correlation(nhanes_m3s: pd.DataFrame, output_path: Path, suffix: str = "") -> None:
    chart: alt.Chart

    missingness_crosstab: pd.DataFrame = pd.DataFrame(
        {
            column: nhanes_m3s.loc[~nhanes_m3s.iloc[:, i].isna(), :].isna().sum()
            for i, column in enumerate(nhanes_m3s.columns)
        }
    )

    chart = (
        alt.Chart(pd.melt(missingness_crosstab, ignore_index=False).reset_index())
        .mark_rect()
        .encode(alt.X("index:N", title=""), alt.Y("variable:N", title=""))
        .properties(width=WIDTH, height=HEIGHT)
    )
    print(f" • Saving plot at {str(output_path)}/missingness{suffix}.png.")
    save(
        chart.encode(color="value:Q") + chart.mark_text(fontSize=8).encode(text="value"),
        str(output_path / f"missingness{suffix}.png"),
        vega_cli_options=["-s 4"],
    )


def missingness_diagnostics(d: pd.DataFrame, output_path: Path) -> None:
    chart: alt.Chart
    m: pd.DataFrame
    m2: pd.DataFrame
    threshold_counts: pd.DataFrame

    m = pd.DataFrame({"missingness": d.isna().sum(axis=1)}).sort_values("missingness")
    m["total_missingness"] = np.cumsum(m.missingness)
    m["average_missingness"] = m.total_missingness / np.cumsum(np.ones(m.shape[0]))
    m["index"] = np.arange(m.shape[0])

    chart = (
        alt.Chart(m.loc[:, ["missingness"]])
        .mark_bar()
        .encode(
            alt.X("missingness", bin=alt.Bin(maxbins=20), title="Missing count (binned)"),
            alt.Y("count()"),
        )
        .properties(width=WIDTH, height=HEIGHT)
    )
    print(f" •  Saving plot at {str(output_path / 'missingness_dist.png')}.")
    save(chart, str(output_path / "missingness_dist.png"), vega_cli_options=["-s 4"])
    chart = (
        alt.Chart(m.loc[:, ["missingness"]])
        .mark_bar()
        .encode(x="missingness", y="count()")
        .properties(width=WIDTH, height=HEIGHT)
    )
    print(f" • Saving plot at {str(output_path / 'missingness_dist_2.png')}.")
    save(chart, str(output_path / "missingness_dist_2.png"), vega_cli_options=["-s 4"])
    chart = (
        alt.Chart(m.loc[:, ["index", "average_missingness"]])
        .mark_line()
        .encode(x="index", y="average_missingness")
        .properties(width=WIDTH, height=HEIGHT)
    )
    print(f" • Saving plot at {str(output_path / 'missingness_cum.png')}.")
    save(chart, str(output_path / "missingness_cum.png"), vega_cli_options=["-s 4"])

    m2 = (d.isna() * 1).loc[m.index, :]
    m2["index"] = np.arange(m2.shape[0])
    melted: pd.DataFrame = m2.melt(id_vars="index")
    chart = (
        alt.Chart(melted)
        .mark_rect()
        .encode(x=alt.X("index", bin=alt.Bin(maxbins=100)), y="variable", color="sum(value)")
        .properties(width=WIDTH, height=HEIGHT)
    )
    print(f" • Saving plot at {str(output_path / 'missingness_all.png')}.")
    save(chart, str(output_path / "missingness_all.png"), vega_cli_options=["-s 4"])

    threshold_counts = pd.DataFrame(
        {
            "threshold": range(1, d.shape[1]),
            "rows_remaining": [d.dropna(thresh=i).shape[0] for i in range(1, d.shape[1])],
        }
    )
    chart = (
        alt.Chart(threshold_counts)
        .mark_line()
        .encode(
            alt.X("threshold", title="Required covariate(s) (threshold)"),
            alt.Y("rows_remaining", title="Rows remaining"),
        )
        .properties(width=WIDTH, height=HEIGHT)
    )
    print(f" • Saving plot at {str(output_path / 'missingness_thresh.png')}.")
    save(chart, str(output_path / "missingness_thresh.png"), vega_cli_options=["-s 4"])


def main(
    base_paths: list[Path],
    suffixes: list[str],
    output_path: Path,
    plot_dists: bool = False,
    plot_missingness_diagnostics: bool = False,
    save_unfiltered_fname: Union[None, str] = None,
) -> pd.DataFrame:
    # initialize our wide DF with the demographics:
    output_path.mkdir(parents=True, exist_ok=True)
    bases: list[pd.DataFrame] = [
        pd.read_sas(base_path / f"DEMO_{suffix}.XPT")
        for base_path, suffix in zip(base_paths, suffixes)
        if make_file(base_path / f"DEMO_{suffix}.XPT")
    ]
    # for i, base_path in enumerate(base_paths):
    #     bases[i]["year"] = int(base_path.name.split("-")[0])
    all_nhanes: pd.DataFrame = pd.concat(
        # [base.loc[:, ["SEQN", "year"] + NHANES_VAR_FILES["DEMO"]] for base in bases]
        [base.loc[:, ["SEQN"] + NHANES_VAR_FILES["DEMO"]] for base in bases]
    )
    print(f"Row count from demographic file is {all_nhanes.shape[0]:,d}")
    mortality: pd.DataFrame = pd.read_parquet(
        "data/processed/LMF_Files/LMF__2003-2004__2005-2006__MORT_2019.parquet"
    )
    print(f"Row count from mortality file is {mortality.shape[0]:,d}")
    print(all_nhanes.head())
    print(mortality.head())
    all_nhanes = all_nhanes.merge(mortality, on=["SEQN"], how="left")
    print(f"Row count from demographic file with mortality is {all_nhanes.shape[0]:,d}")

    all_nhanes = load_NHANES_files(
        all_nhanes,
        base_paths,
        suffixes,
        NHANES_VAR_FILES,
        FILE_MAPPING_PRE2005,
        VAR_MAPPING_PRE2005,
    )

    print("=" * 80)
    print("Full dataframe:")
    print(all_nhanes.head())
    print("=" * 80)
    print("Creating variables.")
    print("-" * 80)

    all_nhanes = process_variables(all_nhanes, LOOKUP)

    if save_unfiltered_fname is not None:
        all_nhanes.to_parquet(save_unfiltered_fname)

    if plot_dists:
        plot_all_dists(all_nhanes, LOOKUP, output_path)

    print("=" * 80)
    print("Full dataframe:")
    print(all_nhanes.head())
    print("=" * 80)
    print("Creating missingness plot.")

    all_nhanes_m3s = all_nhanes.loc[:, list(keep)].copy()

    if plot_missingness_diagnostics:
        missingness_correlation(all_nhanes_m3s, output_path)
        missingness_correlation(
            all_nhanes_m3s.loc[all_nhanes.in_death_file == 1, :].copy(),
            output_path,
            suffix="_followup",
        )

    def print_fwp(num, denom):
        """Print Fraction With Percentage"""
        return f"{int(num):,d}/{int(denom):,d} ({num/denom*100:.0f}%)"

    n_total = all_nhanes.shape[0]
    n_lost = (~(all_nhanes.in_death_file == 1)).sum()
    n_remaining = (all_nhanes.in_death_file == 1).sum()
    n_deaths = all_nhanes.is_dead.sum()
    n_deaths_remaining = all_nhanes.loc[all_nhanes.in_death_file == 1, "is_dead"].sum()
    print(
        f"""Missing followup on {print_fwp(n_lost, n_total)} paricipants, {print_fwp(n_remaining, n_total)} remaining with {print_fwp(n_deaths_remaining, n_deaths)} deaths."""
    )
    all_nhanes_req: pd.DataFrame = (
        all_nhanes.loc[all_nhanes.in_death_file == 1, :].dropna(subset=required, how="any").copy()
    )
    n_total = n_remaining
    n_remaining = all_nhanes_req.shape[0]
    n_lost = n_total - n_remaining
    n_deaths_remaining = all_nhanes_req.is_dead.sum()
    print(
        f"""After filtering for all required vars, lost {print_fwp(n_lost, n_total)} particpants with {print_fwp(n_remaining, n_total)} paricipants remain with {print_fwp(n_deaths_remaining, n_deaths)} deaths."""
    )
    if plot_missingness_diagnostics:
        missingness_correlation(
            all_nhanes_req.loc[:, keep].copy(),
            output_path,
            suffix="_req",
        )

    print("=" * 80)
    print("Imputing.")

    all_nhanes_req = impute_defaults(all_nhanes_req, LOOKUP)

    if plot_missingness_diagnostics:
        missingness_correlation(
            all_nhanes_req.loc[:, keep].copy(),
            output_path,
            suffix="_imputed",
        )

    print("=" * 80)
    print("Making more plots on the missingness.")

    if plot_missingness_diagnostics:
        missingness_diagnostics(all_nhanes_req.loc[:, keep].copy(), output_path)

    thresh: int = 40
    all_nhanes_req_thresh: pd.DataFrame = all_nhanes_req.dropna(subset=predictors, thresh=thresh)
    n_total = n_remaining
    n_remaining = all_nhanes_req_thresh.shape[0]
    n_lost = n_total - n_remaining
    n_deaths = n_deaths_remaining
    n_deaths_remaining = all_nhanes_req_thresh.is_dead.sum()
    print(
        f"""After setting threshold of {thresh} vars, lost {print_fwp(n_lost, n_total)} particpants with {print_fwp(n_remaining, n_total)} paricipants remain with {print_fwp(n_deaths_remaining, n_deaths)} deaths."""
    )

    print("=" * 80)
    print("Done.")

    return all_nhanes_req_thresh
