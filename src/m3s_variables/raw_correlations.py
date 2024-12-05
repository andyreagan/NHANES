# coding: utf-8
# coding: utf-8
from pathlib import Path

import altair as alt
import pandas as pd
from altair_saver import save

# alt.data_transformers.enable('json')
alt.data_transformers.enable("default", max_rows=100000)

base_output_path = Path("output", "combined")
years = {("2003-2004", "C"), ("2005-2006", "D")}

d = pd.read_parquet(Path("data") / "processed" / "nhanes.parquet")

chart = (
    alt.Chart(d.loc[:, ["low_density_lipoprotein_ldl_raw", "low_density_lipoprotein_ldl"]])
    .mark_point()
    .encode(
        alt.X("low_density_lipoprotein_ldl_raw", title="LDL from NHANES (raw)"),
        alt.Y("low_density_lipoprotein_ldl", title="Computed LDL"),
    )
)
save(chart, str(base_output_path / "ldl_cor.png"))
d.loc[
    ~d.low_density_lipoprotein_ldl_raw.isna(),
    ["low_density_lipoprotein_ldl_raw", "low_density_lipoprotein_ldl"],
].low_density_lipoprotein_ldl.isna().value_counts()
d.loc[
    ~d.low_density_lipoprotein_ldl.isna(),
    ["low_density_lipoprotein_ldl_raw", "low_density_lipoprotein_ldl"],
].low_density_lipoprotein_ldl_raw.isna().value_counts()

chart = (
    alt.Chart(d.loc[:, ["total_protein_raw", "total_protein"]])
    .mark_point()
    .encode(
        alt.X("total_protein_raw", title="Total protein from NHANES (raw)"),
        alt.Y("total_protein", title="Computed total protein"),
    )
)
save(chart, str(base_output_path / "totalprotein_cor.png"))
d.loc[
    ~d.total_protein_raw.isna(), ["total_protein_raw", "total_protein"]
].total_protein.isna().value_counts(dropna=False)
d.loc[
    ~d.total_protein.isna(), ["total_protein_raw", "total_protein"]
].total_protein_raw.isna().value_counts(dropna=False)
