"""
Apply a simplified version of the `worn` indicator to filter
to valid days from NHANES accelerometry data.
Only considers full hours (noon-1, 1-2, 2-3, etc)
as `worn` or `unworn` in the whole.
This likely isn't as accurate,
but it reduces the complexity by an order of magnitude.
"""

import click
import numpy as np
import pandas as pd
from util import CHAR_LOOKUP, flatten_columns


@click.command()
@click.argument("year")
def main(year: str):
    paxraw_d = pd.read_sas(f"data/raw/{year}/paxraw_{CHAR_LOOKUP[year].lower()}.xpt")

    agg_columns = ["max_intensity", "out_of_calibration", "unreliable"]
    paxraw_d["max_intensity"] = paxraw_d.PAXINTEN == 32767
    paxraw_d["out_of_calibration"] = paxraw_d.PAXCAL == 2
    paxraw_d["unreliable"] = paxraw_d.PAXSTAT == 2

    if "PAXSTEP" in paxraw_d.columns:
        paxraw_d["zero_steps_with_intensity"] = (paxraw_d.PAXINTEN > 250) & (paxraw_d.PAXSTEP == 0)
        paxraw_d["too_many_steps"] = paxraw_d.PAXSTEP > 200
        # add a variable for steps_filtered, summing steps only if we have intensity over 500
        paxraw_d["steps_filtered_500"] = 0
        paxraw_d.loc[paxraw_d.PAXINTEN >= 500, "steps_filtered_500"] = paxraw_d.PAXSTEP
        paxraw_d["steps_filtered_300"] = 0
        paxraw_d.loc[paxraw_d.PAXINTEN >= 300, "steps_filtered_300"] = paxraw_d.PAXSTEP
        agg_columns += [
            "zero_steps_with_intensity",
            "too_many_steps",
            "steps_filtered_500",
            "steps_filtered_300",
        ]

    tudor2009_filters = (
        paxraw_d.groupby(["SEQN", "PAXDAY"])
        .agg({col: [np.sum, "last"] for col in agg_columns})
        .reset_index()
    )
    tudor2009_filters.columns = flatten_columns(tudor2009_filters.columns.values)

    paxraw_hourgroup = (
        paxraw_d.assign(
            nonzero_inten=lambda d: d.PAXINTEN > 0,
            inten_above_thres=lambda d: d.PAXINTEN >= 100,
        )
        .groupby(["SEQN", "PAXDAY", "PAXHOUR"])
        .aggregate(
            {
                "nonzero_inten": lambda d: d.sum() <= 2,
                "inten_above_thres": lambda d: d.sum() == 0,
            }
        )
        .assign(
            unworn=lambda d: d.nonzero_inten & d.inten_above_thres,
            worn=lambda d: ~d.unworn,
        )
        .groupby(["SEQN", "PAXDAY"])
        .agg({"worn": np.sum})
        .reset_index()
    )

    paxraw_hourgroup.columns = flatten_columns(paxraw_hourgroup.columns.values)

    paxraw_hourgroup = paxraw_hourgroup.merge(tudor2009_filters, how="left", on=["SEQN", "PAXDAY"])

    paxraw_hourgroup.to_parquet(
        f"data/processed/{year}/PAXRAW_{CHAR_LOOKUP[year].upper()}_simple.parquet"
    )


if __name__ == "__main__":
    main()
