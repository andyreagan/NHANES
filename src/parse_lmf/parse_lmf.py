from pathlib import Path

import pandas as pd

DATA_DIR = Path("data")


def main():
    all_year_dfs = []
    for year in {"2003_2004", "2005_2006"}:
        d = pd.read_fwf(
            DATA_DIR / "raw/LMF_Files" / f"NHANES_{year}_MORT_2019_PUBLIC.dat",
            colspecs=[
                (0, 14),
                (14, 15),
                (15, 16),
                (16, 19),
                (19, 20),
                (20, 21),
                (21, 22),
                (22, 26),
                (26, 34),
                (34, 42),
                (42, 45),
                (45, 48),
            ],
            skiprows=0,
            header=None,
            na_values=["."],
        )
        d.columns = [
            "SEQN",  # "publicid"
            "eligstat",
            "mortstat",
            "ucod_leading",
            "diabetes",
            "hyperten",
            "dodqtr",
            "dodyear",
            "wgt_new",
            "sa_wgt_new",
            "permth_int",
            "permth_exm",
        ]
        d.head()
        d.dtypes
        d.shape
        d["year"] = int(year.split("_")[0])

        all_year_dfs.append(d)

    all_year_df = pd.concat(all_year_dfs)
    all_year_df.head()
    all_year_df.dtypes
    all_year_df.shape

    all_year_df.to_parquet(
        DATA_DIR / "processed/LMF_Files/LMF__2003-2004__2005-2006__MORT_2019.parquet"
    )


if __name__ == "__main__":
    main()
