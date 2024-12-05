from pathlib import Path

import numpy as np
import pandas as pd
from dotenv import find_dotenv
from numba import njit

CHAR_LOOKUP = {
    "2003-2004": "c",
    "2005-2006": "d",
}

PROJECT_ROOT = Path(find_dotenv(".gitignore")).parent


def get_datadir(year: str) -> Path:
    return PROJECT_ROOT / f"data/raw/{year}"


def flatten_columns(cols: list[list]) -> list:
    return ["_".join(col).strip("_") for col in cols]


def worn_indicator(input: np.array, tol: int = 2, cap: int = 100, m: int = 60) -> np.array:
    """
    Use numpy array loop, 50ms per person.

    Applying this using numpy takes ~30min.
    """
    output = np.ones(input.shape[0])

    if ((input[:m] > 0).sum() <= tol) and ((input[:m] < cap).sum() == m):
        output[:m] = 0

    for i in range(m + 1, output.shape[0] + 1):
        if ((input[(i - m) : i] > 0).sum() <= tol) and ((input[(i - m) : i] < cap).sum() == m):
            output[(i - m) : i] = 0

    return output


def worn_indicator_SEQN_long(
    input: np.array, SEQN: np.array, tol: int = 2, cap: int = 100, m: int = 60
) -> np.array:
    """
    Use numpy array loop across the whole DF.

    Total time ~10 minutes.
    """
    output = np.ones(input.shape[0])

    if ((input[:m] > 0).sum() <= tol) and ((input[:m] < cap).sum() == m):
        output[:m] = 0

    for i in range(m + 1, output.shape[0] + 1):
        if (
            (SEQN[i - m] == SEQN[i - 1])
            and ((input[(i - m) : i] > 0).sum() <= tol)
            and ((input[(i - m) : i] < cap).sum() == m)
        ):
            output[(i - m) : i] = 0

    return output


@njit
def worn_indicator_fast(input: np.array, tol: int = 2, cap: int = 100, m: int = 60) -> np.array:
    """Compile with numba, 1ms per person.

    Total time around 7s."""
    output = np.ones(input.shape[0])

    if ((input[:m] > 0).sum() <= tol) and ((input[:m] < cap).sum() == m):
        output[:m] = 0

    for i in range(m + 1, output.shape[0] + 1):
        if ((input[(i - m) : i] > 0).sum() <= tol) and ((input[(i - m) : i] < cap).sum() == m):
            output[(i - m) : i] = 0

    return output


@njit
def worn_indicator_SEQN_long_fast(
    input: np.array, SEQN: np.array, tol: int = 2, cap: int = 100, m: int = 60
) -> np.array:
    output = np.ones(input.shape[0])

    if ((input[:m] > 0).sum() <= tol) and ((input[:m] < cap).sum() == m):
        output[:m] = 0

    for i in range(m + 1, output.shape[0] + 1):
        if (
            (SEQN[i - m] == SEQN[i - 1])
            and ((input[(i - m) : i] > 0).sum() <= tol)
            and ((input[(i - m) : i] < cap).sum() == m)
        ):
            output[(i - m) : i] = 0

    return output


@njit
def bout_classifier_SEQN_long(
    input: np.array,
    SEQN: np.array,
    worn: np.array,
    classification: np.array,
    upper: int,
    lower: int,
    # these are extended ranges, with specified tolerances per extra range
    # if outside the specified range, allow up to the tolerance count of values
    # in this extended range
    tol_upper_soft: int,
    tol_lower_soft: int,
    m: int = 60,
    upper_soft: int = 100,
    lower_soft: int = 0,
    # we also allow the tolerance counts to be shared
    combined_soft_tolerances: bool = False,
    check_worn: bool = True,
    check_already_classified: bool = True,
) -> np.array:
    """Generic algorithm to compute metrics on top of PAXRAW data.

    We check that continuous block of `m` minutes are between the ranges of `(lower, upper]`.
    Beyond this,
    we allow for N=`tol_lower_soft` to fall in `(lower_soft, lower]` and
    N=`tol_upper_soft` to fall in `(upper, upper_soft]`.
    More flexibly, if setting `combined_soft_tolerances=T` then we allow a total of up to
    N=`tol_lower_soft+tol_upper_soft` to fall in the ranges of `(lower_soft, lower]` and `(upper, upper_soft]`.

    We also check that
    - Blocks are classified within SEQN values.
    - `worn` is True for the whole block (controlled by check_worn).
    - Any additional indicator is false, given by `classification` and controlled by flag check_already_classified.


    To use this to mimic the worn indicator itself:
    1. Set m=60, upper=0, upper_soft=100, tol_upper_soft=2, worn=np.zeros(...), check_worn=False.
    2. Set all lower to 0, since this is riding 0 in the LHS.
    3. Invert the output (since we'd actually be computed unworn).
    """

    output = np.zeros(input.shape[0])
    combined_tol = tol_upper_soft + tol_lower_soft

    # the first part, the worn check, is this:
    # ((not check_worn) or (worn[:m].sum() == m))
    # which may look weird but consider the possibilites:
    # check_worn worn output
    # T T T
    # T F F
    # F T T
    # F F T
    valid_worn_status = (not check_worn) or (worn[:m].sum() == m)
    not_already_classified = (not check_already_classified) or (classification[:m].sum() == 0)
    same_person = (SEQN[:m] == SEQN[0]).sum() == m
    too_low = (input[:m] <= lower).sum()
    too_high = (input[:m] > upper).sum()
    # compute these on demand
    # too_low_soft = (input[:m] < lower_soft).sum()
    # too_high_soft = (input[:m] > upper_soft).sum()
    # first, we can skip further checks if we are outside of tolerances (don't need to consider soft bounds)
    if not valid_worn_status or not same_person or not not_already_classified:
        pass
    elif (
        not combined_soft_tolerances and (too_low > tol_lower_soft or too_high > tol_upper_soft)
    ) or (combined_soft_tolerances and (too_low + too_high > combined_tol)):
        pass
    elif (input[:m] <= lower_soft).sum() > 0 or (input[:m] > upper_soft).sum() > 0:
        pass
    else:
        output[:m] = 1

    for i in range(m + 1, output.shape[0] + 1):
        valid_worn_status = (not check_worn) or (worn[(i - m) : i].sum() == m)
        not_already_classified = (not check_already_classified) or (
            classification[(i - m) : i].sum() == 0
        )
        same_person = SEQN[i - m] == SEQN[i - 1]
        too_low = (input[(i - m) : i] <= lower).sum()
        too_high = (input[(i - m) : i] > upper).sum()
        # compute these on demand, since they may not be needed
        # too_low_soft = (input[(i - m) : i] < lower_soft).sum()
        # too_high_soft = (input[(i - m) : i] > upper_soft).sum()
        # first, we can skip further checks if we are outside of tolerances (don't need to consider soft bounds)
        # if not valid_worn_status or not same_person or not not_already_classified:
        #     continue
        # elif (not combined_soft_tolerances and (too_low > tol_lower_soft or too_high > tol_upper_soft)) or (combined_soft_tolerances and (too_low + too_high > combined_tol)):
        #     continue
        # elif (input[:m] < lower_soft).sum() > 0 or (input[:m] > upper_soft).sum() > 0:
        #     continue
        # output[(i - m) : i] = 1

        if (
            valid_worn_status
            and same_person
            and not_already_classified
            and (
                (too_low <= tol_lower_soft and too_high <= tol_upper_soft)
                or (combined_soft_tolerances and (too_low + too_high <= combined_tol))
            )
            and (
                (not combined_soft_tolerances and tol_lower_soft == 0)
                or (combined_tol == 0)
                or ((input[(i - m) : i] <= lower_soft).sum() == 0)
            )
            and (
                (not combined_soft_tolerances and tol_upper_soft == 0)
                or (combined_tol == 0)
                or ((input[(i - m) : i] > upper_soft).sum() == 0)
            )
        ):
            output[(i - m) : i] = 1

    return output


def get_person_active_count(
    d: pd.DataFrame,
    min_worn_hours_threshold: int = 10,
    max_nonzero_count_per_unworn_hour: int = 2,
    max_of_nonzero_in_unworn_hour: int = 100,
    m: int = 60,
) -> float:
    """
    Process out active minutes akin to Fishman (2016)

    1. Compute worn/nonworn indicator on each minute, defined as intervals at least 60 minutes of count = 0, with up to two count < 100.
    2. Sum worn time per day.
    3. Discard days with wear time < 10h.
    4. Sum up total count per day.
    5. Measure average total count per day on valid days, per individual.

    Only max_nonzero_count_per_unworn_hour allowed > 0, and all 60 are strictly less than max_of_nonzero_in_unworn_hour.

    This takes approximately 3s per run (projected 7 hours total)
    and adding the pandas overhead of updates on the vector,
    actually takes around 20 hours.
    """

    # need:
    # - PAXDAY
    # - PAXINTEN
    # - worn
    #
    # positions:
    #     PAXINTEN_i = 7
    #     worn_i = -1
    #
    # check positions:
    # assert d.columns[PAXINTEN_i] == 'PAXINTEN'
    # assert d.columns[worn_i] == 'worn'

    # make a copy of a subset the passed data
    d = d.loc[:, ["PAXDAY", "PAXINTEN"]].copy()
    PAXINTEN_i = 1
    worn_i = 2
    # set the indicator to True to start
    d.loc[:, "worn"] = True
    # take the first hour
    # assert d.iloc[:m, :].shape[0] == m
    if ((d.iloc[:m, PAXINTEN_i] > 0).sum() <= max_nonzero_count_per_unworn_hour) and (
        (d.iloc[:m, PAXINTEN_i] < max_of_nonzero_in_unworn_hour).sum() == m
    ):
        d.iloc[:m, worn_i] = False

    for i in range(m + 1, d.shape[0]):
        # assert test_user.iloc[(i-60):i, :].shape[0] == m
        if ((d.iloc[(i - m) : i, PAXINTEN_i] > 0).sum() <= max_nonzero_count_per_unworn_hour) and (
            (d.iloc[(i - m) : i, PAXINTEN_i] < max_of_nonzero_in_unworn_hour).sum() == m
        ):
            d.iloc[(i - m) : i, worn_i] = False

    # sum minutes of wear and activity counts per day
    worn_minutes = d.groupby("PAXDAY").agg({"worn": [sum], "PAXINTEN": [sum]})

    # compute valid days
    worn_minutes["valid_day"] = worn_minutes["worn"]["sum"] > (min_worn_hours_threshold * m)
    # return here to give more insight from this long process
    return worn_minutes

    # filter to valid days
    # worn_minutes = worn_minutes.loc[worn_minutes.valid_day, :]

    # return the average total daily counts for valid days:
    # return np.mean(worn_minutes['PAXINTEN']['sum'])


def get_person_active_count_rolling(
    d: pd.DataFrame,
    min_worn_hours_threshold: int = 10,
    max_nonzero_count_per_unworn_hour: int = 2,
    max_of_nonzero_in_unworn_hour: int = 100,
    m: int = 60,
) -> float:
    """
    Use pandas .rolling() functionality.

    I recall this being a speedup, but I'm not sure how significant (not very).
    Would need to add a SEQN check to the lambda.
    """

    under_max_nonzero_count_per_unworn_hour = (
        d.PAXINTEN.rolling(window=m)
        .apply(lambda x: (x > 0).sum() <= max_nonzero_count_per_unworn_hour)[m - 1 :]
        .values
    )
    under_max_of_nonzero_in_unworn_hour = (
        d.PAXINTEN.rolling(window=m)
        .apply(lambda x: (x < max_of_nonzero_in_unworn_hour).all())[m - 1 :]
        .values
    )
    unworn = under_max_nonzero_count_per_unworn_hour & under_max_of_nonzero_in_unworn_hour

    base = np.zeros(d.shape[0])
    base[: len(unworn)] = 1 - unworn
    base[len(unworn) :] = 1 - unworn[-1]

    d.loc[:, "worn"] = base

    # sum minutes of wear and activity counts per day
    worn_minutes = d.groupby("PAXDAY").agg({"worn": np.sum, "PAXINTEN": np.sum, "PAXSTEP": np.sum})

    # compute valid days
    worn_minutes["valid_day"] = worn_minutes["worn"]["sum"] > (min_worn_hours_threshold * m)
    # return here to give more insight from this long process
    return worn_minutes
