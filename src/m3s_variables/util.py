from collections.abc import Callable
from typing import Any

from numpy import nanmean, nanmin  # noqa: F401
from numpy import round as np_round
from pandas import isna, isnull

# so we can run R expressions too
pmin = min
round = np_round
# and use the ifelse ternary to avoid python specific ternary:
ifelse = lambda boolable, iftrue, iffalse: iftrue if boolable else iffalse  # noqa: E731


def safe_div(val: Any, divisor: float) -> Any:
    """None/NaN-safe division. Returns None for None, NaN for NaN."""
    if val is None:
        return None
    return val / divisor


def safe_mul(val: Any, multiplier: float) -> Any:
    """None/NaN-safe multiplication. Returns None for None, NaN for NaN."""
    if val is None:
        return None
    return val * multiplier


def c(*args) -> list:
    """Behaves like R's c() function"""
    return list(args)


def coalesce_while(*args) -> Any:
    """Return the first non-None value or None if all values are None"""
    i = 0
    while i < len(args):
        if not isna(args[i]):
            return args[i]
        i += 1
    return args[-1]


def coalesce(*values) -> Any:
    """Return the first non-None value or None if all values are None"""
    return next((v for v in values if not isnull(v)), None)


def mylambda(lstr: str, single: bool = False) -> Callable:
    """Create a function that accepts 1 argument, a list,
    and evaluates lstr as an expression with the argument as _."""

    def lambdafun(_, eval_: bool = True, axis: Any = None) -> Any:
        """Evaluate the expr given in lstr with arguments _.

        Pass eval_=False to return the expression string.
        Accept axis kwarg as pd.Series.apply passthrough but doesn't use it."""
        if not eval_:
            return lstr
        # Convert Series to numpy array for consistent integer indexing
        # This avoids the FutureWarning about Series.__getitem__ with integer keys
        if hasattr(_, "values") and hasattr(_, "iloc"):
            _ = _.values
        try:
            return eval(lstr)
        except TypeError:
            # Comparison with None/NaN (e.g., None <= 66666) — treat as NA
            return None

    return lambdafun


def bounds_filter(bounds: list) -> Callable:
    """Create a filter from a list of upper and lower bound.
    Includes both boundaries of the range as valid.
    Returns False for None/NaN values."""
    assert len(bounds) == 2
    assert not (bounds[0] is None and bounds[1] is None)
    if bounds[0] is not None and bounds[1] is not None:
        return lambda _: False if isnull(_) else (_ >= bounds[0] and _ <= bounds[1])
    elif bounds[0] is not None:
        return lambda _: False if isnull(_) else _ >= bounds[0]
    elif bounds[1] is not None:
        return lambda _: False if isnull(_) else _ <= bounds[1]
    # same as the final condition
    return lambda _: False if isnull(_) else _ <= bounds[1]
