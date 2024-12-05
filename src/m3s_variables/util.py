from collections.abc import Callable
from typing import Any

from numpy import nanmean, nanmin  # noqa: F401
from numpy import round as np_round
from pandas import isna, isnull

# so we can run R expressions too
pmin = min
round = np_round
# and use the ifelse ternary to avoid python specific ternary:
ifelse = lambda boolable, iftrue, iffalse: [iffalse, iftrue][boolable]  # noqa: E731


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

    def lambdafun(_: list, eval_: bool = True, axis: Any = None) -> Any:
        """Evaluate the expr given in lstr with arguments _.

        Pass eval_=False to return the expression string.
        Accept axis kwarg as pd.Series.apply passthrough but doesn't use it."""
        if not eval_:
            return lstr
        return eval(lstr)

    return lambdafun


def bounds_filter(bounds: list) -> Callable:
    """Create a filter from a list of uppver and lower bound.
    Includes both boundaries of the range as valid."""
    assert len(bounds) == 2
    assert not (bounds[0] is None and bounds[1] is None)
    if bounds[0] is not None and bounds[1] is not None:
        return lambda _: _ >= bounds[0] and _ <= bounds[1]
    elif bounds[0] is not None:
        return lambda _: _ >= bounds[0]
    elif bounds[1] is not None:
        return lambda _: _ <= bounds[1]
    # same as the final condition
    return lambda _: _ <= bounds[1]
