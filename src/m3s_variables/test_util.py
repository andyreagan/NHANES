from numpy import NaN
from pandas import isna

from .util import c, coalesce, coalesce_while, ifelse, round


def test_ifelse() -> None:
    assert ifelse("a" == "a", "a", "b") == "a"
    assert ifelse(True, "a", "b") == "a"
    assert ifelse("a" == "b", "a", "b") == "b"
    assert ifelse(False, "a", "b") == "b"
    # this doesn't work:
    # assert ifelse(NaN, "a", "b") == "b"


def test_round() -> None:
    # make sure that we handle NaN appropriately
    assert round(0.4) == 0.0
    assert isna(round(NaN))


def test_coalesce_while() -> None:
    # make sure that we handle NaN appropriately
    assert coalesce_while(None, NaN, 0) == 0


def test_coalesce() -> None:
    # make sure that we handle NaN appropriately
    assert coalesce(None, NaN, 0) == 0
    assert coalesce(None, NaN, None) is None
    assert coalesce(None, None) is None


def test_c() -> None:
    assert len(c(1, 2, 3)) == 3
    assert c(1, 2, 3)[0] == 1
    assert c(1, 2, 3)[1] == 2
    assert c(1, 2, 3)[2] == 3


def test_all_helpers() -> None:
    test_ifelse()
    test_coalesce()
    test_coalesce_while()
    test_round()
    test_c()


if __name__ == "__main__":
    test_all_helpers()
