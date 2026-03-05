import numpy as np
import pytest

from pandas.compat import CHAINED_WARNING_DISABLED
from pandas.errors import ChainedAssignmentError

from pandas import DataFrame
import pandas._testing as tm


@pytest.mark.parametrize(
    "indexer", [0, [0, 1], slice(0, 2), np.array([True, False, True])]
)
def test_series_setitem(indexer):
    # ensure we only get a single warning for those typical cases of chained
    # assignment
    df = DataFrame({"a": [1, 2, 3], "b": 1})

    # using custom check instead of tm.assert_produces_warning because that doesn't
    # fail if multiple warnings are raised
    if CHAINED_WARNING_DISABLED:
        return
    with pytest.warns() as record:  # noqa: TID251
        df["a"][indexer] = 0
    assert len(record) == 1
    assert record[0].category == ChainedAssignmentError


@pytest.mark.parametrize(
    "indexer", ["a", ["a", "b"], slice(0, 2), np.array([True, False, True])]
)
def test_frame_setitem(indexer):
    df = DataFrame({"a": [1, 2, 3, 4, 5], "b": 1})

    with tm.raises_chained_assignment_error():
        df[0:3][indexer] = 10


@pytest.mark.parametrize(
    "indexer", [0, [0, 1], slice(0, 2), np.array([True, False, True])]
)
def test_series_iloc_setitem(indexer):
    df = DataFrame({"a": [1, 2, 3], "b": 1})

    with tm.raises_chained_assignment_error():
        df["a"].iloc[indexer] = 0


@pytest.mark.parametrize(
    "indexer", [0, [0, 1], slice(0, 2), np.array([True, False, True])]
)
def test_frame_iloc_setitem(indexer):
    df = DataFrame({"a": [1, 2, 3, 4, 5], "b": 1})

    with tm.raises_chained_assignment_error():
        df[0:3].iloc[indexer] = 10


@pytest.mark.parametrize(
    "indexer", [0, [0, 1], slice(0, 2), np.array([True, False, True])]
)
def test_series_loc_setitem(indexer):
    df = DataFrame({"a": [1, 2, 3], "b": 1})

    with tm.raises_chained_assignment_error():
        df["a"].loc[indexer] = 0


@pytest.mark.parametrize(
    "indexer", [0, [0, 1], (0, "a"), slice(0, 2), np.array([True, False, True])]
)
def test_frame_loc_setitem(indexer):
    df = DataFrame({"a": [1, 2, 3, 4, 5], "b": 1})

    with tm.raises_chained_assignment_error():
        df[0:3].loc[indexer] = 10


def test_series_at_setitem():
    df = DataFrame({"a": [1, 2, 3], "b": 1})

    with tm.raises_chained_assignment_error():
        df["a"].at[0] = 0


def test_frame_at_setitem():
    df = DataFrame({"a": [1, 2, 3, 4, 5], "b": 1})

    with tm.raises_chained_assignment_error():
        df[0:3].at[0, "a"] = 10


def test_series_iat_setitem():
    df = DataFrame({"a": [1, 2, 3], "b": 1})

    with tm.raises_chained_assignment_error():
        df["a"].iat[0] = 0


def test_frame_iat_setitem():
    df = DataFrame({"a": [1, 2, 3, 4, 5], "b": 1})

    with tm.raises_chained_assignment_error():
        df[0:3].iat[0, 0] = 10
