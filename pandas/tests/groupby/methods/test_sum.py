import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    Series,
)
import pandas._testing as tm


@pytest.mark.parametrize("min_count", [0, 10])
def test_groupby_sum_mincount_boolean(min_count):
    b = True
    a = False
    na = np.nan
    dfg = pd.array([b, b, na, na, a, a, b], dtype="boolean")

    df = DataFrame({"A": [1, 1, 2, 2, 3, 3, 1], "B": dfg})
    result = df.groupby("A").sum(min_count=min_count)
    if min_count == 0:
        expected = DataFrame(
            {"B": pd.array([3, 0, 0], dtype="Int64")},
            index=Index([1, 2, 3], name="A"),
        )
        tm.assert_frame_equal(result, expected)
    else:
        expected = DataFrame(
            {"B": pd.array([pd.NA] * 3, dtype="Int64")},
            index=Index([1, 2, 3], name="A"),
        )
        tm.assert_frame_equal(result, expected)


def test_groupby_sum_below_mincount_nullable_integer():
    # https://github.com/pandas-dev/pandas/issues/32861
    df = DataFrame({"a": [0, 1, 2], "b": [0, 1, 2], "c": [0, 1, 2]}, dtype="Int64")
    grouped = df.groupby("a")
    idx = Index([0, 1, 2], name="a", dtype="Int64")

    result = grouped["b"].sum(min_count=2)
    expected = Series([pd.NA] * 3, dtype="Int64", index=idx, name="b")
    tm.assert_series_equal(result, expected)

    result = grouped.sum(min_count=2)
    expected = DataFrame({"b": [pd.NA] * 3, "c": [pd.NA] * 3}, dtype="Int64", index=idx)
    tm.assert_frame_equal(result, expected)


def test_groupby_sum_timedelta_with_nat():
    # GH#42659
    df = DataFrame(
        {
            "a": [1, 1, 2, 2],
            "b": [pd.Timedelta("1d"), pd.Timedelta("2d"), pd.Timedelta("3d"), pd.NaT],
        }
    )
    td3 = pd.Timedelta(days=3)

    gb = df.groupby("a")

    res = gb.sum()
    expected = DataFrame({"b": [td3, td3]}, index=Index([1, 2], name="a"))
    tm.assert_frame_equal(res, expected)

    res = gb["b"].sum()
    tm.assert_series_equal(res, expected["b"])

    res = gb["b"].sum(min_count=2)
    expected = Series([td3, pd.NaT], dtype="m8[ns]", name="b", index=expected.index)
    tm.assert_series_equal(res, expected)
