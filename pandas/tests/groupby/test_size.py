import numpy as np
import pytest

from pandas import DataFrame, Index, PeriodIndex, Series
import pandas._testing as tm


@pytest.mark.parametrize("by", ["A", "B", ["A", "B"]])
def test_size(df, by):
    grouped = df.groupby(by=by)
    result = grouped.size()
    for key, group in grouped:
        assert result[key] == len(group)


@pytest.mark.parametrize("by", ["A", "B", ["A", "B"]])
@pytest.mark.parametrize("sort", [True, False])
def test_size_sort(df, sort, by):
    df = DataFrame(np.random.choice(20, (1000, 3)), columns=list("ABC"))
    left = df.groupby(by=by, sort=sort).size()
    right = df.groupby(by=by, sort=sort)["C"].apply(lambda a: a.shape[0])
    tm.assert_series_equal(left, right, check_names=False)


def test_size_series_dataframe():
    # https://github.com/pandas-dev/pandas/issues/11699
    df = DataFrame(columns=["A", "B"])
    out = Series(dtype="int64", index=Index([], name="A"))
    tm.assert_series_equal(df.groupby("A").size(), out)


def test_size_groupby_all_null():
    # https://github.com/pandas-dev/pandas/issues/23050
    # Assert no 'Value Error : Length of passed values is 2, index implies 0'
    df = DataFrame({"A": [None, None]})  # all-null groups
    result = df.groupby("A").size()
    expected = Series(dtype="int64", index=Index([], name="A"))
    tm.assert_series_equal(result, expected)


def test_size_period_index():
    # https://github.com/pandas-dev/pandas/issues/34010
    ser = Series([1], index=PeriodIndex(["2000"], name="A", freq="D"))
    grp = ser.groupby(level="A")
    result = grp.size()
    tm.assert_series_equal(result, ser)
