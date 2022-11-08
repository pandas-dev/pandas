import io

import pytest

from pandas._config import option_context

import pandas as pd
from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.core.indexes.api import (
    Float64Index,
    Index,
    Int64Index,
    NoIndex,
    RangeIndex,
)

# aliases to make some tests easier to read
NI = NoIndex
RI = RangeIndex
I64 = Int64Index
F64 = Float64Index
OI = Index


# class TestNoIndex(TestRangeIndex):
#     _index_cls = NoIndex


@pytest.fixture
def ser1():
    with option_context("mode.no_default_index", True):
        res = Series([1, 2, 3])
    return res


@pytest.fixture
def ser2():
    with option_context("mode.no_default_index", True):
        res = Series([4, 5, 6])
    return res


@pytest.fixture
def df1():
    with option_context("mode.no_default_index", True):
        res = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    return res


@pytest.fixture
def df2():
    with option_context("mode.no_default_index", True):
        res = DataFrame({"a": [6, 5, 4], "b": [1, 4, 2]})
    return res


def test_boolean_mask(ser1):
    mask = ser1 > 1
    result = ser1[mask]
    expected = Series([2, 3], index=NoIndex(2))
    tm.assert_series_equal(result, expected)

    ser1[mask] = 4
    expected = Series([1, 4, 4], index=NoIndex(3))
    tm.assert_series_equal(ser1, expected)


def test_join(ser1, ser2, df1, df2):
    result = df1.join(df2, lsuffix="_df1")
    expected = DataFrame(
        {
            "a_df1": [1, 2, 3],
            "b_df1": [4, 5, 6],
            "a": [6, 5, 4],
            "b": [1, 4, 2],
        },
        index=NoIndex(3),
    )
    tm.assert_frame_equal(result, expected)


def test_loc(df1):
    with pytest.raises(IndexError, match="Cannot use label-based indexing on NoIndex!"):
        df1.loc[0, "a"]
    with pytest.raises(IndexError, match="Cannot use label-based indexing on NoIndex!"):
        df1.loc[0]

    result = df1.loc[:, "a"]
    expected = Series([1, 2, 3], index=NoIndex(3), name="a")
    tm.assert_series_equal(result, expected)

    mask = df1["a"] > 2
    result = df1.loc[mask]
    expected = DataFrame({"a": [3], "b": [6]}, index=NoIndex(1))
    tm.assert_frame_equal(result, expected)

    result = df1.loc[df1["a"] > 2, "a"]
    expected = Series([3], index=NoIndex(1), name="a")
    tm.assert_series_equal(result, expected)

    result = df1.iloc[1:]
    expected = DataFrame(
        {
            "a": [
                2,
                3,
            ],
            "b": [5, 6],
        }
    )
    tm.assert_frame_equal(result, expected)


def test_alignment(df1):
    with pytest.raises(TypeError, match="Can't join NoIndex of different lengths"):
        result = df1 + df1.iloc[1:]
    result = df1 + df1
    expected = DataFrame({"a": [2, 4, 6], "b": [8, 10, 12]}, index=NoIndex(3))
    tm.assert_frame_equal(result, expected)


def test_reader():
    with option_context("mode.no_default_index", True):
        result = pd.read_csv(io.StringIO("data\n1\n"))
    expected = DataFrame({"data": [1]}, index=NoIndex(1))
    tm.assert_frame_equal(result, expected)


def test_repr():
    with option_context("mode.no_default_index", True):
        df = DataFrame({"a": [1, 2, 3] * 50})
    result = repr(df)
    expected = (
        " a\n"
        " 1\n"
        " 2\n"
        " 3\n"
        " 1\n"
        " 2\n"
        "..\n"
        " 2\n"
        " 3\n"
        " 1\n"
        " 2\n"
        " 3\n"
        "\n[150 rows x 1 columns]"
    )
    assert result == expected

    result = repr(df["a"])
    expected = (
        "1 \n"
        "2 \n"
        "3 \n"
        "1 \n"
        "2 \n"
        "..\n"
        "2 \n"
        "3 \n"
        "1 \n"
        "2 \n"
        "3 \n"
        "Name: a, Length: 150, dtype: int64"
    )
    assert result == expected


def test_concat(df1):
    result = pd.concat([df1, df1])
    expected = DataFrame(
        {
            "a": [1, 2, 3, 1, 2, 3],
            "b": [4, 5, 6, 4, 5, 6],
        },
        index=NoIndex(6),
    )
    tm.assert_frame_equal(result, expected)


def test_merge(df1):
    result = df1.merge(df1, on="a")
    expected = DataFrame(
        {
            "a": [1, 2, 3],
            "b_x": [4, 5, 6],
            "b_y": [4, 5, 6],
        },
        index=NoIndex(3),
    )
    tm.assert_frame_equal(result, expected)
