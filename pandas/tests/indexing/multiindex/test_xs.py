import numpy as np
import pytest

from pandas import DataFrame, IndexSlice, MultiIndex, Series
import pandas._testing as tm


def test_xs_level_series(multiindex_dataframe_random_data):
    # this test is not explicitly testing .xs functionality
    # TODO: move to another module or refactor
    df = multiindex_dataframe_random_data
    s = df["A"]
    result = s[:, "two"]
    expected = df.xs("two", level=1)["A"]
    tm.assert_series_equal(result, expected)


def test_xs_level_series_ymd(multiindex_year_month_day_dataframe_random_data):
    # this test is not explicitly testing .xs functionality
    # TODO: move to another module or refactor
    df = multiindex_year_month_day_dataframe_random_data
    s = df["A"]
    result = s[2000, 5]
    expected = df.loc[2000, 5]["A"]
    tm.assert_series_equal(result, expected)


def test_xs_level_series_slice_not_implemented(
    multiindex_year_month_day_dataframe_random_data,
):
    # this test is not explicitly testing .xs functionality
    # TODO: move to another module or refactor
    # not implementing this for now
    df = multiindex_year_month_day_dataframe_random_data
    s = df["A"]

    msg = r"\(2000, slice\(3, 4, None\)\)"
    with pytest.raises(TypeError, match=msg):
        s[2000, 3:4]


def test_xs_IndexSlice_argument_not_implemented():
    # GH 35301

    index = MultiIndex(
        levels=[[("foo", "bar", 0), ("foo", "baz", 0), ("foo", "qux", 0)], [0, 1]],
        codes=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]],
    )

    series = Series(np.random.randn(6), index=index)
    frame = DataFrame(np.random.randn(6, 4), index=index)

    msg = (
        "Expected label or tuple of labels, got "
        r"\(\('foo', 'qux', 0\), slice\(None, None, None\)\)"
    )
    with pytest.raises(TypeError, match=msg):
        frame.xs(IndexSlice[("foo", "qux", 0), :])
    with pytest.raises(TypeError, match=msg):
        series.xs(IndexSlice[("foo", "qux", 0), :])


def test_xs_levels_raises():
    df = DataFrame({"A": [1, 2, 3]})

    msg = "Index must be a MultiIndex"
    with pytest.raises(TypeError, match=msg):
        df.xs(0, level="as")

    s = df.A
    with pytest.raises(TypeError, match=msg):
        s.xs(0, level="as")
