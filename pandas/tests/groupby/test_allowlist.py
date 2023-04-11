"""
TODO: Existing tests should be moved or deduplicated
Do not add tests here!
"""

from string import ascii_lowercase

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
    date_range,
)
import pandas._testing as tm

AGG_FUNCTIONS = [
    "sum",
    "prod",
    "min",
    "max",
    "median",
    "mean",
    "skew",
    "std",
    "var",
    "sem",
]
AGG_FUNCTIONS_WITH_SKIPNA = ["skew"]


@pytest.fixture
def df():
    return DataFrame(
        {
            "A": ["foo", "bar", "foo", "bar", "foo", "bar", "foo", "foo"],
            "B": ["one", "one", "two", "three", "two", "two", "one", "three"],
            "C": np.random.randn(8),
            "D": np.random.randn(8),
        }
    )


@pytest.fixture
def df_letters():
    letters = np.array(list(ascii_lowercase))
    N = 10
    random_letters = letters.take(np.random.randint(0, 26, N))
    df = DataFrame(
        {
            "floats": N / 10 * Series(np.random.random(N)),
            "letters": Series(random_letters),
        }
    )
    return df


@pytest.fixture
def raw_frame():
    return DataFrame([0])


@pytest.mark.parametrize("op", AGG_FUNCTIONS)
@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize("sort", [True, False])
def test_regression_allowlist_methods(raw_frame, op, axis, skipna, sort):
    # GH6944
    # GH 17537
    # explicitly test the allowlist methods
    if axis == 0:
        frame = raw_frame
        msg = "The 'axis' keyword in DataFrame.groupby is deprecated and will be"
    else:
        frame = raw_frame.T
        msg = "DataFrame.groupby with axis=1 is deprecated"

    with tm.assert_produces_warning(FutureWarning, match=msg):
        grouped = frame.groupby(level=0, axis=axis, sort=sort)

    if op in AGG_FUNCTIONS_WITH_SKIPNA:
        result = getattr(grouped, op)(skipna=skipna)
        expected = frame.groupby(level=0).apply(
            lambda h: getattr(h, op)(axis=axis, skipna=skipna)
        )
        if sort:
            expected = expected.sort_index(axis=axis)
        tm.assert_frame_equal(result, expected)
    else:
        result = getattr(grouped, op)()
        expected = frame.groupby(level=0).apply(lambda h: getattr(h, op)(axis=axis))
        if sort:
            expected = expected.sort_index(axis=axis)
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "method",
    [
        "count",
        "corr",
        "cummax",
        "cummin",
        "cumprod",
        "describe",
        "rank",
        "quantile",
        "diff",
        "shift",
        "all",
        "any",
        "idxmin",
        "idxmax",
        "ffill",
        "bfill",
        "pct_change",
    ],
)
def test_groupby_selection_with_methods(df, method):
    # some methods which require DatetimeIndex
    rng = date_range("2014", periods=len(df))
    df.index = rng

    g = df.groupby(["A"])[["C"]]
    g_exp = df[["C"]].groupby(df["A"])
    # TODO check groupby with > 1 col ?

    res = getattr(g, method)()
    exp = getattr(g_exp, method)()

    # should always be frames!
    tm.assert_frame_equal(res, exp)


def test_groupby_selection_other_methods(df):
    # some methods which require DatetimeIndex
    rng = date_range("2014", periods=len(df))
    df.columns.name = "foo"
    df.index = rng

    g = df.groupby(["A"])[["C"]]
    g_exp = df[["C"]].groupby(df["A"])

    # methods which aren't just .foo()
    tm.assert_frame_equal(g.fillna(0), g_exp.fillna(0))
    msg = "DataFrameGroupBy.dtypes is deprecated"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        tm.assert_frame_equal(g.dtypes, g_exp.dtypes)
    tm.assert_frame_equal(g.apply(lambda x: x.sum()), g_exp.apply(lambda x: x.sum()))

    tm.assert_frame_equal(g.resample("D").mean(), g_exp.resample("D").mean())
    tm.assert_frame_equal(g.resample("D").ohlc(), g_exp.resample("D").ohlc())

    tm.assert_frame_equal(
        g.filter(lambda x: len(x) == 3), g_exp.filter(lambda x: len(x) == 3)
    )
