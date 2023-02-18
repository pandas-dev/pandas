"""
test methods relating to generic function evaluation
the so-called white/black lists
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
from pandas.core.groupby.base import (
    groupby_other_methods,
    reduction_kernels,
    transformation_kernels,
)

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
    else:
        frame = raw_frame.T

    if op in AGG_FUNCTIONS_WITH_SKIPNA:
        grouped = frame.groupby(level=0, axis=axis, sort=sort)
        result = getattr(grouped, op)(skipna=skipna)
        expected = frame.groupby(level=0).apply(
            lambda h: getattr(h, op)(axis=axis, skipna=skipna)
        )
        if sort:
            expected = expected.sort_index(axis=axis)
        tm.assert_frame_equal(result, expected)
    else:
        grouped = frame.groupby(level=0, axis=axis, sort=sort)
        result = getattr(grouped, op)()
        expected = frame.groupby(level=0).apply(lambda h: getattr(h, op)(axis=axis))
        if sort:
            expected = expected.sort_index(axis=axis)
        tm.assert_frame_equal(result, expected)


def test_groupby_blocklist(df_letters):
    df = df_letters
    s = df_letters.floats

    blocklist = [
        "eval",
        "query",
        "abs",
        "where",
        "mask",
        "align",
        "groupby",
        "clip",
        "astype",
        "at",
        "combine",
        "consolidate",
        "convert_objects",
    ]
    to_methods = [method for method in dir(df) if method.startswith("to_")]

    blocklist.extend(to_methods)

    for bl in blocklist:
        for obj in (df, s):
            gb = obj.groupby(df.letters)

            # e.g., to_csv
            defined_but_not_allowed = (
                f"(?:^Cannot.+{repr(bl)}.+'{type(gb).__name__}'.+try "
                f"using the 'apply' method$)"
            )

            # e.g., query, eval
            not_defined = (
                f"(?:^'{type(gb).__name__}' object has no attribute {repr(bl)}$)"
            )

            msg = f"{defined_but_not_allowed}|{not_defined}"

            with pytest.raises(AttributeError, match=msg):
                getattr(gb, bl)


def test_tab_completion(mframe):
    grp = mframe.groupby(level="second")
    results = {v for v in dir(grp) if not v.startswith("_")}
    expected = {
        "A",
        "B",
        "C",
        "agg",
        "aggregate",
        "apply",
        "boxplot",
        "filter",
        "first",
        "get_group",
        "groups",
        "hist",
        "indices",
        "last",
        "max",
        "mean",
        "median",
        "min",
        "ngroups",
        "nth",
        "ohlc",
        "plot",
        "prod",
        "size",
        "std",
        "sum",
        "transform",
        "var",
        "sem",
        "count",
        "nunique",
        "head",
        "describe",
        "cummax",
        "quantile",
        "rank",
        "cumprod",
        "tail",
        "resample",
        "cummin",
        "fillna",
        "cumsum",
        "cumcount",
        "ngroup",
        "all",
        "shift",
        "skew",
        "take",
        "pct_change",
        "any",
        "corr",
        "corrwith",
        "cov",
        "dtypes",
        "ndim",
        "diff",
        "idxmax",
        "idxmin",
        "ffill",
        "bfill",
        "rolling",
        "expanding",
        "pipe",
        "sample",
        "ewm",
        "value_counts",
    }
    assert results == expected


def test_groupby_function_rename(mframe):
    grp = mframe.groupby(level="second")
    for name in ["sum", "prod", "min", "max", "first", "last"]:
        f = getattr(grp, name)
        assert f.__name__ == name


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
    tm.assert_frame_equal(g.dtypes, g_exp.dtypes)
    tm.assert_frame_equal(g.apply(lambda x: x.sum()), g_exp.apply(lambda x: x.sum()))

    tm.assert_frame_equal(g.resample("D").mean(), g_exp.resample("D").mean())
    tm.assert_frame_equal(g.resample("D").ohlc(), g_exp.resample("D").ohlc())

    tm.assert_frame_equal(
        g.filter(lambda x: len(x) == 3), g_exp.filter(lambda x: len(x) == 3)
    )


def test_all_methods_categorized(mframe):
    grp = mframe.groupby(mframe.iloc[:, 0])
    names = {_ for _ in dir(grp) if not _.startswith("_")} - set(mframe.columns)
    new_names = set(names)
    new_names -= reduction_kernels
    new_names -= transformation_kernels
    new_names -= groupby_other_methods

    assert not reduction_kernels & transformation_kernels
    assert not reduction_kernels & groupby_other_methods
    assert not transformation_kernels & groupby_other_methods

    # new public method?
    if new_names:
        msg = f"""
There are uncategorized methods defined on the Grouper class:
{new_names}.

Was a new method recently added?

Every public method On Grouper must appear in exactly one the
following three lists defined in pandas.core.groupby.base:
- `reduction_kernels`
- `transformation_kernels`
- `groupby_other_methods`
see the comments in pandas/core/groupby/base.py for guidance on
how to fix this test.
        """
        raise AssertionError(msg)

    # removed a public method?
    all_categorized = reduction_kernels | transformation_kernels | groupby_other_methods
    if names != all_categorized:
        msg = f"""
Some methods which are supposed to be on the Grouper class
are missing:
{all_categorized - names}.

They're still defined in one of the lists that live in pandas/core/groupby/base.py.
If you removed a method, you should update them
"""
        raise AssertionError(msg)
