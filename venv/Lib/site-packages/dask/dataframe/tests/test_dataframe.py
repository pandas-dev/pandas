from __future__ import annotations

import contextlib
import decimal
import sys
import warnings
import weakref
from datetime import datetime, timedelta
from itertools import product
from operator import add

import numpy as np
import pandas as pd
import pytest
from pandas.errors import PerformanceWarning

import dask
import dask.array as da
import dask.dataframe as dd
import dask.dataframe.groupby
from dask import delayed
from dask._compatibility import WINDOWS
from dask.base import collections_to_dsk
from dask.dataframe import _compat, methods
from dask.dataframe._compat import PANDAS_GE_210, PANDAS_GE_220, PANDAS_GE_300, tm
from dask.dataframe._pyarrow import to_pyarrow_string
from dask.dataframe.core import (
    _concat,
    _map_freq_to_period_start,
    has_parallel_type,
    total_mem_usage,
)
from dask.dataframe.dispatch import meta_nonempty
from dask.dataframe.utils import (
    assert_eq,
    assert_eq_dtypes,
    assert_max_deps,
    get_string_dtype,
    make_meta,
    pyarrow_strings_enabled,
)
from dask.utils import M, is_dataframe_like, is_series_like

try:
    import crick
except ImportError:
    crick = None

try:
    from pyarrow.lib import ArrowNotImplementedError
except ImportError:
    ArrowNotImplementedError = RuntimeError  # some unrelated error to make pytest pass

dsk = {
    ("x", 0): pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=[0, 1, 3]),
    ("x", 1): pd.DataFrame({"a": [4, 5, 6], "b": [3, 2, 1]}, index=[5, 6, 8]),
    ("x", 2): pd.DataFrame({"a": [7, 8, 9], "b": [0, 0, 0]}, index=[9, 9, 9]),
}
meta = make_meta(
    {"a": "i8", "b": "i8"}, index=pd.Index([], "i8"), parent_meta=pd.DataFrame()
)
d = dd.repartition(pd.concat(dsk.values()), divisions=[0, 5, 9, 9])
full = d.compute()


def _drop_mean(df, col=None):
    """TODO: In pandas 2.0, mean is implemented for datetimes, but Dask returns None."""
    if isinstance(df, pd.DataFrame):
        df.at["mean", col] = np.nan
        df.dropna(how="all", inplace=True)
    elif isinstance(df, pd.Series):
        df.drop(labels=["mean"], inplace=True, errors="ignore")
    else:
        raise NotImplementedError("Expected Series or DataFrame with mean")
    return df


def test_Dataframe():
    expected = pd.Series(
        [2, 3, 4, 5, 6, 7, 8, 9, 10], index=[0, 1, 3, 5, 6, 8, 9, 9, 9], name="a"
    )

    assert_eq(d["a"] + 1, expected)

    tm.assert_index_equal(d.columns, pd.Index(["a", "b"]))

    assert_eq(d[d["b"] > 2], full[full["b"] > 2])
    assert_eq(d[["a", "b"]], full[["a", "b"]])
    assert_eq(d.a, full.a)
    assert d.b.mean().compute() == full.b.mean()
    assert np.allclose(d.b.var().compute(), full.b.var())
    assert np.allclose(d.b.std().compute(), full.b.std())

    assert d.index._name == d.index._name  # this is deterministic

    assert repr(d)


def test_head_tail():
    assert_eq(d.head(2), full.head(2))
    assert_eq(d.head(3), full.head(3))
    assert_eq(d.head(2), dsk[("x", 0)].head(2))
    assert_eq(d["a"].head(2), full["a"].head(2))
    assert_eq(d["a"].head(3), full["a"].head(3))
    assert_eq(d["a"].head(2), dsk[("x", 0)]["a"].head(2))
    assert sorted(d.head(2, compute=False).dask) == sorted(
        d.head(2, compute=False).dask
    )
    assert sorted(d.head(2, compute=False).dask) != sorted(
        d.head(3, compute=False).dask
    )

    assert_eq(d.tail(2), full.tail(2))
    assert_eq(d.tail(3), full.tail(3))
    assert_eq(d.tail(2), dsk[("x", 2)].tail(2))
    assert_eq(d["a"].tail(2), full["a"].tail(2))
    assert_eq(d["a"].tail(3), full["a"].tail(3))
    assert_eq(d["a"].tail(2), dsk[("x", 2)]["a"].tail(2))
    assert sorted(d.tail(2, compute=False).dask) == sorted(
        d.tail(2, compute=False).dask
    )
    assert sorted(d.tail(2, compute=False).dask) != sorted(
        d.tail(3, compute=False).dask
    )


def test_head_npartitions():
    assert_eq(d.head(5, npartitions=2), full.head(5))
    assert_eq(d.head(5, npartitions=2, compute=False), full.head(5))
    assert_eq(d.head(5, npartitions=-1), full.head(5))
    assert_eq(d.head(7, npartitions=-1), full.head(7))
    assert_eq(d.head(2, npartitions=-1), full.head(2))
    with pytest.raises(ValueError):
        d.head(2, npartitions=5)


def test_head_npartitions_warn():
    match = "5 elements requested, only 3 elements"
    with pytest.warns(UserWarning, match=match):
        d.head(5)

    match = "Insufficient elements"
    with pytest.warns(UserWarning, match=match):
        d.head(100)

    with pytest.warns(UserWarning, match=match):
        d.head(7)

    with pytest.warns(UserWarning, match=match):
        d.head(7, npartitions=2)

    # No warn if all partitions are inspected
    for n in [3, -1]:
        with warnings.catch_warnings(record=True) as record:
            d.head(10, npartitions=n)
        assert not record

    # With default args, this means that a 1 partition dataframe won't warn
    d2 = dd.from_pandas(pd.DataFrame({"x": [1, 2, 3]}), npartitions=1)
    with warnings.catch_warnings(record=True) as record:
        d2.head()
    assert not record


def test_index_head():
    assert_eq(d.index.head(2), full.index[:2])
    assert_eq(d.index.head(3), full.index[:3])


def test_Series():
    assert isinstance(d.a, dd.Series)
    assert isinstance(d.a + 1, dd.Series)
    assert_eq((d + 1), full + 1)


def test_Index():
    for case in [
        pd.DataFrame(np.random.randn(10, 5), index=list("abcdefghij")),
        pd.DataFrame(
            np.random.randn(10, 5),
            index=pd.date_range("2011-01-01", freq="D", periods=10),
        ),
    ]:
        ddf = dd.from_pandas(case, 3)
        assert_eq(ddf.index, case.index)
        pytest.raises(AttributeError, lambda ddf=ddf: ddf.index.index)


def test_axes():
    pdf = pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})
    df = dd.from_pandas(pdf, npartitions=2)
    assert len(df.axes) == len(pdf.axes)
    assert all(assert_eq(d, p) for d, p in zip(df.axes, pdf.axes))


def test_series_axes():
    ps = pd.Series(["abcde"])
    ds = dd.from_pandas(ps, npartitions=2)
    assert len(ds.axes) == len(ps.axes)
    assert all(assert_eq(d, p) for d, p in zip(ds.axes, ps.axes))


def test_attributes():
    assert "a" in dir(d)
    assert "foo" not in dir(d)
    pytest.raises(AttributeError, lambda: d.foo)

    df = dd.from_pandas(pd.DataFrame({"a b c": [1, 2, 3]}), npartitions=2)
    assert "a b c" not in dir(df)
    df = dd.from_pandas(pd.DataFrame({"a": [1, 2], 5: [1, 2]}), npartitions=2)
    assert "a" in dir(df)
    assert 5 not in dir(df)

    df = dd.from_pandas(_compat.makeTimeDataFrame(), npartitions=3)
    pytest.raises(AttributeError, lambda: df.foo)


def test_column_names():
    tm.assert_index_equal(d.columns, pd.Index(["a", "b"]))
    tm.assert_index_equal(d[["b", "a"]].columns, pd.Index(["b", "a"]))
    assert d["a"].name == "a"
    assert (d["a"] + 1).name == "a"
    assert (d["a"] + d["b"]).name is None


def test_columns_named_divisions_and_meta():
    # https://github.com/dask/dask/issues/7599
    df = pd.DataFrame(
        {"_meta": [1, 2, 3, 4], "divisions": ["a", "b", "c", "d"]},
        index=[0, 1, 3, 5],
    )
    ddf = dd.from_pandas(df, 2)

    assert ddf.divisions == (0, 3, 5)
    assert_eq(ddf["divisions"], df.divisions)
    assert all(ddf._meta.columns == ["_meta", "divisions"])
    assert_eq(ddf["_meta"], df._meta)


def test_index_names():
    assert d.index.name is None

    idx = pd.Index([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], name="x")
    df = pd.DataFrame(np.random.randn(10, 5), idx)
    ddf = dd.from_pandas(df, 3)
    assert ddf.index.name == "x"
    assert ddf.index.compute().name == "x"


def test_rename_columns():
    # GH 819
    df = pd.DataFrame({"a": [1, 2, 3, 4, 5, 6, 7], "b": [7, 6, 5, 4, 3, 2, 1]})
    ddf = dd.from_pandas(df, 2)

    ddf.columns = ["x", "y"]
    df.columns = ["x", "y"]
    tm.assert_index_equal(ddf.columns, pd.Index(["x", "y"]))
    tm.assert_index_equal(ddf._meta.columns, pd.Index(["x", "y"]))
    assert_eq(ddf, df)

    msg = r"Length mismatch: Expected axis has 2 elements, new values have 4 elements"
    with pytest.raises(ValueError) as err:
        ddf.columns = [1, 2, 3, 4]
    assert msg in str(err.value)

    # Multi-index columns
    df = pd.DataFrame({("A", "0"): [1, 2, 2, 3], ("B", 1): [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)

    df.columns = ["x", "y"]
    ddf.columns = ["x", "y"]
    tm.assert_index_equal(ddf.columns, pd.Index(["x", "y"]))
    tm.assert_index_equal(ddf._meta.columns, pd.Index(["x", "y"]))
    assert_eq(ddf, df)


def test_rename_series():
    # GH 819
    s = pd.Series([1, 2, 3, 4, 5, 6, 7], name="x")
    ds = dd.from_pandas(s, 2)

    s.name = "renamed"
    ds.name = "renamed"
    assert s.name == "renamed"
    assert_eq(ds, s)

    ind = s.index
    dind = ds.index
    ind.name = "renamed"
    dind.name = "renamed"
    assert ind.name == "renamed"
    assert_eq(dind, ind)


def test_rename_series_method():
    # Series name
    s = pd.Series([1, 2, 3, 4, 5, 6, 7], name="x")
    ds = dd.from_pandas(s, 2)

    assert_eq(ds.rename("y"), s.rename("y"))
    assert ds.name == "x"  # no mutation
    with pytest.raises(TypeError):
        ds.rename()

    assert_eq(ds, s)


def test_rename_series_method_2():
    # Series index
    s = pd.Series(["a", "b", "c", "d", "e", "f", "g"], name="x")
    ds = dd.from_pandas(s, 2)

    for is_sorted in [True, False]:
        res = ds.rename(lambda x: x**2, sorted_index=is_sorted)
        assert_eq(res, s.rename(lambda x: x**2))
        assert res.known_divisions == is_sorted

        res = ds.rename(s, sorted_index=is_sorted)
        assert_eq(res, s.rename(s))
        assert res.known_divisions == is_sorted

    with pytest.raises(ValueError):
        ds.rename(lambda x: -x, sorted_index=True).divisions
    assert_eq(ds.rename(lambda x: -x), s.rename(lambda x: -x))

    res = ds.rename(ds)
    assert_eq(res, s.rename(s))
    assert not res.known_divisions

    ds2 = ds.clear_divisions()
    res = ds2.rename(lambda x: x**2, sorted_index=True)
    assert_eq(res, s.rename(lambda x: x**2))
    assert not res.known_divisions


@pytest.mark.parametrize(
    "method,test_values",
    [
        pytest.param(
            "tdigest",
            (6, 10),
            marks=pytest.mark.skipif(not crick, reason="Requires crick"),
        ),
        ("dask", (4, 20)),
    ],
)
def test_describe_numeric(method, test_values):
    # prepare test case which approx quantiles will be the same as actuals
    s = pd.Series(list(range(test_values[1])) * test_values[0])
    df = pd.DataFrame(
        {
            "a": list(range(test_values[1])) * test_values[0],
            "b": list(range(test_values[0])) * test_values[1],
        }
    )

    ds = dd.from_pandas(s, test_values[0])
    ddf = dd.from_pandas(df, test_values[0])

    test_quantiles = [0.25, 0.75]

    assert_eq(df.describe(), ddf.describe(percentiles_method=method))
    assert_eq(s.describe(), ds.describe(percentiles_method=method))

    assert_eq(
        df.describe(percentiles=test_quantiles),
        ddf.describe(percentiles=test_quantiles, percentiles_method=method),
    )
    assert_eq(s.describe(), ds.describe(split_every=2, percentiles_method=method))
    assert_eq(df.describe(), ddf.describe(split_every=2, percentiles_method=method))

    # remove string columns
    df = pd.DataFrame(
        {
            "a": list(range(test_values[1])) * test_values[0],
            "b": list(range(test_values[0])) * test_values[1],
            "c": list("abcdef"[: test_values[0]]) * test_values[1],
        }
    )
    ddf = dd.from_pandas(df, test_values[0])
    assert_eq(df.describe(), ddf.describe(percentiles_method=method))
    assert_eq(df.describe(), ddf.describe(split_every=2, percentiles_method=method))


@pytest.mark.parametrize(
    "include,exclude,percentiles,subset",
    [
        (None, None, None, ["c", "d"]),  # numeric
        (None, None, None, ["c", "d", "f"]),  # numeric + timedelta
        (None, None, None, ["c", "d", "g"]),  # numeric + bool
        (None, None, None, ["c", "d", "f", "g"]),  # numeric + bool + timedelta
        (None, None, None, ["f", "g"]),  # bool + timedelta
        ("all", None, None, None),
        (["number"], None, [0.25, 0.5], None),
        ([np.timedelta64], None, None, None),
        (["number", get_string_dtype()], None, [0.25, 0.75], None),
        (None, ["number", get_string_dtype()], None, None),
        ([get_string_dtype(), "datetime", "bool"], None, None, None),
    ],
)
def test_describe(include, exclude, percentiles, subset):
    data = {
        "a": ["aaa", "bbb", "bbb", None, None, "zzz"] * 2,
        "c": [None, 0, 1, 2, 3, 4] * 2,
        "d": [None, 0, 1] * 4,
        "e": [
            pd.Timestamp("2017-05-09 00:00:00.006000"),
            pd.Timestamp("2017-05-09 00:00:00.006000"),
            pd.Timestamp("2017-05-09 07:56:23.858694"),
            pd.Timestamp("2017-05-09 05:59:58.938999"),
            None,
            None,
        ]
        * 2,
        "f": [
            np.timedelta64(3, "D"),
            np.timedelta64(1, "D"),
            None,
            None,
            np.timedelta64(3, "D"),
            np.timedelta64(1, "D"),
        ]
        * 2,
        "g": [True, False, True] * 4,
    }

    # Arrange
    df = pd.DataFrame(data)
    df["a"] = df["a"].astype(get_string_dtype())

    if subset is not None:
        df = df.loc[:, subset]

    ddf = dd.from_pandas(df, 2)

    # Act
    actual = ddf.describe(
        include=include,
        exclude=exclude,
        percentiles=percentiles,
    )
    expected = df.describe(
        include=include,
        exclude=exclude,
        percentiles=percentiles,
    )

    if "e" in expected:
        expected = _drop_mean(expected, "e")

    assert_eq(actual, expected)

    # Check series
    if subset is None:
        for col in ["a", "c", "e", "g"]:
            expected = df[col].describe(include=include, exclude=exclude)
            if col == "e":
                expected = _drop_mean(expected)
            actual = ddf[col].describe(include=include, exclude=exclude)
            assert_eq(expected, actual)


def test_describe_without_datetime_is_numeric():
    data = {
        "a": ["aaa", "bbb", "bbb", None, None, "zzz"] * 2,
        "c": [None, 0, 1, 2, 3, 4] * 2,
        "d": [None, 0, 1] * 4,
        "e": [
            pd.Timestamp("2017-05-09 00:00:00.006000"),
            pd.Timestamp("2017-05-09 00:00:00.006000"),
            pd.Timestamp("2017-05-09 07:56:23.858694"),
            pd.Timestamp("2017-05-09 05:59:58.938999"),
            None,
            None,
        ]
        * 2,
    }
    # Arrange
    df = pd.DataFrame(data)
    ddf = dd.from_pandas(df, 2)

    # Assert
    expected = df.describe()
    expected = _drop_mean(expected, "e")

    assert_eq(ddf.describe(), expected)

    # Check series
    for col in ["a", "c"]:
        assert_eq(df[col].describe(), ddf[col].describe())

    expected = _drop_mean(df.e.describe())
    assert_eq(expected, ddf.e.describe())


# Note: this warning is not always raised on Windows
@pytest.mark.filterwarnings(
    "ignore:invalid value encountered in (true_)?divide:RuntimeWarning"
)
def test_describe_empty():
    df_none = pd.DataFrame({"A": [None, None]})
    ddf_none = dd.from_pandas(df_none, 2)
    df_len0 = pd.DataFrame({"A": [], "B": []})
    ddf_len0 = dd.from_pandas(df_len0, 2)
    ddf_nocols = dd.from_pandas(pd.DataFrame({}), 2)

    assert_eq(df_none.describe(), ddf_none.describe(percentiles_method="dask"))
    assert_eq(df_len0.describe(), ddf_len0.describe(percentiles_method="dask"))
    with pytest.raises(ValueError):
        ddf_nocols.describe(percentiles_method="dask").compute()


def test_describe_empty_tdigest():
    pytest.importorskip("crick")

    df_none = pd.DataFrame({"A": [None, None]})
    ddf_none = dd.from_pandas(df_none, 2)
    df_len0 = pd.DataFrame({"A": []})
    ddf_len0 = dd.from_pandas(df_len0, 2)
    ddf_nocols = dd.from_pandas(pd.DataFrame({}), 2)

    # Pandas have different dtypes for resulting describe dataframe if there are only
    # None-values, pre-compute dask df to bypass _meta check
    assert_eq(
        df_none.describe(), ddf_none.describe(percentiles_method="tdigest").compute()
    )
    with warnings.catch_warnings():
        # dask.dataframe should probably filter this, to match pandas, but
        # it seems quite difficult.
        warnings.simplefilter("ignore", RuntimeWarning)
        assert_eq(df_len0.describe(), ddf_len0.describe(percentiles_method="tdigest"))
        assert_eq(df_len0.describe(), ddf_len0.describe(percentiles_method="tdigest"))

    with pytest.raises(ValueError):
        ddf_nocols.describe(percentiles_method="tdigest").compute()


def test_describe_for_possibly_unsorted_q():
    """make sure describe is sorting percentiles parameter, q, properly and can
    handle lists, tuples and ndarrays.

    See https://github.com/dask/dask/issues/4642.
    """
    # prepare test case where quantiles should equal values
    A = da.arange(0, 101)
    ds = dd.from_dask_array(A)

    for q in [None, [0.25, 0.50, 0.75], [0.25, 0.50, 0.75, 0.99], [0.75, 0.5, 0.25]]:
        for f_convert in [list, tuple, np.array]:
            if q is None:
                r = ds.describe(percentiles=q).compute()
            else:
                r = ds.describe(percentiles=f_convert(q)).compute()

            assert_eq(r["25%"], 25.0)
            assert_eq(r["50%"], 50.0)
            assert_eq(r["75%"], 75.0)


@pytest.mark.skipif(
    WINDOWS and sys.version_info < (3, 11),
    reason="https://github.com/dask/dask/pull/11320#issuecomment-2293798597",
)
def test_cumulative():
    index = [f"row{i:03d}" for i in range(100)]
    df = pd.DataFrame(np.random.randn(100, 5), columns=list("abcde"), index=index)

    ddf = dd.from_pandas(df, 5)

    assert_eq(ddf.cumsum(), df.cumsum())
    assert_eq(ddf.cumprod(), df.cumprod())
    assert_eq(ddf.cummin(), df.cummin())
    assert_eq(ddf.cummax(), df.cummax())

    assert_eq(ddf.cumsum(axis=1), df.cumsum(axis=1))
    assert_eq(ddf.cumprod(axis=1), df.cumprod(axis=1))
    assert_eq(ddf.cummin(axis=1), df.cummin(axis=1))
    assert_eq(ddf.cummax(axis=1), df.cummax(axis=1))

    assert_eq(ddf.a.cumsum(), df.a.cumsum())
    assert_eq(ddf.a.cumprod(), df.a.cumprod())
    assert_eq(ddf.a.cummin(), df.a.cummin())
    assert_eq(ddf.a.cummax(), df.a.cummax())

    assert_eq(np.cumsum(ddf), np.cumsum(df))
    assert_eq(np.cumprod(ddf), np.cumprod(df))
    assert_eq(np.cumsum(ddf, axis=1), np.cumsum(df, axis=1))
    assert_eq(np.cumprod(ddf, axis=1), np.cumprod(df, axis=1))
    assert_eq(np.cumsum(ddf.a), np.cumsum(df.a))
    assert_eq(np.cumprod(ddf.a), np.cumprod(df.a))


def test_cumulative_with_nans():
    df = pd.DataFrame(
        {
            "a": [1, 2, np.nan, 4, 5, 6, 7, 8],
            "b": [1, 2, np.nan, np.nan, np.nan, 5, np.nan, np.nan],
            "c": [np.nan] * 8,
        }
    )
    ddf = dd.from_pandas(df, 3)

    assert_eq(df.cumsum(), ddf.cumsum())
    assert_eq(df.cummin(), ddf.cummin())
    assert_eq(df.cummax(), ddf.cummax())
    assert_eq(df.cumprod(), ddf.cumprod())

    assert_eq(df.cumsum(skipna=False), ddf.cumsum(skipna=False))
    assert_eq(df.cummin(skipna=False), ddf.cummin(skipna=False))
    assert_eq(df.cummax(skipna=False), ddf.cummax(skipna=False))
    assert_eq(df.cumprod(skipna=False), ddf.cumprod(skipna=False))

    assert_eq(df.cumsum(axis=1), ddf.cumsum(axis=1))
    assert_eq(df.cummin(axis=1), ddf.cummin(axis=1))
    assert_eq(df.cummax(axis=1), ddf.cummax(axis=1))
    assert_eq(df.cumprod(axis=1), ddf.cumprod(axis=1))

    assert_eq(df.cumsum(axis=1, skipna=False), ddf.cumsum(axis=1, skipna=False))
    assert_eq(df.cummin(axis=1, skipna=False), ddf.cummin(axis=1, skipna=False))
    assert_eq(df.cummax(axis=1, skipna=False), ddf.cummax(axis=1, skipna=False))
    assert_eq(df.cumprod(axis=1, skipna=False), ddf.cumprod(axis=1, skipna=False))


def test_cumulative_with_duplicate_columns():
    df = pd.DataFrame(np.random.randn(100, 3), columns=list("abb"))
    ddf = dd.from_pandas(df, 3)

    assert_eq(df.cumsum(), ddf.cumsum())
    assert_eq(df.cummin(), ddf.cummin())
    assert_eq(df.cummax(), ddf.cummax())
    assert_eq(df.cumprod(), ddf.cumprod())

    assert_eq(df.cumsum(skipna=False), ddf.cumsum(skipna=False))
    assert_eq(df.cummin(skipna=False), ddf.cummin(skipna=False))
    assert_eq(df.cummax(skipna=False), ddf.cummax(skipna=False))
    assert_eq(df.cumprod(skipna=False), ddf.cumprod(skipna=False))

    assert_eq(df.cumsum(axis=1), ddf.cumsum(axis=1))
    assert_eq(df.cummin(axis=1), ddf.cummin(axis=1))
    assert_eq(df.cummax(axis=1), ddf.cummax(axis=1))
    assert_eq(df.cumprod(axis=1), ddf.cumprod(axis=1))

    assert_eq(df.cumsum(axis=1, skipna=False), ddf.cumsum(axis=1, skipna=False))
    assert_eq(df.cummin(axis=1, skipna=False), ddf.cummin(axis=1, skipna=False))
    assert_eq(df.cummax(axis=1, skipna=False), ddf.cummax(axis=1, skipna=False))
    assert_eq(df.cumprod(axis=1, skipna=False), ddf.cumprod(axis=1, skipna=False))


@pytest.mark.parametrize(
    "func",
    [
        M.cumsum,
        M.cumprod,
        pytest.param(
            M.cummin,
            marks=[
                pytest.mark.xfail(
                    reason="ValueError: Can only compare identically-labeled Series objects"
                )
            ],
        ),
        pytest.param(
            M.cummax,
            marks=[
                pytest.mark.xfail(
                    reason="ValueError: Can only compare identically-labeled Series objects"
                )
            ],
        ),
    ],
)
def test_cumulative_empty_partitions(func):
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=4)
    assert_eq(func(df[df.x < 5]), func(ddf[ddf.x < 5]))

    df = pd.DataFrame({"x": [1, 2, 3, 4, None, 5, 6, None, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=5)
    assert_eq(func(df[df.x < 5]), func(ddf[ddf.x < 5]))


def test_dropna():
    df = pd.DataFrame(
        {
            "x": [np.nan, 2, 3, 4, np.nan, 6],
            "y": [1, 2, np.nan, 4, np.nan, np.nan],
            "z": [1, 2, 3, 4, np.nan, 6],
        },
        index=[10, 20, 30, 40, 50, 60],
    )
    ddf = dd.from_pandas(df, 3)

    assert_eq(ddf.x.dropna(), df.x.dropna())
    assert_eq(ddf.y.dropna(), df.y.dropna())
    assert_eq(ddf.z.dropna(), df.z.dropna())

    assert_eq(ddf.dropna(), df.dropna())
    assert_eq(ddf.dropna(how="all"), df.dropna(how="all"))
    assert_eq(ddf.dropna(subset=["x"]), df.dropna(subset=["x"]))
    assert_eq(ddf.dropna(subset=["y", "z"]), df.dropna(subset=["y", "z"]))
    assert_eq(
        ddf.dropna(subset=["y", "z"], how="all"),
        df.dropna(subset=["y", "z"], how="all"),
    )

    # threshold
    assert_eq(ddf.dropna(thresh=None), df.dropna(thresh=None))
    assert_eq(ddf.dropna(thresh=0), df.dropna(thresh=0))
    assert_eq(ddf.dropna(thresh=1), df.dropna(thresh=1))
    assert_eq(ddf.dropna(thresh=2), df.dropna(thresh=2))
    assert_eq(ddf.dropna(thresh=3), df.dropna(thresh=3))

    # fail when how and thresh are both provided
    # see https://github.com/dask/dask/issues/9365
    with pytest.raises(TypeError, match="cannot set both the how and thresh arguments"):
        ddf.dropna(how="all", thresh=0)

    # Regression test for https://github.com/dask/dask/issues/6540
    df = pd.DataFrame({"_0": [0, 0, np.nan], "_1": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.dropna(subset=["_0"]), df.dropna(subset=["_0"]))


@pytest.mark.parametrize("lower, upper", [(2, 5), (2.5, 3.5)])
def test_clip(lower, upper):
    df = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6, 7, 8, 9], "b": [3, 5, 2, 5, 7, 2, 4, 2, 4]}
    )
    ddf = dd.from_pandas(df, 3)

    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9])
    ds = dd.from_pandas(s, 3)

    assert_eq(ddf.clip(lower=lower, upper=upper), df.clip(lower=lower, upper=upper))
    assert_eq(ddf.clip(lower=lower), df.clip(lower=lower))
    assert_eq(ddf.clip(upper=upper), df.clip(upper=upper))

    assert_eq(ds.clip(lower=lower, upper=upper), s.clip(lower=lower, upper=upper))
    assert_eq(ds.clip(lower=lower), s.clip(lower=lower))
    assert_eq(ds.clip(upper=upper), s.clip(upper=upper))


def test_clip_axis_0():
    df = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5, 6, 7, 8, 9],
            "b": [3, 5, 2, 5, 7, 2, 4, 2, 4],
        }
    )
    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9])
    l = pd.Series([3] * len(df))
    u = pd.Series([7] * len(df))

    ddf = dd.from_pandas(df, 3)
    ds = dd.from_pandas(s, 3)
    dl = dd.from_pandas(l, 3)
    du = dd.from_pandas(u, 3)

    assert_eq(ddf.clip(lower=dl, upper=du, axis=0), df.clip(lower=l, upper=u, axis=0))
    assert_eq(ddf.clip(lower=dl, axis=0), df.clip(lower=l, axis=0))
    assert_eq(ddf.clip(upper=du, axis=0), df.clip(upper=u, axis=0))

    assert_eq(ds.clip(lower=dl, upper=du, axis=0), s.clip(lower=l, upper=u, axis=0))
    assert_eq(ds.clip(lower=dl, axis=0), s.clip(lower=l, axis=0))
    assert_eq(ds.clip(upper=du, axis=0), s.clip(upper=u, axis=0))


def test_clip_axis_1():
    df = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5, 6, 7, 8, 9],
            "b": [3, 5, 2, 5, 7, 2, 4, 2, 4],
        }
    )
    ddf = dd.from_pandas(df, 3)
    l = pd.Series({"a": 2, "b": 3})
    u = pd.Series({"a": 7, "b": 5})

    assert_eq(ddf.clip(lower=l, upper=u, axis=1), df.clip(lower=l, upper=u, axis=1))
    assert_eq(ddf.clip(lower=l, axis=1), df.clip(lower=l, axis=1))
    assert_eq(ddf.clip(upper=u, axis=1), df.clip(upper=u, axis=1))

    with pytest.raises(ValueError, match="No axis named 1 for Series"):
        ddf.a.clip(lower=l, upper=u, axis=1)


def test_squeeze():
    df = pd.DataFrame({"x": [1, 3, 6]})
    df2 = pd.DataFrame({"x": [0]})
    s = pd.Series({"test": 0, "b": 100})

    ddf = dd.from_pandas(df, 3)
    ddf2 = dd.from_pandas(df2, 3)
    ds = dd.from_pandas(s, 2)

    assert_eq(df.squeeze(), ddf.squeeze())
    assert_eq(pd.Series([0], name="x"), ddf2.squeeze())
    assert_eq(ds.squeeze(), s.squeeze())

    with pytest.raises(NotImplementedError) as info:
        ddf.squeeze(axis=0)
    msg = f"{type(ddf)} does not support squeeze along axis 0"
    assert msg in str(info.value)

    with pytest.raises(ValueError) as info:
        ddf.squeeze(axis=2)
    msg = f"No axis {2} for object type {type(ddf)}"
    assert msg in str(info.value)

    with pytest.raises(ValueError) as info:
        ddf.squeeze(axis="test")
    msg = f"No axis test for object type {type(ddf)}"
    assert msg in str(info.value)


def test_where_mask():
    pdf1 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6, 7, 8, 9], "b": [3, 5, 2, 5, 7, 2, 4, 2, 4]}
    )
    ddf1 = dd.from_pandas(pdf1, 2)
    pdf2 = pd.DataFrame({"a": [True, False, True] * 3, "b": [False, False, True] * 3})
    ddf2 = dd.from_pandas(pdf2, 2)

    # different index
    pdf3 = pd.DataFrame(
        {"a": [1, 2, 3, 4, 5, 6, 7, 8, 9], "b": [3, 5, 2, 5, 7, 2, 4, 2, 4]},
        index=[0, 1, 2, 3, 4, 5, 6, 7, 8],
    )
    ddf3 = dd.from_pandas(pdf3, 2)
    pdf4 = pd.DataFrame(
        {"a": [True, False, True] * 3, "b": [False, False, True] * 3},
        index=[5, 6, 7, 8, 9, 10, 11, 12, 13],
    )
    ddf4 = dd.from_pandas(pdf4, 2)

    # different columns
    pdf5 = pd.DataFrame(
        {
            "a": [1, 2, 3, 4, 5, 6, 7, 8, 9],
            "b": [9, 4, 2, 6, 2, 3, 1, 6, 2],
            "c": [5, 6, 7, 8, 9, 10, 11, 12, 13],
        },
        index=[0, 1, 2, 3, 4, 5, 6, 7, 8],
    )
    ddf5 = dd.from_pandas(pdf5, 2)
    pdf6 = pd.DataFrame(
        {
            "a": [True, False, True] * 3,
            "b": [False, False, True] * 3,
            "c": [False] * 9,
            "d": [True] * 9,
        },
        index=[5, 6, 7, 8, 9, 10, 11, 12, 13],
    )
    ddf6 = dd.from_pandas(pdf6, 2)

    cases = [
        (ddf1, ddf2, pdf1, pdf2),
        (ddf1.repartition([0, 3, 6, 8]), ddf2, pdf1, pdf2),
        (ddf1, ddf4, pdf3, pdf4),
        (ddf3.repartition([0, 4, 6, 8]), ddf4.repartition([5, 9, 10, 13]), pdf3, pdf4),
        (ddf5, ddf6, pdf5, pdf6),
        (ddf5.repartition([0, 4, 7, 8]), ddf6, pdf5, pdf6),
        # use pd.DataFrame as cond
        (ddf1, pdf2, pdf1, pdf2),
        (ddf1, pdf4, pdf3, pdf4),
        (ddf5, pdf6, pdf5, pdf6),
    ]

    for ddf, ddcond, pdf, pdcond in cases:
        assert isinstance(ddf, dd.DataFrame)
        assert isinstance(ddcond, (dd.DataFrame, pd.DataFrame))
        assert isinstance(pdf, pd.DataFrame)
        assert isinstance(pdcond, pd.DataFrame)

        assert_eq(ddf.where(ddcond), pdf.where(pdcond))
        assert_eq(ddf.mask(ddcond), pdf.mask(pdcond))
        assert_eq(ddf.where(ddcond, -ddf), pdf.where(pdcond, -pdf))
        assert_eq(ddf.mask(ddcond, -ddf), pdf.mask(pdcond, -pdf))

        assert_eq(ddf.where(ddcond.a, -ddf), pdf.where(pdcond.a, -pdf))
        assert_eq(ddf.mask(ddcond.a, -ddf), pdf.mask(pdcond.a, -pdf))
        assert_eq(ddf.a.where(ddcond.a), pdf.a.where(pdcond.a))
        assert_eq(ddf.a.mask(ddcond.a), pdf.a.mask(pdcond.a))
        assert_eq(ddf.a.where(ddcond.a, -ddf.a), pdf.a.where(pdcond.a, -pdf.a))
        assert_eq(ddf.a.mask(ddcond.a, -ddf.a), pdf.a.mask(pdcond.a, -pdf.a))


def test_map_partitions_multi_argument():
    assert_eq(dd.map_partitions(lambda a, b: a + b, d.a, d.b), full.a + full.b)
    assert_eq(
        dd.map_partitions(lambda a, b, c: a + b + c, d.a, d.b, 1), full.a + full.b + 1
    )


def test_map_partitions():
    assert_eq(d.map_partitions(lambda df: df, meta=d), full)
    assert_eq(d.map_partitions(lambda df: df), full)
    result = d.map_partitions(lambda df: df.sum(axis=1))
    assert_eq(result, full.sum(axis=1))

    assert_eq(
        d.map_partitions(lambda df: 1),
        pd.Series([1, 1, 1], dtype=np.int64),
        check_divisions=False,
        check_index=False,
    )


def test_map_partitions_type():
    result = d.map_partitions(type).compute(scheduler="single-threaded")
    assert isinstance(result, pd.Series)
    assert all(x == pd.DataFrame for x in result)


def test_map_partitions_partition_info():
    def f(df, partition_info=None):
        assert partition_info is not None
        assert "number" in partition_info
        assert "division" in partition_info
        assert dsk[("x", partition_info["number"])].equals(df)
        assert dsk[("x", d.divisions.index(partition_info["division"]))].equals(df)
        return df

    df = d.map_partitions(f, meta=d)
    result = df.compute(scheduler="single-threaded")
    assert type(result) == pd.DataFrame


def test_map_partitions_names():
    func = lambda x: x
    assert sorted(dd.map_partitions(func, d, meta=d).dask) == sorted(
        dd.map_partitions(func, d, meta=d).dask
    )

    func = lambda x, y: x
    assert sorted(dd.map_partitions(func, d, d, meta=d).dask) == sorted(
        dd.map_partitions(func, d, d, meta=d).dask
    )


def test_map_partitions_column_info():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    b = dd.map_partitions(lambda x: x, a, meta=a)
    tm.assert_index_equal(b.columns, a.columns)
    assert_eq(df, b)

    b = dd.map_partitions(lambda x: x, a.x, meta=a.x)
    assert b.name == a.x.name
    assert_eq(df.x, b)

    b = dd.map_partitions(lambda x: x, a.x, meta=a.x)
    assert b.name == a.x.name
    assert_eq(df.x, b)

    b = dd.map_partitions(lambda df: df.x + df.y, a)
    assert isinstance(b, dd.Series)
    assert b.dtype == "i8"

    b = dd.map_partitions(lambda df: df.x + 1, a, meta=("x", "i8"))
    assert isinstance(b, dd.Series)
    assert b.name == "x"
    assert b.dtype == "i8"


def test_map_partitions_method_names():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    b = a.map_partitions(lambda x: x)
    assert isinstance(b, dd.DataFrame)
    tm.assert_index_equal(b.columns, a.columns)

    b = a.map_partitions(lambda df: df.x + 1)
    assert isinstance(b, dd.Series)
    assert b.dtype == "i8"

    b = a.map_partitions(lambda df: df.x + 1, meta=("x", "i8"))
    assert isinstance(b, dd.Series)
    assert b.name == "x"
    assert b.dtype == "i8"


def test_map_partitions_propagates_index_metadata():
    index = pd.Series(list("abcde"), name="myindex")
    df = pd.DataFrame(
        {"A": np.arange(5, dtype=np.int32), "B": np.arange(10, 15, dtype=np.int32)},
        index=index,
    )
    ddf = dd.from_pandas(df, npartitions=2)
    res = ddf.map_partitions(
        lambda df: df.assign(C=df.A + df.B),
        meta=[("A", "i4"), ("B", "i4"), ("C", "i4")],
    )
    sol = df.assign(C=df.A + df.B)
    assert_eq(res, sol)

    res = ddf.map_partitions(lambda df: df.rename_axis("newindex"))
    sol = df.rename_axis("newindex")
    assert_eq(res, sol)


@pytest.mark.xfail(reason="now we use SubgraphCallables")
def test_map_partitions_keeps_kwargs_readable():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    def f(s, x=1):
        return s + x

    b = a.x.map_partitions(f, x=5)

    # NOTE: we'd like to ensure that we keep the keyword arguments readable
    # in the dask graph
    assert "['x', 5]" in str(dict(b.dask)) or "{'x': 5}" in str(dict(b.dask))
    assert_eq(df.x + 5, b)

    assert a.x.map_partitions(f, x=5)._name != a.x.map_partitions(f, x=6)._name


def test_map_partitions_with_delayed_collection():
    # https://github.com/dask/dask/issues/5854
    df = pd.DataFrame(columns=list("abcdefghijk"))
    ddf = dd.from_pandas(df, 2)
    ddf.dropna(subset=list("abcdefghijk")).compute()
    # no error!


def test_metadata_inference_single_partition_aligned_args():
    # https://github.com/dask/dask/issues/3034
    # Previously broadcastable series functionality broke this

    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=1)

    def check(df, df_x):
        assert len(df) == len(df_x)
        assert len(df) > 0
        return df

    res = dd.map_partitions(check, ddf, ddf.x)
    assert_eq(res, ddf)


def test_align_dataframes():
    df1 = pd.DataFrame({"A": [1, 2, 3, 3, 2, 3], "B": [1, 2, 3, 4, 5, 6]})
    df2 = pd.DataFrame({"A": [3, 1, 2], "C": [1, 2, 3]})

    ddf1 = dd.from_pandas(df1, npartitions=2)
    ddf2 = dd.from_pandas(df2, npartitions=1)

    actual = ddf1.map_partitions(
        pd.merge, df2, align_dataframes=False, left_on="A", right_on="A", how="left"
    )
    expected = pd.merge(df1, df2, left_on="A", right_on="A", how="left")
    assert_eq(actual, expected, check_index=False, check_divisions=False)

    actual = ddf2.map_partitions(
        pd.merge, ddf1, align_dataframes=False, left_on="A", right_on="A", how="right"
    )
    expected = pd.merge(df2, df1, left_on="A", right_on="A", how="right")
    assert_eq(actual, expected, check_index=False, check_divisions=False)


def test_drop_duplicates(shuffle_method):
    shuffle_method = None
    res = d.drop_duplicates()
    res2 = d.drop_duplicates(split_every=2, shuffle_method=shuffle_method)
    sol = full.drop_duplicates()
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert res._name != res2._name

    res = d.a.drop_duplicates()
    res2 = d.a.drop_duplicates(split_every=2, shuffle_method=shuffle_method)
    sol = full.a.drop_duplicates()
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert res._name != res2._name

    res = d.index.drop_duplicates()
    res2 = d.index.drop_duplicates(split_every=2, shuffle_method=shuffle_method)
    sol = full.index.drop_duplicates()
    # we shuffle in dask-expr and assert_eq doesn't sort indexes
    assert_eq(res.compute().sort_values(), sol)
    assert_eq(res2.compute().sort_values(), sol)

    _d = d.clear_divisions()
    res = _d.index.drop_duplicates()
    res2 = _d.index.drop_duplicates(split_every=2, shuffle_method=shuffle_method)
    sol = full.index.drop_duplicates()
    assert_eq(res.compute().sort_values(), sol)
    assert_eq(res2.compute().sort_values(), sol)
    assert res._name != res2._name

    with pytest.raises(NotImplementedError):
        d.drop_duplicates(keep=False)


def test_drop_duplicates_subset():
    df = pd.DataFrame({"x": [1, 2, 3, 1, 2, 3], "y": ["a", "a", "b", "b", "c", "c"]})
    ddf = dd.from_pandas(df, npartitions=2)

    for kwarg in [{"keep": "first"}, {"keep": "last"}]:
        assert_eq(df.x.drop_duplicates(**kwarg), ddf.x.drop_duplicates(**kwarg))
        for ss in [["x"], "y", ["x", "y"]]:
            assert_eq(
                df.drop_duplicates(subset=ss, **kwarg),
                ddf.drop_duplicates(subset=ss, **kwarg),
            )
            assert_eq(df.drop_duplicates(ss, **kwarg), ddf.drop_duplicates(ss, **kwarg))


def test_get_partition():
    pdf = pd.DataFrame(np.random.randn(10, 5), columns=list("abcde"))
    ddf = dd.from_pandas(pdf, chunksize=4)
    assert ddf.divisions == (0, 4, 8, 9)

    # DataFrame
    div1 = ddf.get_partition(0)
    assert isinstance(div1, dd.DataFrame)
    assert_eq(div1, pdf.loc[0:3])
    div2 = ddf.get_partition(1)
    assert_eq(div2, pdf.loc[4:7])
    div3 = ddf.get_partition(2)
    assert_eq(div3, pdf.loc[8:9])
    assert len(div1) + len(div2) + len(div3) == len(pdf)

    # Series
    div1 = ddf.a.get_partition(0)
    assert isinstance(div1, dd.Series)
    assert_eq(div1, pdf.a.loc[0:3])
    div2 = ddf.a.get_partition(1)
    assert_eq(div2, pdf.a.loc[4:7])
    div3 = ddf.a.get_partition(2)
    assert_eq(div3, pdf.a.loc[8:9])
    assert len(div1) + len(div2) + len(div3) == len(pdf.a)

    with pytest.raises(ValueError):
        ddf.get_partition(-1)

    with pytest.raises(ValueError):
        ddf.get_partition(3)


def test_ndim():
    assert d.ndim == 2
    assert d.a.ndim == 1
    assert d.index.ndim == 1


def test_dtype():
    assert (d.dtypes == full.dtypes).all()


def test_value_counts():
    df = pd.DataFrame({"x": [1, 2, 1, 3, 3, 1, 4]})
    ddf = dd.from_pandas(df, npartitions=3)
    result = ddf.x.value_counts()
    expected = df.x.value_counts()
    assert_eq(result, expected)
    result2 = ddf.x.value_counts(split_every=2)
    assert_eq(result2, expected)
    assert result._name != result2._name


def test_value_counts_not_sorted():
    df = pd.DataFrame({"x": [1, 2, 1, 3, 3, 1, 4]})
    ddf = dd.from_pandas(df, npartitions=3)
    result = ddf.x.value_counts(sort=False)
    expected = df.x.value_counts(sort=False)
    assert_eq(result, expected)
    result2 = ddf.x.value_counts(split_every=2)
    assert_eq(result2, expected)
    assert result._name != result2._name


def test_value_counts_with_dropna():
    df = pd.DataFrame({"x": [1, 2, 1, 3, np.nan, 1, 4]})
    ddf = dd.from_pandas(df, npartitions=3)
    result = ddf.x.value_counts(dropna=False)
    expected = df.x.value_counts(dropna=False)
    assert_eq(result, expected)
    result2 = ddf.x.value_counts(split_every=2, dropna=False)
    assert_eq(result2, expected)
    assert result._name != result2._name


def test_value_counts_with_normalize():
    df = pd.DataFrame({"x": [1, 2, 1, 3, 3, 1, 4]})
    ddf = dd.from_pandas(df, npartitions=3)
    result = ddf.x.value_counts(normalize=True)
    expected = df.x.value_counts(normalize=True)
    assert_eq(result, expected)

    result2 = ddf.x.value_counts(split_every=2, normalize=True)
    assert_eq(result2, expected)
    assert result._name != result2._name

    result3 = ddf.x.value_counts(split_out=2, normalize=True)
    assert_eq(result3, expected)
    assert result._name != result3._name


@pytest.mark.parametrize("normalize", [True, False])
def test_value_counts_with_normalize_and_dropna(normalize):
    df = pd.DataFrame({"x": [1, 2, 1, 3, np.nan, 1, 4]})
    ddf = dd.from_pandas(df, npartitions=3)

    result = ddf.x.value_counts(dropna=False, normalize=normalize)
    expected = df.x.value_counts(dropna=False, normalize=normalize)
    assert_eq(result, expected)

    result2 = ddf.x.value_counts(split_every=2, dropna=False, normalize=normalize)
    assert_eq(result2, expected)
    assert result._name != result2._name

    result3 = ddf.x.value_counts(split_out=2, dropna=False, normalize=normalize)
    assert_eq(result3, expected)
    assert result._name != result3._name

    result4 = ddf.x.value_counts(dropna=True, normalize=normalize, split_out=2)
    expected4 = df.x.value_counts(dropna=True, normalize=normalize)
    assert_eq(result4, expected4)


def test_unique():
    pdf = pd.DataFrame(
        {
            "x": [1, 2, 1, 3, 3, 1, 4, 2, 3, 1],
            "y": ["a", "c", "b", np.nan, "c", "b", "a", "d", np.nan, "a"],
        }
    )
    ddf = dd.from_pandas(pdf, npartitions=3)
    assert_eq(ddf.x.unique(split_out=1), pd.Series(pdf.x.unique(), name="x"))
    assert_eq(ddf.y.unique(split_out=1), pd.Series(pdf.y.unique(), name="y"))

    assert_eq(
        ddf.x.unique(split_every=2, split_out=1), pd.Series(pdf.x.unique(), name="x")
    )
    assert_eq(
        ddf.y.unique(split_every=2, split_out=1), pd.Series(pdf.y.unique(), name="y")
    )
    assert_eq(ddf.index.unique(split_out=1), pdf.index.unique())

    assert ddf.x.unique(split_every=2)._name != ddf.x.unique()._name


def test_isin():
    f_list = [1, 2, 3]
    f_series = pd.Series(f_list)
    f_dict = {"a": [0, 3], "b": [1, 2]}
    f_list2 = [1, "2"]
    f_list_delayed = [delayed(1), delayed(2), delayed(3)]
    f_list_of_lists = [[1, 2, 3], [2, 3, 4], [3, 4, 5]]

    # Series
    assert_eq(d.a.isin(f_list), full.a.isin(f_list))
    assert_eq(d.a.isin(f_series), full.a.isin(f_series))
    assert_eq(d.a.isin(f_list2), full.a.isin(f_list2))
    assert_eq(
        d.a.isin(f_list_delayed),
        full.a.isin(f_list),
    )
    assert_eq(
        d.a.isin(f_list_of_lists),
        full.a.isin(f_list_of_lists),
    )
    with pytest.raises(NotImplementedError):
        d.a.isin(d.a)

    # DataFrame test
    assert_eq(d.isin(f_list), full.isin(f_list))
    assert_eq(d.isin(f_dict), full.isin(f_dict))
    assert_eq(d.isin(f_list2), full.isin(f_list2))
    assert_eq(
        d.isin(f_list_delayed),
        full.isin(f_list),
    )
    assert_eq(
        d.isin(f_list_of_lists),
        full.isin(f_list_of_lists),
    )
    for obj in [d, f_series, full]:
        with pytest.raises(NotImplementedError):
            d.isin(obj)


def test_contains_frame():
    df = dd.from_pandas(pd.DataFrame({"A": [1, 2], 0: [3, 4]}), 1)
    assert "A" in df
    assert 0 in df
    assert "B" not in df
    assert 1 not in df


def test_len():
    assert len(d) == len(full)
    assert len(d.a) == len(full.a)
    assert len(dd.from_pandas(pd.DataFrame(), npartitions=1)) == 0
    assert len(dd.from_pandas(pd.DataFrame(columns=[1, 2]), npartitions=1)) == 0
    # Regression test for https://github.com/dask/dask/issues/6110
    assert len(dd.from_pandas(pd.DataFrame(columns=["foo", "foo"]), npartitions=1)) == 0


def test_size():
    assert_eq(d.size, full.size)
    assert_eq(d.a.size, full.a.size)
    assert_eq(d.index.size, full.index.size)


def test_shape():
    result = d.shape
    assert_eq((result[0].compute(), result[1]), (len(full), len(full.columns)))
    assert_eq(dd.compute(result)[0], (len(full), len(full.columns)))

    result = d.a.shape
    assert_eq(result[0].compute(), len(full.a))
    assert_eq(dd.compute(result)[0], (len(full.a),))

    sh = dd.from_pandas(pd.DataFrame(index=[1, 2, 3]), npartitions=2).shape
    assert (sh[0].compute(), sh[1]) == (3, 0)
    sh = dd.from_pandas(pd.DataFrame({"a": [], "b": []}, index=[]), npartitions=1).shape
    assert (sh[0].compute(), sh[1]) == (0, 2)


def test_nbytes():
    assert_eq(d.a.nbytes, full.a.nbytes)
    assert_eq(d.index.nbytes, full.index.nbytes)


@pytest.mark.parametrize(
    "method",
    [
        pytest.param(
            "tdigest",
            marks=pytest.mark.skipif(not crick, reason="Requires crick"),
        ),
        "dask",
    ],
)
@pytest.mark.parametrize("quantile", (0.3, 0.5, 0.9))
def test_quantile(method, quantile):
    # https://en.wikipedia.org/wiki/Exponential_distribution
    array = da.random.exponential(1, 10_000, chunks=(10_000 // 2))
    df = dd.from_dask_array(array, columns=["x"])
    exp = -np.log(1 - quantile) * 1

    # dataframe
    result = df.x.quantile([quantile], method=method)
    assert len(result) == 1
    assert result.divisions == (quantile, quantile)
    assert isinstance(result, dd.Series)
    result = result.compute()
    assert isinstance(result, pd.Series)
    assert result.iloc[0] == pytest.approx(exp, rel=0.15)

    # series / single
    result = df.x.quantile(quantile, method=method)
    assert result.ndim == 0
    result = result.compute()
    assert result == pytest.approx(exp, rel=0.15)


@pytest.mark.parametrize(
    "method",
    [
        pytest.param(
            "tdigest", marks=pytest.mark.skipif(not crick, reason="Requires crick")
        ),
        "dask",
    ],
)
def test_quantile_missing(method):
    df = pd.DataFrame({"A": [0, np.nan, 2]})
    # TODO: Test npartitions=2
    # (see https://github.com/dask/dask/issues/9227)
    ddf = dd.from_pandas(df, npartitions=1)
    expected = df.quantile()
    result = ddf.quantile(method=method)
    assert_eq(result, expected)

    expected = df.A.quantile()
    result = ddf.A.quantile(method=method)
    assert_eq(result, expected)


@pytest.mark.parametrize(
    "method",
    [
        pytest.param(
            "tdigest", marks=pytest.mark.skipif(not crick, reason="Requires crick")
        ),
        "dask",
    ],
)
def test_empty_quantile(method):
    with pytest.raises(AssertionError, match="must provide non-"):
        d.b.quantile([], method=method).compute()


# TODO: un-filter once https://github.com/dask/dask/issues/8960 is resolved.
@pytest.mark.filterwarnings(
    "ignore:In future versions of pandas, numeric_only will be set to False:FutureWarning"
)
@pytest.mark.parametrize(
    "method,expected",
    [
        pytest.param(
            "tdigest",
            (
                pd.Series([9.5, 29.5, 19.5], index=["A", "X", "B"]),
                pd.DataFrame(
                    [[4.5, 24.5, 14.5], [14.5, 34.5, 24.5]],
                    index=[0.25, 0.75],
                    columns=["A", "X", "B"],
                ),
            ),
            marks=pytest.mark.skipif(not crick, reason="Requires crick"),
        ),
        (
            "dask",
            (
                pd.Series([7.0, 27.0, 17.0], index=["A", "X", "B"]),
                pd.DataFrame(
                    [[1.50, 21.50, 11.50], [14.0, 34.0, 24.0]],
                    index=[0.25, 0.75],
                    columns=["A", "X", "B"],
                ),
            ),
        ),
    ],
)
@pytest.mark.parametrize("numeric_only", [None, True, False])
def test_dataframe_quantile(method, expected, numeric_only):
    # column X is for test column order and result division
    df = pd.DataFrame(
        {
            "A": np.arange(20),
            "X": np.arange(20, 40),
            "B": np.arange(10, 30),
            "C": ["a", "b", "c", "d"] * 5,
        },
        columns=["A", "X", "B", "C"],
    )
    ddf = dd.from_pandas(df, 3)

    numeric_only_kwarg = {}
    if numeric_only is not None:
        numeric_only_kwarg = {"numeric_only": numeric_only}

    if numeric_only is False or numeric_only is None:
        with pytest.raises(TypeError):
            df.quantile(**numeric_only_kwarg)
        with pytest.raises(
            (TypeError, ArrowNotImplementedError, ValueError),
            match="unsupported operand|no kernel|non-numeric|not supported",
        ):
            ddf.quantile(**numeric_only_kwarg)
    else:
        result = ddf.quantile(method=method, **numeric_only_kwarg)
        assert result.npartitions == 1
        assert result.divisions == ("A", "X")

        result = result.compute()
        assert isinstance(result, pd.Series)
        assert result.name == 0.5
        assert_eq(result, expected[0], check_names=False)

        result = ddf.quantile([0.25, 0.75], method=method, **numeric_only_kwarg)
        assert result.npartitions == 1
        assert result.divisions == (0.25, 0.75)

        result = result.compute()
        assert isinstance(result, pd.DataFrame)
        tm.assert_index_equal(result.index, pd.Index([0.25, 0.75]))
        tm.assert_index_equal(result.columns, pd.Index(["A", "X", "B"]))

        assert (result == expected[1]).all().all()

        with warnings.catch_warnings(record=True):
            warnings.filterwarnings("ignore", category=FutureWarning)
            # pandas issues a warning with 1.5, but not 1.3
            expected = df.quantile(axis=1, **numeric_only_kwarg)

        result = ddf.quantile(axis=1, method=method, **numeric_only_kwarg)

        assert_eq(result, expected)

        with pytest.raises(ValueError):
            ddf.quantile([0.25, 0.75], axis=1, method=method, **numeric_only_kwarg)


def test_quantile_datetime_numeric_only_false():
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "dt": [pd.NaT] + [datetime(2011, i, 1) for i in range(1, 8)],
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
        }
    )
    ddf = dd.from_pandas(df, 1)

    assert_eq(ddf.quantile(numeric_only=False), df.quantile(numeric_only=False))


def test_quantile_for_possibly_unsorted_q():
    """check that quantile is giving correct answers even when quantile parameter, q, may be unsorted.

    See https://github.com/dask/dask/issues/4642.
    """
    # prepare test case where percentiles should equal values
    A = da.arange(0, 101)
    ds = dd.from_dask_array(A)

    for q in [
        [0.25, 0.50, 0.75],
        [0.25, 0.50, 0.75, 0.99],
        [0.75, 0.5, 0.25],
        [0.25, 0.99, 0.75, 0.50],
    ]:
        r = ds.quantile(q).compute()
        assert_eq(r.loc[0.25], 25.0)
        assert_eq(r.loc[0.50], 50.0)
        assert_eq(r.loc[0.75], 75.0)

    r = ds.quantile([0.25]).compute()
    assert_eq(r.loc[0.25], 25.0)

    r = ds.quantile(0.25).compute()
    assert_eq(r, 25.0)


def test_quantile_tiny_partitions():
    """See https://github.com/dask/dask/issues/6551"""
    df = pd.DataFrame({"a": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=3)
    r = ddf["a"].quantile(0.5).compute()
    assert r == 2


def test_quantile_trivial_partitions():
    """See https://github.com/dask/dask/issues/2792"""
    df = pd.DataFrame({"A": []})
    ddf = dd.from_pandas(df, npartitions=2)
    expected = df.quantile(0.5)
    assert_eq(ddf.quantile(0.5), expected)

    df = pd.DataFrame({"A": [np.nan, np.nan, np.nan, np.nan]})
    ddf = dd.from_pandas(df, npartitions=2)
    expected = df.quantile(0.5)
    assert_eq(ddf.quantile(0.5), expected)


def test_index():
    assert_eq(d.index, full.index)


def test_assign():
    df = pd.DataFrame(
        {"a": range(8), "b": [float(i) for i in range(10, 18)]},
        index=pd.Index(list("abcdefgh")),
    )
    ddf = dd.from_pandas(df, npartitions=3)
    ddf_unknown = dd.from_pandas(df, npartitions=3, sort=False).clear_divisions()
    assert not ddf_unknown.known_divisions

    res = ddf.assign(
        c=1,
        d="string",
        e=ddf.a.sum(),
        f=ddf.a + ddf.b,
        g=lambda x: x.a + x.c,
        dt=pd.Timestamp(2018, 2, 13),
    )
    res_unknown = ddf_unknown.assign(
        c=1,
        d="string",
        e=ddf_unknown.a.sum(),
        f=ddf_unknown.a + ddf_unknown.b,
        g=lambda x: x.a + x.c,
        dt=pd.Timestamp(2018, 2, 13),
    )
    sol = df.assign(
        c=1,
        d="string",
        e=df.a.sum(),
        f=df.a + df.b,
        g=lambda x: x.a + x.c,
        dt=pd.Timestamp(2018, 2, 13),
    )
    assert_eq(res, sol)
    assert_eq(res_unknown, sol)

    res = ddf.assign(c=df.a + 1)
    assert_eq(res, df.assign(c=df.a + 1))

    res = ddf.assign(c=ddf.index)
    assert_eq(res, df.assign(c=df.index))

    assert_eq(ddf_unknown.assign(c=df.a + 1), df.assign(c=df.a + 1))

    # unsupported type
    with pytest.raises(TypeError):
        ddf.assign(c=list(range(9)))

    assert_eq(ddf_unknown.assign(foo=ddf.a), df.assign(foo=df.a))
    assert_eq(ddf.assign(foo=ddf_unknown.a), df.assign(foo=df.a))

    df = pd.DataFrame({"A": [1, 2]})
    df.assign(B=lambda df: df["A"], C=lambda df: df.A + df.B)

    ddf = dd.from_pandas(pd.DataFrame({"A": [1, 2]}), npartitions=2)
    ddf.assign(B=lambda df: df["A"], C=lambda df: df.A + df.B)

    assert_eq(df, ddf)


def test_assign_callable():
    df = dd.from_pandas(pd.DataFrame({"A": range(10)}), npartitions=2)
    a = df.assign(B=df.A.shift())
    b = df.assign(B=lambda x: x.A.shift())
    assert_eq(a, b)


def test_assign_dtypes():
    ddf = dd.from_pandas(
        pd.DataFrame(
            data={"col1": ["a", "b"], "col2": [1, 2]}, columns=["col1", "col2"]
        ),
        npartitions=2,
    )

    new_col = {"col3": pd.Series(["0", "1"])}
    res = ddf.assign(**new_col)

    string_dtype = get_string_dtype()
    assert_eq(
        res.dtypes,
        pd.Series(
            data=[string_dtype, "int64", string_dtype],
            index=["col1", "col2", "col3"],
        ),
    )


def test_assign_pandas_series():
    # Make sure we handle when `df.columns.min()` raises TypeError
    df = pd.DataFrame({"a": [1, 2], 1: [5, 6]})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf = ddf.assign(c=df["a"])
    assert_eq(ddf, df.assign(c=df["a"]))


def test_map():
    df = pd.DataFrame(
        {"a": range(9), "b": [4, 5, 6, 1, 2, 3, 0, 0, 0]},
        index=pd.Index([0, 1, 3, 5, 6, 8, 9, 9, 9], name="myindex"),
    )
    ddf = dd.from_pandas(df, npartitions=3)
    with pytest.warns(UserWarning, match="meta"):
        assert_eq(ddf.a.map(lambda x: x + 1), df.a.map(lambda x: x + 1))
        lk = {v: v + 1 for v in df.a.values}
        assert_eq(ddf.a.map(lk), df.a.map(lk))
        assert_eq(ddf.b.map(lk), df.b.map(lk))
        lk = pd.Series(lk)
        assert_eq(ddf.a.map(lk), df.a.map(lk))
        assert_eq(ddf.b.map(lk), df.b.map(lk))
    assert_eq(ddf.b.map(lk, meta=ddf.b), df.b.map(lk))
    assert_eq(ddf.b.map(lk, meta=("b", "i8")), df.b.map(lk))


def test_concat():
    x = _concat([pd.DataFrame(columns=["a", "b"]), pd.DataFrame(columns=["a", "b"])])
    assert list(x.columns) == ["a", "b"]
    assert len(x) == 0


def test_known_divisions():
    assert d.known_divisions
    df = d.clear_divisions()
    assert not df.known_divisions


def test_unknown_divisions():
    dsk = {
        ("x", 0): pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}),
        ("x", 1): pd.DataFrame({"a": [4, 5, 6], "b": [3, 2, 1]}),
        ("x", 2): pd.DataFrame({"a": [7, 8, 9], "b": [0, 0, 0]}),
    }
    d = dd.repartition(
        pd.concat(dsk.values()), divisions=[0, 1, 2, 10]
    ).clear_divisions()
    full = d.compute(scheduler="sync")

    assert_eq(d.a.sum(), full.a.sum())
    assert_eq(d.a + d.b + 1, full.a + full.b + 1)


def test_with_min_count():
    dfs = [
        pd.DataFrame([[None, 2, 3], [None, 5, 6], [5, 4, 9]]),
        pd.DataFrame([[2, None, None], [None, 5, 6], [5, 4, 9]]),
    ]
    ddfs = [dd.from_pandas(df, npartitions=4) for df in dfs]
    axes = [0, 1]

    for df, ddf in zip(dfs, ddfs):
        for axis in axes:
            for min_count in [0, 1, 2, 3]:
                assert_eq(
                    df.sum(min_count=min_count, axis=axis),
                    ddf.sum(min_count=min_count, axis=axis),
                )
                assert_eq(
                    df.prod(min_count=min_count, axis=axis),
                    ddf.prod(min_count=min_count, axis=axis),
                )


@pytest.mark.parametrize("join", ["inner", "outer", "left", "right"])
def test_align(join):
    df1a = pd.DataFrame(
        {"A": np.random.randn(10), "B": np.random.randn(10)},
        index=[1, 12, 5, 6, 3, 9, 10, 4, 13, 11],
    )

    df1b = pd.DataFrame(
        {"A": np.random.randn(10), "B": np.random.randn(10)},
        index=[0, 3, 2, 10, 5, 6, 7, 8, 12, 13],
    )
    ddf1a = dd.from_pandas(df1a, 3)
    ddf1b = dd.from_pandas(df1b, 3)

    # DataFrame
    res1, res2 = ddf1a.align(ddf1b, join=join)
    exp1, exp2 = df1a.align(df1b, join=join)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # Series
    res1, res2 = ddf1a["A"].align(ddf1b["B"], join=join)
    exp1, exp2 = df1a["A"].align(df1b["B"], join=join)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # DataFrame with fill_value
    res1, res2 = ddf1a.align(ddf1b, join=join, fill_value=1)
    exp1, exp2 = df1a.align(df1b, join=join, fill_value=1)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # Series
    res1, res2 = ddf1a["A"].align(ddf1b["B"], join=join, fill_value=1)
    exp1, exp2 = df1a["A"].align(df1b["B"], join=join, fill_value=1)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)


@pytest.mark.parametrize("join", ["inner", "outer", "left", "right"])
def test_align_axis(join):
    df1a = pd.DataFrame(
        {"A": np.random.randn(10), "B": np.random.randn(10), "C": np.random.randn(10)},
        index=[1, 12, 5, 6, 3, 9, 10, 4, 13, 11],
    )

    df1b = pd.DataFrame(
        {"B": np.random.randn(10), "C": np.random.randn(10), "D": np.random.randn(10)},
        index=[0, 3, 2, 10, 5, 6, 7, 8, 12, 13],
    )
    ddf1a = dd.from_pandas(df1a, 3)
    ddf1b = dd.from_pandas(df1b, 3)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis=0)
    exp1, exp2 = df1a.align(df1b, join=join, axis=0)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis=1)
    exp1, exp2 = df1a.align(df1b, join=join, axis=1)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis="index")
    exp1, exp2 = df1a.align(df1b, join=join, axis="index")
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis="columns")
    exp1, exp2 = df1a.align(df1b, join=join, axis="columns")
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # invalid
    with pytest.raises(ValueError):
        ddf1a.align(ddf1b, join=join, axis="XXX")

    with pytest.raises(ValueError):
        ddf1a["A"].align(ddf1b["B"], join=join, axis=1)


def test_combine():
    df1 = pd.DataFrame(
        {
            "A": np.random.choice([1, 2, np.nan], 100),
            "B": np.random.choice(["a", "b", "nan"], 100),
        }
    )

    df2 = pd.DataFrame(
        {
            "A": np.random.choice([1, 2, 3], 100),
            "B": np.random.choice(["a", "b", "c"], 100),
        }
    )
    ddf1 = dd.from_pandas(df1, 4)
    ddf2 = dd.from_pandas(df2, 5)

    first = lambda a, b: a

    # You can add series with strings and nans but you can't add scalars 'a' + np.nan
    str_add = lambda a, b: a + b if a is not np.nan else a

    # DataFrame
    for dda, ddb, a, b, runs in [
        (ddf1, ddf2, df1, df2, [(add, None), (first, None)]),
        (ddf1.A, ddf2.A, df1.A, df2.A, [(add, None), (add, 100), (first, None)]),
        (
            ddf1.B,
            ddf2.B,
            df1.B,
            df2.B,
            [(str_add, None), (str_add, "d"), (first, None)],
        ),
    ]:
        for func, fill_value in runs:
            sol = a.combine(b, func, fill_value=fill_value)
            assert_eq(dda.combine(ddb, func, fill_value=fill_value), sol)
            assert_eq(dda.combine(b, func, fill_value=fill_value), sol)

    assert_eq(
        ddf1.combine(ddf2, add, overwrite=False), df1.combine(df2, add, overwrite=False)
    )
    assert dda.combine(ddb, add)._name == dda.combine(ddb, add)._name


def test_combine_first():
    df1 = pd.DataFrame(
        {
            "A": np.random.choice([1.0, 2.0, np.nan], 100),
            "B": np.random.choice(["a", "b", "nan"], 100),
        }
    )

    df2 = pd.DataFrame(
        {
            "A": np.random.choice([1.0, 2.0, 3.0], 100),
            "B": np.random.choice(["a", "b", "c"], 100),
        }
    )
    ddf1 = dd.from_pandas(df1, 4)
    ddf2 = dd.from_pandas(df2, 5)

    # DataFrame
    assert_eq(ddf1.combine_first(ddf2), df1.combine_first(df2))
    assert_eq(ddf1.combine_first(df2), df1.combine_first(df2))

    # Series
    assert_eq(ddf1.A.combine_first(ddf2.A), df1.A.combine_first(df2.A))
    assert_eq(ddf1.A.combine_first(df2.A), df1.A.combine_first(df2.A))

    assert_eq(ddf1.B.combine_first(ddf2.B), df1.B.combine_first(df2.B))
    assert_eq(ddf1.B.combine_first(df2.B), df1.B.combine_first(df2.B))


# This will likely break when the deprecation is enacted in Pandas 3.0
@pytest.mark.xfail(PANDAS_GE_210, reason="https://github.com/dask/dask/issues/10931")
@pytest.mark.parametrize(
    "dtype_lhs,dtype_rhs",
    [("f8", "i8"), ("f8", "f4"), ("datetime64[s]", "datetime64[ns]")],
)
def test_combine_first_all_nans(dtype_lhs, dtype_rhs):
    """If you call s1.combine_first(s2), where s1 is pandas.Series of all NaNs and s2 is
    a pandas.Series of non-floats, the dtype becomes that of s2. Starting with pandas
    2.1, this comes with a deprecation warning.
    Test behaviour when either a whole dask series or just a chunk is full of NaNs.
    """
    if PANDAS_GE_210:
        ctx = pytest.warns(
            FutureWarning,
            match="The behavior of array concatenation with empty entries is deprecated",
        )
    else:
        ctx = contextlib.nullcontext()

    s1 = pd.Series([np.nan, np.nan], dtype=dtype_lhs)
    s2 = pd.Series([np.nan, 1.0], dtype=dtype_lhs)
    s3 = pd.Series([1, 2], dtype=dtype_rhs)
    ds1 = dd.from_pandas(s1, npartitions=2)
    ds2 = dd.from_pandas(s2, npartitions=2)
    ds3 = dd.from_pandas(s3, npartitions=2)

    # Both dask and pandas use output dtype from rhs, with a FutureWarning
    with ctx:
        s13 = s1.combine_first(s3)
    with ctx:
        assert_eq(ds1.combine_first(ds3), s13)

    # No deprecation; output dtype is from lhs
    assert_eq(ds2.combine_first(ds3), s2.combine_first(s3))


def test_dataframe_picklable():
    from pickle import dumps, loads

    from cloudpickle import dumps as cp_dumps
    from cloudpickle import loads as cp_loads

    d = _compat.makeTimeDataFrame()
    df = dd.from_pandas(d, npartitions=3)
    df = df + 2

    # dataframe
    df2 = loads(dumps(df))
    assert_eq(df, df2)
    df2 = cp_loads(cp_dumps(df))
    assert_eq(df, df2)

    # series
    a2 = loads(dumps(df.A))
    assert_eq(df.A, a2)
    a2 = cp_loads(cp_dumps(df.A))
    assert_eq(df.A, a2)

    # index
    i2 = loads(dumps(df.index))
    assert_eq(df.index, i2)
    i2 = cp_loads(cp_dumps(df.index))
    assert_eq(df.index, i2)

    # scalar
    # lambdas are present, so only test cloudpickle
    s = df.A.sum()
    s2 = cp_loads(cp_dumps(s))
    assert_eq(s, s2)


def test_random_partitions():
    a, b = d.random_split([0.5, 0.5], 42)
    assert isinstance(a, dd.DataFrame)
    assert isinstance(b, dd.DataFrame)
    assert a._name != b._name
    np.testing.assert_array_equal(a.index, sorted(a.index))

    assert len(a.compute()) + len(b.compute()) == len(full)
    a2, b2 = d.random_split([0.5, 0.5], 42)
    assert a2._name == a._name
    assert b2._name == b._name

    a, b = d.random_split([0.5, 0.5], 42, True)
    a2, b2 = d.random_split([0.5, 0.5], 42, True)
    assert_eq(a, a2)
    assert_eq(b, b2)
    with pytest.raises(AssertionError):
        np.testing.assert_array_equal(a.index, sorted(a.index))

    parts = d.random_split([0.4, 0.5, 0.1], 42)
    names = {p._name for p in parts}
    names.update([a._name, b._name])
    assert len(names) == 5

    with pytest.raises(ValueError):
        d.random_split([0.4, 0.5], 42)


def test_series_round():
    ps = pd.Series([1.123, 2.123, 3.123, 1.234, 2.234, 3.234], name="a")
    s = dd.from_pandas(ps, npartitions=3)
    assert_eq(s.round(), ps.round())


def test_repartition():

    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5, 6], "y": list("abdabd")}, index=[10, 20, 30, 40, 50, 60]
    )
    a = dd.from_pandas(df, 2)

    b = a.repartition(divisions=[10, 20, 50, 60])
    assert b.divisions == (10, 20, 50, 60)
    assert_eq(a, b)

    for div in [
        [20, 60],
        [10, 50],
        [1],  # first / last element mismatch
        [0, 60],
        [10, 70],  # do not allow to expand divisions by default
        [10, 50, 20, 60],  # not sorted
        [10, 10, 20, 60],
    ]:  # not unique (last element can be duplicated)
        pytest.raises(
            ValueError, lambda div=div: a.repartition(divisions=div).compute()
        )

    pdf = pd.DataFrame(np.random.randn(7, 5), columns=list("abxyz"))
    ps = pdf.x
    for p in range(1, 7):
        ddf = dd.from_pandas(pdf, p)
        ds = ddf.x
        assert_eq(ddf, pdf)
        assert_eq(ps, ds)
        for div in [
            [0, 6],
            [0, 6, 6],
            [0, 5, 6],
            [0, 4, 6, 6],
            [0, 2, 6],
            [0, 2, 6, 6],
            [0, 2, 3, 6, 6],
            [0, 1, 2, 3, 4, 5, 6, 6],
        ]:
            rddf = ddf.repartition(divisions=div)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ds.repartition(divisions=div)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)

        # expand divisions
        for div in [[-5, 10], [-2, 3, 5, 6], [0, 4, 5, 9, 10]]:
            rddf = ddf.repartition(divisions=div, force=True)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ds.repartition(divisions=div, force=True)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)

    pdf = pd.DataFrame(
        {"x": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], "y": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]},
        index=list("abcdefghij"),
    )
    ps = pdf.x
    for p in range(1, 7):
        ddf = dd.from_pandas(pdf, p)
        ds = ddf.x
        assert_eq(ddf, pdf)
        assert_eq(ps, ds)
        for div in [
            list("aj"),
            list("ajj"),
            list("adj"),
            list("abfj"),
            list("ahjj"),
            list("acdj"),
            list("adfij"),
            list("abdefgij"),
            list("abcdefghij"),
        ]:
            rddf = ddf.repartition(divisions=div)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ds.repartition(divisions=div)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)

        # expand divisions
        for div in [list("Yadijm"), list("acmrxz"), list("Yajz")]:
            rddf = ddf.repartition(divisions=div, force=True)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ds.repartition(divisions=div, force=True)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)


def test_repartition_on_pandas_dataframe():
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5, 6], "y": list("abdabd")}, index=[10, 20, 30, 40, 50, 60]
    )
    ddf = dd.repartition(df, divisions=[10, 20, 50, 60])
    assert isinstance(ddf, dd.DataFrame)
    assert ddf.divisions == (10, 20, 50, 60)
    assert_eq(ddf, df)

    ddf = dd.repartition(df.y, divisions=[10, 20, 50, 60])
    assert isinstance(ddf, dd.Series)
    assert ddf.divisions == (10, 20, 50, 60)
    assert_eq(ddf, df.y)


@pytest.mark.parametrize("use_index", [True, False])
@pytest.mark.parametrize("n", [1, 2, 4, 5])
@pytest.mark.parametrize("k", [1, 2, 4, 5])
@pytest.mark.parametrize("dtype", [float, "M8[ns]"])
@pytest.mark.parametrize("transform", [lambda df: df, lambda df: df.x])
def test_repartition_npartitions(use_index, n, k, dtype, transform):
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5, 6] * 10, "y": list("abdabd") * 10},
        index=pd.Series(list(range(0, 30)) * 2, dtype=dtype),
    )
    df = transform(df)
    a = dd.from_pandas(df, npartitions=n, sort=use_index)
    b = a.repartition(npartitions=k)
    assert_eq(a, b)
    assert b.npartitions == k
    assert all(map(len, b.partitions))


@pytest.mark.parametrize("use_index", [True, False])
@pytest.mark.parametrize("n", [2, 5])
@pytest.mark.parametrize("partition_size", ["1kiB", 379])
@pytest.mark.parametrize("transform", [lambda df: df, lambda df: df.x])
def test_repartition_partition_size(use_index, n, partition_size, transform):
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5, 6] * 10, "y": list("abdabd") * 10},
        index=pd.Series([10, 20, 30, 40, 50, 60] * 10),
    )
    df = transform(df)
    a = dd.from_pandas(df, npartitions=n, sort=use_index)
    b = a.repartition(partition_size=partition_size)
    assert_eq(a, b, check_divisions=False)
    assert np.all(b.map_partitions(total_mem_usage, deep=True).compute() <= 1024)
    parts = dask.get(b.dask, b.__dask_keys__())
    assert all(map(len, parts))


def test_repartition_partition_size_arg():
    df = pd.DataFrame({"x": range(10)})
    a = dd.from_pandas(df, npartitions=2)
    b = a.repartition(partition_size="1 MiB")
    assert b.npartitions == 1


def test_repartition_npartitions_same_limits():
    df = pd.DataFrame(
        {"x": [1, 2, 3]},
        index=[
            pd.Timestamp("2017-05-09 00:00:00.006000"),
            pd.Timestamp("2017-05-09 02:45:00.017999"),
            pd.Timestamp("2017-05-09 05:59:58.938999"),
        ],
    )

    ddf = dd.from_pandas(df, npartitions=2)

    ddf.repartition(npartitions=10)


def test_repartition_npartitions_numeric_edge_case():
    """
    Test that we cover numeric edge cases when
    int(ddf.npartitions / npartitions) * npartitions) != ddf.npartitions
    """
    df = pd.DataFrame({"x": range(100)})
    a = dd.from_pandas(df, npartitions=15)
    assert a.npartitions == 15
    b = a.repartition(npartitions=11)
    assert_eq(a, b)


def test_repartition_object_index():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6] * 10}, index=list("abdabd") * 10)
    a = dd.from_pandas(df, npartitions=5)
    b = a.repartition(npartitions=2)
    assert b.npartitions == 2
    assert_eq(b, df)

    b = a.repartition(npartitions=10)
    assert b.npartitions == 10
    assert_eq(b, df)
    assert not b.known_divisions


def test_repartition_datetime_tz_index():
    # Regression test for https://github.com/dask/dask/issues/8788
    # Use TZ-aware datetime index
    s = pd.Series(range(10))
    s.index = pd.to_datetime(
        [datetime(2020, 1, 1, 12, 0) + timedelta(minutes=x) for x in s], utc=True
    )
    ds = dd.from_pandas(s, npartitions=2)
    assert ds.npartitions == 2
    assert_eq(s, ds)
    result = ds.repartition(npartitions=5)
    assert result.npartitions == 5
    assert_eq(s, result)


@pytest.mark.slow
@pytest.mark.parametrize("npartitions", [1, 20, 243])
@pytest.mark.parametrize("freq", ["1D", "7D", "28h", "1h"])
@pytest.mark.parametrize(
    "end", ["2000-04-15", "2000-04-15 12:37:01", "2000-01-01 12:37:00"]
)
@pytest.mark.parametrize(
    "start", ["2000-01-01", "2000-01-01 12:30:00", "2000-01-01 12:30:00"]
)
def test_repartition_freq(npartitions, freq, start, end):
    start = pd.Timestamp(start)
    end = pd.Timestamp(end)
    ind = pd.date_range(start=start, end=end, freq="60s")
    df = pd.DataFrame({"x": np.arange(len(ind))}, index=ind)
    ddf = dd.from_pandas(df, npartitions=npartitions)

    ddf2 = ddf.repartition(freq=freq)
    assert_eq(ddf2, df)


def test_repartition_freq_divisions():
    df = pd.DataFrame(
        {"x": np.random.random(10)},
        index=pd.DatetimeIndex(np.random.random(10) * 100e9),
    )
    ddf = dd.from_pandas(df, npartitions=3)

    ddf2 = ddf.repartition(freq="15s")
    for div in ddf2.divisions[1:-1]:
        assert div == div.round("15s")
    assert ddf2.divisions[0] == df.index.min()
    assert ddf2.divisions[-1] == df.index.max()
    assert_eq(ddf2, df)


def test_repartition_freq_errors():
    df = pd.DataFrame({"x": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=1)
    with pytest.raises(TypeError) as info:
        ddf.repartition(freq="1s")

    assert "only" in str(info.value)
    assert "timeseries" in str(info.value)


def test_repartition_freq_month():
    ts = pd.date_range("2015-01-01 00:00", "2015-05-01 23:50", freq="10min")
    df = pd.DataFrame(
        np.random.randint(0, 100, size=(len(ts), 4)), columns=list("ABCD"), index=ts
    )
    ddf = dd.from_pandas(df, npartitions=1).repartition(freq="MS")

    assert_eq(df, ddf)

    assert ddf.divisions == (
        pd.Timestamp("2015-1-1 00:00:00"),
        pd.Timestamp("2015-2-1 00:00:00"),
        pd.Timestamp("2015-3-1 00:00:00"),
        pd.Timestamp("2015-4-1 00:00:00"),
        pd.Timestamp("2015-5-1 00:00:00"),
        pd.Timestamp("2015-5-1 23:50:00"),
    )

    assert ddf.npartitions == 5


def test_repartition_freq_day():
    index = [
        pd.Timestamp("2020-1-1"),
        pd.Timestamp("2020-1-1"),
        pd.Timestamp("2020-1-2"),
        pd.Timestamp("2020-1-2"),
    ]
    pdf = pd.DataFrame(index=index, data={"foo": "foo"})
    ddf = dd.from_pandas(pdf, npartitions=1).repartition(freq="D")
    assert_eq(ddf, pdf)
    assert ddf.npartitions == 2
    assert ddf.divisions == (
        pd.Timestamp("2020-1-1"),
        pd.Timestamp("2020-1-2"),
        pd.Timestamp("2020-1-2"),
    )


@pytest.mark.parametrize(
    "freq, expected_freq",
    [
        ("M", "MS"),
        ("ME" if PANDAS_GE_220 else "M", "MS"),
        ("MS", "MS"),
        ("2M", "2MS"),
        ("Q", "QS"),
        ("Q-FEB", "QS-FEB"),
        ("2Q", "2QS"),
        ("2Q-FEB", "2QS-FEB"),
        ("2QS-FEB", "2QS-FEB"),
        ("BQ", "BQS"),
        ("2BQ", "2BQS"),
        ("SM", "SMS"),
        ("A", "YS" if PANDAS_GE_220 else "AS"),
        ("Y", "YS" if PANDAS_GE_220 else "AS"),
        ("A-JUN", "YS-JUN" if PANDAS_GE_220 else "AS-JUN"),
        (
            "Y-JUN" if PANDAS_GE_220 else "A-JUN",
            "YS-JUN" if PANDAS_GE_220 else "AS-JUN",
        ),
        ("BA", "BYS" if PANDAS_GE_220 else "BAS"),
        ("2BA", "2BYS" if PANDAS_GE_220 else "2BAS"),
        ("BY", "BYS" if PANDAS_GE_220 else "BAS"),
        ("Y", "YS" if PANDAS_GE_220 else "AS"),
        (pd.Timedelta(seconds=1), pd.Timedelta(seconds=1)),
    ],
)
def test_map_freq_to_period_start(freq, expected_freq):
    # Handle warnnigs/errors from deprecation cycle in pandas
    ctx = contextlib.nullcontext()
    if PANDAS_GE_300 and (
        freq in ("A", "A-JUN", "BA", "2BA")
        or freq not in ("ME", "MS", pd.Timedelta(seconds=1), "2QS-FEB")
    ):
        ctx = pytest.raises(ValueError, match="Invalid frequency")
    elif PANDAS_GE_220 and freq not in ("ME", "MS", pd.Timedelta(seconds=1), "2QS-FEB"):
        ctx = pytest.warns(
            FutureWarning, match="is deprecated and will be removed in a future version"
        )
    with ctx:
        new_freq = _map_freq_to_period_start(freq)
        assert new_freq == expected_freq


def test_repartition_input_errors():
    df = pd.DataFrame({"x": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=1)
    with pytest.raises(ValueError):
        ddf.repartition(npartitions=5, divisions=[None, None])
    with pytest.raises(ValueError):
        ddf.repartition(npartitions=5, partition_size="5MiB")


def test_embarrassingly_parallel_operations():
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, None, 6], "y": list("abdabd")},
        index=[10, 20, 30, 40, 50, 60],
    )
    a = dd.from_pandas(df, 2)

    assert_eq(a.x.astype("float32"), df.x.astype("float32"))
    assert a.x.astype("float32").compute().dtype == "float32"

    assert_eq(a.x.dropna(), df.x.dropna())

    assert_eq(a.x.between(2, 4), df.x.between(2, 4))

    assert_eq(a.x.clip(2, 4), df.x.clip(2, 4))

    assert_eq(a.x.notnull(), df.x.notnull())
    assert_eq(a.x.isnull(), df.x.isnull())
    assert_eq(a.notnull(), df.notnull())
    assert_eq(a.isnull(), df.isnull())

    assert len(a.sample(frac=0.5).compute()) < len(df)


def test_fillna():
    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.fillna(100), df.fillna(100))
    assert_eq(ddf.A.fillna(100), df.A.fillna(100))
    assert_eq(ddf.A.fillna(ddf["A"].mean()), df.A.fillna(df["A"].mean()))
    assert_eq(ddf.fillna(100, axis=1), df.fillna(100, axis=1))

    pytest.raises(ValueError, lambda: ddf.A.fillna(0, axis=1))


def test_ffill():
    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.ffill(), df.ffill())
    assert_eq(ddf.A.ffill(), df.A.ffill())
    assert_eq(ddf.ffill(limit=2), df.ffill(limit=2))
    assert_eq(ddf.A.ffill(limit=2), df.A.ffill(limit=2))
    assert_eq(ddf.ffill(axis=1), df.ffill(axis=1))
    assert_eq(ddf.ffill(limit=2, axis=1), df.ffill(limit=2, axis=1))

    df = _compat.makeMissingDataframe()
    df.iloc[:15, 0] = np.nan  # all NaN partition
    ddf = dd.from_pandas(df, npartitions=5, sort=False)
    pytest.raises(ValueError, lambda: ddf.ffill().compute())
    assert_eq(df.ffill(limit=3), ddf.ffill(limit=3))


def test_bfill():
    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.bfill(), df.bfill())
    assert_eq(ddf.A.bfill(), df.A.bfill())

    assert_eq(ddf.bfill(limit=2), df.bfill(limit=2))
    assert_eq(ddf.A.bfill(limit=2), df.A.bfill(limit=2))

    df = _compat.makeMissingDataframe()
    df.iloc[:15, 0] = np.nan  # all NaN partition
    ddf = dd.from_pandas(df, npartitions=5, sort=False)
    pytest.raises(ValueError, lambda: ddf.bfill().compute())
    assert_eq(df.bfill(limit=3), ddf.bfill(limit=3))


@pytest.mark.parametrize("optimize", [True, False])
def test_delayed_roundtrip(optimize):
    df1 = d + 1 + 1
    delayed = df1.to_delayed(optimize_graph=optimize)

    assert len(delayed) == df1.npartitions
    dm = d.a.mean().to_delayed(optimize_graph=optimize)

    delayed2 = [x * 2 - dm for x in delayed]

    for x in delayed2:
        x.dask.validate()

    df3 = dd.from_delayed(delayed2, meta=df1, divisions=df1.divisions)
    df4 = df3 - 1 - 1
    assert_eq(df4, (full + 2) * 2 - full.a.mean() - 2)


def test_from_delayed_lazy_if_meta_provided():
    """Ensure that the graph is 100% lazily evaluated if meta is provided"""

    @dask.delayed
    def raise_exception():
        raise RuntimeError()

    tasks = [raise_exception()]
    ddf = dd.from_delayed(tasks, meta=dict(a=float))

    with pytest.raises(RuntimeError):
        ddf.compute()


def test_fillna_duplicate_index():
    @dask.delayed
    def f():
        return pd.DataFrame(dict(a=[1.0], b=[np.nan]))

    ddf = dd.from_delayed([f(), f()], meta=dict(a=float, b=float))
    ddf.b = ddf.b.fillna(ddf.a)
    ddf.compute()


def test_fillna_multi_dataframe():
    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.A.fillna(ddf.B), df.A.fillna(df.B))
    assert_eq(ddf.B.fillna(ddf.A), df.B.fillna(df.A))


def test_fillna_dask_dataframe_input():
    df = _compat.makeMissingDataframe()
    df1 = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=5)
    ddf1 = dd.from_pandas(df1, npartitions=3)

    assert_eq(ddf.fillna(ddf1), df.fillna(df1))

    ddf_unknown = dd.from_pandas(df, npartitions=5, sort=False)
    assert_eq(ddf_unknown.fillna(ddf1), df.fillna(df1))


def test_ffill_bfill():
    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.ffill(), df.ffill())
    assert_eq(ddf.bfill(), df.bfill())
    assert_eq(ddf.ffill(axis=1), df.ffill(axis=1))
    assert_eq(ddf.bfill(axis=1), df.bfill(axis=1))


def test_fillna_series_types():
    # https://github.com/dask/dask/issues/2809
    df = pd.DataFrame({"A": [1, np.nan, 3], "B": [1, np.nan, 3]})
    ddf = dd.from_pandas(df, npartitions=2)
    fill_value = pd.Series([1, 10], index=["A", "C"])
    assert_eq(ddf.fillna(fill_value), df.fillna(fill_value))


def test_sample():
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, None, 6], "y": list("abdabd")},
        index=[10, 20, 30, 40, 50, 60],
    )
    a = dd.from_pandas(df, 2)

    b = a.sample(frac=0.5)

    assert_eq(b, b)

    c = a.sample(frac=0.5, random_state=1234)
    d = a.sample(frac=0.5, random_state=1234)
    assert_eq(c, d)

    assert a.sample(frac=0.5)._name != a.sample(frac=0.5)._name


def test_sample_without_replacement():
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, None, 6], "y": list("abdabd")},
        index=[10, 20, 30, 40, 50, 60],
    )
    a = dd.from_pandas(df, 2)
    b = a.sample(frac=0.7, replace=False)
    bb = b.index.compute()
    assert len(bb) == len(set(bb))


def test_sample_raises():
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, None, 6], "y": list("abdabd")},
        index=[10, 20, 30, 40, 50, 60],
    )
    a = dd.from_pandas(df, 2)

    # Make sure frac is replaced with n when 0 <= n <= 1
    # This is so existing code (i.e. ddf.sample(0.5)) won't break
    with pytest.warns(UserWarning):
        b = a.sample(0.5, random_state=1234)
    c = a.sample(frac=0.5, random_state=1234)
    assert_eq(b, c)

    with pytest.raises(ValueError):
        a.sample(n=10)

    # Make sure frac is provided
    with pytest.raises(ValueError):
        a.sample(frac=None)


def test_query():
    pytest.importorskip("numexpr")

    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.query("x**2 > y"), df.query("x**2 > y"))
    assert_eq(
        ddf.query("x**2 > @value", local_dict={"value": 4}),
        df.query("x**2 > @value", local_dict={"value": 4}),
    )


@pytest.mark.parametrize(
    "include, exclude",
    [
        ([int], None),
        (None, [int]),
        ([np.number, get_string_dtype()], [float]),
        (["datetime"], None),
    ],
)
def test_select_dtypes(include, exclude):
    n = 10
    df = pd.DataFrame(
        {
            "cint": [1] * n,
            "cstr": ["a"] * n,
            "clfoat": [1.0] * n,
            "cdt": pd.date_range("2016-01-01", periods=n),
        }
    )
    df["cstr"] = df["cstr"].astype(get_string_dtype())
    a = dd.from_pandas(df, npartitions=2)
    result = a.select_dtypes(include=include, exclude=exclude)
    expected = df.select_dtypes(include=include, exclude=exclude)
    assert_eq(result, expected)

    # count dtypes
    assert_eq_dtypes(a, df)
    assert_eq_dtypes(result, expected)


def test_deterministic_apply_concat_apply_names():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    assert sorted(a.x.nlargest(2).dask) == sorted(a.x.nlargest(2).dask)
    assert sorted(a.x.nlargest(2).dask) != sorted(a.x.nlargest(3).dask)
    assert sorted(a.x.drop_duplicates().dask) == sorted(a.x.drop_duplicates().dask)
    assert sorted(a.groupby("x").y.mean().dask) == sorted(a.groupby("x").y.mean().dask)


def test_reduction_method():
    df = pd.DataFrame({"x": range(50), "y": range(50, 100)})
    ddf = dd.from_pandas(df, npartitions=4)

    chunk = lambda x, val=0: (x >= val).sum()
    agg = lambda x: x.sum()

    # Output of chunk is a scalar
    res = ddf.x.reduction(chunk, aggregate=agg)
    assert_eq(res, df.x.count())

    # Output of chunk is a series
    res = ddf.reduction(chunk, aggregate=agg)
    assert res._name == ddf.reduction(chunk, aggregate=agg)._name
    assert_eq(res, df.count())

    # Test with keywords
    res2 = ddf.reduction(chunk, aggregate=agg, chunk_kwargs={"val": 25})
    assert (
        res2._name
        == ddf.reduction(chunk, aggregate=agg, chunk_kwargs={"val": 25})._name
    )
    assert res2._name != res._name
    assert_eq(res2, (df >= 25).sum())

    # Output of chunk is a dataframe
    def sum_and_count(x):
        return pd.DataFrame({"sum": x.sum(), "count": x.count()})

    res = ddf.reduction(sum_and_count, aggregate=lambda x: x.groupby(level=0).sum())

    assert_eq(res, pd.DataFrame({"sum": df.sum(), "count": df.count()}))


def test_reduction_method_split_every():
    df = pd.Series([1] * 60)
    ddf = dd.from_pandas(df, npartitions=15)

    def chunk(x, constant=0):
        return x.sum() + constant

    def combine(x, constant=0):
        return x.sum() + constant + 1

    def agg(x, constant=0):
        return x.sum() + constant + 2

    f = lambda n: ddf.reduction(
        chunk,
        aggregate=agg,
        combine=combine,
        chunk_kwargs=dict(constant=1.0),
        combine_kwargs=dict(constant=2.0),
        aggregate_kwargs=dict(constant=3.0),
        split_every=n,
    )

    assert_max_deps(f(3), 3)
    assert_max_deps(f(4), 4, False)
    assert_max_deps(f(5), 5)
    assert f(15).dask.keys() == f(ddf.npartitions).dask.keys()

    r3 = f(3)
    r4 = f(4)
    assert r3._name != r4._name

    # Keywords are different for each step
    assert f(3).compute() == 60 + 15 + 7 * (2 + 1) + (3 + 2)
    # Keywords are same for each step
    res = ddf.reduction(
        chunk, aggregate=agg, combine=combine, constant=3.0, split_every=3
    )
    assert res.compute() == 60 + 15 * 3 + 7 * (3 + 1) + (3 + 2)
    # No combine provided, combine is agg
    res = ddf.reduction(chunk, aggregate=agg, constant=3.0, split_every=3)
    assert res.compute() == 60 + 15 * 3 + 8 * (3 + 2)

    # split_every must be >= 2
    with pytest.raises(ValueError):
        f(1)

    # combine_kwargs with no combine provided
    with pytest.raises(ValueError):
        ddf.reduction(
            chunk,
            aggregate=agg,
            split_every=3,
            chunk_kwargs=dict(constant=1.0),
            combine_kwargs=dict(constant=2.0),
            aggregate_kwargs=dict(constant=3.0),
        )


def test_pipe():
    df = pd.DataFrame({"x": range(50), "y": range(50, 100)})
    ddf = dd.from_pandas(df, npartitions=4)

    def f(x, y, z=0):
        return x + y + z

    assert_eq(ddf.pipe(f, 1, z=2), f(ddf, 1, z=2))
    assert_eq(ddf.x.pipe(f, 1, z=2), f(ddf.x, 1, z=2))


def test_gh_517():
    arr = np.random.randn(100, 2)
    df = pd.DataFrame(arr, columns=["a", "b"])
    ddf = dd.from_pandas(df, 2)
    assert ddf.index.nunique().compute() == 100

    ddf2 = dd.from_pandas(pd.concat([df, df]), 5)
    assert ddf2.index.nunique().compute() == 100


def test_drop_axis_1():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [5, 6, 7, 8], "z": [9, 10, 11, 12]})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(ddf.drop("y", axis=1), df.drop("y", axis=1))
    assert_eq(ddf.drop(["y", "z"], axis=1), df.drop(["y", "z"], axis=1))
    with pytest.raises((ValueError, KeyError)):
        ddf.drop(["a", "x"], axis=1)
    assert_eq(
        ddf.drop(["a", "x"], axis=1, errors="ignore"),
        df.drop(["a", "x"], axis=1, errors="ignore"),
    )
    assert_eq(ddf.drop(columns=["y", "z"]), df.drop(columns=["y", "z"]))


@pytest.mark.parametrize("columns", [["b"], []])
def test_drop_columns(columns):
    # Check both populated and empty list argument
    # https://github.com/dask/dask/issues/6870

    df = pd.DataFrame(
        {
            "a": [2, 4, 6, 8],
            "b": ["1a", "2b", "3c", "4d"],
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = ddf.drop(columns=columns)
    ddf["new"] = ddf["a"] + 1  # Check that ddf2 is not modified

    assert_eq(df.drop(columns=columns), ddf2)


def test_gh580():
    df = pd.DataFrame({"x": np.arange(10, dtype=float)})
    ddf = dd.from_pandas(df, 2)
    assert_eq(np.cos(df["x"]), np.cos(ddf["x"]))
    assert_eq(np.cos(df["x"]), np.cos(ddf["x"]))


def test_rename_dict():
    renamer = {"a": "A", "b": "B"}
    assert_eq(d.rename(columns=renamer), full.rename(columns=renamer))


def test_rename_function():
    renamer = lambda x: x.upper()
    assert_eq(d.rename(columns=renamer), full.rename(columns=renamer))


def test_rename_index():
    renamer = {0: 1}
    pytest.raises(ValueError, lambda: d.rename(index=renamer))


def test_to_timestamp():
    index = pd.period_range(freq="Y", start="1/1/2001", end="12/1/2004")
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]}, index=index)
    ddf = dd.from_pandas(df, npartitions=3)
    assert_eq(ddf.to_timestamp(), df.to_timestamp(), check_freq=False)
    assert_eq(
        ddf.to_timestamp(freq="M", how="s").compute(),
        df.to_timestamp(freq="M", how="s"),
        check_freq=False,
    )
    assert_eq(ddf.x.to_timestamp(), df.x.to_timestamp(), check_freq=False)
    assert_eq(
        ddf.x.to_timestamp(freq="M", how="s").compute(),
        df.x.to_timestamp(freq="M", how="s"),
        check_freq=False,
    )

    ddf = dd.from_pandas(df, npartitions=4)
    assert_eq(ddf.to_timestamp(how="end"), df.to_timestamp(how="end"))
    assert_eq(ddf.x.to_timestamp(how="end"), df.x.to_timestamp(how="end"))


def test_to_frame():
    s = pd.Series([1, 2, 3], name="foo")
    a = dd.from_pandas(s, npartitions=2)

    assert_eq(s.to_frame(), a.to_frame())
    assert_eq(s.to_frame("bar"), a.to_frame("bar"))


@pytest.mark.parametrize("as_frame", [False, False])
def test_to_dask_array_raises(as_frame):
    s = pd.Series([1, 2, 3, 4, 5, 6], name="foo")
    a = dd.from_pandas(s, npartitions=2)

    if as_frame:
        a = a.to_frame()

    with pytest.raises(ValueError, match="4 != 2"):
        a.to_dask_array((1, 2, 3, 4))

    with pytest.raises(ValueError, match="Unexpected value"):
        a.to_dask_array(5)


@pytest.mark.parametrize("as_frame", [False, True])
def test_to_dask_array_unknown(as_frame):
    s = pd.Series([1, 2, 3, 4, 5], name="foo")
    a = dd.from_pandas(s, npartitions=2)

    if as_frame:
        a = a.to_frame()

    result = a.to_dask_array()
    assert isinstance(result, da.Array)
    result = result.chunks

    if as_frame:
        assert len(result) == 2
        assert result[1] == (1,)
    else:
        assert len(result) == 1

    result = result[0]
    assert len(result) == 2
    assert all(np.isnan(x) for x in result)


@pytest.mark.parametrize(
    "lengths,as_frame,meta",
    [
        ([2, 2, 1], False, None),
        (True, False, None),
        (True, False, np.array([], dtype="f4")),
    ],
)
def test_to_dask_array(meta, as_frame, lengths):
    s = pd.Series([1, 2, 3, 4, 5], name="foo", dtype="i4")
    a = dd.from_pandas(s, chunksize=2)

    if as_frame:
        a = a.to_frame()

    result = a.to_dask_array(lengths=lengths, meta=meta)
    assert isinstance(result, da.Array)

    expected_chunks = ((2, 2, 1),)

    if as_frame:
        expected_chunks = expected_chunks + ((1,),)

    assert result.chunks == expected_chunks


def test_apply():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(
        ddf.x.apply(lambda x: x + 1, meta=("x", int)), df.x.apply(lambda x: x + 1)
    )

    # specify meta
    assert_eq(
        ddf.apply(lambda xy: xy.iloc[0] + xy.iloc[1], axis=1, meta=(None, int)),
        df.apply(lambda xy: xy.iloc[0] + xy.iloc[1], axis=1),
    )
    assert_eq(
        ddf.apply(lambda xy: xy.iloc[0] + xy.iloc[1], axis="columns", meta=(None, int)),
        df.apply(lambda xy: xy.iloc[0] + xy.iloc[1], axis="columns"),
    )

    # inference
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        assert_eq(
            ddf.apply(lambda xy: xy.iloc[0] + xy.iloc[1], axis=1),
            df.apply(lambda xy: xy.iloc[0] + xy.iloc[1], axis=1),
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        assert_eq(ddf.apply(lambda xy: xy, axis=1), df.apply(lambda xy: xy, axis=1))

    # specify meta
    func = lambda x: pd.Series([x, x])
    assert_eq(ddf.x.apply(func, meta=[(0, int), (1, int)]), df.x.apply(func))
    # inference
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        assert_eq(ddf.x.apply(func), df.x.apply(func))

    # axis=0
    with pytest.raises(NotImplementedError):
        ddf.apply(lambda xy: xy, axis=0)

    with pytest.raises(NotImplementedError):
        ddf.apply(lambda xy: xy, axis="index")


def test_apply_warns():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    func = lambda row: row["x"] + row["y"]

    with pytest.warns(UserWarning) as w:
        ddf.apply(func, axis=1)
    assert len(w) == 1

    with warnings.catch_warnings(record=True) as record:
        ddf.apply(func, axis=1, meta=(None, int))
    assert not record

    with pytest.warns(UserWarning) as w:
        ddf.apply(lambda x: x, axis=1)
    assert len(w) == 1
    assert "'x'" in str(w[0].message)
    assert "int64" in str(w[0].message)


@pytest.mark.skipif(not PANDAS_GE_210, reason="Not available before")
@pytest.mark.parametrize("na_action", [None, "ignore"])
def test_dataframe_map(na_action):
    df = pd.DataFrame({"x": [1, 2, 3, np.nan], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)
    with pytest.warns(UserWarning, match="meta"):
        assert_eq(
            ddf.map(lambda x: x + 1, na_action=na_action),
            df.map(lambda x: x + 1, na_action=na_action),
        )
        assert_eq(ddf.map(lambda x: (x, x)), df.map(lambda x: (x, x)))


@pytest.mark.skipif(PANDAS_GE_210, reason="Available at 2.1")
def test_dataframe_map_raises():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)
    with pytest.raises(NotImplementedError, match="DataFrame.map requires pandas"):
        ddf.map(lambda x: x + 1)


def test_add_prefix():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5], "y": [4, 5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.add_prefix("abc"), df.add_prefix("abc"))
    assert_eq(ddf.x.add_prefix("abc"), df.x.add_prefix("abc"))


def test_add_suffix():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5], "y": [4, 5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.add_suffix("abc"), df.add_suffix("abc"))
    assert_eq(ddf.x.add_suffix("abc"), df.x.add_suffix("abc"))


def test_abs():
    df = pd.DataFrame(
        {
            "A": [1, -2, 3, -4, 5],
            "B": [-6.0, -7, -8, -9, 10],
            "C": ["a", "b", "c", "d", "e"],
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.A.abs(), df.A.abs())
    assert_eq(ddf[["A", "B"]].abs(), df[["A", "B"]].abs())
    # raises TypeError with object dtype, but NotImplementedError with string[pyarrow]
    with pytest.raises((TypeError, NotImplementedError, ValueError)):
        ddf.C.abs()
    with pytest.raises((TypeError, NotImplementedError)):
        ddf.abs()


def test_round():
    df = pd.DataFrame({"col1": [1.123, 2.123, 3.123], "col2": [1.234, 2.234, 3.234]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.round(), df.round())
    assert_eq(ddf.round(2), df.round(2))


@pytest.mark.parametrize(
    "numeric_only",
    [
        None,
        True,
        False,
    ],
)
def test_cov_dataframe(numeric_only):
    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=6)

    numeric_only_kwarg = {}
    if numeric_only is not None:
        numeric_only_kwarg = {"numeric_only": numeric_only}

    res = ddf.cov(**numeric_only_kwarg)
    res2 = ddf.cov(**numeric_only_kwarg, split_every=2)
    res3 = ddf.cov(10, **numeric_only_kwarg)
    res4 = ddf.cov(10, **numeric_only_kwarg, split_every=2)
    sol = df.cov(**numeric_only_kwarg)
    sol2 = df.cov(10, **numeric_only_kwarg)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == ddf.cov(**numeric_only_kwarg)._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name


def test_cov_series():
    df = _compat.makeMissingDataframe()
    a = df.A
    b = df.B
    da = dd.from_pandas(a, npartitions=6)
    db = dd.from_pandas(b, npartitions=6)

    res = da.cov(db)
    res2 = da.cov(db, split_every=2)
    res3 = da.cov(db, 10)
    res4 = da.cov(db, 10, split_every=2)
    sol = a.cov(b)
    sol2 = a.cov(b, 10)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == da.cov(db)._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name


@pytest.mark.gpu
@pytest.mark.skip(reason="https://github.com/rapidsai/cudf/issues/16560")
@pytest.mark.parametrize(
    "numeric_only",
    [None, True, False],
)
def test_cov_gpu(numeric_only):
    cudf = pytest.importorskip("cudf")

    # cudf DataFrame
    df = cudf.from_pandas(_compat.makeDataFrame())
    ddf = dd.from_pandas(df, npartitions=6)

    numeric_only_kwarg = {}
    if numeric_only is not None:
        numeric_only_kwarg = {"numeric_only": numeric_only}

    res = ddf.cov(**numeric_only_kwarg)
    res2 = ddf.cov(**numeric_only_kwarg, split_every=2)
    sol = df.cov(**numeric_only_kwarg)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert res._name == ddf.cov(**numeric_only_kwarg)._name
    assert res._name != res2._name


def test_corr():
    # DataFrame
    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=6)

    res = ddf.corr()
    res2 = ddf.corr(split_every=2)
    res3 = ddf.corr(min_periods=10)
    res4 = ddf.corr(min_periods=10, split_every=2)
    sol = df.corr()
    sol2 = df.corr(min_periods=10)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == ddf.corr()._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name

    pytest.raises(NotImplementedError, lambda: ddf.corr(method="spearman"))

    # Series
    a = df.A
    b = df.B
    da = dd.from_pandas(a, npartitions=6)
    db = dd.from_pandas(b, npartitions=6)

    res = da.corr(db)
    res2 = da.corr(db, split_every=2)
    res3 = da.corr(db, min_periods=10)
    res4 = da.corr(db, min_periods=10, split_every=2)
    sol = da.corr(db)
    sol2 = da.corr(db, min_periods=10)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == da.corr(db)._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name

    pytest.raises(NotImplementedError, lambda: da.corr(db, method="spearman"))
    pytest.raises(TypeError, lambda: da.corr(ddf))


@pytest.mark.gpu
@pytest.mark.skip(reason="https://github.com/rapidsai/cudf/issues/16560")
def test_corr_gpu():
    cudf = pytest.importorskip("cudf")

    # cudf DataFrame
    df = cudf.from_pandas(_compat.makeDataFrame())
    ddf = dd.from_pandas(df, npartitions=6)

    res = ddf.corr()
    res2 = ddf.corr(split_every=2)
    sol = df.corr()
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert res._name == ddf.corr()._name
    assert res._name != res2._name


def test_corr_same_name():
    # Series with same names (see https://github.com/dask/dask/issues/4906)

    df = _compat.makeMissingDataframe()
    ddf = dd.from_pandas(df, npartitions=6)

    result = ddf.A.corr(ddf.B.rename("A"))
    expected = ddf.A.corr(ddf.B)
    assert_eq(result, expected)

    # test with split_every
    result2 = ddf.A.corr(ddf.B.rename("A"), split_every=2)
    assert_eq(result2, expected)


@pytest.mark.parametrize("chunksize", [1, 2])
def test_cov_corr_meta(chunksize):
    df = pd.DataFrame(
        {
            "a": np.array([1, 2, 3, 4]),
            "b": np.array([1.0, 2.0, 3.0, 4.0], dtype="f4"),
            "c": np.array([1.0, 2.0, 3.0, 4.0]),
        },
        index=pd.Index([1, 2, 3, 4], name="myindex"),
    )
    ddf = dd.from_pandas(df, chunksize=chunksize)
    assert_eq(ddf.corr(), df.corr())
    assert_eq(ddf.cov(), df.cov())
    assert ddf.a.cov(ddf.b)._meta.dtype == "f8"
    assert ddf.a.corr(ddf.b)._meta.dtype == "f8"


@pytest.mark.slow
def test_cov_corr_stable():
    df = pd.DataFrame(np.random.uniform(-1, 1, (20000000, 2)), columns=["a", "b"])
    ddf = dd.from_pandas(df, npartitions=50)
    assert_eq(ddf.cov(split_every=8), df.cov())
    assert_eq(ddf.corr(split_every=8), df.corr())


@pytest.mark.parametrize(
    "numeric_only",
    [
        pytest.param(
            None,
            marks=pytest.mark.xfail(reason="fails with non-numeric data"),
        ),
        True,
        pytest.param(
            False,
            marks=[
                pytest.mark.xfail(reason="fails with non-numeric data"),
            ],
        ),
    ],
)
def test_cov_corr_mixed(numeric_only):
    size = 1000
    d = {
        "dates": pd.date_range("2015-01-01", periods=size, freq="1min"),
        "unique_id": np.arange(0, size),
        "ints": np.random.randint(0, size, size=size),
        "floats": np.random.randn(size),
        "bools": np.random.choice([0, 1], size=size),
        "int_nans": np.random.choice([0, 1, np.nan], size=size),
        "float_nans": np.random.choice([0.0, 1.0, np.nan], size=size),
        "constant": 1,
        "int_categorical": np.random.choice([10, 20, 30, 40, 50], size=size),
        "categorical_binary": np.random.choice(["a", "b"], size=size),
        "categorical_nans": np.random.choice(["a", "b", "c"], size=size),
    }
    df = pd.DataFrame(d)
    df["hardbools"] = df["bools"] == 1
    df["categorical_nans"] = df["categorical_nans"].replace("c", np.nan)
    df["categorical_binary"] = df["categorical_binary"].astype("category")
    df["unique_id"] = df["unique_id"].astype(str)

    ddf = dd.from_pandas(df, npartitions=20)

    numeric_only_kwarg = {}
    if numeric_only is not None:
        numeric_only_kwarg = {"numeric_only": numeric_only}

    # Corr
    expected = df.corr(**numeric_only_kwarg)
    result = ddf.corr(split_every=4, **numeric_only_kwarg)
    assert_eq(result, expected, check_divisions=False)

    # Cov
    expected = df.cov(**numeric_only_kwarg)
    result = ddf.cov(split_every=4, **numeric_only_kwarg)
    assert_eq(result, expected, check_divisions=False)


def test_autocorr():
    x = pd.Series(np.random.random(100))
    dx = dd.from_pandas(x, npartitions=10)
    assert_eq(dx.autocorr(2), x.autocorr(2))
    assert_eq(dx.autocorr(0), x.autocorr(0))
    assert_eq(dx.autocorr(-2), x.autocorr(-2))
    assert_eq(dx.autocorr(2, split_every=3), x.autocorr(2))
    pytest.raises(TypeError, lambda: dx.autocorr(1.5))


def test_apply_infer_columns():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    def return_df(x):
        # will create new DataFrame which columns is ['sum', 'mean']
        return pd.Series([x.sum(), x.mean()], index=["sum", "mean"])

    # DataFrame to completely different DataFrame
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = ddf.apply(return_df, axis=1)
    assert isinstance(result, dd.DataFrame)
    tm.assert_index_equal(result.columns, pd.Index(["sum", "mean"]))
    assert_eq(result, df.apply(return_df, axis=1))

    # DataFrame to Series
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = ddf.apply(lambda x: 1, axis=1)
    assert isinstance(result, dd.Series)
    assert result.name is None
    assert_eq(result, df.apply(lambda x: 1, axis=1))

    def return_df2(x):
        return pd.Series([x * 2, x * 3], index=["x2", "x3"])

    # Series to completely different DataFrame
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = ddf.x.apply(return_df2)
    assert isinstance(result, dd.DataFrame)
    tm.assert_index_equal(result.columns, pd.Index(["x2", "x3"]))
    assert_eq(result, df.x.apply(return_df2))

    # Series to Series
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        result = ddf.x.apply(lambda x: 1)
    assert isinstance(result, dd.Series)
    assert result.name == "x"
    assert_eq(result, df.x.apply(lambda x: 1))


def test_index_time_properties():
    i = _compat.makeTimeSeries()
    a = dd.from_pandas(i, npartitions=3)

    assert "day" in dir(a.index)
    # returns a numpy array in pandas, but an Index in dask
    assert_eq(a.index.day, pd.Index(i.index.day))
    assert_eq(a.index.month, pd.Index(i.index.month))


def test_nlargest_nsmallest():
    from string import ascii_lowercase

    df = pd.DataFrame(
        {
            "a": np.random.permutation(20),
            "b": list(ascii_lowercase[:20]),
            "c": np.random.permutation(20).astype("float64"),
        }
    )
    ddf = dd.from_pandas(df, npartitions=3)

    for m in ["nlargest", "nsmallest"]:
        f = lambda df=df, m=m, *args, **kwargs: getattr(df, m)(*args, **kwargs)

        res = f(ddf, m, 5, "a")
        res2 = f(ddf, m, 5, "a", split_every=2)
        sol = f(df, m, 5, "a")
        assert_eq(res, sol)
        assert_eq(res2, sol)
        assert res._name != res2._name

        res = f(ddf, m, 5, ["a", "c"])
        res2 = f(ddf, m, 5, ["a", "c"], split_every=2)
        sol = f(df, m, 5, ["a", "c"])
        assert_eq(res, sol)
        assert_eq(res2, sol)
        assert res._name != res2._name

        res = f(ddf.a, m, 5)
        res2 = f(ddf.a, m, 5, split_every=2)
        sol = f(df.a, m, 5)
        assert_eq(res, sol)
        assert_eq(res2, sol)
        assert res._name != res2._name


def test_nlargest_nsmallest_raises():
    df = pd.DataFrame({"a": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=2)
    with pytest.raises(TypeError, match="required positional"):
        ddf.nlargest()
    with pytest.raises(TypeError, match="required positional"):
        ddf.nsmallest()


def test_reset_index():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    sol = df.reset_index()
    res = ddf.reset_index()
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)

    sol = df.reset_index(drop=True)
    res = ddf.reset_index(drop=True)
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)

    sol = df.x.reset_index()
    res = ddf.x.reset_index()
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)

    sol = df.x.reset_index(drop=True)
    res = ddf.x.reset_index(drop=True)
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)


def test_dataframe_compute_forward_kwargs():
    x = dd.from_pandas(pd.DataFrame({"a": range(10)}), npartitions=2).a.sum()
    x.compute(bogus_keyword=10)


def test_series_iter():
    s = pd.DataFrame({"x": [1, 2, 3, 4]})
    ds = dd.from_pandas(s, npartitions=2)
    for a, b in zip(s["x"], ds["x"]):
        assert a == b


def test_dataframe_iterrows():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    for a, b in zip(df.iterrows(), ddf.iterrows()):
        tm.assert_series_equal(a[1], b[1])


def test_dataframe_itertuples():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    for a, b in zip(df.itertuples(), ddf.itertuples()):
        assert a == b


@pytest.mark.parametrize(
    "columns",
    [
        ("x", "y"),
    ],
)
def test_dataframe_items(columns):
    df = pd.DataFrame([[1, 10], [2, 20], [3, 30], [4, 40]], columns=columns)
    ddf = dd.from_pandas(df, npartitions=2)
    for a, b in zip(df.items(), ddf.items()):
        assert a[0] == b[0]  # column name
        assert_eq(a[1], b[1].compute())  # column values


def test_dataframe_itertuples_with_index_false():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    for a, b in zip(df.itertuples(index=False), ddf.itertuples(index=False)):
        assert a == b


def test_dataframe_itertuples_with_name_none():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    for a, b in zip(df.itertuples(name=None), ddf.itertuples(name=None)):
        assert a == b
        assert type(a) is type(b)


def test_astype():
    df = pd.DataFrame(
        {"x": [1, 2, 3, None], "y": [10, 20, 30, 40]}, index=[10, 20, 30, 40]
    )
    a = dd.from_pandas(df, 2)

    assert_eq(a.astype(float), df.astype(float))
    assert_eq(a.x.astype(float), df.x.astype(float))


def test_astype_categoricals():
    df = pd.DataFrame(
        {
            "x": ["a", "b", "c", "b", "c"],
            "y": ["x", "y", "z", "x", "y"],
            "z": [1, 2, 3, 4, 5],
        }
    )
    df = df.astype({"y": "category"})
    ddf = dd.from_pandas(df, 2)
    assert ddf.y.cat.known

    ddf2 = ddf.astype({"x": "category"})
    assert not ddf2.x.cat.known
    assert ddf2.y.cat.known
    assert ddf2.x.dtype == "category"
    assert ddf2.compute().x.dtype == "category"

    dx = ddf.x.astype("category")
    assert not dx.cat.known
    assert dx.dtype == "category"
    assert dx.compute().dtype == "category"


def test_astype_categoricals_known():
    df = pd.DataFrame(
        {
            "x": ["a", "b", "c", "b", "c"],
            "y": ["x", "y", "z", "y", "z"],
            "z": ["b", "b", "b", "c", "b"],
            "other": [1, 2, 3, 4, 5],
        }
    )
    ddf = dd.from_pandas(df, 2)

    abc = pd.api.types.CategoricalDtype(["a", "b", "c"], ordered=False)
    category = pd.api.types.CategoricalDtype(ordered=False)

    # DataFrame
    ddf2 = ddf.astype({"x": abc, "y": category, "z": "category", "other": "f8"})

    for col, known in [("x", True), ("y", False), ("z", False)]:
        x = getattr(ddf2, col)
        assert isinstance(x.dtype, pd.CategoricalDtype)
        assert x.cat.known == known

    # Series
    for dtype, known in [("category", False), (category, False), (abc, True)]:
        dx2 = ddf.x.astype(dtype)
        assert isinstance(dx2.dtype, pd.CategoricalDtype)
        assert dx2.cat.known == known


def test_groupby_callable():
    a = pd.DataFrame({"x": [1, 2, 3, None], "y": [10, 20, 30, 40]}, index=[1, 2, 3, 4])
    b = dd.from_pandas(a, 2)

    def iseven(x):
        return x % 2 == 0

    assert_eq(a.groupby(iseven).y.sum(), b.groupby(iseven).y.sum())
    assert_eq(a.y.groupby(iseven).sum(), b.y.groupby(iseven).sum())


def test_methods_tokenize_differently():
    df = pd.DataFrame({"x": [1, 2, 3, 4]})
    df = dd.from_pandas(df, npartitions=1)
    assert (
        df.x.map_partitions(lambda x: pd.Series(x.min()))._name
        != df.x.map_partitions(lambda x: pd.Series(x.max()))._name
    )


def _assert_info(df, ddf, memory_usage=True):
    from io import StringIO

    assert isinstance(df, pd.DataFrame)
    assert isinstance(ddf, dd.DataFrame)

    buf_pd, buf_da = StringIO(), StringIO()

    df.info(buf=buf_pd, memory_usage=memory_usage)
    ddf.info(buf=buf_da, verbose=True, memory_usage=memory_usage)

    stdout_pd = buf_pd.getvalue()
    stdout_da = buf_da.getvalue()
    stdout_da = stdout_da.replace(str(type(ddf)), str(type(df)))
    # TODO
    assert stdout_pd == stdout_da


def test_gh_1301():
    df = pd.DataFrame([["1", "2"], ["3", "4"]])
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = ddf.assign(y=ddf[1].astype(int))
    assert_eq(ddf2, df.assign(y=df[1].astype(int)))

    assert ddf2.dtypes["y"] == np.dtype(int)


def test_timeseries_sorted():
    df = _compat.makeTimeDataFrame()
    ddf = dd.from_pandas(df.reset_index(), npartitions=2)
    df.index.name = "index"
    assert_eq(ddf.set_index("index", sorted=True, drop=True), df)


def test_index_errors():
    df = _compat.makeTimeDataFrame()
    ddf = dd.from_pandas(df.reset_index(), npartitions=2)
    with pytest.raises(NotImplementedError, match="You tried to index with this index"):
        ddf.set_index([["A"]])  # should index with ["A"] instead of [["A"]]
    with pytest.raises(NotImplementedError, match="You tried to index with a frame"):
        ddf.set_index(ddf[["A"]])  # should index with ddf["A"] instead of ddf[["A"]]
    with pytest.raises(KeyError, match="has no column|None of"):
        ddf.set_index("foo")  # a column that doesn't exist
    with pytest.raises(KeyError, match="has no column|None of"):
        # column name doesn't need to be a string, but anyhow a KeyError should be raised if not found
        ddf.set_index(0)


@pytest.mark.parametrize("null_value", [None, pd.NaT, pd.NA])
def test_index_nulls(null_value):
    "Setting the index with some non-numeric null raises error"
    df = pd.DataFrame(
        {"numeric": [1, 2, 3, 4], "non_numeric": ["foo", "bar", "foo", "bar"]}
    )
    # an object column with only some nulls fails
    ddf = dd.from_pandas(df, npartitions=2)
    with pytest.raises(NotImplementedError, match="presence of nulls"):
        with pytest.warns(UserWarning, match="meta"):
            ddf.set_index(
                ddf["non_numeric"].map({"foo": "foo", "bar": null_value})
            ).compute()


def test_set_index_with_index():
    "Setting the index with the existing index is a no-op"
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [1, 0, 1, 0]}).set_index("x")
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = ddf.set_index(ddf.index)
    assert ddf2 is ddf
    ddf = ddf.set_index("x", drop=False)
    assert ddf2 is ddf


def test_column_assignment():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [1, 0, 1, 0]})
    ddf = dd.from_pandas(df, npartitions=2)
    orig = ddf.copy()
    ddf["z"] = ddf.x + ddf.y
    df["z"] = df.x + df.y

    assert_eq(df, ddf)
    assert "z" not in orig.columns


def test_array_assignment():
    df = pd.DataFrame({"x": np.random.normal(size=50), "y": np.random.normal(size=50)})
    ddf = dd.from_pandas(df, npartitions=2)
    orig = ddf.copy()

    arr = np.array(np.random.normal(size=50))
    darr = da.from_array(arr, chunks=25)

    df["z"] = arr
    ddf["z"] = darr
    assert_eq(df, ddf)
    assert "z" not in orig.columns

    arr = np.array(np.random.normal(size=(50, 50)))
    darr = da.from_array(arr, chunks=25)
    msg = "Array assignment only supports 1-D arrays"
    with pytest.raises(ValueError, match=msg):
        ddf["z"] = darr

    arr = np.array(np.random.normal(size=50))
    darr = da.from_array(arr, chunks=10)
    msg = "Number of partitions do not match"
    with pytest.raises(ValueError, match=msg):
        ddf["z"] = darr


def test_columns_assignment():
    df = pd.DataFrame({"x": [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)

    df2 = df.assign(y=df.x + 1, z=df.x - 1)
    df[["a", "b"]] = df2[["y", "z"]]

    ddf2 = ddf.assign(y=ddf.x + 1, z=ddf.x - 1)
    ddf[["a", "b"]] = ddf2[["y", "z"]]

    assert_eq(df, ddf)


def test_attribute_assignment():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5], "y": [1.0, 2.0, 3.0, 4.0, 5.0]})
    ddf = dd.from_pandas(df, npartitions=2)

    ddf.y = ddf.x + ddf.y
    assert_eq(ddf, df.assign(y=df.x + df.y))


def test_setitem_triggering_realign():
    a = dd.from_pandas(pd.DataFrame({"A": range(12)}), npartitions=3)
    b = dd.from_pandas(pd.Series(range(12), name="B"), npartitions=4)
    a["C"] = b
    assert len(a) == 12


def test_inplace_operators():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5], "y": [1.0, 2.0, 3.0, 4.0, 5.0]})
    ddf = dd.from_pandas(df, npartitions=2)

    ddf["y"] **= 0.5

    assert_eq(ddf.y, df.y**0.5)
    assert_eq(ddf, df.assign(y=df.y**0.5))


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize(
    "idx",
    [
        np.arange(100),
        sorted(np.random.random(size=100)),
        pd.date_range("20150101", periods=100),
    ],
)
def test_idxmaxmin(idx, skipna):
    df = pd.DataFrame(np.random.randn(100, 5), columns=list("abcde"), index=idx)
    df.iloc[31, 1] = np.nan
    df.iloc[78, 3] = np.nan
    ddf = dd.from_pandas(df, npartitions=3)

    # https://github.com/pandas-dev/pandas/issues/43587
    check_dtype = not all((skipna is False, isinstance(idx, pd.DatetimeIndex)))

    ctx = contextlib.nullcontext()
    if PANDAS_GE_300 and not skipna:
        ctx = pytest.raises(ValueError, match="Encountered an NA")

    with ctx:
        with warnings.catch_warnings(record=True):
            if not skipna and PANDAS_GE_210:
                warnings.simplefilter("ignore", category=FutureWarning)
            assert_eq(
                df.idxmax(axis=1, skipna=skipna), ddf.idxmax(axis=1, skipna=skipna)
            )
            assert_eq(
                df.idxmin(axis=1, skipna=skipna), ddf.idxmin(axis=1, skipna=skipna)
            )

            assert_eq(
                df.idxmax(skipna=skipna),
                ddf.idxmax(skipna=skipna),
                check_dtype=check_dtype,
            )
            assert_eq(
                df.idxmax(skipna=skipna),
                ddf.idxmax(skipna=skipna, split_every=2),
                check_dtype=check_dtype,
            )
            assert (
                ddf.idxmax(skipna=skipna)._name
                != ddf.idxmax(skipna=skipna, split_every=2)._name
            )

            assert_eq(
                df.idxmin(skipna=skipna),
                ddf.idxmin(skipna=skipna),
                check_dtype=check_dtype,
            )
            assert_eq(
                df.idxmin(skipna=skipna),
                ddf.idxmin(skipna=skipna, split_every=2),
                check_dtype=check_dtype,
            )
            assert (
                ddf.idxmin(skipna=skipna)._name
                != ddf.idxmin(skipna=skipna, split_every=2)._name
            )

            assert_eq(df.a.idxmax(skipna=skipna), ddf.a.idxmax(skipna=skipna))
            assert_eq(
                df.a.idxmax(skipna=skipna), ddf.a.idxmax(skipna=skipna, split_every=2)
            )
            assert (
                ddf.a.idxmax(skipna=skipna)._name
                != ddf.a.idxmax(skipna=skipna, split_every=2)._name
            )

            assert_eq(df.a.idxmin(skipna=skipna), ddf.a.idxmin(skipna=skipna))
            assert_eq(
                df.a.idxmin(skipna=skipna), ddf.a.idxmin(skipna=skipna, split_every=2)
            )
            assert (
                ddf.a.idxmin(skipna=skipna)._name
                != ddf.a.idxmin(skipna=skipna, split_every=2)._name
            )


@pytest.mark.parametrize("func", ["idxmin", "idxmax"])
def test_idxmaxmin_numeric_only(func):
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2010, i, 1) for i in range(1, 8)],
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
            "bool": [True, False] * 4,
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(
        getattr(ddf, func)(numeric_only=False),
        getattr(df, func)(numeric_only=False).sort_index(),
    )
    assert_eq(
        getattr(ddf, func)(numeric_only=True),
        getattr(df, func)(numeric_only=True).sort_index(),
    )

    assert_eq(
        getattr(ddf.drop(columns="bool"), func)(numeric_only=True, axis=1),
        getattr(df.drop(columns="bool"), func)(numeric_only=True, axis=1).sort_index(),
    )


def test_idxmaxmin_empty_partitions():
    df = pd.DataFrame(
        {"a": [1, 2, 3], "b": [1.5, 2, 3], "c": [np.nan] * 3, "d": [1, 2, np.nan]}
    )
    empty = df.iloc[:0]

    ddf = dd.concat(
        [dd.from_pandas(df, npartitions=1)]
        + [dd.from_pandas(empty, npartitions=1)] * 10
    )

    if PANDAS_GE_300:
        ctx = pytest.raises(
            ValueError, match="Encountered all NA values|Encountered an NA value with"
        )
    elif PANDAS_GE_210:
        ctx = pytest.warns(FutureWarning, match="all-NA values")
    else:
        ctx = contextlib.nullcontext()

    for skipna in [True, False]:
        with ctx:
            expected = df.idxmin(skipna=skipna)
        if not PANDAS_GE_300:
            # No warning at graph construction time because we don't know
            # about empty partitions prior to computing
            result = ddf.idxmin(skipna=skipna, split_every=3)
            with ctx:
                assert_eq(result, expected)

    assert_eq(
        ddf[["a", "b", "d"]].idxmin(skipna=True, split_every=3),
        df[["a", "b", "d"]].idxmin(skipna=True),
    )

    assert_eq(ddf.b.idxmax(split_every=3), df.b.idxmax())

    # Completely empty raises
    ddf = dd.concat([dd.from_pandas(empty, npartitions=1)] * 10)
    with pytest.raises(ValueError):
        ddf.idxmax().compute()
    with pytest.raises(ValueError):
        ddf.b.idxmax().compute()


def test_mode_numeric_only():
    df = pd.DataFrame(
        {
            "int": [1, 2, 3, 4, 5, 6, 7, 8],
            "float": [1.0, 2.0, 3.0, 4.0, np.nan, 6.0, 7.0, 8.0],
            "dt": [pd.NaT] + [datetime(2010, i, 1) for i in range(1, 8)],
            "timedelta": pd.to_timedelta([1, 2, 3, 4, 5, 6, 7, np.nan]),
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(ddf.mode(numeric_only=False), df.mode(numeric_only=False))
    assert_eq(ddf.mode(), df.mode())
    assert_eq(ddf.mode(numeric_only=True), df.mode(numeric_only=True))


def test_getitem_meta():
    data = {"col1": ["a", "a", "b"], "col2": [0, 1, 0]}

    df = pd.DataFrame(data=data, columns=["col1", "col2"])
    ddf = dd.from_pandas(df, npartitions=1)

    assert_eq(df.col2[df.col1 == "a"], ddf.col2[ddf.col1 == "a"])


def test_getitem_string_subclass():
    df = pd.DataFrame({"column_1": list(range(10))})
    ddf = dd.from_pandas(df, npartitions=3)

    class string_subclass(str):
        pass

    column_1 = string_subclass("column_1")

    assert_eq(df[column_1], ddf[column_1])


@pytest.mark.parametrize("col_type", [list, np.array, pd.Series, pd.Index])
def test_getitem_column_types(col_type):
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    ddf = dd.from_pandas(df, 2)
    cols = col_type(["C", "A", "B"])

    assert_eq(df[cols], ddf[cols])


def test_getitem_with_bool_dataframe_as_key():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    ddf = dd.from_pandas(df, 2)
    assert_eq(df[df > 3], ddf[ddf > 3])


def test_getitem_with_non_series():
    s = pd.Series(list(range(10)), index=list("abcdefghij"))
    ds = dd.from_pandas(s, npartitions=3)
    assert_eq(s[["a", "b"]], ds[["a", "b"]])


def test_ipython_completion():
    df = pd.DataFrame({"a": [1], "b": [2]})
    ddf = dd.from_pandas(df, npartitions=1)

    completions = ddf._ipython_key_completions_()
    assert "a" in completions
    assert "b" in completions
    assert "c" not in completions


def test_diff():
    df = pd.DataFrame(np.random.randn(100, 5), columns=list("abcde"))
    ddf = dd.from_pandas(df, 5)

    assert_eq(ddf.diff(), df.diff())
    assert_eq(ddf.diff(0), df.diff(0))
    assert_eq(ddf.diff(2), df.diff(2))
    assert_eq(ddf.diff(-2), df.diff(-2))

    assert_eq(ddf.diff(2, axis=1), df.diff(2, axis=1))

    assert_eq(ddf.a.diff(), df.a.diff())
    assert_eq(ddf.a.diff(0), df.a.diff(0))
    assert_eq(ddf.a.diff(2), df.a.diff(2))
    assert_eq(ddf.a.diff(-2), df.a.diff(-2))

    assert ddf.diff(2)._name == ddf.diff(2)._name
    assert ddf.diff(2)._name != ddf.diff(3)._name
    with pytest.raises((TypeError, ValueError)):
        ddf.diff(1.5)


def test_shift():
    df = _compat.makeTimeDataFrame()
    ddf = dd.from_pandas(df, npartitions=4)

    # DataFrame
    assert_eq(ddf.shift(), df.shift())
    assert_eq(ddf.shift(0), df.shift(0))
    assert_eq(ddf.shift(2), df.shift(2))
    assert_eq(ddf.shift(-2), df.shift(-2))

    assert_eq(ddf.shift(2, axis=1), df.shift(2, axis=1))

    # Series
    assert_eq(ddf.A.shift(), df.A.shift())
    assert_eq(ddf.A.shift(0), df.A.shift(0))
    assert_eq(ddf.A.shift(2), df.A.shift(2))
    assert_eq(ddf.A.shift(-2), df.A.shift(-2))

    with pytest.raises(TypeError):
        ddf.shift(1.5)


@pytest.mark.parametrize("data_freq,divs1", [("B", False), ("D", True), ("h", True)])
def test_shift_with_freq_DatetimeIndex(data_freq, divs1):
    df = _compat.makeTimeDataFrame()
    df = df.set_index(_compat.makeDateIndex(30, freq=data_freq))
    ddf = dd.from_pandas(df, npartitions=4)
    for freq, divs2 in [("s", True), ("W", False), (pd.Timedelta(10, unit="h"), True)]:
        for d, p in [(ddf, df), (ddf.A, df.A), (ddf.index, df.index)]:
            res = d.shift(2, freq=freq)
            assert_eq(res, p.shift(2, freq=freq))
            assert res.known_divisions == divs2
    # Index shifts also work with freq=None
    res = ddf.index.shift(2)
    assert_eq(res, df.index.shift(2))
    assert res.known_divisions == divs1


@pytest.mark.parametrize("data_freq,divs", [("D", True), ("h", True)])
def test_shift_with_freq_PeriodIndex(data_freq, divs):
    df = _compat.makeTimeDataFrame()
    # PeriodIndex
    df = df.set_index(pd.period_range("2000-01-01", periods=30, freq=data_freq))
    ddf = dd.from_pandas(df, npartitions=4)
    for d, p in [(ddf, df), (ddf.A, df.A)]:
        res = d.shift(2, freq=data_freq)
        assert_eq(res, p.shift(2, freq=data_freq))
        assert res.known_divisions == divs
    # PeriodIndex.shift doesn't have `freq` parameter
    res = ddf.index.shift(2)
    assert_eq(res, df.index.shift(2))
    assert res.known_divisions == divs

    with pytest.raises((ValueError, TypeError)):
        ddf.index.shift(2, freq="D")  # freq keyword not supported


def test_shift_with_freq_TimedeltaIndex():
    df = _compat.makeTimeDataFrame()
    # TimedeltaIndex
    for data_freq in ["min", "D", "h"]:
        df = df.set_index(_compat.makeTimedeltaIndex(30, freq=data_freq))
        ddf = dd.from_pandas(df, npartitions=4)
        for freq in ["s", pd.Timedelta(10, unit="h")]:
            for d, p in [(ddf, df), (ddf.A, df.A), (ddf.index, df.index)]:
                res = d.shift(2, freq=freq)
                assert_eq(res, p.shift(2, freq=freq))
                assert res.known_divisions
        # Index shifts also work with freq=None
        res = ddf.index.shift(2)
        assert_eq(res, df.index.shift(2))
        assert res.known_divisions


def test_shift_with_freq_errors():
    # Other index types error
    df = _compat.makeDataFrame()
    ddf = dd.from_pandas(df, npartitions=4)
    pytest.raises(NotImplementedError, lambda: ddf.shift(2, freq="s"))
    pytest.raises(NotImplementedError, lambda: ddf.A.shift(2, freq="s"))
    pytest.raises(NotImplementedError, lambda: ddf.index.shift(2))


@pytest.mark.parametrize("npartitions", [1, 4, 20])
@pytest.mark.parametrize("split_every", [2, 5])
@pytest.mark.parametrize("split_out", [1, 5, 20])
def test_hash_split_unique(npartitions, split_every, split_out):
    from string import ascii_lowercase

    s = pd.Series(np.random.choice(list(ascii_lowercase), 1000, replace=True))
    ds = dd.from_pandas(s, npartitions=npartitions)

    dropped = ds.unique(split_every=split_every, split_out=split_out)

    dsk = dropped.__dask_optimize__(dropped.dask, dropped.__dask_keys__())
    from dask.core import get_deps

    dependencies, dependents = get_deps(dsk)

    assert dropped.npartitions == (split_out or 1)
    assert sorted(dropped.compute(scheduler="sync")) == sorted(s.unique())


@pytest.mark.parametrize("split_every", [None, 2])
def test_split_out_drop_duplicates(split_every):
    x = np.concatenate([np.arange(10)] * 100)[:, None]
    y = x.copy()
    z = np.concatenate([np.arange(20)] * 50)[:, None]
    rs = np.random.RandomState(1)
    rs.shuffle(x)
    rs.shuffle(y)
    rs.shuffle(z)
    df = pd.DataFrame(np.concatenate([x, y, z], axis=1), columns=["x", "y", "z"])
    ddf = dd.from_pandas(df, npartitions=20)

    for subset, keep in product([None, ["x", "z"]], ["first", "last"]):
        sol = df.drop_duplicates(subset=subset, keep=keep)
        res = ddf.drop_duplicates(
            subset=subset, keep=keep, split_every=split_every, split_out=10
        )
        assert res.npartitions == 10
        assert_eq(sol, res)


@pytest.mark.parametrize("split_every", [None, 2])
def test_split_out_value_counts(split_every):
    df = pd.DataFrame({"x": [1, 2, 3] * 100})
    ddf = dd.from_pandas(df, npartitions=5)

    assert ddf.x.value_counts(split_out=10, split_every=split_every).npartitions == 10
    assert_eq(
        ddf.x.value_counts(split_out=10, split_every=split_every), df.x.value_counts()
    )


def test_values():
    from dask.array.utils import assert_eq

    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
        index=pd.Index([1.0, 2.0, 3.0, 4.0], name="ind"),
    )

    ddf = dd.from_pandas(df, 2)

    assert_eq(df.values, ddf.values)
    # When using pyarrow strings, we emit a warning about converting
    # pandas extension dtypes to object. Same as `test_values_extension_dtypes`.
    ctx = contextlib.nullcontext()
    if pyarrow_strings_enabled():
        ctx = pytest.warns(UserWarning, match="object dtype")
    with ctx:
        result = ddf.x.values
    assert_eq(df.x.values, result)
    assert_eq(df.y.values, ddf.y.values)
    assert_eq(df.index.values, ddf.index.values)


def test_values_extension_dtypes():
    from dask.array.utils import assert_eq

    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
        index=pd.Index([1.0, 2.0, 3.0, 4.0], dtype="Float64", name="ind"),
    )
    df = df.astype({"x": "string[python]", "y": "Int64"})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(df.values, ddf.values)

    with pytest.warns(UserWarning, match="object dtype"):
        result = ddf.x.values
    assert_eq(result, df.x.values.astype(object))

    with pytest.warns(UserWarning, match="object dtype"):
        result = ddf.y.values
    assert_eq(result, df.y.values.astype(object))

    with pytest.warns(UserWarning, match="object dtype"):
        result = ddf.index.values
    assert_eq(result, df.index.values.astype(object))


def test_copy():
    df = pd.DataFrame({"x": [1, 2, 3]})

    a = dd.from_pandas(df, npartitions=2)
    b = a.copy()
    c = a.copy(deep=False)

    a["y"] = a.x * 2

    assert_eq(b, df)
    assert_eq(c, df)

    deep_err = (
        "The `deep` value must be False. This is strictly a shallow copy "
        "of the underlying computational graph."
    )
    for deep in [True, None, ""]:
        with pytest.raises(ValueError, match=deep_err):
            a.copy(deep=deep)


def test_del():
    df = pd.DataFrame(
        {"x": ["a", "b", "c", "d"], "y": [2, 3, 4, 5]},
        index=pd.Index([1.0, 2.0, 3.0, 4.0], name="ind"),
    )
    a = dd.from_pandas(df, 2)
    b = a.copy()

    del a["x"]
    assert_eq(b, df)

    del df["x"]
    assert_eq(a, df)


@pytest.mark.parametrize("index", [True, False])
@pytest.mark.parametrize("deep", [True, False])
def test_memory_usage_dataframe(index, deep):
    df = pd.DataFrame(
        {"x": [1, 2, 3], "y": [1.0, 2.0, 3.0], "z": ["a", "b", "c"]},
        # Dask will use more memory for a multi-partition
        # RangeIndex, so we must set an index explicitly
        index=[1, 2, 3],
    )
    if pyarrow_strings_enabled():
        # pandas should measure memory usage of pyarrow strings
        df = to_pyarrow_string(df)
    ddf = dd.from_pandas(df, npartitions=2)
    expected = df.memory_usage(index=index, deep=deep)
    result = ddf.memory_usage(index=index, deep=deep)
    assert_eq(expected, result)


@pytest.mark.parametrize("index", [True, False])
@pytest.mark.parametrize("deep", [True, False])
def test_memory_usage_series(index, deep):
    s = pd.Series([1, 2, 3, 4], index=["a", "b", "c", "d"])
    if pyarrow_strings_enabled():
        # pandas should measure memory usage of pyarrow strings
        s = to_pyarrow_string(s)
    expected = s.memory_usage(index=index, deep=deep)
    ds = dd.from_pandas(s, npartitions=2)
    result = ds.memory_usage(index=index, deep=deep)
    assert_eq(expected, result)


@pytest.mark.parametrize("deep", [True, False])
def test_memory_usage_index(deep):
    s = pd.Series([1, 2, 3, 4], index=["a", "b", "c", "d"])
    if pyarrow_strings_enabled():
        # pandas should measure memory usage of pyarrow strings
        s = to_pyarrow_string(s)
    expected = s.index.memory_usage(deep=deep)
    ds = dd.from_pandas(s, npartitions=2)
    result = ds.index.memory_usage(deep=deep)
    assert_eq(expected, result)


@pytest.mark.parametrize("index", [True, False])
@pytest.mark.parametrize("deep", [True, False])
def test_memory_usage_per_partition(index, deep):
    df = pd.DataFrame(
        {
            "x": [1, 2, 3, 4, 5],
            "y": [1.0, 2.0, 3.0, 4.0, 5.0],
            "z": ["a", "b", "c", "d", "e"],
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)

    # DataFrame.memory_usage_per_partition
    expected = pd.Series(
        part.compute().memory_usage(index=index, deep=deep).sum()
        for part in ddf.partitions
    )
    result = ddf.memory_usage_per_partition(index=index, deep=deep)
    assert_eq(expected, result, check_index=False)

    # Series.memory_usage_per_partition
    expected = pd.Series(
        part.x.compute().memory_usage(index=index, deep=deep) for part in ddf.partitions
    )
    result = ddf.x.memory_usage_per_partition(index=index, deep=deep)
    assert_eq(expected, result, check_index=False)


@pytest.mark.parametrize(
    "reduction",
    [
        "sum",
        "mean",
        "std",
        "var",
        "count",
        "min",
        "max",
        "idxmin",
        "idxmax",
        "prod",
        "all",
        "sem",
    ],
)
def test_dataframe_reductions_arithmetic(reduction):
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5], "y": [1.1, 2.2, 3.3, 4.4, 5.5]})
    ddf = dd.from_pandas(df, npartitions=3)

    assert_eq(
        ddf - (getattr(ddf, reduction)() + 1), df - (getattr(df, reduction)() + 1)
    )


def test_median():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5], "y": [1.1, 2.2, 3.3, 4.4, 5.5]})
    ddf = dd.from_pandas(df, npartitions=3)

    # Exact medians work when `axis=1` or when there's only have a single partition
    assert_eq(ddf.median(axis=1), df.median(axis=1))
    ddf_single = dd.from_pandas(df, npartitions=1)
    assert_eq(ddf_single.median(axis=1), df.median(axis=1))
    assert_eq(ddf_single.median(axis=0), df.median(axis=0))
    assert_eq(ddf_single.x.median(), df.x.median())

    # Ensure `median` redirects to `median_approximate` appropriately
    for axis in [None, 0, "rows"]:
        with pytest.raises(
            NotImplementedError, match="See the `median_approximate` method instead"
        ):
            ddf.median(axis=axis)

    with pytest.raises(
        NotImplementedError, match="See the `median_approximate` method instead"
    ):
        ddf.x.median()


@pytest.mark.parametrize(
    "method",
    [
        "dask",
        pytest.param(
            "tdigest", marks=pytest.mark.skipif(not crick, reason="Requires crick")
        ),
    ],
)
def test_median_approximate(method):
    df = pd.DataFrame({"x": range(100), "y": range(100, 200)})
    ddf = dd.from_pandas(df, npartitions=10)
    assert_eq(
        ddf.median_approximate(method=method),
        df.median(),
        atol=1,
    )


def test_datetime_loc_open_slicing():
    dtRange = pd.date_range("01.01.2015", "05.05.2015")
    df = pd.DataFrame(np.random.random((len(dtRange), 2)), index=dtRange)
    ddf = dd.from_pandas(df, npartitions=5)
    assert_eq(df.loc[:"02.02.2015"], ddf.loc[:"02.02.2015"])
    assert_eq(df.loc["02.02.2015":], ddf.loc["02.02.2015":])
    assert_eq(df[0].loc[:"02.02.2015"], ddf[0].loc[:"02.02.2015"])
    assert_eq(df[0].loc["02.02.2015":], ddf[0].loc["02.02.2015":])


@pytest.mark.parametrize("gpu", [False, pytest.param(True, marks=pytest.mark.gpu)])
def test_to_datetime(gpu):
    xd = pd if not gpu else pytest.importorskip("cudf")

    # meta dtype is inconsistent for cuDF-backed frames
    check_dtype = not gpu

    df = xd.DataFrame({"year": [2015, 2016], "month": ["2", "3"], "day": [4, 5]})
    df.index.name = "ix"
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(xd.to_datetime(df), dd.to_datetime(ddf), check_dtype=check_dtype)
    assert_eq(xd.to_datetime(df), dd.to_datetime(df), check_dtype=check_dtype)

    s = xd.Series(
        ["3/11/2000", "3/12/2000", "3/13/2000"] * 100,
        index=["3/11/2000", "3/12/2000", "3/13/2000"] * 100,
    )
    ds = dd.from_pandas(s, npartitions=10, sort=False)

    # cuDF does not yet support timezone-aware datetimes
    if not gpu:
        assert_eq(
            xd.to_datetime(s, utc=True),
            dd.to_datetime(ds, utc=True),
        )
        assert_eq(
            xd.to_datetime(s, utc=True),
            dd.to_datetime(s, utc=True),
        )

    for arg in ("2021-08-03", 2021, s.index):
        with pytest.raises(NotImplementedError, match="non-index-able arguments"):
            dd.to_datetime(arg)


def test_to_timedelta():
    s = pd.Series(range(10))
    ds = dd.from_pandas(s, npartitions=2)

    assert_eq(pd.to_timedelta(s), dd.to_timedelta(ds))
    assert_eq(pd.to_timedelta(s, unit="h"), dd.to_timedelta(ds, unit="h"))

    s = pd.Series([1, 2, "this will error"])
    ds = dd.from_pandas(s, npartitions=2)
    assert_eq(pd.to_timedelta(s, errors="coerce"), dd.to_timedelta(ds, errors="coerce"))

    s = pd.Series(["1", 2, "1 day 2 hours"])
    ds = dd.from_pandas(s, npartitions=2)
    assert_eq(pd.to_timedelta(s), dd.to_timedelta(ds))

    with pytest.raises(
        ValueError, match="unit must not be specified if the input contains a str"
    ):
        dd.to_timedelta(ds, unit="s").compute()


@pytest.mark.parametrize("values", [[np.nan, 0], [1, 1]])
def test_isna(values):
    s = pd.Series(values)
    ds = dd.from_pandas(s, npartitions=2)

    assert_eq(pd.isna(s), dd.isna(ds))


@pytest.mark.parametrize("drop", [0, 9])
def test_slice_on_filtered_boundary(drop):
    # https://github.com/dask/dask/issues/2211
    x = np.arange(10)
    x[[5, 6]] -= 2
    df = pd.DataFrame({"A": x, "B": np.arange(len(x))})
    pdf = df.set_index("A").query(f"B != {drop}")
    ddf = dd.from_pandas(df, 1).set_index("A").query(f"B != {drop}")

    result = dd.concat([ddf, ddf.rename(columns={"B": "C"})], axis=1)
    expected = pd.concat([pdf, pdf.rename(columns={"B": "C"})], axis=1)
    assert_eq(result, expected)


def test_boundary_slice_nonmonotonic():
    x = np.array([-1, -2, 2, 4, 3])
    df = pd.DataFrame({"B": range(len(x))}, index=x)
    result = methods.boundary_slice(df, 0, 4)
    expected = df.iloc[2:]
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -1, 4)
    expected = df.drop(-2)
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -2, 3)
    expected = df.drop(4)
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -2, 3.5)
    expected = df.drop(4)
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -2, 4)
    expected = df
    tm.assert_frame_equal(result, expected)


def test_boundary_slice_empty():
    df = pd.DataFrame()
    result = methods.boundary_slice(df, 1, 4)
    expected = pd.DataFrame()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "start, stop, right_boundary, left_boundary, drop",
    [
        (-1, None, False, False, [-1, -2]),
        (-1, None, False, True, [-2]),
        (None, 3, False, False, [3, 4]),
        (None, 3, True, False, [4]),
        # Missing keys
        (-0.5, None, False, False, [-1, -2]),
        (-0.5, None, False, True, [-1, -2]),
        (-1.5, None, False, True, [-2]),
        (None, 3.5, False, False, [4]),
        (None, 3.5, True, False, [4]),
        (None, 2.5, False, False, [3, 4]),
    ],
)
def test_with_boundary(start, stop, right_boundary, left_boundary, drop):
    x = np.array([-1, -2, 2, 4, 3])
    df = pd.DataFrame({"B": range(len(x))}, index=x)
    result = methods.boundary_slice(df, start, stop, right_boundary, left_boundary)
    expected = df.drop(drop)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "index, left, right",
    [
        (range(10), 0, 9),
        (range(10), -1, None),
        (range(10), None, 10),
        ([-1, 0, 2, 1], None, None),
        ([-1, 0, 2, 1], -1, None),
        ([-1, 0, 2, 1], None, 2),
        ([-1, 0, 2, 1], -2, 3),
        (pd.date_range("2017", periods=10), None, None),
        (pd.date_range("2017", periods=10), pd.Timestamp("2017"), None),
        (pd.date_range("2017", periods=10), None, pd.Timestamp("2017-01-10")),
        (pd.date_range("2017", periods=10), pd.Timestamp("2016"), None),
        (pd.date_range("2017", periods=10), None, pd.Timestamp("2018")),
    ],
)
def test_boundary_slice_same(index, left, right):
    df = pd.DataFrame({"A": range(len(index))}, index=index)
    result = methods.boundary_slice(df, left, right)
    tm.assert_frame_equal(result, df)


def test_better_errors_object_reductions():
    # GH2452
    s = pd.Series(["a", "b", "c", "d"])
    ds = dd.from_pandas(s, npartitions=2)
    with pytest.raises(ValueError) as err:
        ds.mean()
    assert str(err.value) == f"`mean` not supported with {ds.dtype} series"


def test_sample_empty_partitions():
    @dask.delayed
    def make_df(n):
        return pd.DataFrame(np.zeros((n, 4)), columns=list("abcd"))

    ddf = dd.from_delayed([make_df(0), make_df(100), make_df(0)])
    ddf2 = ddf.sample(frac=0.2)
    # smoke test sample on empty partitions
    res = ddf2.compute()
    assert res.dtypes.equals(ddf2.dtypes)


def test_coerce():
    df = pd.DataFrame(np.arange(100).reshape((10, 10)))
    ddf = dd.from_pandas(df, npartitions=2)
    funcs = (int, float, complex)
    for d, t in product(funcs, (ddf, ddf[0])):
        pytest.raises(TypeError, lambda t=t, d=d: t(d))


def test_bool():
    df = pd.DataFrame(np.arange(100).reshape((10, 10)))
    ddf = dd.from_pandas(df, npartitions=2)
    conditions = [ddf, ddf[0], ddf == ddf, ddf[0] == ddf[0]]
    for cond in conditions:
        with pytest.raises(ValueError):
            bool(cond)


def test_cumulative_multiple_columns():
    # GH 3037
    df = pd.DataFrame(np.random.randn(100, 5), columns=list("abcde"))
    ddf = dd.from_pandas(df, 5)

    for d in [ddf, df]:
        for c in df.columns:
            d[c + "cs"] = d[c].cumsum()
            d[c + "cmin"] = d[c].cummin()
            d[c + "cmax"] = d[c].cummax()
            d[c + "cp"] = d[c].cumprod()

    assert_eq(ddf, df)


@pytest.mark.parametrize("func", [np.asarray, M.to_records])
def test_map_partition_array(func):
    from dask.array.utils import assert_eq

    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5], "y": [6.0, 7.0, 8.0, 9.0, 10.0]},
        index=["a", "b", "c", "d", "e"],
    )
    ddf = dd.from_pandas(df, npartitions=2)

    for pre in [lambda a: a, lambda a: a.x, lambda a: a.y, lambda a: a.index]:
        try:
            expected = func(pre(df))
        except Exception:
            continue
        x = pre(ddf).map_partitions(func)
        assert_eq(x, expected, check_type=False)  # TODO: make check_type pass

        assert isinstance(x, da.Array)
        assert x.chunks[0] == (np.nan, np.nan)


def test_map_partition_sparse():
    sparse = pytest.importorskip("sparse")
    # Avoid searchsorted failure.
    pytest.importorskip("numba", minversion="0.40.0")

    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5], "y": [6.0, 7.0, 8.0, 9.0, 10.0]},
        index=["a", "b", "c", "d", "e"],
    )
    ddf = dd.from_pandas(df, npartitions=2)

    def f(d):
        return sparse.COO(np.array(d))

    for pre in [lambda a: a, lambda a: a.x]:
        expected = f(pre(df))
        result = pre(ddf).map_partitions(f)
        assert isinstance(result, da.Array)
        computed = result.compute()
        assert (computed.data == expected.data).all()
        assert (computed.coords == expected.coords).all()


def test_mixed_dask_array_operations():
    df = pd.DataFrame({"x": [1, 2, 3]}, index=[4, 5, 6])
    ddf = dd.from_pandas(df, npartitions=2)
    assert (ddf.x + ddf.x.values).npartitions >= (
        ddf.x + ddf.x.values
    ).optimize().npartitions

    assert_eq(df.x + df.x.values, ddf.x + ddf.x.values)
    assert_eq(df.x.values + df.x, ddf.x.values + ddf.x)

    assert_eq(df.x + df.index.values, ddf.x + ddf.index.values)
    assert_eq(df.index.values + df.x, ddf.index.values + ddf.x)


def test_mixed_dask_array_operations_errors():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]}, index=[4, 5, 6, 7, 8])
    ddf = dd.from_pandas(df, npartitions=2)

    x = da.arange(5, chunks=((1, 4),))
    x._chunks = ((np.nan, np.nan),)

    with pytest.raises(ValueError):
        (ddf.x + x).compute()

    x = da.arange(5, chunks=((2, 2, 1),))
    with pytest.raises(ValueError) as info:
        ddf.x + x

    assert "add" in str(info.value) or "different numbers" in str(info.value)


def test_mixed_dask_array_multi_dimensional():
    df = pd.DataFrame(
        {"x": [1, 2, 3, 4, 5], "y": [5.0, 6.0, 7.0, 8.0, 9.0]}, columns=["x", "y"]
    )
    ddf = dd.from_pandas(df, npartitions=2)

    x = (df.values + 1).astype(float)
    dx = (ddf.values + 1).astype(float)

    assert_eq(ddf + dx + 1, df + x + 1)
    assert_eq(ddf + dx.rechunk((None, 1)) + 1, df + x + 1)
    assert_eq(ddf[["y", "x"]] + dx + 1, df[["y", "x"]] + x + 1)


def test_meta_raises():
    # Raise when we use a user defined function
    s = pd.Series(["abcd", "abcd"])
    ds = dd.from_pandas(s, npartitions=2)
    try:
        ds.map(lambda x: x[3])
    except ValueError as e:
        assert "meta=" in str(e)


@pytest.mark.skip_with_pyarrow_strings  # DateOffset has to be an object
def test_meta_nonempty_uses_meta_value_if_provided():
    # https://github.com/dask/dask/issues/6958
    base = pd.Series([1, 2, 3], dtype="datetime64[ns]")
    offsets = pd.Series([pd.offsets.DateOffset(years=o) for o in range(3)])
    dask_base = dd.from_pandas(base, npartitions=1)
    dask_offsets = dd.from_pandas(offsets, npartitions=1)

    with warnings.catch_warnings():  # not vectorized performance warning
        warnings.simplefilter("ignore", PerformanceWarning)
        warnings.simplefilter("ignore", UserWarning)
        expected = base + offsets
        actual = dask_base + dask_offsets
        assert_eq(expected, actual)


@pytest.mark.parametrize("container", ["array", "matrix"])
def test_dask_dataframe_holds_scipy_sparse_containers(container):
    sparse = pytest.importorskip("scipy.sparse")
    da = pytest.importorskip("dask.array")
    if container == "array" and not hasattr(sparse, "csr_array"):
        pytest.skip("scipy<1.11 has no sparray")
    cls = sparse.csr_matrix if container == "matrix" else sparse.csr_array

    x = da.random.random((1000, 10), chunks=(100, 10))
    x[x < 0.9] = 0
    df = dd.from_dask_array(x)
    y = df.map_partitions(cls)

    assert isinstance(y, da.Array)

    vs = y.to_delayed().flatten().tolist()
    values = dask.compute(*vs, scheduler="single-threaded")
    assert all(isinstance(v, cls) for v in values)


@pytest.mark.xfail(reason="we can't do this yet")
def test_map_partitions_delays_large_inputs():
    df = pd.DataFrame({"x": [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)

    big = np.ones(1000000)

    b = ddf.map_partitions(lambda x, y: x, y=big)
    assert any(big is v for v in b.dask.values())

    a = ddf.map_partitions(lambda x, y: x, big)
    assert any(big is v for v in a.dask.values())


def test_partitions_indexer():
    df = pd.DataFrame({"x": range(10)})
    ddf = dd.from_pandas(df, npartitions=5)

    assert_eq(ddf.partitions[0], ddf.get_partition(0))
    assert_eq(ddf.partitions[3], ddf.get_partition(3))
    assert_eq(ddf.partitions[-1], ddf.get_partition(4))

    assert ddf.partitions[:3].npartitions == 3
    assert ddf.x.partitions[:3].npartitions == 3

    assert ddf.x.partitions[::2].compute().tolist() == [0, 1, 4, 5, 8, 9]


def test_mod_eq():
    df = pd.DataFrame({"a": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=1)
    assert_eq(df, ddf)
    assert_eq(df.a, ddf.a)
    assert_eq(df.a + 2, ddf.a + 2)
    assert_eq(df.a + 2 == 0, ddf.a + 2 == 0)


def test_setitem():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    ddf = dd.from_pandas(df.copy(), 2)
    df[df.columns] = 1
    ddf[ddf.columns] = 1
    assert_eq(df, ddf)


def test_setitem_with_bool_dataframe_as_key():
    df = pd.DataFrame({"A": [1, 4], "B": [3, 2]})
    ddf = dd.from_pandas(df.copy(), 2)
    df[df > 2] = 5
    ddf[ddf > 2] = 5
    assert_eq(df, ddf)


def test_setitem_with_bool_series_as_key():
    df = pd.DataFrame({"A": [1, 4], "B": [3, 2]})
    ddf = dd.from_pandas(df.copy(), 2)
    df[df["A"] > 2] = 5
    ddf[ddf["A"] > 2] = 5
    assert_eq(df, ddf)


def test_setitem_with_numeric_column_name_raises_not_implemented():
    df = pd.DataFrame({0: [1, 4], 1: [3, 2]})
    ddf = dd.from_pandas(df.copy(), 2)
    # works for pandas
    df[0] = 5
    # raises error for dask
    with pytest.raises(NotImplementedError, match="not supported"):
        ddf[0] = 5


def test_broadcast():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf - (ddf.sum() + 1), df - (df.sum() + 1))


def test_scalar_with_array():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)

    da.utils.assert_eq(df.x.values + df.x.mean(), ddf.x.values + ddf.x.mean())


def test_has_parallel_type():
    assert has_parallel_type(pd.DataFrame())
    assert has_parallel_type(pd.Series(dtype=float))
    assert not has_parallel_type(123)


def test_map_index():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert ddf.known_divisions is True
    with pytest.warns(UserWarning, match="meta"):
        cleared = ddf.index.map(lambda x: x * 10)
    assert cleared.known_divisions is False

    with pytest.warns(UserWarning, match="meta"):
        applied = ddf.index.map(lambda x: x * 10, is_monotonic=True)
    assert applied.known_divisions is True
    assert applied.divisions == tuple(x * 10 for x in ddf.divisions)


def test_assign_index():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)

    ddf_copy = ddf.copy()

    ddf.index = ddf.index * 10

    expected = df.copy()
    expected.index = expected.index * 10

    assert_eq(ddf, expected)
    assert_eq(ddf_copy, df)


def test_index_divisions():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(ddf.index + 1, df.index + 1)
    assert_eq(10 * ddf.index, 10 * df.index)
    assert_eq(-ddf.index, -df.index)


def test_replace():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(df.replace(1, 10), ddf.replace(1, 10))
    assert_eq(df.replace({1: 10, 2: 20}), ddf.replace({1: 10, 2: 20}))
    assert_eq(df.x.replace(1, 10), ddf.x.replace(1, 10))
    assert_eq(df.x.replace({1: 10, 2: 20}), ddf.x.replace({1: 10, 2: 20}))


@pytest.mark.xfail(reason="we can't do this yet")
def test_map_partitions_delays_lists():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)

    L = list(range(100))
    out = ddf.map_partitions(lambda x, y: x + sum(y), y=L)
    assert any(str(L) == str(v) for v in out.__dask_graph__().values())

    out = ddf.map_partitions(lambda x, y: x + sum(y), L)
    assert any(str(L) == str(v) for v in out.__dask_graph__().values())


def test_dtype_cast():
    df = pd.DataFrame(
        {
            "A": np.arange(10, dtype=np.int32),
            "B": np.arange(10, dtype=np.int64),
            "C": np.arange(10, dtype=np.float32),
        }
    )
    ddf = dd.from_pandas(df, npartitions=2)
    assert ddf.A.dtype == np.int32
    assert ddf.B.dtype == np.int64
    assert ddf.C.dtype == np.float32

    col = pd.Series(np.arange(10, dtype=np.float32)) / 2
    assert col.dtype == np.float32

    ddf = ddf.assign(D=col)
    assert ddf.D.dtype == np.float32
    assert ddf.C.dtype == np.float32
    # fails
    assert ddf.B.dtype == np.int64
    # fails
    assert ddf.A.dtype == np.int32


@pytest.mark.parametrize("base_npart", [1, 4])
@pytest.mark.parametrize("map_npart", [1, 3])
@pytest.mark.parametrize("sorted_index", [False, True])
@pytest.mark.parametrize("sorted_map_index", [False, True])
def test_series_map(base_npart, map_npart, sorted_index, sorted_map_index):
    if map_npart != base_npart:
        pytest.xfail(reason="not yet implemented")
    base = pd.Series(
        ["".join(np.random.choice(["a", "b", "c"], size=3)) for x in range(100)]
    )
    if not sorted_index:
        index = np.arange(100)
        np.random.shuffle(index)
        base.index = index
    map_index = ["".join(x) for x in product("abc", repeat=3)]
    mapper = pd.Series(np.random.randint(50, size=len(map_index)), index=map_index)
    if not sorted_map_index:
        map_index = np.array(map_index)
        np.random.shuffle(map_index)
        mapper.index = map_index
    expected = base.map(mapper)
    dask_base = dd.from_pandas(base, npartitions=base_npart, sort=False)
    dask_map = dd.from_pandas(mapper, npartitions=map_npart, sort=False)
    with pytest.warns(UserWarning, match="meta"):
        result = dask_base.map(dask_map)
    assert_eq(expected, result)


@pytest.mark.skip_with_pyarrow_strings  # has to be array to explode
def test_dataframe_explode():
    df = pd.DataFrame({"A": [[1, 2, 3], "foo", [3, 4]], "B": 1})
    exploded_df = df.explode("A")

    ddf = dd.from_pandas(df, npartitions=2)

    exploded_ddf = ddf.explode("A")
    assert ddf.divisions == exploded_ddf.divisions
    assert_eq(exploded_ddf.compute(), exploded_df)


@pytest.mark.skip_with_pyarrow_strings  # has to be array to explode
def test_series_explode():
    s = pd.Series([[1, 2, 3], "foo", [3, 4]])
    exploded_s = s.explode()
    ds = dd.from_pandas(s, npartitions=2)
    exploded_ds = ds.explode()
    assert_eq(exploded_ds, exploded_s)
    assert ds.divisions == exploded_ds.divisions


def test_pop():
    df = pd.DataFrame({"x": range(10), "y": range(10)})

    ddf = dd.from_pandas(df, npartitions=2)

    s = ddf.pop("y")
    assert s.name == "y"
    assert ddf.columns == ["x"]
    assert_eq(ddf, df[["x"]])


@pytest.mark.parametrize("dropna", [True, False])
@pytest.mark.parametrize("axis", [0, 1])
def test_nunique(dropna, axis):
    df = pd.DataFrame(
        {"x": ["a", "a", "c"], "y": [None, 1, 2], "c": np.arange(0, 1, 0.4)}
    )
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf["y"].nunique(dropna=dropna), df["y"].nunique(dropna=dropna))
    assert_eq(
        ddf.nunique(dropna=dropna, axis=axis), df.nunique(dropna=dropna, axis=axis)
    )


def test_iter():
    df = pd.DataFrame({"A": [1, 2, 3, 4], "B": [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, 2)

    assert list(df) == list(ddf)
    for col, expected in zip(ddf, ["A", "B"]):
        assert col == expected


def test_dataframe_groupby_cumsum_agg_empty_partitions():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=4)
    assert_eq(ddf[ddf.x < 5].x.cumsum(), df[df.x < 5].x.cumsum())
    assert_eq(ddf[ddf.x > 5].x.cumsum(), df[df.x > 5].x.cumsum())


def test_dataframe_groupby_cumprod_agg_empty_partitions():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=4)
    assert_eq(ddf[ddf.x < 5].x.cumprod(), df[df.x < 5].x.cumprod())
    assert_eq(ddf[ddf.x > 5].x.cumprod(), df[df.x > 5].x.cumprod())


@pytest.mark.skipif(reason="not important now")
def test_attrs_dataframe():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4], "C": [5, 6]})
    df.attrs = {"date": "2020-10-16"}
    ddf = dd.from_pandas(df, 2)

    assert df.attrs == ddf.attrs
    assert df.abs().attrs == ddf.abs().attrs


@pytest.mark.skipif(reason="not important now")
def test_attrs_series():
    s = pd.Series([1, 2], name="A")
    s.attrs["unit"] = "kg"
    ds = dd.from_pandas(s, 2)

    assert s.attrs == ds.attrs
    assert s.fillna(1).attrs == ds.fillna(1).attrs


def test_join_series():
    df = pd.DataFrame({"x": [1, 2, 3, 4, 5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=1)
    expected_df = dd.from_pandas(df.join(df["x"], lsuffix="_"), npartitions=1)
    actual_df = ddf.join(ddf["x"], lsuffix="_")
    assert_eq(actual_df, expected_df)


def test_assign_na_float_columns():
    # See https://github.com/dask/dask/issues/7156
    df_pandas = pd.DataFrame({"a": [1.1]}, dtype="Float64")
    df = dd.from_pandas(df_pandas, npartitions=1)

    df = df.assign(new_col=df["a"])

    assert df.compute()["a"].dtypes == "Float64"
    assert df.compute()["new_col"].dtypes == "Float64"


def test_assign_no_warning_fragmented():
    df = pd.DataFrame({"a": [1, 2, 3, 4, 5] * 10})
    df = dd.from_pandas(df, npartitions=50)
    with warnings.catch_warnings(record=True) as w:
        for i in range(105):
            df[str(i)] = 5
        assert len(w) == 0


def test_dot():
    s1 = pd.Series([1, 2, 3, 4])
    s2 = pd.Series([4, 5, 6, 6])
    df = pd.DataFrame({"one": s1, "two": s2})

    dask_s1 = dd.from_pandas(s1, npartitions=1)
    dask_df = dd.from_pandas(df, npartitions=1)
    dask_s2 = dd.from_pandas(s2, npartitions=1)

    assert_eq(s1.dot(s2), dask_s1.dot(dask_s2))
    assert_eq(s1.dot(df), dask_s1.dot(dask_df))

    # With partitions
    partitioned_s1 = dd.from_pandas(s1, npartitions=2)
    partitioned_df = dd.from_pandas(df, npartitions=2)
    partitioned_s2 = dd.from_pandas(s2, npartitions=2)

    assert_eq(s1.dot(s2), partitioned_s1.dot(partitioned_s2))
    assert_eq(s1.dot(df), partitioned_s1.dot(partitioned_df))

    # Test passing meta kwarg
    res = dask_s1.dot(dask_df, meta=pd.Series([1], name="test_series")).compute()
    assert res.name == "test_series"

    # Test validation of second operand
    with pytest.raises(TypeError):
        dask_s1.dot(da.array([1, 2, 3, 4]))


def test_dot_nan():
    # Test that nan inputs match pandas' behavior
    s1 = pd.Series([1, 2, 3, 4])
    dask_s1 = dd.from_pandas(s1, npartitions=1)

    s2 = pd.Series([np.nan, np.nan, np.nan, np.nan])
    dask_s2 = dd.from_pandas(s2, npartitions=1)

    df = pd.DataFrame({"one": s1, "two": s2})
    dask_df = dd.from_pandas(df, npartitions=1)

    assert_eq(s1.dot(s2), dask_s1.dot(dask_s2))
    assert_eq(s2.dot(df), dask_s2.dot(dask_df))


def test_use_of_weakref_proxy():
    """Testing wrapping frames in proxy wrappers"""
    df = pd.DataFrame({"data": [1, 2, 3]})
    df_pxy = weakref.proxy(df)
    ser = pd.Series({"data": [1, 2, 3]})
    ser_pxy = weakref.proxy(ser)

    assert is_dataframe_like(df_pxy)
    assert is_series_like(ser_pxy)

    assert dask.dataframe.groupby._cov_chunk(df_pxy, "data")
    assert isinstance(
        dask.dataframe.groupby._groupby_apply_funcs(df_pxy, "data", funcs=[]),
        pd.DataFrame,
    )

    # Test wrapping each Dask dataframe chunk in a proxy
    l = []

    def f(x):
        l.append(x)  # Keep `x` alive
        return weakref.proxy(x)

    d = pd.DataFrame({"g": [0, 0, 1] * 3, "b": [1, 2, 3] * 3})
    a = dd.from_pandas(d, npartitions=1)
    a = a.map_partitions(f, meta=a._meta)
    pxy = weakref.proxy(a)
    res = pxy["b"].groupby(pxy["g"]).sum()
    isinstance(res.compute(), pd.Series)


@pytest.mark.parametrize(
    "series",
    [
        [0, 1, 0, 0, 1, 0],  # not monotonic
        [0, 1, 1, 2, 2, 3],  # monotonic increasing
        [0, 1, 2, 3, 4, 5],  # strictly monotonic increasing
        [0, 0, 0, 0, 0, 0],  # both monotonic increasing and monotonic decreasing
        [0, 1, 2, 0, 1, 2],  # Partitions are individually monotonic; whole series isn't
    ],
)
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("cls", ["Series", "Index"])
def test_is_monotonic_numeric(series, reverse, cls):
    if reverse:
        series = series[::-1]
    pds = pd.Series(series, index=series)
    ds = dd.from_pandas(pds, chunksize=3, sort=False)
    if cls == "Index":
        pds = pds.index
        ds = ds.index

    assert ds.is_monotonic_increasing.compute() == pds.is_monotonic_increasing
    assert ds.is_monotonic_decreasing.compute() == pds.is_monotonic_decreasing


def test_is_monotonic_dt64():
    s = pd.Series(pd.date_range("20130101", periods=10))
    ds = dd.from_pandas(s, npartitions=5)
    assert_eq(s.is_monotonic_increasing, ds.is_monotonic_increasing)
    assert_eq(s.is_monotonic_decreasing, ds.is_monotonic_decreasing)

    s_2 = pd.Series(list(reversed(s)))
    ds_2 = dd.from_pandas(s_2, npartitions=5)
    assert_eq(s_2.is_monotonic_increasing, ds_2.is_monotonic_increasing)
    assert_eq(s_2.is_monotonic_decreasing, ds_2.is_monotonic_decreasing)


def test_index_is_monotonic_dt64():
    s = pd.Series(1, index=pd.date_range("20130101", periods=10))
    ds = dd.from_pandas(s, npartitions=5, sort=False)
    assert_eq(s.index.is_monotonic_increasing, ds.index.is_monotonic_increasing)
    assert_eq(s.index.is_monotonic_decreasing, ds.index.is_monotonic_decreasing)

    s_2 = s[::-1]
    ds_2 = dd.from_pandas(s_2, npartitions=5, sort=False)
    assert_eq(s_2.index.is_monotonic_increasing, ds_2.index.is_monotonic_increasing)
    assert_eq(s_2.index.is_monotonic_decreasing, ds_2.index.is_monotonic_decreasing)


def test_is_monotonic_empty_partitions():
    df = pd.DataFrame({"a": [1, 2, 3, 4], "b": [4, 3, 2, 1]})
    ddf = dd.from_pandas(df, npartitions=2)

    # slice it to get empty partitions
    df = df[df["a"] >= 3]
    ddf = ddf[ddf["a"] >= 3]
    assert_eq(df["a"].is_monotonic_increasing, ddf["a"].is_monotonic_increasing)
    assert_eq(df.index.is_monotonic_increasing, ddf.index.is_monotonic_increasing)
    assert_eq(df["a"].is_monotonic_decreasing, ddf["a"].is_monotonic_decreasing)
    assert_eq(df.index.is_monotonic_decreasing, ddf.index.is_monotonic_decreasing)


@pytest.mark.parametrize("dtype", [int, float])
@pytest.mark.parametrize("orient", ["columns", "index"])
@pytest.mark.parametrize("npartitions", [2, 5])
def test_from_dict(dtype, orient, npartitions):
    data = {"a": range(10), "b": range(10)}
    expected = pd.DataFrame.from_dict(data, dtype=dtype, orient=orient)
    result = dd.DataFrame.from_dict(
        data, npartitions=npartitions, dtype=dtype, orient=orient
    )
    if orient == "index":
        # DataFrame only has two rows with this orientation
        assert result.npartitions == min(npartitions, 2)
    else:
        assert result.npartitions == npartitions
    assert_eq(result, expected)


def test_from_dict_raises():
    s = pd.Series(range(10))
    ds = dd.from_pandas(s, npartitions=2)
    with pytest.raises(NotImplementedError, match="Dask collections as inputs"):
        dd.DataFrame.from_dict({"a": ds}, npartitions=2)


def test_empty():
    with pytest.raises(NotImplementedError, match="may be expensive"):
        d.empty
    with pytest.raises(AttributeError, match="may be expensive"):
        d.empty


@pytest.mark.parametrize(
    "dtype",
    [
        "int64[pyarrow]",
        "int32[pyarrow]",
        "float64[pyarrow]",
        "float32[pyarrow]",
        "uint8[pyarrow]",
    ],
)
def test_pyarrow_extension_dtype(dtype):
    # Ensure simple Dask DataFrame operations work with pyarrow extension dtypes
    pytest.importorskip("pyarrow")
    df = pd.DataFrame({"x": range(10)}, dtype=dtype)
    ddf = dd.from_pandas(df, npartitions=3)
    expected = (df.x + df.x) * 2
    result = (ddf.x + ddf.x) * 2
    assert_eq(expected, result)


def test_pyarrow_decimal_extension_dtype():
    # Similar to `test_pyarrow_extension_dtype` but for pyarrow decimal dtypes
    pa = pytest.importorskip("pyarrow")
    pa_dtype = pa.decimal128(precision=7, scale=3)

    data = pa.array(
        [
            decimal.Decimal("8093.234"),
            decimal.Decimal("8094.234"),
            decimal.Decimal("8095.234"),
            decimal.Decimal("8096.234"),
            decimal.Decimal("8097.234"),
            decimal.Decimal("8098.234"),
        ],
        type=pa_dtype,
    )
    df = pd.DataFrame({"x": data}, dtype=pd.ArrowDtype(pa_dtype))
    ddf = dd.from_pandas(df, npartitions=3)
    expected = (df.x + df.x) * 2
    result = (ddf.x + ddf.x) * 2
    assert_eq(expected, result)


def test_to_backend():
    # Test that `DataFrame.to_backend` works as expected
    with dask.config.set({"dataframe.backend": "pandas"}):
        # Start with pandas-backed data
        df = dd.from_dict({"a": range(10)}, npartitions=2)
        assert isinstance(df._meta, pd.DataFrame)

        # Default `to_backend` shouldn't change data
        assert_eq(df, df.to_backend())

        # Moving to a "missing" backend should raise an error
        with pytest.raises(ValueError, match="No backend dispatch registered"):
            df.to_backend("missing")


@pytest.mark.parametrize("func", ["max", "sum"])
def test_transform_getitem_works(func):
    df = pd.DataFrame({"ints": [1, 2, 3], "grouper": [0, 1, 0]})

    ddf = dd.from_pandas(df, npartitions=2)

    # what happens here is not exactly a transform, but a reduction
    # the result of which is broadcasted back to the original shape.
    # Broadcasting an aggregation in this manner is very performant in pandas
    meta = df.groupby("grouper").transform(func)
    df["new"] = df.groupby("grouper").transform(func)["ints"]
    ddf["new"] = ddf.groupby("grouper").transform(func, meta=meta)["ints"]

    assert_eq(df, ddf)


@pytest.mark.parametrize(
    "df,cond",
    [
        (pd.DataFrame({"x": [1, 2]}, index=[1, 2]), [[True], [False]]),
        (pd.DataFrame({"x": [1, 2], "y": [3, 4]}), [[True, False], [True, False]]),
        (pd.DataFrame({"x": [1, 2], "y": [3, 4]}), [[True, True], [False, False]]),
        (
            pd.DataFrame({"x": [1, 2, 3, 4], "y": [3, 4, 5, 6]}),
            [[True, True], [True, True], [False, False], [False, False]],
        ),
        (
            pd.DataFrame({"x": [1, 2, 3, 4], "y": [3, 4, 5, 6]}),
            [[True, False], [True, False], [True, False], [True, False]],
        ),
    ],
)
def test_mask_where_array_like(df, cond):
    """DataFrame.mask fails for single-row partitions
    https://github.com/dask/dask/issues/9848
    """
    ddf = dd.from_pandas(df, npartitions=2)

    # ensure raises when list is provided
    with pytest.raises(ValueError, match="can be aligned|shape"):
        ddf.mask(cond=cond, other=5)

    with pytest.raises(ValueError, match="can be aligned|shape"):
        ddf.where(cond=cond, other=5)

    # but works when DataFrame is provided, with matching index
    dd_cond = pd.DataFrame(cond, index=df.index, columns=df.columns)
    expected = df.mask(cond=cond, other=5)
    result = ddf.mask(cond=dd_cond, other=5)
    assert_eq(expected, result)

    expected = df.where(cond=cond, other=5)
    result = ddf.where(cond=dd_cond, other=5)
    assert_eq(expected, result)


def test_mask_where_callable():
    """https://github.com/dask/dask/issues/10282"""
    pdf = pd.DataFrame({"x": [1, None]})
    ddf = dd.from_pandas(pdf, npartitions=1)

    # dataframe
    assert_eq(pdf.where(lambda d: d.isna(), 3), ddf.where(lambda d: d.isna(), 3))

    #  series
    assert_eq(pdf.x.where(lambda d: d == 1, 2), ddf.x.where(lambda d: d == 1, 2))


def test_pyarrow_schema_dispatch():
    from dask.dataframe.dispatch import (
        pyarrow_schema_dispatch,
        to_pyarrow_table_dispatch,
    )

    pytest.importorskip("pyarrow")

    df = pd.DataFrame(np.random.randn(10, 3), columns=list("abc"))
    df["d"] = pd.Series(["cat", "dog"] * 5, dtype="string[pyarrow]")
    table = to_pyarrow_table_dispatch(df)
    schema = pyarrow_schema_dispatch(df)

    assert schema.equals(table.schema)


@pytest.mark.parametrize("preserve_index", [True, False])
def test_pyarrow_schema_dispatch_preserves_index(preserve_index):
    from dask.dataframe.dispatch import (
        pyarrow_schema_dispatch,
        to_pyarrow_table_dispatch,
    )

    pytest.importorskip("pyarrow")

    df = pd.DataFrame(np.random.randn(10, 3), columns=list("abc"))
    df["d"] = pd.Series(["cat", "dog"] * 5, dtype="string[pyarrow]")
    table = to_pyarrow_table_dispatch(df, preserve_index=preserve_index)
    schema = pyarrow_schema_dispatch(df, preserve_index=preserve_index)

    assert schema.equals(table.schema)


@pytest.mark.parametrize("self_destruct", [True, False])
def test_pyarrow_conversion_dispatch(self_destruct):
    from dask.dataframe.dispatch import (
        from_pyarrow_table_dispatch,
        to_pyarrow_table_dispatch,
    )

    pytest.importorskip("pyarrow")

    df1 = pd.DataFrame(np.random.randn(10, 3), columns=list("abc"))
    df1["d"] = pd.Series(["cat", "dog"] * 5, dtype="string[pyarrow]")
    df2 = from_pyarrow_table_dispatch(
        df1,
        to_pyarrow_table_dispatch(df1),
        self_destruct=self_destruct,
    )

    assert type(df1) == type(df2)
    assert_eq(df1, df2)


@pytest.mark.gpu
def test_pyarrow_conversion_dispatch_cudf():
    from dask.dataframe.dispatch import (
        from_pyarrow_table_dispatch,
        to_pyarrow_table_dispatch,
    )

    cudf = pytest.importorskip("cudf")

    df1 = cudf.DataFrame(np.random.randn(10, 3), columns=list("abc"))
    df2 = from_pyarrow_table_dispatch(df1, to_pyarrow_table_dispatch(df1))

    assert type(df1) == type(df2)
    assert_eq(df1, df2)


def test_enforce_runtime_divisions():
    pdf = pd.DataFrame({"x": range(50)})
    ddf = dd.from_pandas(pdf, 5)
    # Default divisions should be correct
    assert_eq(pdf, ddf.enforce_runtime_divisions())


def test_preserve_ts_unit_in_meta_creation():
    pdf = pd.DataFrame(
        {
            "a": [1],
            "timestamp": pd.Series(
                [pd.Timestamp.now("UTC")], dtype="datetime64[us, UTC]"
            ),
        }
    )
    df = dd.from_pandas(pdf, npartitions=1)
    assert_eq(meta_nonempty(df._meta).dtypes, pdf.dtypes)
    assert_eq(df, pdf)


def test_from_xarray():
    xr = pytest.importorskip("xarray")

    a = da.zeros(shape=(6000, 45000), chunks=(30, 6000))
    a = xr.DataArray(
        a, coords={"c1": list(range(6000)), "c2": list(range(45000))}, dims=["c1", "c2"]
    )
    sel_coords = {
        "c1": pd.date_range("2024-11-01", "2024-11-10").to_numpy(),
        "c2": [29094],
    }

    result = a.reindex(**sel_coords).to_dask_dataframe()
    assert len(collections_to_dsk([result])) < 30  # previously 3000


def test_from_xarray_string_conversion():
    xr = pytest.importorskip("xarray")

    x = np.random.randn(10)
    y = np.arange(10, dtype="uint8")
    t = list("abcdefghij")
    ds = xr.Dataset(
        {"a": ("t", da.from_array(x, chunks=4)), "b": ("t", y), "t": ("t", t)}
    )
    expected_pd = pd.DataFrame({"a": x, "b": y}, index=pd.Index(t, name="t"))
    expected = dd.from_pandas(expected_pd, chunksize=4)
    actual = ds.to_dask_dataframe(set_index=True)
    assert isinstance(actual, dd.DataFrame)
    pd.testing.assert_frame_equal(actual.compute(), expected.compute())


def test_dataframe_into_delayed():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    df = dd.from_pandas(pdf, npartitions=2)

    def test_func(df):
        return df.sum().sum()

    def delayed_func(i):
        # sanity check
        assert i.sum() == 6

    df = df[["a"]].map_partitions(test_func, meta=(None, int))
    result = delayed(delayed_func)(df)
    assert sum(map(len, result.dask.layers.values())) == 6
    result.compute()


def test_array_to_df_conversion():
    def foo(arr):
        return pd.DataFrame(arr)

    arr = da.random.random((100, 10), chunks=(10, 10))
    result = arr.map_blocks(foo, meta=pd.DataFrame(np.random.random((1, 10))))
    expected = pd.DataFrame(arr.compute())
    assert_eq(result, expected, check_index=False)


@pytest.mark.parametrize("npartitions", [1, 5])
def test_map_partitions_always_series(npartitions):
    pdf = pd.DataFrame({"x": range(50)})
    ddf = dd.from_pandas(pdf, npartitions)

    def assert_series(x):
        assert isinstance(x, pd.Series)
        return x

    res = ddf.x.map_partitions(M.min).map_partitions(assert_series).compute()
    assert len(res) == ddf.npartitions
    assert min(res) == min(pdf.x)

    assert ddf.x.map_partitions(M.min).count().compute() == ddf.npartitions


def test_partitions_are_plain_scalars():
    # https://github.com/dask/dask/issues/11765
    df = pd.DataFrame(
        {"a": range(8), "index": [1, 5, 10, 11, 12, 100, 200, 300]}
    ).set_index("index")
    ddf = dd.from_pandas(df, npartitions=3)

    result = ddf.divisions[0]
    assert type(result) is int


def test_datetime_partitions_are_plain_scalars():
    df = pd.DataFrame(
        {"a": range(8), "index": pd.date_range("2000", periods=8)}
    ).set_index("index")
    ddf = dd.from_pandas(df, npartitions=3)

    result = ddf.divisions[0]
    assert type(result) is pd.Timestamp


def test_loc_partitions_are_plain_scalars():
    indexer = pd.Series([1, 100, 300])
    df = pd.DataFrame(
        {"a": range(8), "index": [1, 5, 10, 11, 12, 100, 200, 300]}
    ).set_index("index")
    ddf = dd.from_pandas(df, npartitions=3)
    result = ddf.loc[indexer]
    assert type(result.divisions[0]) is int

    result = ddf.loc[indexer.tolist()]
    assert type(result.divisions[0]) is int
