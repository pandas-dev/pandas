from __future__ import annotations

import contextlib
import warnings

import numpy as np
import pandas as pd
import pytest

import dask.dataframe as dd
from dask.dataframe._compat import tm
from dask.dataframe.utils import assert_eq


@pytest.mark.parametrize(
    "data",
    [
        pd.Series([1, 1, 1, 2, 2, 1, 3, 4], dtype="category"),
        pd.Series(pd.Categorical([1, 1, 1, 2, 2, 1, 3, 4], categories=[4, 3, 2, 1])),
        pd.DataFrame(
            {"a": [1, 2, 3, 4, 4, 3, 2, 1], "b": pd.Categorical(list("abcdabcd"))}
        ),
    ],
)
def test_get_dummies(data):
    exp = pd.get_dummies(data)

    ddata = dd.from_pandas(data, 2)
    res = dd.get_dummies(ddata)
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, exp.columns)


def test_get_dummies_categories_order():
    df = pd.DataFrame({"a": [0.0, 0.0, 1.0, 1.0, 0.0], "b": [1.0, 0.0, 1.0, 0.0, 1.0]})
    ddf = dd.from_pandas(df, npartitions=1)
    ddf = ddf.categorize(columns=["a", "b"])
    res_p = pd.get_dummies(df.astype("category"))
    res_d = dd.get_dummies(ddf)
    assert_eq(res_d, res_p)


def test_get_dummies_object():
    df = pd.DataFrame(
        {
            "a": pd.Categorical([1, 2, 3, 4, 4, 3, 2, 1]),
            "b": list("abcdabcd"),
            "c": pd.Categorical(list("abcdabcd")),
        }
    )
    ddf = dd.from_pandas(df, 2)

    # Explicitly exclude object columns
    exp = pd.get_dummies(df, columns=["a", "c"])
    res = dd.get_dummies(ddf, columns=["a", "c"])
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, exp.columns)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf.b)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf, columns=["b"])


def test_get_dummies_kwargs():
    s = pd.Series([1, 1, 1, 2, 2, 1, 3, 4], dtype="category")
    exp = pd.get_dummies(s, prefix="X", prefix_sep="-")

    ds = dd.from_pandas(s, 2)
    res = dd.get_dummies(ds, prefix="X", prefix_sep="-")
    assert_eq(res, exp)

    exp = pd.get_dummies(s, drop_first=True)
    res = dd.get_dummies(ds, drop_first=True)
    assert_eq(res, exp)

    # nan
    s = pd.Series([1, 1, 1, 2, np.nan, 3, np.nan, 5], dtype="category")
    exp = pd.get_dummies(s)

    ds = dd.from_pandas(s, 2)
    res = dd.get_dummies(ds)
    assert_eq(res, exp)

    # dummy_na
    exp = pd.get_dummies(s, dummy_na=True)
    res = dd.get_dummies(ds, dummy_na=True)
    assert_eq(res, exp)


@contextlib.contextmanager
def ignore_numpy_bool8_deprecation():
    # This warning comes from inside `pandas`. We can't do anything about it, so we ignore the warning.
    # Note it's been fixed upstream in `pandas` https://github.com/pandas-dev/pandas/pull/49886.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=DeprecationWarning,
            message="`np.bool8` is a deprecated alias for `np.bool_`",
        )
        yield


def test_get_dummies_sparse():
    s = pd.Series(pd.Categorical(["a", "b", "a"], categories=["a", "b", "c"]))
    ds = dd.from_pandas(s, 2)

    exp = pd.get_dummies(s, sparse=True)
    res = dd.get_dummies(ds, sparse=True)
    with ignore_numpy_bool8_deprecation():
        assert_eq(exp, res)

    dtype = res.compute().a.dtype
    assert dtype.fill_value == bool(0)
    assert dtype.subtype == bool
    assert isinstance(res.a.compute().dtype, pd.SparseDtype)

    exp = pd.get_dummies(s.to_frame(name="a"), sparse=True)
    res = dd.get_dummies(ds.to_frame(name="a"), sparse=True)
    with ignore_numpy_bool8_deprecation():
        assert_eq(exp, res)
    assert isinstance(res.a_a.compute().dtype, pd.SparseDtype)


def test_get_dummies_sparse_mix():
    df = pd.DataFrame(
        {
            "A": pd.Categorical(["a", "b", "a"], categories=["a", "b", "c"]),
            "B": [0, 0, 1],
        }
    )
    ddf = dd.from_pandas(df, 2)

    exp = pd.get_dummies(df, sparse=True)
    res = dd.get_dummies(ddf, sparse=True)
    with ignore_numpy_bool8_deprecation():
        assert_eq(exp, res)

    dtype = res.compute().A_a.dtype
    assert dtype.fill_value == bool(0)
    assert dtype.subtype == bool
    assert isinstance(res.A_a.compute().dtype, pd.SparseDtype)


def test_get_dummies_dtype():
    df = pd.DataFrame(
        {
            "A": pd.Categorical(["a", "b", "a"], categories=["a", "b", "c"]),
            "B": [0, 0, 1],
        }
    )
    ddf = dd.from_pandas(df, 2)

    exp = pd.get_dummies(df, dtype="float64")
    res = dd.get_dummies(ddf, dtype="float64")
    assert_eq(exp, res)
    assert res.compute().A_a.dtype == "float64"

    # dask's get_dummies on a pandas dataframe.
    assert_eq(dd.get_dummies(df, dtype="float64"), exp)
    assert res.compute().A_a.dtype == "float64"


def test_get_dummies_errors():
    with pytest.raises(NotImplementedError):
        # not Categorical
        s = pd.Series([1, 1, 1, 2, 2, 1, 3, 4])
        ds = dd.from_pandas(s, 2)
        dd.get_dummies(ds)

    # unknown categories
    df = pd.DataFrame({"x": list("abcbc"), "y": list("bcbcb")})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf = ddf.astype("category")

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf, columns=["x", "y"])

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf.x)


@pytest.mark.parametrize("values", ["B", ["B"], ["B", "D"]])
@pytest.mark.parametrize("aggfunc", ["mean", "sum", "count", "first", "last"])
def test_pivot_table(values, aggfunc):
    df = pd.DataFrame(
        {
            "A": np.random.choice(list("XYZ"), size=100),
            "B": np.random.randn(100),
            "C": pd.Categorical(np.random.choice(list("abc"), size=100)),
            "D": np.random.randn(100),
        }
    )
    ddf = dd.from_pandas(df, 5).repartition((0, 20, 40, 60, 80, 98, 99))

    res = dd.pivot_table(ddf, index="A", columns="C", values=values, aggfunc=aggfunc)
    exp = pd.pivot_table(
        df, index="A", columns="C", values=values, aggfunc=aggfunc, observed=False
    )
    if aggfunc == "count":
        # dask result cannot be int64 dtype depending on divisions because of NaN
        exp = exp.astype(np.float64)

    assert_eq(res, exp)

    # method
    res = ddf.pivot_table(index="A", columns="C", values=values, aggfunc=aggfunc)
    exp = df.pivot_table(
        index="A", columns="C", values=values, aggfunc=aggfunc, observed=False
    )
    if aggfunc == "count":
        # dask result cannot be int64 dtype depending on divisions because of NaN
        exp = exp.astype(np.float64)
    assert_eq(res, exp)


@pytest.mark.parametrize("values", ["B", ["D"], ["B", "D"]])
@pytest.mark.parametrize("aggfunc", ["first", "last"])
def test_pivot_table_firstlast(values, aggfunc):
    df = pd.DataFrame(
        {
            "A": np.random.choice(list("XYZ"), size=100),
            "B": np.random.randn(100),
            "C": pd.Categorical(np.random.choice(list("abc"), size=100)),
            "D": np.random.choice(list("abc"), size=100),
        }
    )
    ddf = dd.from_pandas(df, 5).repartition((0, 20, 40, 60, 80, 98, 99))

    res = dd.pivot_table(ddf, index="A", columns="C", values=values, aggfunc=aggfunc)
    exp = pd.pivot_table(
        df, index="A", columns="C", values=values, aggfunc=aggfunc, observed=False
    )

    assert_eq(exp, res)

    # method
    res = ddf.pivot_table(index="A", columns="C", values=values, aggfunc=aggfunc)
    exp = df.pivot_table(
        index="A", columns="C", values=values, aggfunc=aggfunc, observed=False
    )

    assert_eq(exp, res)


def test_pivot_table_dtype():
    df = pd.DataFrame(
        {"A": list("AABB"), "B": pd.Categorical(list("ABAB")), "C": [1, 2, 3, 4]}
    )
    ddf = dd.from_pandas(df, 2)
    res = dd.pivot_table(ddf, index="A", columns="B", values="C", aggfunc="count")

    exp_index = pd.CategoricalIndex(["A", "B"], name="B")
    exp = pd.Series([np.float64] * 2, index=exp_index)
    tm.assert_series_equal(res.dtypes, exp)

    exp = pd.pivot_table(
        df, index="A", columns="B", values="C", aggfunc="count", observed=False
    ).astype(np.float64)

    assert_eq(res, exp)


def test_pivot_table_index_dtype():
    df = pd.DataFrame(
        {
            "A": pd.date_range(start="2019-08-01", periods=3, freq="1D"),
            "B": pd.Categorical(list("abc")),
            "C": [1, 2, 3],
        }
    )
    ddf = dd.from_pandas(df, 2)
    res = dd.pivot_table(ddf, index="A", columns="B", values="C", aggfunc="count")

    assert res.index.dtype == np.dtype("datetime64[ns]")


def test_pivot_table_errors():
    df = pd.DataFrame(
        {
            "A": np.random.choice(list("abc"), size=10),
            "B": np.random.randn(10),
            "C": pd.Categorical(np.random.choice(list("abc"), size=10)),
        }
    )
    ddf = dd.from_pandas(df, 2)

    msg = "'index' must be the name of an existing column"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index=["A"], columns="C", values="B")
    assert msg in str(err.value)
    msg = "'columns' must be the name of an existing column"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index="A", columns=["C"], values="B")
    assert msg in str(err.value)
    msg = "'values' must refer to an existing column or columns"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index="A", columns="C", values=[["B"]])
    assert msg in str(err.value)

    msg = "aggfunc must be either 'mean', 'sum', 'count', 'first', 'last'"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index="A", columns="C", values="B", aggfunc=["sum"])
    assert msg in str(err.value)

    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index="A", columns="C", values="B", aggfunc="xx")
    assert msg in str(err.value)

    # unknown categories
    ddf["C"] = ddf.C.cat.as_unknown()
    msg = "'columns' must have known categories"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index="A", columns="C", values=["B"])
    assert msg in str(err.value)

    df = pd.DataFrame(
        {
            "A": np.random.choice(list("abc"), size=10),
            "B": np.random.randn(10),
            "C": np.random.choice(list("abc"), size=10),
        }
    )
    ddf = dd.from_pandas(df, 2)
    msg = "'columns' must be category dtype"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index="A", columns="C", values="B")
    assert msg in str(err.value)
