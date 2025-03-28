from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    idx = pd.date_range("2000-01-01", periods=12, freq="min")
    pdf = pd.DataFrame({"foo": range(len(idx))}, index=idx)
    pdf["bar"] = 1
    yield pdf


@pytest.fixture
def df(pdf, request):
    npartitions = getattr(request, "param", 2)
    yield from_pandas(pdf, npartitions=npartitions)


@pytest.mark.parametrize(
    "api,how_args",
    [
        ("count", ()),
        ("mean", ()),
        ("sum", ()),
        ("min", ()),
        ("max", ()),
        ("var", ()),
        ("std", ()),
        ("median", ()),
        ("skew", ()),
        ("quantile", (0.5,)),
        ("kurt", ()),
    ],
)
@pytest.mark.parametrize("window,min_periods", ((1, None), (3, 2), (3, 3)))
@pytest.mark.parametrize("center", (True, False))
@pytest.mark.parametrize("df", (1, 2), indirect=True)
def test_rolling_apis(df, pdf, window, api, how_args, min_periods, center):
    args = (window,)
    kwargs = dict(min_periods=min_periods, center=center)

    result = getattr(df.rolling(*args, **kwargs), api)(*how_args)
    expected = getattr(pdf.rolling(*args, **kwargs), api)(*how_args)
    assert_eq(result, expected)

    result = getattr(df.rolling(*args, **kwargs), api)(*how_args)["foo"]
    expected = getattr(pdf.rolling(*args, **kwargs), api)(*how_args)["foo"]
    assert_eq(result, expected)

    q = result.simplify()
    eq = getattr(df["foo"].rolling(*args, **kwargs), api)(*how_args).simplify()
    assert q._name == eq._name


@pytest.mark.parametrize("window", (1, 2))
@pytest.mark.parametrize("df", (1, 2), indirect=True)
def test_rolling_agg(df, pdf, window):
    def my_sum(vals, foo=None, *, bar=None):
        return vals.sum()

    result = df.rolling(window).agg(my_sum, "foo", bar="bar")
    expected = pdf.rolling(window).agg(my_sum, "foo", bar="bar")
    assert_eq(result, expected)

    result = df.rolling(window).agg(my_sum)["foo"]
    expected = pdf.rolling(window).agg(my_sum)["foo"]
    assert_eq(result, expected)

    # simplify up disabled for `agg`, function may access other columns
    q = df.rolling(window).agg(my_sum)["foo"].simplify()
    eq = df["foo"].rolling(window).agg(my_sum).simplify()
    assert q._name != eq._name


@pytest.mark.parametrize("window", (1, 2))
@pytest.mark.parametrize("df", (1, 2), indirect=True)
@pytest.mark.parametrize("raw", (True, False))
@pytest.mark.parametrize("foo", (1, None))
@pytest.mark.parametrize("bar", (2, None))
def test_rolling_apply(df, pdf, window, raw, foo, bar):
    def my_sum(vals, foo_=None, *, bar_=None):
        assert foo_ == foo
        assert bar_ == bar
        if raw:
            assert isinstance(vals, np.ndarray)
        else:
            assert isinstance(vals, pd.Series)
        return vals.sum()

    kwargs = dict(raw=raw, args=(foo,), kwargs=dict(bar_=bar))

    result = df.rolling(window).apply(my_sum, **kwargs)
    expected = pdf.rolling(window).apply(my_sum, **kwargs)
    assert_eq(result, expected)

    result = df.rolling(window).apply(my_sum, **kwargs)["foo"]
    expected = pdf.rolling(window).apply(my_sum, **kwargs)["foo"]
    assert_eq(result, expected)

    # simplify up disabled for `apply`, function may access other columns
    q = df.rolling(window).apply(my_sum, **kwargs)["foo"].simplify()
    eq = df["foo"].rolling(window).apply(my_sum, **kwargs).simplify()
    assert q._name == eq._name


def test_rolling_one_element_window(df, pdf):
    pdf.index = pd.date_range("2000-01-01", periods=12, freq="2s")
    df = from_pandas(pdf, npartitions=3)
    result = pdf.foo.rolling("1s").count()
    expected = df.foo.rolling("1s").count()
    assert_eq(result, expected)


@pytest.mark.parametrize("window", ["2s", "5s", "20s", "10h"])
def test_time_rolling_large_window_variable_chunks(window):
    df = pd.DataFrame(
        {
            "a": pd.date_range("2016-01-01 00:00:00", periods=100, freq="1s"),
            "b": np.random.randint(100, size=(100,)),
        }
    )
    ddf = from_pandas(df, 5)
    ddf = ddf.repartition(divisions=[0, 5, 20, 28, 33, 54, 79, 80, 82, 99])
    df = df.set_index("a")
    ddf = ddf.set_index("a")
    assert_eq(ddf.rolling(window).sum(), df.rolling(window).sum())
    assert_eq(ddf.rolling(window).count(), df.rolling(window).count())
    assert_eq(ddf.rolling(window).mean(), df.rolling(window).mean())


def test_rolling_one_element_window_empty_after(df, pdf):
    pdf.index = pd.date_range("2000-01-01", periods=12, freq="2s")
    df = from_pandas(pdf, npartitions=3)
    result = df.map_overlap(lambda x: x.rolling("1s").count(), before="1s", after="1s")
    expected = pdf.rolling("1s").count()
    assert_eq(result, expected)


@pytest.mark.parametrize("window", [1, 2, 4, 5])
@pytest.mark.parametrize("center", [True, False])
def test_rolling_cov(df, pdf, window, center):
    # DataFrame
    prolling = pdf.drop("foo", axis=1).rolling(window, center=center)
    drolling = df.drop("foo", axis=1).rolling(window, center=center)
    assert_eq(prolling.cov(), drolling.cov())

    # Series
    prolling = pdf.bar.rolling(window, center=center)
    drolling = df.bar.rolling(window, center=center)
    assert_eq(prolling.cov(), drolling.cov())

    # Projection
    actual = df.rolling(window, center=center).cov()[["foo", "bar"]].simplify()
    expected = df[["foo", "bar"]].rolling(window, center=center).cov().simplify()
    assert actual._name == expected._name


def test_rolling_raises():
    df = pd.DataFrame(
        {"a": np.random.randn(25).cumsum(), "b": np.random.randint(100, size=(25,))}
    )
    ddf = from_pandas(df, npartitions=2)

    pytest.raises(ValueError, lambda: ddf.rolling(1.5))
    pytest.raises(ValueError, lambda: ddf.rolling(-1))
    pytest.raises(ValueError, lambda: ddf.rolling(3, min_periods=1.2))
    pytest.raises(ValueError, lambda: ddf.rolling(3, min_periods=-2))
    pytest.raises(NotImplementedError, lambda: ddf.rolling(100).mean().compute())


def test_time_rolling_constructor(df):
    result = df.rolling("4s")
    assert result.window == "4s"
    assert result.min_periods is None
    assert result.win_type is None
