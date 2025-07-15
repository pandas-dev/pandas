from __future__ import annotations

from itertools import product

import pytest

from dask.dataframe._compat import PANDAS_GE_220
from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


def resample(df, freq, how="mean", **kwargs):
    return getattr(df.resample(freq, **kwargs), how)()


@pytest.fixture
def pdf():
    idx = pd.date_range("2000-01-01", periods=12, freq="min")
    pdf = pd.DataFrame({"foo": range(len(idx))}, index=idx)
    pdf["bar"] = 1
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=4)


@pytest.mark.parametrize("kwargs", [{}, {"closed": "left"}])
@pytest.mark.parametrize(
    "api",
    [
        "count",
        "prod",
        "mean",
        "sum",
        "min",
        "max",
        "first",
        "last",
        "var",
        "std",
        "size",
        "nunique",
        "median",
        "quantile",
        "ohlc",
        "sem",
    ],
)
def test_resample_apis(df, pdf, api, kwargs):
    result = getattr(df.resample("2min", **kwargs), api)()
    expected = getattr(pdf.resample("2min", **kwargs), api)()
    assert_eq(result, expected)

    # No column output
    if api not in ("size",):
        result = getattr(df.resample("2min"), api)()["foo"]
        expected = getattr(pdf.resample("2min"), api)()["foo"]
        assert_eq(result, expected)

        if api != "ohlc":
            # ohlc actually gives back a DataFrame, so this doesn't work
            q = result.simplify()
            eq = getattr(df["foo"].resample("2min"), api)().simplify()
            assert q._name == eq._name


@pytest.mark.parametrize(
    ["obj", "method", "npartitions", "freq", "closed", "label"],
    list(
        product(
            ["series", "frame"],
            ["count", "mean", "ohlc"],
            [2, 5],
            ["30min", "h", "D", "W"],
            ["right", "left"],
            ["right", "left"],
        )
    ),
)
def test_series_resample(obj, method, npartitions, freq, closed, label):
    index = pd.date_range("1-1-2000", "2-15-2000", freq="h")
    index = index.union(pd.date_range("4-15-2000", "5-15-2000", freq="h"))
    if obj == "series":
        ps = pd.Series(range(len(index)), index=index)
    elif obj == "frame":
        ps = pd.DataFrame({"a": range(len(index))}, index=index)
    ds = from_pandas(ps, npartitions=npartitions)
    # Series output

    result = resample(ds, freq, how=method, closed=closed, label=label)
    expected = resample(ps, freq, how=method, closed=closed, label=label)

    assert_eq(result, expected, check_dtype=False)

    divisions = result.divisions

    assert expected.index[0] == divisions[0]
    assert expected.index[-1] == divisions[-1]


def test_resample_agg(df, pdf):
    def my_sum(vals, foo=None, *, bar=None):
        return vals.sum()

    result = df.resample("2min").agg(my_sum, "foo", bar="bar")
    expected = pdf.resample("2min").agg(my_sum, "foo", bar="bar")
    assert_eq(result, expected)

    result = df.resample("2min").agg(my_sum)["foo"]
    expected = pdf.resample("2min").agg(my_sum)["foo"]
    assert_eq(result, expected)

    # simplify up disabled for `agg`, function may access other columns
    q = df.resample("2min").agg(my_sum)["foo"].simplify()
    eq = df["foo"].resample("2min").agg(my_sum).simplify()
    assert q._name != eq._name


@pytest.mark.parametrize("method", ["count", "nunique", "size", "sum"])
def test_resample_has_correct_fill_value(method):
    index = pd.date_range("2000-01-01", "2000-02-15", freq="h")
    index = index.union(pd.date_range("4-15-2000", "5-15-2000", freq="h"))
    ps = pd.Series(range(len(index)), index=index)
    ds = from_pandas(ps, npartitions=2)

    assert_eq(
        getattr(ds.resample("30min"), method)(), getattr(ps.resample("30min"), method)()
    )


@pytest.mark.skipif(not PANDAS_GE_220, reason="freq not supported")
def test_resample_divisions_propagation():
    idx = pd.date_range(start="10:00:00.873821", end="10:05:00", freq="0.002s")
    pdf = pd.DataFrame({"data": 1}, index=idx)
    df = from_pandas(pdf, npartitions=10)
    result = df.resample("0.03s").mean()
    result = result.repartition(freq="1D")
    expected = pdf.resample("0.03s").mean()
    assert_eq(result, expected)

    result = df.resample("0.03s").mean().partitions[1]
    expected = pdf.resample("0.03s").mean()[997 : 2 * 997]
    assert_eq(result, expected)
