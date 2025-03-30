from __future__ import annotations

import contextlib
from itertools import product

import pandas as pd
import pytest

import dask.dataframe as dd
from dask.dataframe._compat import PANDAS_GE_220, PANDAS_GE_300
from dask.dataframe.utils import assert_eq


def resample(df, freq, how="mean", **kwargs):
    return getattr(df.resample(freq, **kwargs), how)()


ME = "ME" if PANDAS_GE_220 else "M"


@pytest.mark.parametrize(
    ["obj", "method", "npartitions", "freq", "closed", "label"],
    list(
        product(
            ["series", "frame"],
            ["count", "mean", "ohlc"],
            [2, 5],
            ["30min", "h", "D", "W", ME],
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
    ds = dd.from_pandas(ps, npartitions=npartitions)
    # Series output

    result = resample(ds, freq, how=method, closed=closed, label=label)
    expected = resample(ps, freq, how=method, closed=closed, label=label)

    assert_eq(result, expected, check_dtype=False)

    divisions = result.divisions

    assert expected.index[0] == divisions[0]
    assert expected.index[-1] == divisions[-1]


@pytest.mark.parametrize("method", ["count", "nunique", "size", "sum"])
def test_resample_has_correct_fill_value(method):
    index = pd.date_range("2000-01-01", "2000-02-15", freq="h")
    index = index.union(pd.date_range("4-15-2000", "5-15-2000", freq="h"))
    ps = pd.Series(range(len(index)), index=index)
    ds = dd.from_pandas(ps, npartitions=2)

    assert_eq(
        getattr(ds.resample("30min"), method)(), getattr(ps.resample("30min"), method)()
    )


def test_resample_agg():
    index = pd.date_range("2000-01-01", "2000-02-15", freq="h")
    ps = pd.Series(range(len(index)), index=index)
    ds = dd.from_pandas(ps, npartitions=2)

    assert_eq(ds.resample("10min").agg("mean"), ps.resample("10min").agg("mean"))
    assert_eq(
        ds.resample("10min").agg(["mean", "min"]),
        ps.resample("10min").agg(["mean", "min"]),
    )


def test_resample_agg_passes_kwargs():
    index = pd.date_range("2000-01-01", "2000-02-15", freq="h")
    ps = pd.Series(range(len(index)), index=index)
    ds = dd.from_pandas(ps, npartitions=2)

    def foo(series, bar=1, *args, **kwargs):
        return bar

    assert_eq(ds.resample("2h").agg(foo, bar=2), ps.resample("2h").agg(foo, bar=2))
    assert (ds.resample("2h").agg(foo, bar=2) == 2).compute().all()


def test_resample_throws_error_when_parition_index_does_not_match_index():
    index = pd.date_range("1-1-2000", "2-15-2000", freq="D")
    index = index.union(pd.date_range("4-15-2000", "5-15-2000", freq="D"))
    ps = pd.Series(range(len(index)), index=index)
    ds = dd.from_pandas(ps, npartitions=5)
    with pytest.raises(ValueError, match="Index is not contained within new index."):
        ds.resample(f"2{ME}").count().compute()


def test_resample_pads_last_division_to_avoid_off_by_one():
    # https://github.com/dask/dask/issues/6230
    times = [
        1545362463409128000,
        1545362504369352000,
        1545362545326966000,
        1545363118769636000,
        1545363159726490000,
        1545363200687178000,
        1545363241648824000,
        1573318190393973000,
        1573318231353350000,
        1573318272313774000,
        1573318313275299000,
        1573318354233962000,
        1573318395195456000,
        1573318436154609000,
        1580687544437145000,
        1580687585394881000,
        1580687667316809000,
        1580687708275414000,
        1580687790195742000,
        1580687831154951000,
        1580687872115363000,
        1580687954035133000,
        1559127673402811000,
    ]

    freq = "1QE" if PANDAS_GE_220 else "1Q"
    df = pd.DataFrame({"Time": times, "Counts": range(len(times))})
    df["Time"] = pd.to_datetime(df["Time"], utc=True)
    expected = df.set_index("Time").resample(freq).size()

    ddf = dd.from_pandas(df, npartitions=2).set_index("Time")
    actual = ddf.resample(freq).size().compute()
    assert_eq(actual, expected)


def test_resample_does_not_evenly_divide_day():
    import numpy as np

    index = pd.date_range("2012-01-02", "2012-02-02", freq="h")
    index = index.union(pd.date_range("2012-03-02", "2012-04-02", freq="h"))
    df = pd.DataFrame({"p": np.random.random(len(index))}, index=index)
    ddf = dd.from_pandas(df, npartitions=5)
    # Frequency doesn't evenly divide day
    expected = df.resample("2D").count()
    result = ddf.resample("2D").count().compute()

    assert_eq(result, expected)


def test_series_resample_does_not_evenly_divide_day():
    index = pd.date_range("2012-01-02 00:00:00", "2012-01-02 01:00:00", freq="min")
    index = index.union(
        pd.date_range("2012-01-02 06:00:00", "2012-01-02 08:00:00", freq="min")
    )
    s = pd.Series(range(len(index)), index=index)
    ds = dd.from_pandas(s, npartitions=5)
    # Frequency doesn't evenly divide day
    expected = s.resample("57min").mean()
    result = ds.resample("57min").mean().compute()

    assert_eq(result, expected)


def test_unknown_divisions_error():
    df = pd.DataFrame({"x": [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=2, sort=False).clear_divisions()
    try:
        ddf.x.resample("1m").mean()
        assert False
    except ValueError as e:
        assert "divisions" in str(e)


def test_resample_index_name():
    from datetime import datetime, timedelta

    import numpy as np

    date_today = datetime.now()
    days = pd.date_range(date_today, date_today + timedelta(20), freq="D")
    data = np.random.randint(1, high=100, size=len(days))

    df = pd.DataFrame({"date": days, "values": data})
    df = df.set_index("date")

    ddf = dd.from_pandas(df, npartitions=4)

    assert ddf.resample("D").mean().head().index.name == "date"


def test_series_resample_non_existent_datetime():
    index = [
        pd.Timestamp("2016-10-15 00:00:00"),
        pd.Timestamp("2016-10-16 10:00:00"),
        pd.Timestamp("2016-10-17 00:00:00"),
    ]
    df = pd.DataFrame([[1], [2], [3]], index=index)
    df.index = df.index.tz_localize("America/Sao_Paulo")
    ddf = dd.from_pandas(df, npartitions=1)
    result = ddf.resample("1D").mean()
    expected = df.resample("1D").mean()

    assert_eq(result, expected, check_freq=False)


@pytest.mark.parametrize("agg", ["nunique", "mean", "count", "size", "quantile"])
def test_common_aggs(agg):
    index = pd.date_range("2000-01-01", "2000-02-15", freq="h")
    ps = pd.Series(range(len(index)), index=index)
    ds = dd.from_pandas(ps, npartitions=2)

    f = lambda df: getattr(df, agg)()

    res = f(ps.resample("1D"))
    expected = f(ds.resample("1D"))

    assert_eq(res, expected, check_dtype=False)


def test_rule_deprecated():
    index = pd.date_range("2000-01-01", "2000-02-15", freq="h")
    s = pd.Series(range(len(index)), index=index)
    ds = dd.from_pandas(s, npartitions=2)

    if PANDAS_GE_300:
        ctx = pytest.warns(FutureWarning, match="'d' is deprecated")
    else:
        ctx = contextlib.nullcontext()

    with ctx:
        res = s.resample("1d").count()
    with ctx:
        expected = ds.resample("1d").count()

    with ctx:
        assert_eq(res, expected, check_dtype=False)
