from datetime import datetime
from operator import methodcaller

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
from pandas.core.groupby.grouper import Grouper
from pandas.core.indexes.datetimes import date_range
import pandas.util.testing as tm
from pandas.util.testing import assert_frame_equal, assert_series_equal

test_series = Series(np.random.randn(1000), index=date_range("1/1/2000", periods=1000))


def test_apply():
    grouper = Grouper(freq="A", label="right", closed="right")

    grouped = test_series.groupby(grouper)

    def f(x):
        return x.sort_values()[-3:]

    applied = grouped.apply(f)
    expected = test_series.groupby(lambda x: x.year).apply(f)

    applied.index = applied.index.droplevel(0)
    expected.index = expected.index.droplevel(0)
    assert_series_equal(applied, expected)


def test_count():
    test_series[::3] = np.nan

    expected = test_series.groupby(lambda x: x.year).count()

    grouper = Grouper(freq="A", label="right", closed="right")
    result = test_series.groupby(grouper).count()
    expected.index = result.index
    assert_series_equal(result, expected)

    result = test_series.resample("A").count()
    expected.index = result.index
    assert_series_equal(result, expected)


def test_numpy_reduction():
    result = test_series.resample("A", closed="right").prod()

    expected = test_series.groupby(lambda x: x.year).agg(np.prod)
    expected.index = result.index

    assert_series_equal(result, expected)


def test_apply_iteration():
    # #2300
    N = 1000
    ind = pd.date_range(start="2000-01-01", freq="D", periods=N)
    df = DataFrame({"open": 1, "close": 2}, index=ind)
    tg = Grouper(freq="M")

    _, grouper, _ = tg._get_grouper(df)

    # Errors
    grouped = df.groupby(grouper, group_keys=False)

    def f(df):
        return df["close"] / df["open"]

    # it works!
    result = grouped.apply(f)
    tm.assert_index_equal(result.index, df.index)


@pytest.mark.parametrize(
    "name, func",
    [
        ("Int64Index", tm.makeIntIndex),
        ("Index", tm.makeUnicodeIndex),
        ("Float64Index", tm.makeFloatIndex),
        ("MultiIndex", lambda m: tm.makeCustomIndex(m, 2)),
    ],
)
def test_fails_on_no_datetime_index(name, func):
    n = 2
    index = func(n)
    df = DataFrame({"a": np.random.randn(n)}, index=index)

    msg = (
        "Only valid with DatetimeIndex, TimedeltaIndex "
        "or PeriodIndex, but got an instance of '{}'".format(name)
    )
    with pytest.raises(TypeError, match=msg):
        df.groupby(Grouper(freq="D"))


def test_aaa_group_order():
    # GH 12840
    # check TimeGrouper perform stable sorts
    n = 20
    data = np.random.randn(n, 4)
    df = DataFrame(data, columns=["A", "B", "C", "D"])
    df["key"] = [
        datetime(2013, 1, 1),
        datetime(2013, 1, 2),
        datetime(2013, 1, 3),
        datetime(2013, 1, 4),
        datetime(2013, 1, 5),
    ] * 4
    grouped = df.groupby(Grouper(key="key", freq="D"))

    tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 1)), df[::5])
    tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 2)), df[1::5])
    tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 3)), df[2::5])
    tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 4)), df[3::5])
    tm.assert_frame_equal(grouped.get_group(datetime(2013, 1, 5)), df[4::5])


def test_aggregate_normal(resample_method):
    """Check TimeGrouper's aggregation is identical as normal groupby."""

    if resample_method == "ohlc":
        pytest.xfail(reason="DataError: No numeric types to aggregate")

    data = np.random.randn(20, 4)
    normal_df = DataFrame(data, columns=["A", "B", "C", "D"])
    normal_df["key"] = [1, 2, 3, 4, 5] * 4

    dt_df = DataFrame(data, columns=["A", "B", "C", "D"])
    dt_df["key"] = [
        datetime(2013, 1, 1),
        datetime(2013, 1, 2),
        datetime(2013, 1, 3),
        datetime(2013, 1, 4),
        datetime(2013, 1, 5),
    ] * 4

    normal_grouped = normal_df.groupby("key")
    dt_grouped = dt_df.groupby(Grouper(key="key", freq="D"))

    expected = getattr(normal_grouped, resample_method)()
    dt_result = getattr(dt_grouped, resample_method)()
    expected.index = date_range(start="2013-01-01", freq="D", periods=5, name="key")
    tm.assert_equal(expected, dt_result)

    # if TimeGrouper is used included, 'nth' doesn't work yet

    """
    for func in ['nth']:
        expected = getattr(normal_grouped, func)(3)
        expected.index = date_range(start='2013-01-01',
                                    freq='D', periods=5, name='key')
        dt_result = getattr(dt_grouped, func)(3)
        assert_frame_equal(expected, dt_result)
    """


@pytest.mark.parametrize(
    "method, method_args, unit",
    [
        ("sum", dict(), 0),
        ("sum", dict(min_count=0), 0),
        ("sum", dict(min_count=1), np.nan),
        ("prod", dict(), 1),
        ("prod", dict(min_count=0), 1),
        ("prod", dict(min_count=1), np.nan),
    ],
)
def test_resample_entirly_nat_window(method, method_args, unit):
    s = pd.Series([0] * 2 + [np.nan] * 2, index=pd.date_range("2017", periods=4))
    result = methodcaller(method, **method_args)(s.resample("2d"))
    expected = pd.Series(
        [0.0, unit], index=pd.to_datetime(["2017-01-01", "2017-01-03"])
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "func, fill_value",
    [("min", np.nan), ("max", np.nan), ("sum", 0), ("prod", 1), ("count", 0)],
)
def test_aggregate_with_nat(func, fill_value):
    # check TimeGrouper's aggregation is identical as normal groupby
    # if NaT is included, 'var', 'std', 'mean', 'first','last'
    # and 'nth' doesn't work yet

    n = 20
    data = np.random.randn(n, 4).astype("int64")
    normal_df = DataFrame(data, columns=["A", "B", "C", "D"])
    normal_df["key"] = [1, 2, np.nan, 4, 5] * 4

    dt_df = DataFrame(data, columns=["A", "B", "C", "D"])
    dt_df["key"] = [
        datetime(2013, 1, 1),
        datetime(2013, 1, 2),
        pd.NaT,
        datetime(2013, 1, 4),
        datetime(2013, 1, 5),
    ] * 4

    normal_grouped = normal_df.groupby("key")
    dt_grouped = dt_df.groupby(Grouper(key="key", freq="D"))

    normal_result = getattr(normal_grouped, func)()
    dt_result = getattr(dt_grouped, func)()

    pad = DataFrame([[fill_value] * 4], index=[3], columns=["A", "B", "C", "D"])
    expected = normal_result.append(pad)
    expected = expected.sort_index()
    expected.index = date_range(start="2013-01-01", freq="D", periods=5, name="key")
    assert_frame_equal(expected, dt_result)
    assert dt_result.index.name == "key"


def test_aggregate_with_nat_size():
    # GH 9925
    n = 20
    data = np.random.randn(n, 4).astype("int64")
    normal_df = DataFrame(data, columns=["A", "B", "C", "D"])
    normal_df["key"] = [1, 2, np.nan, 4, 5] * 4

    dt_df = DataFrame(data, columns=["A", "B", "C", "D"])
    dt_df["key"] = [
        datetime(2013, 1, 1),
        datetime(2013, 1, 2),
        pd.NaT,
        datetime(2013, 1, 4),
        datetime(2013, 1, 5),
    ] * 4

    normal_grouped = normal_df.groupby("key")
    dt_grouped = dt_df.groupby(Grouper(key="key", freq="D"))

    normal_result = normal_grouped.size()
    dt_result = dt_grouped.size()

    pad = Series([0], index=[3])
    expected = normal_result.append(pad)
    expected = expected.sort_index()
    expected.index = date_range(start="2013-01-01", freq="D", periods=5, name="key")
    assert_series_equal(expected, dt_result)
    assert dt_result.index.name == "key"


def test_repr():
    # GH18203
    result = repr(Grouper(key="A", freq="H"))
    expected = (
        "TimeGrouper(key='A', freq=<Hour>, axis=0, sort=True, "
        "closed='left', label='left', how='mean', "
        "convention='e', base=0)"
    )
    assert result == expected


@pytest.mark.parametrize(
    "method, method_args, expected_values",
    [
        ("sum", dict(), [1, 0, 1]),
        ("sum", dict(min_count=0), [1, 0, 1]),
        ("sum", dict(min_count=1), [1, np.nan, 1]),
        ("sum", dict(min_count=2), [np.nan, np.nan, np.nan]),
        ("prod", dict(), [1, 1, 1]),
        ("prod", dict(min_count=0), [1, 1, 1]),
        ("prod", dict(min_count=1), [1, np.nan, 1]),
        ("prod", dict(min_count=2), [np.nan, np.nan, np.nan]),
    ],
)
def test_upsample_sum(method, method_args, expected_values):
    s = pd.Series(1, index=pd.date_range("2017", periods=2, freq="H"))
    resampled = s.resample("30T")
    index = pd.to_datetime(
        ["2017-01-01T00:00:00", "2017-01-01T00:30:00", "2017-01-01T01:00:00"]
    )
    result = methodcaller(method, **method_args)(resampled)
    expected = pd.Series(expected_values, index=index)
    tm.assert_series_equal(result, expected)
