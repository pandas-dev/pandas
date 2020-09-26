import numpy as np
import pytest

from pandas import DataFrame, Series, concat, isna, notna
import pandas._testing as tm

import pandas.tseries.offsets as offsets


def test_series(series):
    result = series.rolling(50).sum()
    assert isinstance(result, Series)
    tm.assert_almost_equal(result.iloc[-1], np.nansum(series[-50:]))


def test_frame(raw, frame):
    result = frame.rolling(50).sum()
    assert isinstance(result, DataFrame)
    tm.assert_series_equal(
        result.iloc[-1, :],
        frame.iloc[-50:, :].apply(np.nansum, axis=0, raw=raw),
        check_names=False,
    )


def test_time_rule_series(series):
    win = 25
    minp = 10
    ser = series[::2].resample("B").mean()
    series_result = ser.rolling(window=win, min_periods=minp).sum()
    last_date = series_result.index[-1]
    prev_date = last_date - 24 * offsets.BDay()

    trunc_series = series[::2].truncate(prev_date, last_date)
    tm.assert_almost_equal(series_result[-1], np.nansum(trunc_series))


def test_time_rule_frame(raw, frame):
    win = 25
    minp = 10
    frm = frame[::2].resample("B").mean()
    frame_result = frm.rolling(window=win, min_periods=minp).sum()
    last_date = frame_result.index[-1]
    prev_date = last_date - 24 * offsets.BDay()

    trunc_frame = frame[::2].truncate(prev_date, last_date)
    tm.assert_series_equal(
        frame_result.xs(last_date),
        trunc_frame.apply(np.nansum, raw=raw),
        check_names=False,
    )


def test_nans():
    obj = Series(np.random.randn(50))
    obj[:10] = np.NaN
    obj[-10:] = np.NaN

    result = obj.rolling(50, min_periods=30).sum()
    tm.assert_almost_equal(result.iloc[-1], np.nansum(obj[10:-10]))

    # min_periods is working correctly
    result = obj.rolling(20, min_periods=15).sum()
    assert isna(result.iloc[23])
    assert not isna(result.iloc[24])

    assert not isna(result.iloc[-6])
    assert isna(result.iloc[-5])

    obj2 = Series(np.random.randn(20))
    result = obj2.rolling(10, min_periods=5).sum()
    assert isna(result.iloc[3])
    assert notna(result.iloc[4])


@pytest.mark.parametrize("minp", [0, 99, 100])
def test_min_periods(series, minp):
    result = series.rolling(len(series) + 1, min_periods=minp).sum()
    expected = series.rolling(len(series), min_periods=minp).sum()
    nan_mask = isna(result)
    tm.assert_series_equal(nan_mask, isna(expected))

    nan_mask = ~nan_mask
    tm.assert_almost_equal(result[nan_mask], expected[nan_mask])


def test_center():
    obj = Series(np.random.randn(50))
    obj[:10] = np.NaN
    obj[-10:] = np.NaN

    result = obj.rolling(20, min_periods=15, center=True).sum()
    expected = (
        concat([obj, Series([np.NaN] * 9)])
        .rolling(20, min_periods=15)
        .sum()[9:]
        .reset_index(drop=True)
    )
    tm.assert_series_equal(result, expected)


def test_center_reindex_series(series):
    # shifter index
    s = [f"x{x:d}" for x in range(12)]
    minp = 10

    series_xp = (
        series.reindex(list(series.index) + s)
        .rolling(window=25, min_periods=minp)
        .sum()
        .shift(-12)
        .reindex(series.index)
    )
    series_rs = series.rolling(window=25, min_periods=minp, center=True).sum()
    tm.assert_series_equal(series_xp, series_rs)


def test_center_reindex_frame(frame):
    # shifter index
    s = [f"x{x:d}" for x in range(12)]
    minp = 10

    frame_xp = (
        frame.reindex(list(frame.index) + s)
        .rolling(window=25, min_periods=minp)
        .sum()
        .shift(-12)
        .reindex(frame.index)
    )
    frame_rs = frame.rolling(window=25, min_periods=minp, center=True).sum()
    tm.assert_frame_equal(frame_xp, frame_rs)
