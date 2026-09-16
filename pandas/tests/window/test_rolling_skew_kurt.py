from functools import partial

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

from pandas.tseries import offsets


@pytest.mark.parametrize("sp_func, roll_func", [["kurtosis", "kurt"], ["skew", "skew"]])
def test_series(series, sp_func, roll_func):
    sp_stats = pytest.importorskip("scipy.stats")

    compare_func = partial(getattr(sp_stats, sp_func), bias=False)
    result = getattr(series.rolling(50), roll_func)()
    assert isinstance(result, pd.Series)
    tm.assert_almost_equal(result.iloc[-1], compare_func(series[-50:]))


@pytest.mark.parametrize("sp_func, roll_func", [["kurtosis", "kurt"], ["skew", "skew"]])
def test_frame(raw, frame, sp_func, roll_func):
    sp_stats = pytest.importorskip("scipy.stats")

    compare_func = partial(getattr(sp_stats, sp_func), bias=False)
    result = getattr(frame.rolling(50), roll_func)()
    assert isinstance(result, pd.DataFrame)
    tm.assert_series_equal(
        result.iloc[-1, :],
        frame.iloc[-50:, :].apply(compare_func, axis=0, raw=raw),
        check_names=False,
    )


@pytest.mark.parametrize("sp_func, roll_func", [["kurtosis", "kurt"], ["skew", "skew"]])
def test_time_rule_series(series, sp_func, roll_func):
    sp_stats = pytest.importorskip("scipy.stats")

    compare_func = partial(getattr(sp_stats, sp_func), bias=False)
    win = 25
    ser = series[::2].resample("B").mean()
    series_result = getattr(ser.rolling(window=win, min_periods=10), roll_func)()
    last_date = series_result.index[-1]
    prev_date = last_date - 24 * offsets.BDay()

    trunc_series = series[::2].truncate(prev_date, last_date)
    tm.assert_almost_equal(series_result.iloc[-1], compare_func(trunc_series))


@pytest.mark.parametrize("sp_func, roll_func", [["kurtosis", "kurt"], ["skew", "skew"]])
def test_time_rule_frame(raw, frame, sp_func, roll_func):
    sp_stats = pytest.importorskip("scipy.stats")

    compare_func = partial(getattr(sp_stats, sp_func), bias=False)
    win = 25
    frm = frame[::2].resample("B").mean()
    frame_result = getattr(frm.rolling(window=win, min_periods=10), roll_func)()
    last_date = frame_result.index[-1]
    prev_date = last_date - 24 * offsets.BDay()

    trunc_frame = frame[::2].truncate(prev_date, last_date)
    tm.assert_series_equal(
        frame_result.xs(last_date),
        trunc_frame.apply(compare_func, raw=raw),
        check_names=False,
    )


@pytest.mark.parametrize("sp_func, roll_func", [["kurtosis", "kurt"], ["skew", "skew"]])
def test_nans(sp_func, roll_func):
    sp_stats = pytest.importorskip("scipy.stats")

    compare_func = partial(getattr(sp_stats, sp_func), bias=False)
    obj = pd.Series(np.random.default_rng(2).standard_normal(50))
    obj[:10] = np.nan
    obj[-10:] = np.nan

    result = getattr(obj.rolling(50, min_periods=30), roll_func)()
    tm.assert_almost_equal(result.iloc[-1], compare_func(obj[10:-10]))

    # min_periods is working correctly
    result = getattr(obj.rolling(20, min_periods=15), roll_func)()
    assert pd.isna(result.iloc[23])
    assert not pd.isna(result.iloc[24])

    assert not pd.isna(result.iloc[-6])
    assert pd.isna(result.iloc[-5])

    obj2 = pd.Series(np.random.default_rng(2).standard_normal(20))
    result = getattr(obj2.rolling(10, min_periods=5), roll_func)()
    assert pd.isna(result.iloc[3])
    assert pd.notna(result.iloc[4])

    result0 = getattr(obj.rolling(20, min_periods=0), roll_func)()
    result1 = getattr(obj.rolling(20, min_periods=1), roll_func)()
    tm.assert_almost_equal(result0, result1)


@pytest.mark.parametrize("minp", [0, 99, 100])
@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_min_periods(series, minp, roll_func, step):
    result = getattr(
        series.rolling(len(series) + 1, min_periods=minp, step=step), roll_func
    )()
    expected = getattr(
        series.rolling(len(series), min_periods=minp, step=step), roll_func
    )()
    nan_mask = pd.isna(result)
    tm.assert_series_equal(nan_mask, pd.isna(expected))

    nan_mask = ~nan_mask
    tm.assert_almost_equal(result[nan_mask], expected[nan_mask])


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_center(roll_func):
    obj = pd.Series(np.random.default_rng(2).standard_normal(50))
    obj[:10] = np.nan
    obj[-10:] = np.nan

    result = getattr(obj.rolling(20, center=True), roll_func)()
    expected = (
        getattr(pd.concat([obj, pd.Series([np.nan] * 9)]).rolling(20), roll_func)()
        .iloc[9:]
        .reset_index(drop=True)
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_center_reindex_series(series, roll_func):
    # shifter index
    s = [f"x{x:d}" for x in range(12)]

    series_xp = (
        getattr(
            series.reindex(list(series.index) + s).rolling(window=25),
            roll_func,
        )()
        .shift(-12)
        .reindex(series.index)
    )
    series_rs = getattr(series.rolling(window=25, center=True), roll_func)()
    tm.assert_series_equal(series_xp, series_rs)


@pytest.mark.slow
@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_center_reindex_frame(frame, roll_func):
    # shifter index
    s = [f"x{x:d}" for x in range(12)]

    frame_xp = (
        getattr(
            frame.reindex(list(frame.index) + s).rolling(window=25),
            roll_func,
        )()
        .shift(-12)
        .reindex(frame.index)
    )
    frame_rs = getattr(frame.rolling(window=25, center=True), roll_func)()
    tm.assert_frame_equal(frame_xp, frame_rs)


def test_rolling_skew_edge_cases(step):
    expected = pd.Series([np.nan] * 5)[::step]
    # yields all NaN (0 variance)
    d = pd.Series([1] * 5)
    x = d.rolling(window=5, step=step).skew()
    # index 4 should be NaN as it contains 5 same obs
    tm.assert_series_equal(expected, x)

    expected = pd.Series([np.nan] * 5)[::step]
    # yields all NaN (window too small)
    d = pd.Series(np.random.default_rng(2).standard_normal(5))
    x = d.rolling(window=2, step=step).skew()
    tm.assert_series_equal(expected, x)

    # yields [NaN, NaN, NaN, 0.177994, 1.548824]
    d = pd.Series([-1.50837035, -0.1297039, 0.19501095, 1.73508164, 0.41941401])
    expected = pd.Series([np.nan, np.nan, np.nan, 0.177994, 1.548824])[::step]
    x = d.rolling(window=4, step=step).skew()
    tm.assert_series_equal(expected, x)


def test_rolling_kurt_edge_cases(step):
    expected = pd.Series([np.nan] * 5)[::step]

    # yields all NaN (0 variance)
    d = pd.Series([1] * 5)
    x = d.rolling(window=5, step=step).kurt()
    tm.assert_series_equal(expected, x)

    # yields all NaN (window too small)
    expected = pd.Series([np.nan] * 5)[::step]
    d = pd.Series(np.random.default_rng(2).standard_normal(5))
    x = d.rolling(window=3, step=step).kurt()
    tm.assert_series_equal(expected, x)

    # yields [NaN, NaN, NaN, 1.224307, 2.671499]
    d = pd.Series([-1.50837035, -0.1297039, 0.19501095, 1.73508164, 0.41941401])
    expected = pd.Series([np.nan, np.nan, np.nan, 1.224307, 2.671499])[::step]
    x = d.rolling(window=4, step=step).kurt()
    tm.assert_series_equal(expected, x)


def test_rolling_skew_eq_value_fperr(step):
    # #18804 all rolling skew for all equal values should return Nan
    a = pd.Series([1.1] * 15).rolling(window=10, step=step).skew()
    expected = pd.Series([np.nan] * 15)[::step]
    tm.assert_series_equal(a, expected)


def test_rolling_kurt_eq_value_fperr(step):
    # #18804 all rolling kurt for all equal values should return Nan
    a = pd.Series([1.1] * 15).rolling(window=10, step=step).kurt()
    expected = pd.Series([np.nan] * 15)[::step]
    tm.assert_series_equal(a, expected)


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
@pytest.mark.parametrize("scale_factor", [1e-20, 1e20])
def test_skew_kurt_is_scale_invariant(roll_func, scale_factor):
    # GH-62946
    obj = pd.Series(np.random.default_rng(2).standard_normal(50))
    obj_scaled = obj * scale_factor
    result = getattr(obj.rolling(20), roll_func)()
    result_scaled = getattr(obj_scaled.rolling(20), roll_func)()
    tm.assert_series_equal(result, result_scaled)


def _window_reduction(series, window, roll_func):
    # oracle: the matching whole-array reduction over each window on its own. It
    # accumulates independently of the sliding kernels under test and centres the
    # values first, so it stays exact on the offset data below. It also skips
    # NaN, so windows holding one are blanked to match rolling's default
    # min_periods=window.
    values = series.to_numpy()
    expected = [np.nan] * (window - 1)
    for stop in range(window, len(values) + 1):
        chunk = values[stop - window : stop]
        if np.isnan(chunk).any():
            expected.append(np.nan)
        else:
            expected.append(getattr(pd.Series(chunk), roll_func)())
    return pd.Series(expected)


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
@pytest.mark.parametrize("offset", [1e6, 1e10])
def test_rolling_skew_kurt_shared_offset(roll_func, offset):
    # GH#68920 an offset shared by the whole series left almost no significant
    # digits in the deviations the accumulators are built from, so results were
    # wrong with no outlier anywhere in the data
    window = 5
    series = pd.Series([1, 2, 4, 7, 3, 5, 9, 2, 6, 8], dtype="float64") + offset

    result = getattr(series.rolling(window), roll_func)()

    tm.assert_series_equal(
        result, _window_reduction(series, window, roll_func), rtol=1e-10, atol=0
    )


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_rolling_skew_kurt_low_variance_offset(roll_func):
    # GH#68920 values one float64 ulp apart on a large offset: the incremental
    # mean update is a no-op at that magnitude, so the accumulators drifted into
    # garbage. kurt returned 213 here, well outside the [-6, 4] a 4-point window
    # can attain at all.
    window = 4
    ulp = np.spacing(1e8)
    codes = [0, 1, 0, 2, 1, 0, 3, 1, 2, 0, 1, 4, 0, 2, 1, 0]
    series = pd.Series([1e8 + ulp * code for code in codes])

    result = getattr(series.rolling(window), roll_func)()

    # skew and kurt are unchanged by translation, so the same window without the
    # offset is the answer. atol covers the windows whose true skew is 0, where
    # both sides are round-off and a relative tolerance means nothing.
    expected = getattr(pd.Series(codes, dtype="float64").rolling(window), roll_func)()
    tm.assert_series_equal(result, expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_rolling_skew_kurt_nan_gap_recovery(roll_func):
    # GH#68920 a window whose observation count drops to 0 must not poison
    # subsequent overlapping windows with NaN
    window = 4
    series = pd.Series([1.0, 2.0, 3.0, 5.0] + [np.nan] * 4 + [4.0, 5.0, 6.0, 8.0])

    result = getattr(series.rolling(window), roll_func)()

    assert not result.iloc[-1:].isna().any()
    tm.assert_series_equal(
        result, _window_reduction(series, window, roll_func), rtol=1e-12, atol=0
    )


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_rolling_skew_kurt_extreme_range_recovers(roll_func):
    # GH#68920 a window spanning nearly the whole float64 range must not leave
    # the accumulators holding NaN, which would blank every later window
    window = 5
    series = pd.Series(
        [1e308, 1.0, 2.0, 3.0, -1e308, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 3.0]
    )

    result = getattr(series.rolling(window), roll_func)()

    # once the extreme values leave, the windows are ordinary data and exact
    assert not result.iloc[9:].isna().any()
    tm.assert_series_equal(
        result.iloc[9:],
        _window_reduction(series, window, roll_func).iloc[9:],
        rtol=1e-12,
        atol=0,
    )


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_rolling_skew_kurt_drifting_level(roll_func):
    # GH#68920 anchoring the accumulators to a window's first value only helps
    # while the data stays near it. On a series whose level drifts -- a timestamp
    # column, a counter -- the anchor goes stale and the result was wrong again,
    # here by three orders of magnitude, with nothing in the window itself to
    # show for it: the same window recomputed on its own is exact.
    window = 20
    n = 60_000
    rng = np.random.default_rng(0)
    values = np.sort(1.7e9 + np.arange(n) + rng.normal(size=n) * 0.3)

    result = getattr(pd.Series(values).rolling(window), roll_func)()

    # over a slice, not one index: a single position can land on a good value
    # while its neighbours are wrong
    for i in range(n - 50, n):
        expected = getattr(pd.Series(values[i - window + 1 : i + 1]), roll_func)()
        assert result.iloc[i] == pytest.approx(expected, rel=1e-8)


@pytest.mark.parametrize("roll_func", ["kurt", "skew"])
def test_rolling_skew_kurt_midband_outlier_recovers(roll_func):
    # GH#68920 a deviation much beyond 1e77 overflowed the instability test's own
    # arithmetic, leaving it comparing inf against inf -- which is False, so the
    # accumulators were never recomputed and every later window returned the same
    # frozen garbage. 1e308 does not reach this: there m3/m4 go NaN and the NaN
    # arm fires, which is why an extreme-range test alone misses it. Only the
    # skew half pins the overflow -- kurt's m4 reaches NaN here too -- but both
    # are kept so the pair reads the same as the rest of the file.
    window = 20
    rng = np.random.default_rng(4)
    values = rng.normal(size=120)
    values[40] = 1e90

    result = getattr(pd.Series(values).rolling(window), roll_func)()

    tail = result.iloc[60:]
    assert tail.nunique() == len(tail)
    expected = getattr(pd.Series(values[-window:]), roll_func)()
    assert tail.iloc[-1] == pytest.approx(expected, rel=1e-12)
