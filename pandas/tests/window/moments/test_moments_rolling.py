import copy
import warnings

import numpy as np
from numpy.random import randn
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, Series, isna, notna
import pandas._testing as tm

import pandas.tseries.offsets as offsets


def _check_moment_func(
    static_comp,
    name,
    raw,
    has_min_periods=True,
    has_center=True,
    has_time_rule=True,
    fill_value=None,
    zero_min_periods_equal=True,
    series=None,
    frame=None,
    **kwargs,
):

    # inject raw
    if name == "apply":
        kwargs = copy.copy(kwargs)
        kwargs["raw"] = raw

    def get_result(obj, window, min_periods=None, center=False):
        r = obj.rolling(window=window, min_periods=min_periods, center=center)
        return getattr(r, name)(**kwargs)

    series_result = get_result(series, window=50)
    assert isinstance(series_result, Series)
    tm.assert_almost_equal(series_result.iloc[-1], static_comp(series[-50:]))

    frame_result = get_result(frame, window=50)
    assert isinstance(frame_result, DataFrame)
    tm.assert_series_equal(
        frame_result.iloc[-1, :],
        frame.iloc[-50:, :].apply(static_comp, axis=0, raw=raw),
        check_names=False,
    )

    # check time_rule works
    if has_time_rule:
        win = 25
        minp = 10
        ser = series[::2].resample("B").mean()
        frm = frame[::2].resample("B").mean()

        if has_min_periods:
            series_result = get_result(ser, window=win, min_periods=minp)
            frame_result = get_result(frm, window=win, min_periods=minp)
        else:
            series_result = get_result(ser, window=win, min_periods=0)
            frame_result = get_result(frm, window=win, min_periods=0)

        last_date = series_result.index[-1]
        prev_date = last_date - 24 * offsets.BDay()

        trunc_series = series[::2].truncate(prev_date, last_date)
        trunc_frame = frame[::2].truncate(prev_date, last_date)

        tm.assert_almost_equal(series_result[-1], static_comp(trunc_series))

        tm.assert_series_equal(
            frame_result.xs(last_date),
            trunc_frame.apply(static_comp, raw=raw),
            check_names=False,
        )

    # excluding NaNs correctly
    obj = Series(randn(50))
    obj[:10] = np.NaN
    obj[-10:] = np.NaN
    if has_min_periods:
        result = get_result(obj, 50, min_periods=30)
        tm.assert_almost_equal(result.iloc[-1], static_comp(obj[10:-10]))

        # min_periods is working correctly
        result = get_result(obj, 20, min_periods=15)
        assert isna(result.iloc[23])
        assert not isna(result.iloc[24])

        assert not isna(result.iloc[-6])
        assert isna(result.iloc[-5])

        obj2 = Series(randn(20))
        result = get_result(obj2, 10, min_periods=5)
        assert isna(result.iloc[3])
        assert notna(result.iloc[4])

        if zero_min_periods_equal:
            # min_periods=0 may be equivalent to min_periods=1
            result0 = get_result(obj, 20, min_periods=0)
            result1 = get_result(obj, 20, min_periods=1)
            tm.assert_almost_equal(result0, result1)
    else:
        result = get_result(obj, 50)
        tm.assert_almost_equal(result.iloc[-1], static_comp(obj[10:-10]))

    # window larger than series length (#7297)
    if has_min_periods:
        for minp in (0, len(series) - 1, len(series)):
            result = get_result(series, len(series) + 1, min_periods=minp)
            expected = get_result(series, len(series), min_periods=minp)
            nan_mask = isna(result)
            tm.assert_series_equal(nan_mask, isna(expected))

            nan_mask = ~nan_mask
            tm.assert_almost_equal(result[nan_mask], expected[nan_mask])
    else:
        result = get_result(series, len(series) + 1, min_periods=0)
        expected = get_result(series, len(series), min_periods=0)
        nan_mask = isna(result)
        tm.assert_series_equal(nan_mask, isna(expected))

        nan_mask = ~nan_mask
        tm.assert_almost_equal(result[nan_mask], expected[nan_mask])

    # check center=True
    if has_center:
        if has_min_periods:
            result = get_result(obj, 20, min_periods=15, center=True)
            expected = get_result(
                pd.concat([obj, Series([np.NaN] * 9)]), 20, min_periods=15
            )[9:].reset_index(drop=True)
        else:
            result = get_result(obj, 20, min_periods=0, center=True)
            print(result)
            expected = get_result(
                pd.concat([obj, Series([np.NaN] * 9)]), 20, min_periods=0
            )[9:].reset_index(drop=True)

        tm.assert_series_equal(result, expected)

        # shifter index
        s = [f"x{x:d}" for x in range(12)]

        if has_min_periods:
            minp = 10

            series_xp = (
                get_result(
                    series.reindex(list(series.index) + s), window=25, min_periods=minp,
                )
                .shift(-12)
                .reindex(series.index)
            )
            frame_xp = (
                get_result(
                    frame.reindex(list(frame.index) + s), window=25, min_periods=minp,
                )
                .shift(-12)
                .reindex(frame.index)
            )

            series_rs = get_result(series, window=25, min_periods=minp, center=True)
            frame_rs = get_result(frame, window=25, min_periods=minp, center=True)

        else:
            series_xp = (
                get_result(
                    series.reindex(list(series.index) + s), window=25, min_periods=0,
                )
                .shift(-12)
                .reindex(series.index)
            )
            frame_xp = (
                get_result(
                    frame.reindex(list(frame.index) + s), window=25, min_periods=0,
                )
                .shift(-12)
                .reindex(frame.index)
            )

            series_rs = get_result(series, window=25, min_periods=0, center=True)
            frame_rs = get_result(frame, window=25, min_periods=0, center=True)

        if fill_value is not None:
            series_xp = series_xp.fillna(fill_value)
            frame_xp = frame_xp.fillna(fill_value)
        tm.assert_series_equal(series_xp, series_rs)
        tm.assert_frame_equal(frame_xp, frame_rs)


def test_centered_axis_validation():

    # ok
    Series(np.ones(10)).rolling(window=3, center=True, axis=0).mean()

    # bad axis
    with pytest.raises(ValueError):
        Series(np.ones(10)).rolling(window=3, center=True, axis=1).mean()

    # ok ok
    DataFrame(np.ones((10, 10))).rolling(window=3, center=True, axis=0).mean()
    DataFrame(np.ones((10, 10))).rolling(window=3, center=True, axis=1).mean()

    # bad axis
    with pytest.raises(ValueError):
        (DataFrame(np.ones((10, 10))).rolling(window=3, center=True, axis=2).mean())


def test_rolling_sum(raw, series, frame):
    _check_moment_func(
        np.nansum,
        name="sum",
        zero_min_periods_equal=False,
        raw=raw,
        series=series,
        frame=frame,
    )


def test_rolling_count(raw, series, frame):
    counter = lambda x: np.isfinite(x).astype(float).sum()
    _check_moment_func(
        counter,
        name="count",
        has_min_periods=False,
        fill_value=0,
        raw=raw,
        series=series,
        frame=frame,
    )


def test_rolling_mean(raw, series, frame):
    _check_moment_func(np.mean, name="mean", raw=raw, series=series, frame=frame)


@td.skip_if_no_scipy
def test_cmov_mean():
    # GH 8238
    vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48])
    result = Series(vals).rolling(5, center=True).mean()
    expected_values = [
        np.nan,
        np.nan,
        9.962,
        11.27,
        11.564,
        12.516,
        12.818,
        12.952,
        np.nan,
        np.nan,
    ]
    expected = Series(expected_values)
    tm.assert_series_equal(expected, result)


@td.skip_if_no_scipy
def test_cmov_window():
    # GH 8238
    vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48])
    result = Series(vals).rolling(5, win_type="boxcar", center=True).mean()
    expected_values = [
        np.nan,
        np.nan,
        9.962,
        11.27,
        11.564,
        12.516,
        12.818,
        12.952,
        np.nan,
        np.nan,
    ]
    expected = Series(expected_values)
    tm.assert_series_equal(expected, result)


@td.skip_if_no_scipy
def test_cmov_window_corner():
    # GH 8238
    # all nan
    vals = pd.Series([np.nan] * 10)
    result = vals.rolling(5, center=True, win_type="boxcar").mean()
    assert np.isnan(result).all()

    # empty
    vals = pd.Series([], dtype=object)
    result = vals.rolling(5, center=True, win_type="boxcar").mean()
    assert len(result) == 0

    # shorter than window
    vals = pd.Series(np.random.randn(5))
    result = vals.rolling(10, win_type="boxcar").mean()
    assert np.isnan(result).all()
    assert len(result) == 5


@td.skip_if_no_scipy
@pytest.mark.parametrize(
    "f,xp",
    [
        (
            "mean",
            [
                [np.nan, np.nan],
                [np.nan, np.nan],
                [9.252, 9.392],
                [8.644, 9.906],
                [8.87, 10.208],
                [6.81, 8.588],
                [7.792, 8.644],
                [9.05, 7.824],
                [np.nan, np.nan],
                [np.nan, np.nan],
            ],
        ),
        (
            "std",
            [
                [np.nan, np.nan],
                [np.nan, np.nan],
                [3.789706, 4.068313],
                [3.429232, 3.237411],
                [3.589269, 3.220810],
                [3.405195, 2.380655],
                [3.281839, 2.369869],
                [3.676846, 1.801799],
                [np.nan, np.nan],
                [np.nan, np.nan],
            ],
        ),
        (
            "var",
            [
                [np.nan, np.nan],
                [np.nan, np.nan],
                [14.36187, 16.55117],
                [11.75963, 10.48083],
                [12.88285, 10.37362],
                [11.59535, 5.66752],
                [10.77047, 5.61628],
                [13.51920, 3.24648],
                [np.nan, np.nan],
                [np.nan, np.nan],
            ],
        ),
        (
            "sum",
            [
                [np.nan, np.nan],
                [np.nan, np.nan],
                [46.26, 46.96],
                [43.22, 49.53],
                [44.35, 51.04],
                [34.05, 42.94],
                [38.96, 43.22],
                [45.25, 39.12],
                [np.nan, np.nan],
                [np.nan, np.nan],
            ],
        ),
    ],
)
def test_cmov_window_frame(f, xp):
    # Gh 8238
    df = DataFrame(
        np.array(
            [
                [12.18, 3.64],
                [10.18, 9.16],
                [13.24, 14.61],
                [4.51, 8.11],
                [6.15, 11.44],
                [9.14, 6.21],
                [11.31, 10.67],
                [2.94, 6.51],
                [9.42, 8.39],
                [12.44, 7.34],
            ]
        )
    )
    xp = DataFrame(np.array(xp))

    roll = df.rolling(5, win_type="boxcar", center=True)
    rs = getattr(roll, f)()

    tm.assert_frame_equal(xp, rs)


@td.skip_if_no_scipy
def test_cmov_window_na_min_periods():
    # min_periods
    vals = Series(np.random.randn(10))
    vals[4] = np.nan
    vals[8] = np.nan

    xp = vals.rolling(5, min_periods=4, center=True).mean()
    rs = vals.rolling(5, win_type="boxcar", min_periods=4, center=True).mean()
    tm.assert_series_equal(xp, rs)


@td.skip_if_no_scipy
def test_cmov_window_regular(win_types):
    # GH 8238
    vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48])
    xps = {
        "hamming": [
            np.nan,
            np.nan,
            8.71384,
            9.56348,
            12.38009,
            14.03687,
            13.8567,
            11.81473,
            np.nan,
            np.nan,
        ],
        "triang": [
            np.nan,
            np.nan,
            9.28667,
            10.34667,
            12.00556,
            13.33889,
            13.38,
            12.33667,
            np.nan,
            np.nan,
        ],
        "barthann": [
            np.nan,
            np.nan,
            8.4425,
            9.1925,
            12.5575,
            14.3675,
            14.0825,
            11.5675,
            np.nan,
            np.nan,
        ],
        "bohman": [
            np.nan,
            np.nan,
            7.61599,
            9.1764,
            12.83559,
            14.17267,
            14.65923,
            11.10401,
            np.nan,
            np.nan,
        ],
        "blackmanharris": [
            np.nan,
            np.nan,
            6.97691,
            9.16438,
            13.05052,
            14.02156,
            15.10512,
            10.74574,
            np.nan,
            np.nan,
        ],
        "nuttall": [
            np.nan,
            np.nan,
            7.04618,
            9.16786,
            13.02671,
            14.03559,
            15.05657,
            10.78514,
            np.nan,
            np.nan,
        ],
        "blackman": [
            np.nan,
            np.nan,
            7.73345,
            9.17869,
            12.79607,
            14.20036,
            14.57726,
            11.16988,
            np.nan,
            np.nan,
        ],
        "bartlett": [
            np.nan,
            np.nan,
            8.4425,
            9.1925,
            12.5575,
            14.3675,
            14.0825,
            11.5675,
            np.nan,
            np.nan,
        ],
    }

    xp = Series(xps[win_types])
    rs = Series(vals).rolling(5, win_type=win_types, center=True).mean()
    tm.assert_series_equal(xp, rs)


@td.skip_if_no_scipy
def test_cmov_window_regular_linear_range(win_types):
    # GH 8238
    vals = np.array(range(10), dtype=np.float)
    xp = vals.copy()
    xp[:2] = np.nan
    xp[-2:] = np.nan
    xp = Series(xp)

    rs = Series(vals).rolling(5, win_type=win_types, center=True).mean()
    tm.assert_series_equal(xp, rs)


@td.skip_if_no_scipy
def test_cmov_window_regular_missing_data(win_types):
    # GH 8238
    vals = np.array(
        [6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, np.nan, 10.63, 14.48]
    )
    xps = {
        "bartlett": [
            np.nan,
            np.nan,
            9.70333,
            10.5225,
            8.4425,
            9.1925,
            12.5575,
            14.3675,
            15.61667,
            13.655,
        ],
        "blackman": [
            np.nan,
            np.nan,
            9.04582,
            11.41536,
            7.73345,
            9.17869,
            12.79607,
            14.20036,
            15.8706,
            13.655,
        ],
        "barthann": [
            np.nan,
            np.nan,
            9.70333,
            10.5225,
            8.4425,
            9.1925,
            12.5575,
            14.3675,
            15.61667,
            13.655,
        ],
        "bohman": [
            np.nan,
            np.nan,
            8.9444,
            11.56327,
            7.61599,
            9.1764,
            12.83559,
            14.17267,
            15.90976,
            13.655,
        ],
        "hamming": [
            np.nan,
            np.nan,
            9.59321,
            10.29694,
            8.71384,
            9.56348,
            12.38009,
            14.20565,
            15.24694,
            13.69758,
        ],
        "nuttall": [
            np.nan,
            np.nan,
            8.47693,
            12.2821,
            7.04618,
            9.16786,
            13.02671,
            14.03673,
            16.08759,
            13.65553,
        ],
        "triang": [
            np.nan,
            np.nan,
            9.33167,
            9.76125,
            9.28667,
            10.34667,
            12.00556,
            13.82125,
            14.49429,
            13.765,
        ],
        "blackmanharris": [
            np.nan,
            np.nan,
            8.42526,
            12.36824,
            6.97691,
            9.16438,
            13.05052,
            14.02175,
            16.1098,
            13.65509,
        ],
    }

    xp = Series(xps[win_types])
    rs = Series(vals).rolling(5, win_type=win_types, min_periods=3).mean()
    tm.assert_series_equal(xp, rs)


@td.skip_if_no_scipy
def test_cmov_window_special(win_types_special):
    # GH 8238
    kwds = {
        "kaiser": {"beta": 1.0},
        "gaussian": {"std": 1.0},
        "general_gaussian": {"power": 2.0, "width": 2.0},
        "exponential": {"tau": 10},
    }

    vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48, 10.63, 14.48])

    xps = {
        "gaussian": [
            np.nan,
            np.nan,
            8.97297,
            9.76077,
            12.24763,
            13.89053,
            13.65671,
            12.01002,
            np.nan,
            np.nan,
        ],
        "general_gaussian": [
            np.nan,
            np.nan,
            9.85011,
            10.71589,
            11.73161,
            13.08516,
            12.95111,
            12.74577,
            np.nan,
            np.nan,
        ],
        "kaiser": [
            np.nan,
            np.nan,
            9.86851,
            11.02969,
            11.65161,
            12.75129,
            12.90702,
            12.83757,
            np.nan,
            np.nan,
        ],
        "exponential": [
            np.nan,
            np.nan,
            9.83364,
            11.10472,
            11.64551,
            12.66138,
            12.92379,
            12.83770,
            np.nan,
            np.nan,
        ],
    }

    xp = Series(xps[win_types_special])
    rs = (
        Series(vals)
        .rolling(5, win_type=win_types_special, center=True)
        .mean(**kwds[win_types_special])
    )
    tm.assert_series_equal(xp, rs)


@td.skip_if_no_scipy
def test_cmov_window_special_linear_range(win_types_special):
    # GH 8238
    kwds = {
        "kaiser": {"beta": 1.0},
        "gaussian": {"std": 1.0},
        "general_gaussian": {"power": 2.0, "width": 2.0},
        "slepian": {"width": 0.5},
        "exponential": {"tau": 10},
    }

    vals = np.array(range(10), dtype=np.float)
    xp = vals.copy()
    xp[:2] = np.nan
    xp[-2:] = np.nan
    xp = Series(xp)

    rs = (
        Series(vals)
        .rolling(5, win_type=win_types_special, center=True)
        .mean(**kwds[win_types_special])
    )
    tm.assert_series_equal(xp, rs)


def test_rolling_median(raw, series, frame):
    _check_moment_func(np.median, name="median", raw=raw, series=series, frame=frame)


def test_rolling_min(raw, series, frame):
    _check_moment_func(np.min, name="min", raw=raw, series=series, frame=frame)

    a = pd.Series([1, 2, 3, 4, 5])
    result = a.rolling(window=100, min_periods=1).min()
    expected = pd.Series(np.ones(len(a)))
    tm.assert_series_equal(result, expected)

    with pytest.raises(ValueError):
        pd.Series([1, 2, 3]).rolling(window=3, min_periods=5).min()


def test_rolling_max(raw, series, frame):
    _check_moment_func(np.max, name="max", raw=raw, series=series, frame=frame)

    a = pd.Series([1, 2, 3, 4, 5], dtype=np.float64)
    b = a.rolling(window=100, min_periods=1).max()
    tm.assert_almost_equal(a, b)

    with pytest.raises(ValueError):
        pd.Series([1, 2, 3]).rolling(window=3, min_periods=5).max()


@pytest.mark.parametrize("q", [0.0, 0.1, 0.5, 0.9, 1.0])
def test_rolling_quantile(q, raw, series, frame):
    def scoreatpercentile(a, per):
        values = np.sort(a, axis=0)

        idx = int(per / 1.0 * (values.shape[0] - 1))

        if idx == values.shape[0] - 1:
            retval = values[-1]

        else:
            qlow = float(idx) / float(values.shape[0] - 1)
            qhig = float(idx + 1) / float(values.shape[0] - 1)
            vlow = values[idx]
            vhig = values[idx + 1]
            retval = vlow + (vhig - vlow) * (per - qlow) / (qhig - qlow)

        return retval

    def quantile_func(x):
        return scoreatpercentile(x, q)

    _check_moment_func(
        quantile_func, name="quantile", quantile=q, raw=raw, series=series, frame=frame
    )


def test_rolling_quantile_np_percentile():
    # #9413: Tests that rolling window's quantile default behavior
    # is analogous to Numpy's percentile
    row = 10
    col = 5
    idx = pd.date_range("20100101", periods=row, freq="B")
    df = DataFrame(np.random.rand(row * col).reshape((row, -1)), index=idx)

    df_quantile = df.quantile([0.25, 0.5, 0.75], axis=0)
    np_percentile = np.percentile(df, [25, 50, 75], axis=0)

    tm.assert_almost_equal(df_quantile.values, np.array(np_percentile))


@pytest.mark.parametrize("quantile", [0.0, 0.1, 0.45, 0.5, 1])
@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "nearest", "midpoint"]
)
@pytest.mark.parametrize(
    "data",
    [
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        [8.0, 1.0, 3.0, 4.0, 5.0, 2.0, 6.0, 7.0],
        [0.0, np.nan, 0.2, np.nan, 0.4],
        [np.nan, np.nan, np.nan, np.nan],
        [np.nan, 0.1, np.nan, 0.3, 0.4, 0.5],
        [0.5],
        [np.nan, 0.7, 0.6],
    ],
)
def test_rolling_quantile_interpolation_options(quantile, interpolation, data):
    # Tests that rolling window's quantile behavior is analogous to
    # Series' quantile for each interpolation option
    s = Series(data)

    q1 = s.quantile(quantile, interpolation)
    q2 = s.expanding(min_periods=1).quantile(quantile, interpolation).iloc[-1]

    if np.isnan(q1):
        assert np.isnan(q2)
    else:
        assert q1 == q2


def test_invalid_quantile_value():
    data = np.arange(5)
    s = Series(data)

    msg = "Interpolation 'invalid' is not supported"
    with pytest.raises(ValueError, match=msg):
        s.rolling(len(data), min_periods=1).quantile(0.5, interpolation="invalid")


def test_rolling_quantile_param():
    ser = Series([0.0, 0.1, 0.5, 0.9, 1.0])

    with pytest.raises(ValueError):
        ser.rolling(3).quantile(-0.1)

    with pytest.raises(ValueError):
        ser.rolling(3).quantile(10.0)

    with pytest.raises(TypeError):
        ser.rolling(3).quantile("foo")


def test_rolling_apply(raw, series, frame):
    # suppress warnings about empty slices, as we are deliberately testing
    # with a 0-length Series

    def f(x):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*(empty slice|0 for slice).*",
                category=RuntimeWarning,
            )
            return x[np.isfinite(x)].mean()

    _check_moment_func(
        np.mean, name="apply", func=f, raw=raw, series=series, frame=frame
    )


def test_rolling_std(raw, series, frame):
    _check_moment_func(
        lambda x: np.std(x, ddof=1), name="std", raw=raw, series=series, frame=frame
    )
    _check_moment_func(
        lambda x: np.std(x, ddof=0),
        name="std",
        ddof=0,
        raw=raw,
        series=series,
        frame=frame,
    )


def test_rolling_std_1obs():
    vals = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])

    result = vals.rolling(1, min_periods=1).std()
    expected = pd.Series([np.nan] * 5)
    tm.assert_series_equal(result, expected)

    result = vals.rolling(1, min_periods=1).std(ddof=0)
    expected = pd.Series([0.0] * 5)
    tm.assert_series_equal(result, expected)

    result = pd.Series([np.nan, np.nan, 3, 4, 5]).rolling(3, min_periods=2).std()
    assert np.isnan(result[2])


def test_rolling_std_neg_sqrt():
    # unit test from Bottleneck

    # Test move_nanstd for neg sqrt.

    a = pd.Series(
        [
            0.0011448196318903589,
            0.00028718669878572767,
            0.00028718669878572767,
            0.00028718669878572767,
            0.00028718669878572767,
        ]
    )
    b = a.rolling(window=3).std()
    assert np.isfinite(b[2:]).all()

    b = a.ewm(span=3).std()
    assert np.isfinite(b[2:]).all()


def test_rolling_var(raw, series, frame):
    _check_moment_func(
        lambda x: np.var(x, ddof=1), name="var", raw=raw, series=series, frame=frame
    )
    _check_moment_func(
        lambda x: np.var(x, ddof=0),
        name="var",
        ddof=0,
        raw=raw,
        series=series,
        frame=frame,
    )


@td.skip_if_no_scipy
def test_rolling_skew(raw, series, frame):
    from scipy.stats import skew

    _check_moment_func(
        lambda x: skew(x, bias=False), name="skew", raw=raw, series=series, frame=frame
    )


@td.skip_if_no_scipy
def test_rolling_kurt(raw, series, frame):
    from scipy.stats import kurtosis

    _check_moment_func(
        lambda x: kurtosis(x, bias=False),
        name="kurt",
        raw=raw,
        series=series,
        frame=frame,
    )
