"""
Tests for Series cumulative operations.

See also
--------
tests.frame.test_cumulative
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

methods = {
    "cumsum": np.cumsum,
    "cumprod": np.cumprod,
    "cummin": np.minimum.accumulate,
    "cummax": np.maximum.accumulate,
}


def _check_accum_op(name, series, check_dtype=True):
    func = getattr(np, name)
    tm.assert_numpy_array_equal(
        func(series).values, func(np.array(series)), check_dtype=check_dtype
    )

    # with missing values
    ts = series.copy()
    ts[::2] = np.NaN

    result = func(ts)[1::2]
    expected = func(np.array(ts.dropna()))

    tm.assert_numpy_array_equal(result.values, expected, check_dtype=False)


class TestSeriesCumulativeOps:
    def test_cumsum(self, datetime_series):
        _check_accum_op("cumsum", datetime_series)

    def test_cumprod(self, datetime_series):
        _check_accum_op("cumprod", datetime_series)

    @pytest.mark.parametrize("method", ["cummin", "cummax"])
    def test_cummin_cummax(self, datetime_series, method):
        ufunc = methods[method]

        result = getattr(datetime_series, method)().values
        expected = ufunc(np.array(datetime_series))

        tm.assert_numpy_array_equal(result, expected)
        ts = datetime_series.copy()
        ts[::2] = np.NaN
        result = getattr(ts, method)()[1::2]
        expected = ufunc(ts.dropna())

        result.index = result.index._with_freq(None)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "ts",
        [
            pd.Timedelta(0),
            pd.Timestamp("1999-12-31"),
            pd.Timestamp("1999-12-31").tz_localize("US/Pacific"),
        ],
    )
    def test_cummin_cummax_datetimelike(self, ts):
        # with ts==pd.Timedelta(0), we are testing td64; with naive Timestamp
        #  we are testing datetime64[ns]; with Timestamp[US/Pacific]
        #  we are testing dt64tz
        tdi = pd.to_timedelta(["NaT", "2 days", "NaT", "1 days", "NaT", "3 days"])
        ser = pd.Series(tdi + ts)

        exp_tdi = pd.to_timedelta(["NaT", "2 days", "NaT", "2 days", "NaT", "3 days"])
        expected = pd.Series(exp_tdi + ts)
        result = ser.cummax(skipna=True)
        tm.assert_series_equal(expected, result)

        exp_tdi = pd.to_timedelta(["NaT", "2 days", "NaT", "1 days", "NaT", "1 days"])
        expected = pd.Series(exp_tdi + ts)
        result = ser.cummin(skipna=True)
        tm.assert_series_equal(expected, result)

        exp_tdi = pd.to_timedelta(
            ["NaT", "2 days", "2 days", "2 days", "2 days", "3 days"]
        )
        expected = pd.Series(exp_tdi + ts)
        result = ser.cummax(skipna=False)
        tm.assert_series_equal(expected, result)

        exp_tdi = pd.to_timedelta(
            ["NaT", "2 days", "2 days", "1 days", "1 days", "1 days"]
        )
        expected = pd.Series(exp_tdi + ts)
        result = ser.cummin(skipna=False)
        tm.assert_series_equal(expected, result)

    def test_cummethods_bool(self):
        # GH#6270
        # checking Series method vs the ufunc applied to the values

        a = pd.Series([False, False, False, True, True, False, False])
        c = pd.Series([False] * len(a))

        for method in methods:
            for ser in [a, ~a, c, ~c]:
                ufunc = methods[method]

                exp_vals = ufunc(ser.values)
                expected = pd.Series(exp_vals)

                result = getattr(ser, method)()

                tm.assert_series_equal(result, expected)

    def test_cummethods_bool_in_object_dtype(self):

        ser = pd.Series([False, True, np.nan, False])
        cse = pd.Series([0, 1, np.nan, 1], dtype=object)
        cpe = pd.Series([False, 0, np.nan, 0])
        cmin = pd.Series([False, False, np.nan, False])
        cmax = pd.Series([False, True, np.nan, True])
        expecteds = {"cumsum": cse, "cumprod": cpe, "cummin": cmin, "cummax": cmax}

        for method in methods:
            res = getattr(ser, method)()
            tm.assert_series_equal(res, expecteds[method])
