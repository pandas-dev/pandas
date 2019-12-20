"""
Tests for Series cumulative operations.

See also
--------
tests.frame.test_cumulative
"""
import numpy as np
import pytest

from pandas.compat.numpy import _np_version_under1p18

import pandas as pd
import pandas.util.testing as tm


def _check_accum_op(name, series, check_dtype=True):
    func = getattr(np, name)
    tm.assert_numpy_array_equal(
        func(series).values, func(np.array(series)), check_dtype=check_dtype,
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

    def test_cummin(self, datetime_series):
        tm.assert_numpy_array_equal(
            datetime_series.cummin().values,
            np.minimum.accumulate(np.array(datetime_series)),
        )
        ts = datetime_series.copy()
        ts[::2] = np.NaN
        result = ts.cummin()[1::2]
        expected = np.minimum.accumulate(ts.dropna())

        tm.assert_series_equal(result, expected)

    def test_cummax(self, datetime_series):
        tm.assert_numpy_array_equal(
            datetime_series.cummax().values,
            np.maximum.accumulate(np.array(datetime_series)),
        )
        ts = datetime_series.copy()
        ts[::2] = np.NaN
        result = ts.cummax()[1::2]
        expected = np.maximum.accumulate(ts.dropna())

        tm.assert_series_equal(result, expected)

    @pytest.mark.xfail(
        not _np_version_under1p18, reason="numpy 1.18 changed min/max behavior for NaT"
    )
    def test_cummin_datetime64(self):
        s = pd.Series(
            pd.to_datetime(["NaT", "2000-1-2", "NaT", "2000-1-1", "NaT", "2000-1-3"])
        )

        expected = pd.Series(
            pd.to_datetime(["NaT", "2000-1-2", "NaT", "2000-1-1", "NaT", "2000-1-1"])
        )
        result = s.cummin(skipna=True)
        tm.assert_series_equal(expected, result)

        expected = pd.Series(
            pd.to_datetime(
                ["NaT", "2000-1-2", "2000-1-2", "2000-1-1", "2000-1-1", "2000-1-1"]
            )
        )
        result = s.cummin(skipna=False)
        tm.assert_series_equal(expected, result)

    @pytest.mark.xfail(
        not _np_version_under1p18, reason="numpy 1.18 changed min/max behavior for NaT"
    )
    def test_cummax_datetime64(self):
        s = pd.Series(
            pd.to_datetime(["NaT", "2000-1-2", "NaT", "2000-1-1", "NaT", "2000-1-3"])
        )

        expected = pd.Series(
            pd.to_datetime(["NaT", "2000-1-2", "NaT", "2000-1-2", "NaT", "2000-1-3"])
        )
        result = s.cummax(skipna=True)
        tm.assert_series_equal(expected, result)

        expected = pd.Series(
            pd.to_datetime(
                ["NaT", "2000-1-2", "2000-1-2", "2000-1-2", "2000-1-2", "2000-1-3"]
            )
        )
        result = s.cummax(skipna=False)
        tm.assert_series_equal(expected, result)

    @pytest.mark.xfail(
        not _np_version_under1p18, reason="numpy 1.18 changed min/max behavior for NaT"
    )
    def test_cummin_timedelta64(self):
        s = pd.Series(pd.to_timedelta(["NaT", "2 min", "NaT", "1 min", "NaT", "3 min"]))

        expected = pd.Series(
            pd.to_timedelta(["NaT", "2 min", "NaT", "1 min", "NaT", "1 min"])
        )
        result = s.cummin(skipna=True)
        tm.assert_series_equal(expected, result)

        expected = pd.Series(
            pd.to_timedelta(["NaT", "2 min", "2 min", "1 min", "1 min", "1 min"])
        )
        result = s.cummin(skipna=False)
        tm.assert_series_equal(expected, result)

    @pytest.mark.xfail(
        not _np_version_under1p18, reason="numpy 1.18 changed min/max behavior for NaT"
    )
    def test_cummax_timedelta64(self):
        s = pd.Series(pd.to_timedelta(["NaT", "2 min", "NaT", "1 min", "NaT", "3 min"]))

        expected = pd.Series(
            pd.to_timedelta(["NaT", "2 min", "NaT", "2 min", "NaT", "3 min"])
        )
        result = s.cummax(skipna=True)
        tm.assert_series_equal(expected, result)

        expected = pd.Series(
            pd.to_timedelta(["NaT", "2 min", "2 min", "2 min", "2 min", "3 min"])
        )
        result = s.cummax(skipna=False)
        tm.assert_series_equal(expected, result)
