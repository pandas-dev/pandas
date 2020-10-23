from datetime import time

import numpy as np
import pytest

from pandas._libs.tslibs import timezones

from pandas import DataFrame, Series, date_range
import pandas._testing as tm


class TestAtTime:
    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_localized_at_time(self, tzstr):
        tz = timezones.maybe_get_tz(tzstr)

        rng = date_range("4/16/2012", "5/1/2012", freq="H")
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_local = ts.tz_localize(tzstr)

        result = ts_local.at_time(time(10, 0))
        expected = ts.at_time(time(10, 0)).tz_localize(tzstr)
        tm.assert_series_equal(result, expected)
        assert timezones.tz_compare(result.index.tz, tz)

    def test_at_time(self):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = Series(np.random.randn(len(rng)), index=rng)
        rs = ts.at_time(rng[1])
        assert (rs.index.hour == rng[1].hour).all()
        assert (rs.index.minute == rng[1].minute).all()
        assert (rs.index.second == rng[1].second).all()

        result = ts.at_time("9:30")
        expected = ts.at_time(time(9, 30))
        tm.assert_series_equal(result, expected)

        df = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts[time(9, 30)]
        result_df = df.loc[time(9, 30)]
        expected = ts[(rng.hour == 9) & (rng.minute == 30)]
        exp_df = df[(rng.hour == 9) & (rng.minute == 30)]

        result.index = result.index._with_freq(None)
        tm.assert_series_equal(result, expected)
        tm.assert_frame_equal(result_df, exp_df)

        chunk = df.loc["1/4/2000":]
        result = chunk.loc[time(9, 30)]
        expected = result_df[-1:]

        # Without resetting the freqs, these are 5 min and 1440 min, respectively
        result.index = result.index._with_freq(None)
        expected.index = expected.index._with_freq(None)
        tm.assert_frame_equal(result, expected)

        # midnight, everything
        rng = date_range("1/1/2000", "1/31/2000")
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.at_time(time(0, 0))
        tm.assert_series_equal(result, ts)

        # time doesn't exist
        rng = date_range("1/1/2012", freq="23Min", periods=384)
        ts = Series(np.random.randn(len(rng)), rng)
        rs = ts.at_time("16:00")
        assert len(rs) == 0

    def test_at_time_raises(self):
        # GH20725
        ser = Series("a b c".split())
        msg = "Index must be DatetimeIndex"
        with pytest.raises(TypeError, match=msg):
            ser.at_time("00:00")
