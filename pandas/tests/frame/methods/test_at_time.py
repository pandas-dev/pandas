from datetime import time

import numpy as np
import pytest
import pytz

from pandas import DataFrame, date_range
import pandas._testing as tm


class TestAtTime:
    def test_at_time(self):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        rs = ts.at_time(rng[1])
        assert (rs.index.hour == rng[1].hour).all()
        assert (rs.index.minute == rng[1].minute).all()
        assert (rs.index.second == rng[1].second).all()

        result = ts.at_time("9:30")
        expected = ts.at_time(time(9, 30))
        tm.assert_frame_equal(result, expected)

        result = ts.loc[time(9, 30)]
        expected = ts.loc[(rng.hour == 9) & (rng.minute == 30)]

        tm.assert_frame_equal(result, expected)

        # midnight, everything
        rng = date_range("1/1/2000", "1/31/2000")
        ts = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts.at_time(time(0, 0))
        tm.assert_frame_equal(result, ts)

        # time doesn't exist
        rng = date_range("1/1/2012", freq="23Min", periods=384)
        ts = DataFrame(np.random.randn(len(rng), 2), rng)
        rs = ts.at_time("16:00")
        assert len(rs) == 0

    @pytest.mark.parametrize(
        "hour", ["1:00", "1:00AM", time(1), time(1, tzinfo=pytz.UTC)]
    )
    def test_at_time_errors(self, hour):
        # GH#24043
        dti = date_range("2018", periods=3, freq="H")
        df = DataFrame(list(range(len(dti))), index=dti)
        if getattr(hour, "tzinfo", None) is None:
            result = df.at_time(hour)
            expected = df.iloc[1:2]
            tm.assert_frame_equal(result, expected)
        else:
            with pytest.raises(ValueError, match="Index must be timezone"):
                df.at_time(hour)

    def test_at_time_tz(self):
        # GH#24043
        dti = date_range("2018", periods=3, freq="H", tz="US/Pacific")
        df = DataFrame(list(range(len(dti))), index=dti)
        result = df.at_time(time(4, tzinfo=pytz.timezone("US/Eastern")))
        expected = df.iloc[1:2]
        tm.assert_frame_equal(result, expected)

    def test_at_time_raises(self):
        # GH#20725
        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        msg = "Index must be DatetimeIndex"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            df.at_time("00:00")

    @pytest.mark.parametrize("axis", ["index", "columns", 0, 1])
    def test_at_time_axis(self, axis):
        # issue 8839
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), len(rng)))
        ts.index, ts.columns = rng, rng

        indices = rng[(rng.hour == 9) & (rng.minute == 30) & (rng.second == 0)]

        if axis in ["index", 0]:
            expected = ts.loc[indices, :]
        elif axis in ["columns", 1]:
            expected = ts.loc[:, indices]

        result = ts.at_time("9:30", axis=axis)

        # Without clearing freq, result has freq 1440T and expected 5T
        result.index = result.index._with_freq(None)
        expected.index = expected.index._with_freq(None)
        tm.assert_frame_equal(result, expected)
