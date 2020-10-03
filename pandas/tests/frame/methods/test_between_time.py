from datetime import time

import numpy as np
import pytest

from pandas import DataFrame, date_range
import pandas._testing as tm


class TestBetweenTime:
    def test_between_time(self, close_open_fixture):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(0, 0)
        etime = time(1, 0)
        inc_start, inc_end = close_open_fixture

        filtered = ts.between_time(stime, etime, inc_start, inc_end)
        exp_len = 13 * 4 + 1
        if not inc_start:
            exp_len -= 5
        if not inc_end:
            exp_len -= 4

        assert len(filtered) == exp_len
        for rs in filtered.index:
            t = rs.time()
            if inc_start:
                assert t >= stime
            else:
                assert t > stime

            if inc_end:
                assert t <= etime
            else:
                assert t < etime

        result = ts.between_time("00:00", "01:00")
        expected = ts.between_time(stime, etime)
        tm.assert_frame_equal(result, expected)

        # across midnight
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(22, 0)
        etime = time(9, 0)

        filtered = ts.between_time(stime, etime, inc_start, inc_end)
        exp_len = (12 * 11 + 1) * 4 + 1
        if not inc_start:
            exp_len -= 4
        if not inc_end:
            exp_len -= 4

        assert len(filtered) == exp_len
        for rs in filtered.index:
            t = rs.time()
            if inc_start:
                assert (t >= stime) or (t <= etime)
            else:
                assert (t > stime) or (t <= etime)

            if inc_end:
                assert (t <= etime) or (t >= stime)
            else:
                assert (t < etime) or (t >= stime)

    def test_between_time_raises(self):
        # GH#20725
        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        msg = "Index must be DatetimeIndex"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            df.between_time(start_time="00:00", end_time="12:00")

    def test_between_time_axis(self, axis):
        # GH#8839
        rng = date_range("1/1/2000", periods=100, freq="10min")
        ts = DataFrame(np.random.randn(len(rng), len(rng)))
        stime, etime = ("08:00:00", "09:00:00")
        exp_len = 7

        if axis in ["index", 0]:
            ts.index = rng
            assert len(ts.between_time(stime, etime)) == exp_len
            assert len(ts.between_time(stime, etime, axis=0)) == exp_len

        if axis in ["columns", 1]:
            ts.columns = rng
            selected = ts.between_time(stime, etime, axis=1).columns
            assert len(selected) == exp_len

    def test_between_time_axis_raises(self, axis):
        # issue 8839
        rng = date_range("1/1/2000", periods=100, freq="10min")
        mask = np.arange(0, len(rng))
        rand_data = np.random.randn(len(rng), len(rng))
        ts = DataFrame(rand_data, index=rng, columns=rng)
        stime, etime = ("08:00:00", "09:00:00")

        msg = "Index must be DatetimeIndex"
        if axis in ["columns", 1]:
            ts.index = mask
            with pytest.raises(TypeError, match=msg):
                ts.between_time(stime, etime)
            with pytest.raises(TypeError, match=msg):
                ts.between_time(stime, etime, axis=0)

        if axis in ["index", 0]:
            ts.columns = mask
            with pytest.raises(TypeError, match=msg):
                ts.between_time(stime, etime, axis=1)
