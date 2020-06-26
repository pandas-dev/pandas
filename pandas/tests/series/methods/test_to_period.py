import numpy as np
import pytest

from pandas import (
    DataFrame,
    DatetimeIndex,
    PeriodIndex,
    Series,
    date_range,
    period_range,
)
import pandas._testing as tm


class TestToPeriod:
    def test_to_period(self):
        rng = date_range("1/1/2000", "1/1/2001", freq="D")
        ts = Series(np.random.randn(len(rng)), index=rng)

        pts = ts.to_period()
        exp = ts.copy()
        exp.index = period_range("1/1/2000", "1/1/2001")
        tm.assert_series_equal(pts, exp)

        pts = ts.to_period("M")
        exp.index = exp.index.asfreq("M")
        tm.assert_index_equal(pts.index, exp.index.asfreq("M"))
        tm.assert_series_equal(pts, exp)

        # GH#7606 without freq
        idx = DatetimeIndex(["2011-01-01", "2011-01-02", "2011-01-03", "2011-01-04"])
        exp_idx = PeriodIndex(
            ["2011-01-01", "2011-01-02", "2011-01-03", "2011-01-04"], freq="D"
        )

        s = Series(np.random.randn(4), index=idx)
        expected = s.copy()
        expected.index = exp_idx
        tm.assert_series_equal(s.to_period(), expected)

        df = DataFrame(np.random.randn(4, 4), index=idx, columns=idx)
        expected = df.copy()
        expected.index = exp_idx
        tm.assert_frame_equal(df.to_period(), expected)

        expected = df.copy()
        expected.columns = exp_idx
        tm.assert_frame_equal(df.to_period(axis=1), expected)

    def test_to_period_raises(self, index):
        # https://github.com/pandas-dev/pandas/issues/33327
        ser = Series(index=index, dtype=object)
        if not isinstance(index, DatetimeIndex):
            msg = f"unsupported Type {type(index).__name__}"
            with pytest.raises(TypeError, match=msg):
                ser.to_period()
