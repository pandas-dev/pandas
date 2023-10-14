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
    def test_to_period(self, frame_or_series):
        K = 5

        dr = date_range("1/1/2000", "1/1/2001", freq="D")
        obj = DataFrame(
            np.random.default_rng(2).standard_normal((len(dr), K)),
            index=dr,
            columns=["A", "B", "C", "D", "E"],
        )
        obj["mix"] = "a"
        obj = tm.get_obj(obj, frame_or_series)

        msg = "(Series|DataFrame).to_period is deprecated"
        with tm.assert_produces_warning(FutureWarning, match=msg):
            pts = obj.to_period()
        exp = obj.copy()
        exp.index = period_range("1/1/2000", "1/1/2001")
        tm.assert_equal(pts, exp)

        with tm.assert_produces_warning(FutureWarning, match=msg):
            pts = obj.to_period("M")
        exp.index = exp.index.asfreq("M")
        tm.assert_equal(pts, exp)

    def test_to_period_without_freq(self, frame_or_series):
        # GH#7606 without freq
        idx = DatetimeIndex(["2011-01-01", "2011-01-02", "2011-01-03", "2011-01-04"])
        exp_idx = PeriodIndex(
            ["2011-01-01", "2011-01-02", "2011-01-03", "2011-01-04"], freq="D"
        )

        obj = DataFrame(
            np.random.default_rng(2).standard_normal((4, 4)), index=idx, columns=idx
        )
        obj = tm.get_obj(obj, frame_or_series)
        expected = obj.copy()
        expected.index = exp_idx
        msg = f"{type(expected).__name__}.to_period is deprecated"
        with tm.assert_produces_warning(FutureWarning, match=msg):
            res = obj.to_period()
        tm.assert_equal(res, expected)

        if frame_or_series is DataFrame:
            expected = obj.copy()
            expected.columns = exp_idx
            with tm.assert_produces_warning(FutureWarning, match=msg):
                res = obj.to_period(axis=1)
            tm.assert_frame_equal(res, expected)

    def test_to_period_columns(self):
        dr = date_range("1/1/2000", "1/1/2001")
        df = DataFrame(np.random.default_rng(2).standard_normal((len(dr), 5)), index=dr)
        df["mix"] = "a"

        df = df.T
        msg = "DataFrame.to_period is deprecated"
        with tm.assert_produces_warning(FutureWarning, match=msg):
            pts = df.to_period(axis=1)
        exp = df.copy()
        exp.columns = period_range("1/1/2000", "1/1/2001")
        tm.assert_frame_equal(pts, exp)

        with tm.assert_produces_warning(FutureWarning, match=msg):
            pts = df.to_period("M", axis=1)
        tm.assert_index_equal(pts.columns, exp.columns.asfreq("M"))

    def test_to_period_invalid_axis(self):
        dr = date_range("1/1/2000", "1/1/2001")
        df = DataFrame(np.random.default_rng(2).standard_normal((len(dr), 5)), index=dr)
        df["mix"] = "a"

        msg = "No axis named 2 for object type DataFrame"
        depr_msg = "DataFrame.to_period is deprecated"
        with pytest.raises(ValueError, match=msg):
            with tm.assert_produces_warning(FutureWarning, match=depr_msg):
                df.to_period(axis=2)

    def test_to_period_raises(self, index, frame_or_series):
        # https://github.com/pandas-dev/pandas/issues/33327
        obj = Series(index=index, dtype=object)
        if frame_or_series is DataFrame:
            obj = obj.to_frame()

        depr_msg = rf"{type(obj).__name__}\.to_period is deprecated"
        if not isinstance(index, DatetimeIndex):
            msg = f"unsupported Type {type(index).__name__}"
            with pytest.raises(TypeError, match=msg):
                with tm.assert_produces_warning(FutureWarning, match=depr_msg):
                    obj.to_period()
