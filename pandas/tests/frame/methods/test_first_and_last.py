"""
Note: includes tests for `last`
"""
import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    bdate_range,
)
import pandas._testing as tm


class TestFirst:
    def test_first_subset(self, frame_or_series):
        ts = tm.makeTimeDataFrame(freq="12h")
        ts = tm.get_obj(ts, frame_or_series)
        result = ts.first("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(freq="D")
        ts = tm.get_obj(ts, frame_or_series)
        result = ts.first("10d")
        assert len(result) == 10

        result = ts.first("3M")
        expected = ts[:"3/31/2000"]
        tm.assert_equal(result, expected)

        result = ts.first("21D")
        expected = ts[:21]
        tm.assert_equal(result, expected)

        result = ts[:0].first("3M")
        tm.assert_equal(result, ts[:0])

    def test_first_last_raises(self, frame_or_series):
        # GH#20725
        obj = DataFrame([[1, 2, 3], [4, 5, 6]])
        obj = tm.get_obj(obj, frame_or_series)

        msg = "'first' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            obj.first("1D")

        msg = "'last' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            obj.last("1D")

    def test_last_subset(self, frame_or_series):
        ts = tm.makeTimeDataFrame(freq="12h")
        ts = tm.get_obj(ts, frame_or_series)
        result = ts.last("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(nper=30, freq="D")
        ts = tm.get_obj(ts, frame_or_series)
        result = ts.last("10d")
        assert len(result) == 10

        result = ts.last("21D")
        expected = ts["2000-01-10":]
        tm.assert_equal(result, expected)

        result = ts.last("21D")
        expected = ts[-21:]
        tm.assert_equal(result, expected)

        result = ts[:0].last("3M")
        tm.assert_equal(result, ts[:0])

    @pytest.mark.parametrize("start, periods", [("2010-03-31", 1), ("2010-03-30", 2)])
    def test_first_with_first_day_last_of_month(self, frame_or_series, start, periods):
        # GH#29623
        x = frame_or_series([1] * 100, index=bdate_range(start, periods=100))
        result = x.first("1M")
        expected = frame_or_series(
            [1] * periods, index=bdate_range(start, periods=periods)
        )
        tm.assert_equal(result, expected)

    def test_first_with_first_day_end_of_frq_n_greater_one(self, frame_or_series):
        # GH#29623
        x = frame_or_series([1] * 100, index=bdate_range("2010-03-31", periods=100))
        result = x.first("2M")
        expected = frame_or_series(
            [1] * 23, index=bdate_range("2010-03-31", "2010-04-30")
        )
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("func", ["first", "last"])
    def test_empty_not_input(self, func):
        # GH#51032
        df = DataFrame(index=pd.DatetimeIndex([]))
        result = getattr(df, func)(offset=1)
        tm.assert_frame_equal(df, result)
        assert df is not result

    @pytest.mark.parametrize("start, periods", [("2010-03-31", 1), ("2010-03-30", 2)])
    def test_last_day_of_months_with_date_offset(self, frame_or_series, start, periods):
        x = frame_or_series([1] * 100, index=pd.date_range(start, periods=100))
        result = x.first(pd.DateOffset(days=periods))
        expected = frame_or_series(
            [1] * periods, index=pd.date_range(start, periods=periods)
        )
        tm.assert_equal(result, expected)

    def test_date_offset_multiple_days(self, frame_or_series):
        x = frame_or_series([1] * 100, index=pd.date_range("2010-03-31", periods=100))
        result = x.first(pd.DateOffset(days=2))
        expected = frame_or_series(
            [1] * 2, index=pd.date_range("2010-03-31", "2010-04-01")
        )
        tm.assert_equal(result, expected)

    def test_first_with_date_offset(self):
        # GH#51284
        i = pd.to_datetime(["2018-04-09", "2018-04-10", "2018-04-11", "2018-04-12"])
        x = DataFrame({"A": [1, 2, 3, 4]}, index=i)
        result = x.first(pd.DateOffset(days=2))
        expected = DataFrame(
            {"A": [1, 2]}, index=pd.to_datetime(["2018-04-09", "2018-04-10"])
        )
        tm.assert_equal(result, expected)

    def test_date_offset_15_days(self):
        # GH#45908
        i = pd.date_range("2018-04-09", periods=30, freq="2D")
        x = DataFrame({"A": np.arange(30)}, index=i)
        result = x.first(pd.DateOffset(days=15))
        i2 = pd.date_range("2018-04-09", periods=8, freq="2D")
        expected = DataFrame({"A": np.arange(8)}, index=i2)
        tm.assert_equal(result, expected)

    def test_first_with_date_offset_months(self, frame_or_series):
        periods = 40
        x = frame_or_series(
            [1] * periods, index=pd.date_range("2010-03-31", periods=periods)
        )
        result = x.first(pd.DateOffset(months=1))
        expected = frame_or_series(
            [1] * 30, index=pd.date_range("2010-03-31", periods=30)
        )
        tm.assert_equal(result, expected)
