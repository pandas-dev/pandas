"""
Note: includes tests for `last`
"""
import pytest

from pandas import DataFrame
import pandas._testing as tm


class TestFirst:
    def test_first_subset(self, frame_or_series):
        ts = tm.makeTimeDataFrame(freq="12h")
        if frame_or_series is not DataFrame:
            ts = ts["A"]
        result = ts.first("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(freq="D")
        if frame_or_series is not DataFrame:
            ts = ts["A"]
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
        if frame_or_series is not DataFrame:
            obj = obj[0]

        msg = "'first' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            obj.first("1D")

        msg = "'last' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            obj.last("1D")

    def test_last_subset(self, frame_or_series):
        ts = tm.makeTimeDataFrame(freq="12h")
        if frame_or_series is not DataFrame:
            ts = ts["A"]
        result = ts.last("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(nper=30, freq="D")
        if frame_or_series is not DataFrame:
            ts = ts["A"]
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
