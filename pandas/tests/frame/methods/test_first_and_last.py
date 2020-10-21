"""
Note: includes tests for `last`
"""
import pytest

from pandas import DataFrame
import pandas._testing as tm


class TestFirst:
    def test_first_subset(self):
        ts = tm.makeTimeDataFrame(freq="12h")
        result = ts.first("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(freq="D")
        result = ts.first("10d")
        assert len(result) == 10

        result = ts.first("3M")
        expected = ts[:"3/31/2000"]
        tm.assert_frame_equal(result, expected)

        result = ts.first("21D")
        expected = ts[:21]
        tm.assert_frame_equal(result, expected)

        result = ts[:0].first("3M")
        tm.assert_frame_equal(result, ts[:0])

    def test_first_raises(self):
        # GH#20725
        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        msg = "'first' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            df.first("1D")

    def test_last_subset(self):
        ts = tm.makeTimeDataFrame(freq="12h")
        result = ts.last("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(nper=30, freq="D")
        result = ts.last("10d")
        assert len(result) == 10

        result = ts.last("21D")
        expected = ts["2000-01-10":]
        tm.assert_frame_equal(result, expected)

        result = ts.last("21D")
        expected = ts[-21:]
        tm.assert_frame_equal(result, expected)

        result = ts[:0].last("3M")
        tm.assert_frame_equal(result, ts[:0])

    def test_last_raises(self):
        # GH20725
        df = DataFrame([[1, 2, 3], [4, 5, 6]])
        msg = "'last' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):  # index is not a DatetimeIndex
            df.last("1D")
