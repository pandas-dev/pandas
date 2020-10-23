"""
Note: includes tests for `last`
"""

import numpy as np
import pytest

from pandas import Series, date_range
import pandas._testing as tm


class TestFirst:
    def test_first_subset(self):
        rng = date_range("1/1/2000", "1/1/2010", freq="12h")
        ts = Series(np.random.randn(len(rng)), index=rng)
        result = ts.first("10d")
        assert len(result) == 20

        rng = date_range("1/1/2000", "1/1/2010", freq="D")
        ts = Series(np.random.randn(len(rng)), index=rng)
        result = ts.first("10d")
        assert len(result) == 10

        result = ts.first("3M")
        expected = ts[:"3/31/2000"]
        tm.assert_series_equal(result, expected)

        result = ts.first("21D")
        expected = ts[:21]
        tm.assert_series_equal(result, expected)

        result = ts[:0].first("3M")
        tm.assert_series_equal(result, ts[:0])

    def test_first_raises(self):
        # GH#20725
        ser = Series("a b c".split())
        msg = "'first' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):
            ser.first("1D")

    def test_last_subset(self):
        rng = date_range("1/1/2000", "1/1/2010", freq="12h")
        ts = Series(np.random.randn(len(rng)), index=rng)
        result = ts.last("10d")
        assert len(result) == 20

        rng = date_range("1/1/2000", "1/1/2010", freq="D")
        ts = Series(np.random.randn(len(rng)), index=rng)
        result = ts.last("10d")
        assert len(result) == 10

        result = ts.last("21D")
        expected = ts["12/12/2009":]
        tm.assert_series_equal(result, expected)

        result = ts.last("21D")
        expected = ts[-21:]
        tm.assert_series_equal(result, expected)

        result = ts[:0].last("3M")
        tm.assert_series_equal(result, ts[:0])

    def test_last_raises(self):
        # GH#20725
        ser = Series("a b c".split())
        msg = "'last' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):
            ser.last("1D")
