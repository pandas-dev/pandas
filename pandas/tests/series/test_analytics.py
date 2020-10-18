import operator

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm


class TestSeriesAnalytics:
    def test_matmul(self):
        # matmul test is for GH #10259
        a = Series(np.random.randn(4), index=["p", "q", "r", "s"])
        b = DataFrame(
            np.random.randn(3, 4), index=["1", "2", "3"], columns=["p", "q", "r", "s"]
        ).T

        # Series @ DataFrame -> Series
        result = operator.matmul(a, b)
        expected = Series(np.dot(a.values, b.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        # DataFrame @ Series -> Series
        result = operator.matmul(b.T, a)
        expected = Series(np.dot(b.T.values, a.T.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        # Series @ Series -> scalar
        result = operator.matmul(a, a)
        expected = np.dot(a.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # vector (1D np.array) @ Series (__rmatmul__)
        result = operator.matmul(a.values, a)
        expected = np.dot(a.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # vector (1D list) @ Series (__rmatmul__)
        result = operator.matmul(a.values.tolist(), a)
        expected = np.dot(a.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # matrix (2D np.array) @ Series (__rmatmul__)
        result = operator.matmul(b.T.values, a)
        expected = np.dot(b.T.values, a.values)
        tm.assert_almost_equal(result, expected)

        # GH 21530
        # matrix (2D nested lists) @ Series (__rmatmul__)
        result = operator.matmul(b.T.values.tolist(), a)
        expected = np.dot(b.T.values, a.values)
        tm.assert_almost_equal(result, expected)

        # mixed dtype DataFrame @ Series
        a["p"] = int(a.p)
        result = operator.matmul(b.T, a)
        expected = Series(np.dot(b.T.values, a.T.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        # different dtypes DataFrame @ Series
        a = a.astype(int)
        result = operator.matmul(b.T, a)
        expected = Series(np.dot(b.T.values, a.T.values), index=["1", "2", "3"])
        tm.assert_series_equal(result, expected)

        msg = r"Dot product shape mismatch, \(4,\) vs \(3,\)"
        # exception raised is of type Exception
        with pytest.raises(Exception, match=msg):
            a.dot(a.values[:3])
        msg = "matrices are not aligned"
        with pytest.raises(ValueError, match=msg):
            a.dot(b.T)

    def test_ptp(self):
        # GH21614
        N = 1000
        arr = np.random.randn(N)
        ser = Series(arr)
        assert np.ptp(ser) == np.ptp(arr)

    def test_repeat(self):
        s = Series(np.random.randn(3), index=["a", "b", "c"])

        reps = s.repeat(5)
        exp = Series(s.values.repeat(5), index=s.index.values.repeat(5))
        tm.assert_series_equal(reps, exp)

        to_rep = [2, 3, 4]
        reps = s.repeat(to_rep)
        exp = Series(s.values.repeat(to_rep), index=s.index.values.repeat(to_rep))
        tm.assert_series_equal(reps, exp)

    def test_numpy_repeat(self):
        s = Series(np.arange(3), name="x")
        expected = Series(s.values.repeat(2), name="x", index=s.index.values.repeat(2))
        tm.assert_series_equal(np.repeat(s, 2), expected)

        msg = "the 'axis' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            np.repeat(s, 2, axis=0)

    def test_is_monotonic(self):

        s = Series(np.random.randint(0, 10, size=1000))
        assert not s.is_monotonic
        s = Series(np.arange(1000))
        assert s.is_monotonic is True
        assert s.is_monotonic_increasing is True
        s = Series(np.arange(1000, 0, -1))
        assert s.is_monotonic_decreasing is True

        s = Series(pd.date_range("20130101", periods=10))
        assert s.is_monotonic is True
        assert s.is_monotonic_increasing is True
        s = Series(list(reversed(s.tolist())))
        assert s.is_monotonic is False
        assert s.is_monotonic_decreasing is True
