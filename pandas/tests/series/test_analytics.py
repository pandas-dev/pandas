import operator

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm


class TestSeriesAnalytics:
    def test_prod_numpy16_bug(self):
        s = Series([1.0, 1.0, 1.0], index=range(3))
        result = s.prod()

        assert not isinstance(result, Series)

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

    @pytest.mark.parametrize("func", [np.any, np.all])
    @pytest.mark.parametrize("kwargs", [dict(keepdims=True), dict(out=object())])
    @td.skip_if_np_lt("1.15")
    def test_validate_any_all_out_keepdims_raises(self, kwargs, func):
        s = pd.Series([1, 2])
        param = list(kwargs)[0]
        name = func.__name__

        msg = (
            f"the '{param}' parameter is not "
            "supported in the pandas "
            fr"implementation of {name}\(\)"
        )
        with pytest.raises(ValueError, match=msg):
            func(s, **kwargs)

    @td.skip_if_np_lt("1.15")
    def test_validate_sum_initial(self):
        s = pd.Series([1, 2])
        msg = (
            r"the 'initial' parameter is not "
            r"supported in the pandas "
            r"implementation of sum\(\)"
        )
        with pytest.raises(ValueError, match=msg):
            np.sum(s, initial=10)

    def test_validate_median_initial(self):
        s = pd.Series([1, 2])
        msg = (
            r"the 'overwrite_input' parameter is not "
            r"supported in the pandas "
            r"implementation of median\(\)"
        )
        with pytest.raises(ValueError, match=msg):
            # It seems like np.median doesn't dispatch, so we use the
            # method instead of the ufunc.
            s.median(overwrite_input=True)

    @td.skip_if_np_lt("1.15")
    def test_validate_stat_keepdims(self):
        s = pd.Series([1, 2])
        msg = (
            r"the 'keepdims' parameter is not "
            r"supported in the pandas "
            r"implementation of sum\(\)"
        )
        with pytest.raises(ValueError, match=msg):
            np.sum(s, keepdims=True)

    def test_td64_summation_overflow(self):
        # GH 9442
        s = pd.Series(pd.date_range("20130101", periods=100000, freq="H"))
        s[0] += pd.Timedelta("1s 1ms")

        # mean
        result = (s - s.min()).mean()
        expected = pd.Timedelta((pd.TimedeltaIndex((s - s.min())).asi8 / len(s)).sum())

        # the computation is converted to float so
        # might be some loss of precision
        assert np.allclose(result.value / 1000, expected.value / 1000)

        # sum
        msg = "overflow in timedelta operation"
        with pytest.raises(ValueError, match=msg):
            (s - s.min()).sum()

        s1 = s[0:10000]
        with pytest.raises(ValueError, match=msg):
            (s1 - s1.min()).sum()
        s2 = s[0:1000]
        (s2 - s2.min()).sum()
