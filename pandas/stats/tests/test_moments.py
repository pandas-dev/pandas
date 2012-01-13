import unittest
import nose

from datetime import datetime
from numpy.random import randn
import numpy as np

from pandas.core.api import Series, DataFrame, DateRange
from pandas.util.testing import assert_almost_equal
import pandas.core.datetools as datetools
import pandas.stats.moments as mom
import pandas.util.testing as tm

N, K = 100, 10

class TestMoments(unittest.TestCase):

    _nan_locs = np.arange(20, 40)
    _inf_locs = np.array([])

    def setUp(self):
        arr = randn(N)
        arr[self._nan_locs] = np.NaN

        self.arr = arr
        self.rng = DateRange(datetime(2009, 1, 1), periods=N)

        self.series = Series(arr.copy(), index=self.rng)

        self.frame = DataFrame(randn(N, K), index=self.rng,
                               columns=np.arange(K))

    def test_rolling_sum(self):
        self._check_moment_func(mom.rolling_sum, np.sum)

    def test_rolling_count(self):
        counter = lambda x: np.isfinite(x).astype(float).sum()
        self._check_moment_func(mom.rolling_count, counter,
                                has_min_periods=False,
                                preserve_nan=False)

    def test_rolling_mean(self):
        self._check_moment_func(mom.rolling_mean, np.mean)

    def test_rolling_median(self):
        self._check_moment_func(mom.rolling_median, np.median)

    def test_rolling_min(self):
        self._check_moment_func(mom.rolling_min, np.min)

    def test_rolling_max(self):
        self._check_moment_func(mom.rolling_max, np.max)

    def test_rolling_quantile(self):
        qs = [.1, .5, .9]

        def scoreatpercentile(a, per):
            values = np.sort(a,axis=0)

            idx = per /1. * (values.shape[0] - 1)
            return values[int(idx)]

        for q in qs:
            def f(x, window, min_periods=None, time_rule=None):
                return mom.rolling_quantile(x, window, q,
                                                min_periods=min_periods,
                                                time_rule=time_rule)
            def alt(x):
                return scoreatpercentile(x, q)

            self._check_moment_func(f, alt)

    def test_rolling_apply(self):
        def roll_mean(x, window, min_periods=None, time_rule=None):
            return mom.rolling_apply(x, window,
                                         lambda x: x[np.isfinite(x)].mean(),
                                         min_periods=min_periods,
                                         time_rule=time_rule)
        self._check_moment_func(roll_mean, np.mean)

    def test_rolling_std(self):
        self._check_moment_func(mom.rolling_std,
                                lambda x: np.std(x, ddof=1))

    def test_rolling_var(self):
        self._check_moment_func(mom.rolling_var,
                                lambda x: np.var(x, ddof=1))

    def test_rolling_skew(self):
        try:
            from scipy.stats import skew
        except ImportError:
            raise nose.SkipTest('no scipy')
        self._check_moment_func(mom.rolling_skew,
                                lambda x: skew(x, bias=False))

    def test_rolling_kurt(self):
        try:
            from scipy.stats import kurtosis
        except ImportError:
            raise nose.SkipTest('no scipy')
        self._check_moment_func(mom.rolling_kurt,
                                lambda x: kurtosis(x, bias=False))

    def _check_moment_func(self, func, static_comp, window=50,
                           has_min_periods=True,
                           has_time_rule=True,
                           preserve_nan=True):

        self._check_ndarray(func, static_comp, window=window,
                            has_min_periods=has_min_periods,
                            preserve_nan=preserve_nan)

        self._check_structures(func, static_comp,
                               has_min_periods=has_min_periods,
                               has_time_rule=has_time_rule)

    def _check_ndarray(self, func, static_comp, window=50,
                       has_min_periods=True,
                       preserve_nan=True):

        result = func(self.arr, window)
        assert_almost_equal(result[-1],
                            static_comp(self.arr[-50:]))

        if preserve_nan:
            assert(np.isnan(result[self._nan_locs]).all())

        # excluding NaNs correctly
        arr = randn(50)
        arr[:10] = np.NaN
        arr[-10:] = np.NaN

        if has_min_periods:
            result = func(arr, 50, min_periods=30)
            assert_almost_equal(result[-1], static_comp(arr[10:-10]))

            # min_periods is working correctly
            result = func(arr, 20, min_periods=15)
            self.assert_(np.isnan(result[23]))
            self.assert_(not np.isnan(result[24]))

            self.assert_(not np.isnan(result[-6]))
            self.assert_(np.isnan(result[-5]))

            # min_periods=0
            result0 = func(arr, 20, min_periods=0)
            result1 = func(arr, 20, min_periods=1)
            assert_almost_equal(result0, result1)
        else:
            result = func(arr, 50)
            assert_almost_equal(result[-1], static_comp(arr[10:-10]))

    def _check_structures(self, func, static_comp,
                          has_min_periods=True, has_time_rule=True):

        series_result = func(self.series, 50)
        self.assert_(isinstance(series_result, Series))

        frame_result = func(self.frame, 50)
        self.assertEquals(type(frame_result), DataFrame)

        # check time_rule works
        if has_time_rule:
            win = 25
            minp = 10

            if has_min_periods:
                series_result = func(self.series[::2], win, min_periods=minp,
                                     time_rule='WEEKDAY')
                frame_result = func(self.frame[::2], win, min_periods=minp,
                                    time_rule='WEEKDAY')
            else:
                series_result = func(self.series[::2], win, time_rule='WEEKDAY')
                frame_result = func(self.frame[::2], win, time_rule='WEEKDAY')

            last_date = series_result.index[-1]
            prev_date = last_date - 24 * datetools.bday

            trunc_series = self.series[::2].truncate(prev_date, last_date)
            trunc_frame = self.frame[::2].truncate(prev_date, last_date)

            assert_almost_equal(series_result[-1], static_comp(trunc_series))

            assert_almost_equal(frame_result.xs(last_date),
                                trunc_frame.apply(static_comp))

    def test_ewma(self):
        self._check_ew(mom.ewma)

    def test_ewmvar(self):
        self._check_ew(mom.ewmvar)

    def test_ewmvol(self):
        self._check_ew(mom.ewmvol)

    def test_ewma_span_com_args(self):
        A = mom.ewma(self.arr, com=9.5)
        B = mom.ewma(self.arr, span=20)
        assert_almost_equal(A, B)

        self.assertRaises(Exception, mom.ewma, self.arr, com=9.5, span=20)
        self.assertRaises(Exception, mom.ewma, self.arr)

    def _check_ew(self, func):
        self._check_ew_ndarray(func)
        self._check_ew_structures(func)

    def _check_ew_ndarray(self, func, preserve_nan=False):
        result = func(self.arr, com=10)
        if preserve_nan:
            assert(np.isnan(result[self._nan_locs]).all())

        # excluding NaNs correctly
        arr = randn(50)
        arr[:10] = np.NaN
        arr[-10:] = np.NaN

        # ??? check something

        # pass in ints
        result2 = func(np.arange(50), span=10)
        self.assert_(result.dtype == np.float_)

    def _check_ew_structures(self, func):
        series_result = func(self.series, com=10)
        self.assert_(isinstance(series_result, Series))
        frame_result = func(self.frame, com=10)
        self.assertEquals(type(frame_result), DataFrame)

    # binary moments
    def test_rolling_cov(self):
        A = self.series
        B = A + randn(len(A))

        result = mom.rolling_cov(A, B, 50, min_periods=25)
        assert_almost_equal(result[-1], np.cov(A[-50:], B[-50:])[0, 1])

    def test_rolling_corr(self):
        A = self.series
        B = A + randn(len(A))

        result = mom.rolling_corr(A, B, 50, min_periods=25)
        assert_almost_equal(result[-1], np.corrcoef(A[-50:], B[-50:])[0, 1])

        # test for correct bias correction
        a = tm.makeTimeSeries()
        b = tm.makeTimeSeries()
        a[:5] = np.nan
        b[:10] = np.nan

        result = mom.rolling_corr(a, b, len(a), min_periods=1)
        assert_almost_equal(result[-1], a.corr(b))

    def test_rolling_corr_pairwise(self):
        panel = mom.rolling_corr_pairwise(self.frame, 10, min_periods=5)

        correl = panel.ix[:, 1, 5]
        exp = mom.rolling_corr(self.frame[1], self.frame[5],
                               10, min_periods=5)
        tm.assert_series_equal(correl, exp)

    def test_flex_binary_frame(self):
        def _check(method):
            series = self.frame[1]

            res = method(series, self.frame, 10)
            res2 = method(self.frame, series, 10)
            exp = self.frame.apply(lambda x: method(series, x, 10))

            tm.assert_frame_equal(res, exp)
            tm.assert_frame_equal(res2, exp)

            frame2 = self.frame.copy()
            frame2.values[:] = np.random.randn(*frame2.shape)

            res3 = method(self.frame, frame2, 10)
            exp = DataFrame(dict((k, method(self.frame[k], frame2[k], 10))
                                 for k in self.frame))
            tm.assert_frame_equal(res3, exp)

        methods = [mom.rolling_corr, mom.rolling_cov]
        for meth in methods:
            _check(meth)

    def test_ewmcov(self):
        self._check_binary_ew(mom.ewmcov)

    def test_ewmcorr(self):
        self._check_binary_ew(mom.ewmcorr)

    def _check_binary_ew(self, func):
        A = Series(randn(50), index=np.arange(50))
        B = A[2:] + randn(48)

        A[:10] = np.NaN
        B[-10:] = np.NaN

        result = func(A, B, 20, min_periods=5)

        self.assert_(np.isnan(result.values[:15]).all())
        self.assert_(not np.isnan(result.values[15:]).any())

        self.assertRaises(Exception, func, A, randn(50), 20, min_periods=5)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
