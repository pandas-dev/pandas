import unittest
import nose
import sys
import functools

from datetime import datetime
from numpy.random import randn
import numpy as np

from pandas import Series, DataFrame, bdate_range, isnull, notnull
from pandas.util.testing import (
    assert_almost_equal, assert_series_equal, assert_frame_equal
)
from pandas.util.py3compat import PY3
import pandas.core.datetools as datetools
import pandas.stats.moments as mom
import pandas.util.testing as tm

N, K = 100, 10


class TestMoments(unittest.TestCase):

    _multiprocess_can_split_ = True

    _nan_locs = np.arange(20, 40)
    _inf_locs = np.array([])

    def setUp(self):
        arr = randn(N)
        arr[self._nan_locs] = np.NaN

        self.arr = arr
        self.rng = bdate_range(datetime(2009, 1, 1), periods=N)

        self.series = Series(arr.copy(), index=self.rng)

        self.frame = DataFrame(randn(N, K), index=self.rng,
                               columns=np.arange(K))

    def test_centered_axis_validation(self):
        # ok
        mom.rolling_mean(Series(np.ones(10)),3,center=True ,axis=0)
        # bad axis
        self.assertRaises(ValueError, mom.rolling_mean,Series(np.ones(10)),3,center=True ,axis=1)

        # ok ok
        mom.rolling_mean(DataFrame(np.ones((10,10))),3,center=True ,axis=0)
        mom.rolling_mean(DataFrame(np.ones((10,10))),3,center=True ,axis=1)
        # bad axis
        self.assertRaises(ValueError, mom.rolling_mean,DataFrame(np.ones((10,10))),3,center=True ,axis=2)

    def test_rolling_sum(self):
        self._check_moment_func(mom.rolling_sum, np.sum)

    def test_rolling_count(self):
        counter = lambda x: np.isfinite(x).astype(float).sum()
        self._check_moment_func(mom.rolling_count, counter,
                                has_min_periods=False,
                                preserve_nan=False,
                                fill_value=0)

    def test_rolling_mean(self):
        self._check_moment_func(mom.rolling_mean, np.mean)

    def test_cmov_mean(self):
        try:
            from scikits.timeseries.lib import cmov_mean
        except ImportError:
            raise nose.SkipTest

        vals = np.random.randn(10)
        xp = cmov_mean(vals, 5)

        rs = mom.rolling_mean(vals, 5, center=True)
        assert_almost_equal(xp.compressed(), rs[2:-2])
        assert_almost_equal(xp.mask, np.isnan(rs))

        xp = Series(rs)
        rs = mom.rolling_mean(Series(vals), 5, center=True)
        assert_series_equal(xp, rs)

    def test_cmov_window(self):
        try:
            from scikits.timeseries.lib import cmov_window
        except ImportError:
            raise nose.SkipTest

        vals = np.random.randn(10)
        xp = cmov_window(vals, 5, 'boxcar')

        rs = mom.rolling_window(vals, 5, 'boxcar', center=True)
        assert_almost_equal(xp.compressed(), rs[2:-2])
        assert_almost_equal(xp.mask, np.isnan(rs))

        xp = Series(rs)
        rs = mom.rolling_window(Series(vals), 5, 'boxcar', center=True)
        assert_series_equal(xp, rs)

    def test_cmov_window_corner(self):
        try:
            from scikits.timeseries.lib import cmov_window
        except ImportError:
            raise nose.SkipTest

        # all nan
        vals = np.empty(10, dtype=float)
        vals.fill(np.nan)
        rs = mom.rolling_window(vals, 5, 'boxcar', center=True)
        self.assert_(np.isnan(rs).all())

        # empty
        vals = np.array([])
        rs = mom.rolling_window(vals, 5, 'boxcar', center=True)
        self.assert_(len(rs) == 0)

        # shorter than window
        vals = np.random.randn(5)
        rs = mom.rolling_window(vals, 10, 'boxcar')
        self.assert_(np.isnan(rs).all())
        self.assert_(len(rs) == 5)

    def test_cmov_window_frame(self):
        try:
            from scikits.timeseries.lib import cmov_window
        except ImportError:
            raise nose.SkipTest

        # DataFrame
        vals = np.random.randn(10, 2)
        xp = cmov_window(vals, 5, 'boxcar')
        rs = mom.rolling_window(DataFrame(vals), 5, 'boxcar', center=True)
        assert_frame_equal(DataFrame(xp), rs)

    def test_cmov_window_na_min_periods(self):
        try:
            from scikits.timeseries.lib import cmov_window
        except ImportError:
            raise nose.SkipTest

        # min_periods
        vals = Series(np.random.randn(10))
        vals[4] = np.nan
        vals[8] = np.nan

        xp = mom.rolling_mean(vals, 5, min_periods=4, center=True)
        rs = mom.rolling_window(vals, 5, 'boxcar', min_periods=4, center=True)

        assert_series_equal(xp, rs)

    def test_cmov_window_regular(self):
        try:
            from scikits.timeseries.lib import cmov_window
        except ImportError:
            raise nose.SkipTest

        win_types = ['triang', 'blackman', 'hamming', 'bartlett', 'bohman',
                     'blackmanharris', 'nuttall', 'barthann']
        for wt in win_types:
            vals = np.random.randn(10)
            xp = cmov_window(vals, 5, wt)

            rs = mom.rolling_window(Series(vals), 5, wt, center=True)
            assert_series_equal(Series(xp), rs)

    def test_cmov_window_special(self):
        try:
            from scikits.timeseries.lib import cmov_window
        except ImportError:
            raise nose.SkipTest

        win_types = ['kaiser', 'gaussian', 'general_gaussian', 'slepian']
        kwds = [{'beta': 1.}, {'std': 1.}, {'power': 2., 'width': 2.},
                {'width': 0.5}]

        for wt, k in zip(win_types, kwds):
            vals = np.random.randn(10)
            xp = cmov_window(vals, 5, (wt,) + tuple(k.values()))

            rs = mom.rolling_window(Series(vals), 5, wt, center=True,
                                    **k)
            assert_series_equal(Series(xp), rs)

    def test_rolling_median(self):
        self._check_moment_func(mom.rolling_median, np.median)

    def test_rolling_min(self):
        self._check_moment_func(mom.rolling_min, np.min)

        a = np.array([1, 2, 3, 4, 5])
        b = mom.rolling_min(a, window=100, min_periods=1)
        assert_almost_equal(b, np.ones(len(a)))

        self.assertRaises(ValueError, mom.rolling_min, np.array([1,
                          2, 3]), window=3, min_periods=5)

    def test_rolling_max(self):
        self._check_moment_func(mom.rolling_max, np.max)

        a = np.array([1, 2, 3, 4, 5])
        b = mom.rolling_max(a, window=100, min_periods=1)
        assert_almost_equal(a, b)

        self.assertRaises(ValueError, mom.rolling_max, np.array([1,
                          2, 3]), window=3, min_periods=5)

    def test_rolling_quantile(self):
        qs = [.1, .5, .9]

        def scoreatpercentile(a, per):
            values = np.sort(a, axis=0)

            idx = per / 1. * (values.shape[0] - 1)
            return values[int(idx)]

        for q in qs:
            def f(x, window, min_periods=None, freq=None, center=False):
                return mom.rolling_quantile(x, window, q,
                                            min_periods=min_periods,
                                            freq=freq,
                                            center=center)

            def alt(x):
                return scoreatpercentile(x, q)

            self._check_moment_func(f, alt)

    def test_rolling_apply(self):
        ser = Series([])
        assert_series_equal(
            ser, mom.rolling_apply(ser, 10, lambda x: x.mean()))

        def roll_mean(x, window, min_periods=None, freq=None, center=False):
            return mom.rolling_apply(x, window,
                                     lambda x: x[np.isfinite(x)].mean(),
                                     min_periods=min_periods,
                                     freq=freq,
                                     center=center)
        self._check_moment_func(roll_mean, np.mean)

    def test_rolling_apply_out_of_bounds(self):
        # #1850
        arr = np.arange(4)

        # it works!
        result = mom.rolling_apply(arr, 10, np.sum)
        self.assert_(isnull(result).all())

        result = mom.rolling_apply(arr, 10, np.sum, min_periods=1)
        assert_almost_equal(result, result)

    def test_rolling_std(self):
        self._check_moment_func(mom.rolling_std,
                                lambda x: np.std(x, ddof=1))
        self._check_moment_func(functools.partial(mom.rolling_std, ddof=0),
                                lambda x: np.std(x, ddof=0))

    def test_rolling_std_1obs(self):
        result = mom.rolling_std(np.array([1., 2., 3., 4., 5.]),
                                 1, min_periods=1)
        expected = np.zeros(5)

        assert_almost_equal(result, expected)

        result = mom.rolling_std(np.array([np.nan, np.nan, 3., 4., 5.]),
                                 3, min_periods=2)
        self.assert_(np.isnan(result[2]))

    def test_rolling_std_neg_sqrt(self):
        # unit test from Bottleneck

        # Test move_nanstd for neg sqrt.

        a = np.array([0.0011448196318903589,
                      0.00028718669878572767,
                      0.00028718669878572767,
                      0.00028718669878572767,
                      0.00028718669878572767])
        b = mom.rolling_std(a, window=3)
        self.assert_(np.isfinite(b[2:]).all())

        b = mom.ewmstd(a, span=3)
        self.assert_(np.isfinite(b[2:]).all())

    def test_rolling_var(self):
        self._check_moment_func(mom.rolling_var,
                                lambda x: np.var(x, ddof=1))
        self._check_moment_func(functools.partial(mom.rolling_var, ddof=0),
                                lambda x: np.var(x, ddof=0))

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

    def test_fperr_robustness(self):
        # TODO: remove this once python 2.5 out of picture
        if PY3:
            raise nose.SkipTest

        # #2114
        data = '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x1a@\xaa\xaa\xaa\xaa\xaa\xaa\x02@8\x8e\xe38\x8e\xe3\xe8?z\t\xed%\xb4\x97\xd0?\xa2\x0c<\xdd\x9a\x1f\xb6?\x82\xbb\xfa&y\x7f\x9d?\xac\'\xa7\xc4P\xaa\x83?\x90\xdf\xde\xb0k8j?`\xea\xe9u\xf2zQ?*\xe37\x9d\x98N7?\xe2.\xf5&v\x13\x1f?\xec\xc9\xf8\x19\xa4\xb7\x04?\x90b\xf6w\x85\x9f\xeb>\xb5A\xa4\xfaXj\xd2>F\x02\xdb\xf8\xcb\x8d\xb8>.\xac<\xfb\x87^\xa0>\xe8:\xa6\xf9_\xd3\x85>\xfb?\xe2cUU\xfd?\xfc\x7fA\xed8\x8e\xe3?\xa5\xaa\xac\x91\xf6\x12\xca?n\x1cs\xb6\xf9a\xb1?\xe8%D\xf3L-\x97?5\xddZD\x11\xe7~?#>\xe7\x82\x0b\x9ad?\xd9R4Y\x0fxK?;7x;\nP2?N\xf4JO\xb8j\x18?4\xf81\x8a%G\x00?\x9a\xf5\x97\r2\xb4\xe5>\xcd\x9c\xca\xbcB\xf0\xcc>3\x13\x87(\xd7J\xb3>\x99\x19\xb4\xe0\x1e\xb9\x99>ff\xcd\x95\x14&\x81>\x88\x88\xbc\xc7p\xddf>`\x0b\xa6_\x96|N>@\xb2n\xea\x0eS4>U\x98\x938i\x19\x1b>\x8eeb\xd0\xf0\x10\x02>\xbd\xdc-k\x96\x16\xe8=(\x93\x1e\xf2\x0e\x0f\xd0=\xe0n\xd3Bii\xb5=*\xe9\x19Y\x8c\x8c\x9c=\xc6\xf0\xbb\x90]\x08\x83=]\x96\xfa\xc0|`i=>d\xfc\xd5\xfd\xeaP=R0\xfb\xc7\xa7\x8e6=\xc2\x95\xf9_\x8a\x13\x1e=\xd6c\xa6\xea\x06\r\x04=r\xda\xdd8\t\xbc\xea<\xf6\xe6\x93\xd0\xb0\xd2\xd1<\x9d\xdeok\x96\xc3\xb7<&~\xea9s\xaf\x9f<UUUUUU\x13@q\x1c\xc7q\x1c\xc7\xf9?\xf6\x12\xdaKh/\xe1?\xf2\xc3"e\xe0\xe9\xc6?\xed\xaf\x831+\x8d\xae?\xf3\x1f\xad\xcb\x1c^\x94?\x15\x1e\xdd\xbd>\xb8\x02@\xc6\xd2&\xfd\xa8\xf5\xe8?\xd9\xe1\x19\xfe\xc5\xa3\xd0?v\x82"\xa8\xb2/\xb6?\x9dX\x835\xee\x94\x9d?h\x90W\xce\x9e\xb8\x83?\x8a\xc0th~Kj?\\\x80\xf8\x9a\xa9\x87Q?%\xab\xa0\xce\x8c_7?1\xe4\x80\x13\x11*\x1f? \x98\x00\r\xb6\xc6\x04?\x80u\xabf\x9d\xb3\xeb>UNrD\xbew\xd2>\x1c\x13C[\xa8\x9f\xb8>\x12b\xd7<pj\xa0>m-\x1fQ@\xe3\x85>\xe6\x91)l\x00/m>Da\xc6\xf2\xaatS>\x05\xd7]\xee\xe3\xf09>'

        arr = np.frombuffer(data, dtype='<f8')
        if sys.byteorder != "little":
            arr = arr.byteswap().newbyteorder()

        result = mom.rolling_sum(arr, 2)
        self.assertTrue((result[1:] >= 0).all())

        result = mom.rolling_mean(arr, 2)
        self.assertTrue((result[1:] >= 0).all())

        result = mom.rolling_var(arr, 2)
        self.assertTrue((result[1:] >= 0).all())

        # #2527, ugh
        arr = np.array([0.00012456, 0.0003, 0])
        result = mom.rolling_mean(arr, 1)
        self.assertTrue(result[-1] >= 0)

        result = mom.rolling_mean(-arr, 1)
        self.assertTrue(result[-1] <= 0)

    def _check_moment_func(self, func, static_comp, window=50,
                           has_min_periods=True,
                           has_center=True,
                           has_time_rule=True,
                           preserve_nan=True,
                           fill_value=None):

        self._check_ndarray(func, static_comp, window=window,
                            has_min_periods=has_min_periods,
                            preserve_nan=preserve_nan,
                            has_center=has_center,
                            fill_value=fill_value)

        self._check_structures(func, static_comp,
                               has_min_periods=has_min_periods,
                               has_time_rule=has_time_rule,
                               fill_value=fill_value,
                               has_center=has_center)

    def _check_ndarray(self, func, static_comp, window=50,
                       has_min_periods=True,
                       preserve_nan=True,
                       has_center=True,
                       fill_value=None):

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

            arr2 = randn(20)
            result = func(arr2, 10, min_periods=5)
            self.assert_(isnull(result[3]))
            self.assert_(notnull(result[4]))

            # min_periods=0
            result0 = func(arr, 20, min_periods=0)
            result1 = func(arr, 20, min_periods=1)
            assert_almost_equal(result0, result1)
        else:
            result = func(arr, 50)
            assert_almost_equal(result[-1], static_comp(arr[10:-10]))

        if has_center:
            if has_min_periods:
                result = func(arr, 20, min_periods=15, center=True)
                expected = func(arr, 20, min_periods=15)
            else:
                result = func(arr, 20, center=True)
                expected = func(arr, 20)

            assert_almost_equal(result[1], expected[10])
            if fill_value is None:
                self.assert_(np.isnan(result[-9:]).all())
            else:
                self.assert_((result[-9:] == 0).all())
            if has_min_periods:
                self.assert_(np.isnan(expected[23]))
                self.assert_(np.isnan(result[14]))
                self.assert_(np.isnan(expected[-5]))
                self.assert_(np.isnan(result[-14]))

    def _check_structures(self, func, static_comp,
                          has_min_periods=True, has_time_rule=True,
                          has_center=True,
                          fill_value=None):

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
                                     freq='B')
                frame_result = func(self.frame[::2], win, min_periods=minp,
                                    freq='B')
            else:
                series_result = func(self.series[::2], win, freq='B')
                frame_result = func(self.frame[::2], win, freq='B')

            last_date = series_result.index[-1]
            prev_date = last_date - 24 * datetools.bday

            trunc_series = self.series[::2].truncate(prev_date, last_date)
            trunc_frame = self.frame[::2].truncate(prev_date, last_date)

            assert_almost_equal(series_result[-1], static_comp(trunc_series))

            assert_almost_equal(frame_result.xs(last_date),
                                trunc_frame.apply(static_comp))

        if has_center:
            if has_min_periods:
                minp = 10
                series_xp = func(self.series, 25, min_periods=minp).shift(-12)
                frame_xp = func(self.frame, 25, min_periods=minp).shift(-12)

                series_rs = func(self.series, 25, min_periods=minp,
                                 center=True)
                frame_rs = func(self.frame, 25, min_periods=minp,
                                center=True)

            else:
                series_xp = func(self.series, 25).shift(-12)
                frame_xp = func(self.frame, 25).shift(-12)

                series_rs = func(self.series, 25, center=True)
                frame_rs = func(self.frame, 25, center=True)

            if fill_value is not None:
                series_xp = series_xp.fillna(fill_value)
                frame_xp = frame_xp.fillna(fill_value)
            assert_series_equal(series_xp, series_rs)
            assert_frame_equal(frame_xp, frame_rs)

    def test_legacy_time_rule_arg(self):
        from StringIO import StringIO
        # suppress deprecation warnings
        sys.stderr = StringIO()

        rng = bdate_range('1/1/2000', periods=20)
        ts = Series(np.random.randn(20), index=rng)
        ts = ts.take(np.random.permutation(len(ts))[:12]).sort_index()

        try:
            result = mom.rolling_mean(ts, 1, min_periods=1, freq='B')
            expected = mom.rolling_mean(ts, 1, min_periods=1,
                                        time_rule='WEEKDAY')
            tm.assert_series_equal(result, expected)

            result = mom.ewma(ts, span=5, freq='B')
            expected = mom.ewma(ts, span=5, time_rule='WEEKDAY')
            tm.assert_series_equal(result, expected)

        finally:
            sys.stderr = sys.__stderr__

    def test_ewma(self):
        self._check_ew(mom.ewma)

        arr = np.zeros(1000)
        arr[5] = 1
        result = mom.ewma(arr, span=100, adjust=False).sum()
        self.assert_(np.abs(result - 1) < 1e-2)

    def test_ewma_nan_handling(self):
        s = Series([1.] + [np.nan] * 5 + [1.])

        result = mom.ewma(s, com=5)
        assert_almost_equal(result, [1] * len(s))

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

    def test_ew_empty_arrays(self):
        arr = np.array([], dtype=np.float64)

        funcs = [mom.ewma, mom.ewmvol, mom.ewmvar]
        for f in funcs:
            result = f(arr, 3)
            assert_almost_equal(result, arr)

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
        self.assert_(result2.dtype == np.float_)

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

    def test_flex_binary_moment(self):
        # GH3155
        # don't blow the stack
        self.assertRaises(ValueError, mom._flex_binary_moment,5,6,None)

    def test_corr_sanity(self):
        #GH 3155
        df = DataFrame(
            np.array(
                    [[ 0.87024726,  0.18505595],
                      [ 0.64355431,  0.3091617 ],
                      [ 0.92372966,  0.50552513],
                      [ 0.00203756,  0.04520709],
                      [ 0.84780328,  0.33394331],
                      [ 0.78369152,  0.63919667]])
            )

        res = mom.rolling_corr(df[0],df[1],5,center=True)
        self.assertTrue(all([np.abs(np.nan_to_num(x)) <=1 for x in res]))

        # and some fuzzing
        for i in range(10):
            df = DataFrame(np.random.rand(30,2))
            res = mom.rolling_corr(df[0],df[1],5,center=True)
            print( res)
            self.assertTrue(all([np.abs(np.nan_to_num(x)) <=1 for x in res]))

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

    def test_expanding_apply(self):
        ser = Series([])
        assert_series_equal(ser, mom.expanding_apply(ser, lambda x: x.mean()))

        def expanding_mean(x, min_periods=1, freq=None):
            return mom.expanding_apply(x,
                                       lambda x: x.mean(),
                                       min_periods=min_periods,
                                       freq=freq)
        self._check_expanding(expanding_mean, np.mean)

    def test_expanding_corr(self):
        A = self.series.dropna()
        B = (A + randn(len(A)))[:-5]

        result = mom.expanding_corr(A, B)

        rolling_result = mom.rolling_corr(A, B, len(A), min_periods=1)

        assert_almost_equal(rolling_result, result)

    def test_expanding_count(self):
        result = mom.expanding_count(self.series)
        assert_almost_equal(result, mom.rolling_count(self.series,
                                                      len(self.series)))

    def test_expanding_quantile(self):
        result = mom.expanding_quantile(self.series, 0.5)

        rolling_result = mom.rolling_quantile(self.series,
                                              len(self.series),
                                              0.5, min_periods=1)

        assert_almost_equal(result, rolling_result)

    def test_expanding_cov(self):
        A = self.series
        B = (A + randn(len(A)))[:-5]

        result = mom.expanding_cov(A, B)

        rolling_result = mom.rolling_cov(A, B, len(A), min_periods=1)

        assert_almost_equal(rolling_result, result)

    def test_expanding_max(self):
        self._check_expanding(mom.expanding_max, np.max, preserve_nan=False)

    def test_expanding_corr_pairwise(self):
        result = mom.expanding_corr_pairwise(self.frame)

        rolling_result = mom.rolling_corr_pairwise(self.frame,
                                                   len(self.frame),
                                                   min_periods=1)

        for i in result.items:
            assert_almost_equal(result[i], rolling_result[i])

    def _check_expanding_ndarray(self, func, static_comp, has_min_periods=True,
                                 has_time_rule=True, preserve_nan=True):
        result = func(self.arr)

        assert_almost_equal(result[10],
                            static_comp(self.arr[:11]))

        if preserve_nan:
            assert(np.isnan(result[self._nan_locs]).all())

        arr = randn(50)

        if has_min_periods:
            result = func(arr, min_periods=30)
            assert(np.isnan(result[:29]).all())
            assert_almost_equal(result[-1], static_comp(arr[:50]))

            # min_periods is working correctly
            result = func(arr, min_periods=15)
            self.assert_(np.isnan(result[13]))
            self.assert_(not np.isnan(result[14]))

            arr2 = randn(20)
            result = func(arr2, min_periods=5)
            self.assert_(isnull(result[3]))
            self.assert_(notnull(result[4]))

            # min_periods=0
            result0 = func(arr, min_periods=0)
            result1 = func(arr, min_periods=1)
            assert_almost_equal(result0, result1)
        else:
            result = func(arr)
            assert_almost_equal(result[-1], static_comp(arr[:50]))

    def _check_expanding_structures(self, func):
        series_result = func(self.series)
        self.assert_(isinstance(series_result, Series))
        frame_result = func(self.frame)
        self.assertEquals(type(frame_result), DataFrame)

    def _check_expanding(self, func, static_comp, has_min_periods=True,
                         has_time_rule=True,
                         preserve_nan=True):
        self._check_expanding_ndarray(func, static_comp,
                                      has_min_periods=has_min_periods,
                                      has_time_rule=has_time_rule,
                                      preserve_nan=preserve_nan)
        self._check_expanding_structures(func)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
