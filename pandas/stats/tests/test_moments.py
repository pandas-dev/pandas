import nose
import sys
import functools
import warnings

from datetime import datetime
from numpy.random import randn
from numpy.testing.decorators import slow
import numpy as np
from distutils.version import LooseVersion

from pandas import Series, DataFrame, Panel, bdate_range, isnull, notnull, concat
from pandas.util.testing import (
    assert_almost_equal, assert_series_equal, assert_frame_equal, assert_panel_equal, assert_index_equal
)
import pandas.core.datetools as datetools
import pandas.stats.moments as mom
import pandas.util.testing as tm
from pandas.compat import range, zip, PY3, StringIO

N, K = 100, 10

class Base(tm.TestCase):

    _multiprocess_can_split_ = True

    _nan_locs = np.arange(20, 40)
    _inf_locs = np.array([])

    def _create_data(self):
        arr = randn(N)
        arr[self._nan_locs] = np.NaN

        self.arr = arr
        self.rng = bdate_range(datetime(2009, 1, 1), periods=N)

        self.series = Series(arr.copy(), index=self.rng)

        self.frame = DataFrame(randn(N, K), index=self.rng,
                               columns=np.arange(K))

class TestMoments(Base):

    def setUp(self):
        self._create_data()
        warnings.simplefilter("ignore", category=FutureWarning)

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
        # GH 8238
        tm._skip_if_no_scipy()

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49,
                         16.68, 9.48, 10.63, 14.48])
        xp = np.array([np.nan, np.nan, 9.962, 11.27 , 11.564, 12.516,
                       12.818,  12.952, np.nan, np.nan])

        rs = mom.rolling_mean(vals, 5, center=True)
        assert_almost_equal(xp, rs)

        xp = Series(rs)
        rs = mom.rolling_mean(Series(vals), 5, center=True)
        assert_series_equal(xp, rs)

    def test_cmov_window(self):
        # GH 8238
        tm._skip_if_no_scipy()

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81,
                         13.49, 16.68, 9.48, 10.63, 14.48])
        xp = np.array([np.nan, np.nan, 9.962, 11.27 , 11.564, 12.516,
                       12.818,  12.952, np.nan, np.nan])

        rs = mom.rolling_window(vals, 5, 'boxcar', center=True)
        assert_almost_equal(xp, rs)

        xp = Series(rs)
        rs = mom.rolling_window(Series(vals), 5, 'boxcar', center=True)
        assert_series_equal(xp, rs)

    def test_cmov_window_corner(self):
        # GH 8238
        tm._skip_if_no_scipy()

        # all nan
        vals = np.empty(10, dtype=float)
        vals.fill(np.nan)
        rs = mom.rolling_window(vals, 5, 'boxcar', center=True)
        self.assertTrue(np.isnan(rs).all())

        # empty
        vals = np.array([])
        rs = mom.rolling_window(vals, 5, 'boxcar', center=True)
        self.assertEqual(len(rs), 0)

        # shorter than window
        vals = np.random.randn(5)
        rs = mom.rolling_window(vals, 10, 'boxcar')
        self.assertTrue(np.isnan(rs).all())
        self.assertEqual(len(rs), 5)

    def test_cmov_window_frame(self):
        # Gh 8238
        tm._skip_if_no_scipy()

        vals = np.array([[ 12.18,   3.64],
                         [ 10.18,   9.16],
                         [ 13.24,  14.61],
                         [  4.51,   8.11],
                         [  6.15,  11.44],
                         [  9.14,   6.21],
                         [ 11.31,  10.67],
                         [  2.94,   6.51],
                         [  9.42,   8.39],
                         [ 12.44,   7.34 ]])

        xp = np.array([[ np.nan,  np.nan],
                       [ np.nan,  np.nan],
                       [  9.252,   9.392],
                       [  8.644,   9.906],
                       [  8.87 ,  10.208],
                       [  6.81 ,   8.588],
                       [  7.792,   8.644],
                       [  9.05 ,   7.824],
                       [ np.nan,  np.nan],
                       [ np.nan,  np.nan]])

        # DataFrame
        rs = mom.rolling_window(DataFrame(vals), 5, 'boxcar', center=True)
        assert_frame_equal(DataFrame(xp), rs)

    def test_cmov_window_na_min_periods(self):
        tm._skip_if_no_scipy()

        # min_periods
        vals = Series(np.random.randn(10))
        vals[4] = np.nan
        vals[8] = np.nan

        xp = mom.rolling_mean(vals, 5, min_periods=4, center=True)
        rs = mom.rolling_window(vals, 5, 'boxcar', min_periods=4, center=True)

        assert_series_equal(xp, rs)

    def test_cmov_window_regular(self):
        # GH 8238
        tm._skip_if_no_scipy()

        win_types = ['triang', 'blackman', 'hamming', 'bartlett', 'bohman',
                     'blackmanharris', 'nuttall', 'barthann']

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81,
                         13.49, 16.68, 9.48, 10.63, 14.48])
        xps = {
            'hamming': [np.nan, np.nan, 8.71384, 9.56348, 12.38009,
                        14.03687, 13.8567, 11.81473, np.nan, np.nan],
            'triang': [np.nan, np.nan, 9.28667, 10.34667, 12.00556,
                       13.33889, 13.38, 12.33667, np.nan, np.nan],
            'barthann': [np.nan, np.nan, 8.4425, 9.1925, 12.5575,
                         14.3675, 14.0825, 11.5675, np.nan, np.nan],
            'bohman': [np.nan, np.nan, 7.61599, 9.1764, 12.83559,
                       14.17267, 14.65923, 11.10401, np.nan, np.nan],
            'blackmanharris': [np.nan, np.nan, 6.97691, 9.16438, 13.05052,
                               14.02156, 15.10512, 10.74574, np.nan, np.nan],
            'nuttall': [np.nan, np.nan, 7.04618, 9.16786, 13.02671,
                        14.03559, 15.05657, 10.78514, np.nan, np.nan],
            'blackman': [np.nan, np.nan, 7.73345, 9.17869, 12.79607,
                         14.20036, 14.57726, 11.16988, np.nan, np.nan],
            'bartlett': [np.nan, np.nan, 8.4425, 9.1925, 12.5575,
                         14.3675, 14.0825, 11.5675, np.nan, np.nan]}

        for wt in win_types:
            xp = Series(xps[wt])
            rs = mom.rolling_window(Series(vals), 5, wt, center=True)
            assert_series_equal(xp, rs)

    def test_cmov_window_regular_linear_range(self):
        # GH 8238
        tm._skip_if_no_scipy()

        win_types = ['triang', 'blackman', 'hamming', 'bartlett', 'bohman',
                     'blackmanharris', 'nuttall', 'barthann']

        vals = np.array(range(10), dtype=np.float)
        xp = vals.copy()
        xp[:2] = np.nan
        xp[-2:] = np.nan
        xp = Series(xp)

        for wt in win_types:
            rs = mom.rolling_window(Series(vals), 5, wt, center=True)
            assert_series_equal(xp, rs)

    def test_cmov_window_regular_missing_data(self):
        # GH 8238
        tm._skip_if_no_scipy()

        win_types = ['triang', 'blackman', 'hamming', 'bartlett', 'bohman',
                     'blackmanharris', 'nuttall', 'barthann']

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81,
                         13.49, 16.68, np.nan, 10.63, 14.48])
        xps = {
            'bartlett': [np.nan, np.nan, 9.70333, 10.5225, 8.4425,
                         9.1925, 12.5575, 14.3675, 15.61667, 13.655],
            'blackman': [np.nan, np.nan, 9.04582, 11.41536, 7.73345,
                         9.17869, 12.79607, 14.20036, 15.8706, 13.655],
            'barthann': [np.nan, np.nan, 9.70333, 10.5225, 8.4425,
                         9.1925, 12.5575, 14.3675, 15.61667, 13.655],
            'bohman': [np.nan, np.nan, 8.9444, 11.56327, 7.61599,
                       9.1764, 12.83559, 14.17267, 15.90976, 13.655],
            'hamming': [np.nan, np.nan, 9.59321, 10.29694, 8.71384,
                        9.56348, 12.38009, 14.20565, 15.24694, 13.69758],
            'nuttall': [np.nan, np.nan, 8.47693, 12.2821, 7.04618,
                        9.16786, 13.02671, 14.03673, 16.08759, 13.65553],
            'triang': [np.nan, np.nan, 9.33167, 9.76125, 9.28667,
                       10.34667, 12.00556, 13.82125, 14.49429, 13.765],
            'blackmanharris': [np.nan, np.nan, 8.42526, 12.36824, 6.97691,
                               9.16438, 13.05052, 14.02175, 16.1098,
                               13.65509]
            }

        for wt in win_types:
            xp = Series(xps[wt])
            rs = mom.rolling_window(Series(vals), 5, wt, min_periods=3)
            assert_series_equal(xp, rs)

    def test_cmov_window_special(self):
        # GH 8238
        tm._skip_if_no_scipy()

        win_types = ['kaiser', 'gaussian', 'general_gaussian', 'slepian']
        kwds = [{'beta': 1.}, {'std': 1.}, {'power': 2., 'width': 2.},
                {'width': 0.5}]

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81,
                         13.49, 16.68, 9.48, 10.63, 14.48])

        xps = {
            'gaussian': [np.nan, np.nan, 8.97297, 9.76077, 12.24763,
                         13.89053, 13.65671, 12.01002, np.nan, np.nan],
            'general_gaussian': [np.nan, np.nan, 9.85011, 10.71589,
                                 11.73161, 13.08516, 12.95111, 12.74577,
                                 np.nan, np.nan],
            'slepian': [np.nan, np.nan, 9.81073, 10.89359, 11.70284,
                        12.88331, 12.96079, 12.77008, np.nan, np.nan],
            'kaiser': [np.nan, np.nan, 9.86851, 11.02969, 11.65161,
                       12.75129, 12.90702, 12.83757, np.nan, np.nan]
        }

        for wt, k in zip(win_types, kwds):
            xp = Series(xps[wt])

            rs = mom.rolling_window(Series(vals), 5, wt, center=True,
                                    **k)
            assert_series_equal(xp, rs)

    def test_cmov_window_special_linear_range(self):
        # GH 8238
        tm._skip_if_no_scipy()

        win_types = ['kaiser', 'gaussian', 'general_gaussian', 'slepian']
        kwds = [{'beta': 1.}, {'std': 1.}, {'power': 2., 'width': 2.},
                {'width': 0.5}]

        vals = np.array(range(10), dtype=np.float)
        xp = vals.copy()
        xp[:2] = np.nan
        xp[-2:] = np.nan
        xp = Series(xp)

        for wt, k in zip(win_types, kwds):
            rs = mom.rolling_window(Series(vals), 5, wt, center=True,
                                    **k)
            assert_series_equal(xp, rs)

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
        # suppress warnings about empty slices, as we are deliberately testing with a 0-length Series
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*(empty slice|0 for slice).*", category=RuntimeWarning)

            ser = Series([])
            assert_series_equal(ser, mom.rolling_apply(ser, 10, lambda x: x.mean()))

            def roll_mean(x, window, min_periods=None, freq=None, center=False):
                return mom.rolling_apply(x, window,
                                         lambda x: x[np.isfinite(x)].mean(),
                                         min_periods=min_periods,
                                         freq=freq,
                                         center=center)
            self._check_moment_func(roll_mean, np.mean)

        # GH 8080
        s = Series([None, None, None])
        result = mom.rolling_apply(s, 2, lambda x: len(x), min_periods=0)
        expected = Series([1., 2., 2.])
        assert_series_equal(result, expected)

    def test_rolling_apply_out_of_bounds(self):
        # #1850
        arr = np.arange(4)

        # it works!
        result = mom.rolling_apply(arr, 10, np.sum)
        self.assertTrue(isnull(result).all())

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
        expected = np.array([np.nan] * 5)
        assert_almost_equal(result, expected)

        result = mom.rolling_std(np.array([1., 2., 3., 4., 5.]),
                                 1, min_periods=1, ddof=0)
        expected = np.zeros(5)
        assert_almost_equal(result, expected)

        result = mom.rolling_std(np.array([np.nan, np.nan, 3., 4., 5.]),
                                 3, min_periods=2)
        self.assertTrue(np.isnan(result[2]))

    def test_rolling_std_neg_sqrt(self):
        # unit test from Bottleneck

        # Test move_nanstd for neg sqrt.

        a = np.array([0.0011448196318903589,
                      0.00028718669878572767,
                      0.00028718669878572767,
                      0.00028718669878572767,
                      0.00028718669878572767])
        b = mom.rolling_std(a, window=3)
        self.assertTrue(np.isfinite(b[2:]).all())

        b = mom.ewmstd(a, span=3)
        self.assertTrue(np.isfinite(b[2:]).all())

    def test_rolling_var(self):
        self._check_moment_func(mom.rolling_var,
                                lambda x: np.var(x, ddof=1),
                                test_stable=True)
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
            raise nose.SkipTest("doesn't work on python 3")

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
                           fill_value=None,
                           test_stable=False):

        self._check_ndarray(func, static_comp, window=window,
                            has_min_periods=has_min_periods,
                            preserve_nan=preserve_nan,
                            has_center=has_center,
                            fill_value=fill_value,
                            test_stable=test_stable)

        self._check_structures(func, static_comp,
                               has_min_periods=has_min_periods,
                               has_time_rule=has_time_rule,
                               fill_value=fill_value,
                               has_center=has_center)

    def _check_ndarray(self, func, static_comp, window=50,
                       has_min_periods=True,
                       preserve_nan=True,
                       has_center=True,
                       fill_value=None,
                       test_stable=False,
                       test_window=True):

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
            self.assertTrue(np.isnan(result[23]))
            self.assertFalse(np.isnan(result[24]))

            self.assertFalse(np.isnan(result[-6]))
            self.assertTrue(np.isnan(result[-5]))

            arr2 = randn(20)
            result = func(arr2, 10, min_periods=5)
            self.assertTrue(isnull(result[3]))
            self.assertTrue(notnull(result[4]))

            # min_periods=0
            result0 = func(arr, 20, min_periods=0)
            result1 = func(arr, 20, min_periods=1)
            assert_almost_equal(result0, result1)
        else:
            result = func(arr, 50)
            assert_almost_equal(result[-1], static_comp(arr[10:-10]))

        # GH 7925
        if has_center:
            if has_min_periods:
                result = func(arr, 20, min_periods=15, center=True)
                expected = func(np.concatenate((arr, np.array([np.NaN] * 9))), 20, min_periods=15)[9:]
            else:
                result = func(arr, 20, center=True)
                expected = func(np.concatenate((arr, np.array([np.NaN] * 9))), 20)[9:]

            self.assert_numpy_array_equal(result, expected)

        if test_stable:
            result = func(self.arr + 1e9, window)
            assert_almost_equal(result[-1],
                                static_comp(self.arr[-50:] + 1e9))

        # Test window larger than array, #7297
        if test_window:
            if has_min_periods:
                for minp in (0, len(self.arr)-1, len(self.arr)):
                    result = func(self.arr, len(self.arr)+1, min_periods=minp)
                    expected = func(self.arr, len(self.arr), min_periods=minp)
                    nan_mask = np.isnan(result)
                    self.assertTrue(np.array_equal(nan_mask,
                                                   np.isnan(expected)))
                    nan_mask = ~nan_mask
                    assert_almost_equal(result[nan_mask], expected[nan_mask])
            else:
                result = func(self.arr, len(self.arr)+1)
                expected = func(self.arr, len(self.arr))
                nan_mask = np.isnan(result)
                self.assertTrue(np.array_equal(nan_mask, np.isnan(expected)))
                nan_mask = ~nan_mask
                assert_almost_equal(result[nan_mask], expected[nan_mask])




    def _check_structures(self, func, static_comp,
                          has_min_periods=True, has_time_rule=True,
                          has_center=True,
                          fill_value=None):

        series_result = func(self.series, 50)
        tm.assertIsInstance(series_result, Series)

        frame_result = func(self.frame, 50)
        self.assertEqual(type(frame_result), DataFrame)

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

        # GH 7925
        if has_center:
            if has_min_periods:
                minp = 10
                series_xp = func(self.series.reindex(list(self.series.index)+['x%d'%x for x in range(12)]), 25, min_periods=minp).shift(-12).reindex(self.series.index)
                frame_xp = func(self.frame.reindex(list(self.frame.index)+['x%d'%x for x in range(12)]), 25, min_periods=minp).shift(-12).reindex(self.frame.index)

                series_rs = func(self.series, 25, min_periods=minp,
                                 center=True)
                frame_rs = func(self.frame, 25, min_periods=minp,
                                center=True)

            else:
                series_xp = func(self.series.reindex(list(self.series.index)+['x%d'%x for x in range(12)]), 25).shift(-12).reindex(self.series.index)
                frame_xp = func(self.frame.reindex(list(self.frame.index)+['x%d'%x for x in range(12)]), 25).shift(-12).reindex(self.frame.index)

                series_rs = func(self.series, 25, center=True)
                frame_rs = func(self.frame, 25, center=True)

            if fill_value is not None:
                series_xp = series_xp.fillna(fill_value)
                frame_xp = frame_xp.fillna(fill_value)
            assert_series_equal(series_xp, series_rs)
            assert_frame_equal(frame_xp, frame_rs)

    def test_ewma(self):
        self._check_ew(mom.ewma)

        arr = np.zeros(1000)
        arr[5] = 1
        result = mom.ewma(arr, span=100, adjust=False).sum()
        self.assertTrue(np.abs(result - 1) < 1e-2)

        s = Series([1.0, 2.0, 4.0, 8.0])

        expected = Series([1.0, 1.6, 2.736842, 4.923077])
        for f in [lambda s: mom.ewma(s, com=2.0, adjust=True),
                  lambda s: mom.ewma(s, com=2.0, adjust=True, ignore_na=False),
                  lambda s: mom.ewma(s, com=2.0, adjust=True, ignore_na=True),
                 ]:
            result = f(s)
            assert_series_equal(result, expected)

        expected = Series([1.0, 1.333333, 2.222222, 4.148148])
        for f in [lambda s: mom.ewma(s, com=2.0, adjust=False),
                  lambda s: mom.ewma(s, com=2.0, adjust=False, ignore_na=False),
                  lambda s: mom.ewma(s, com=2.0, adjust=False, ignore_na=True),
                 ]:
            result = f(s)
            assert_series_equal(result, expected)

    def test_ewma_nan_handling(self):
        s = Series([1.] + [np.nan] * 5 + [1.])
        result = mom.ewma(s, com=5)
        assert_almost_equal(result, [1.] * len(s))

        s = Series([np.nan] * 2 + [1.] + [np.nan] * 2 + [1.])
        result = mom.ewma(s, com=5)
        assert_almost_equal(result, [np.nan] * 2 + [1.] * 4)

        # GH 7603
        s0 = Series([np.nan, 1., 101.])
        s1 = Series([1., np.nan, 101.])
        s2 = Series([np.nan, 1., np.nan, np.nan, 101., np.nan])
        s3 = Series([1., np.nan, 101., 50.])
        com = 2.
        alpha = 1. / (1. + com)

        def simple_wma(s, w):
            return (s.multiply(w).cumsum() / w.cumsum()).fillna(method='ffill')

        for (s, adjust, ignore_na, w) in [
                (s0, True, False, [np.nan, (1. - alpha), 1.]),
                (s0, True, True, [np.nan, (1. - alpha), 1.]),
                (s0, False, False, [np.nan, (1. - alpha), alpha]),
                (s0, False, True, [np.nan, (1. - alpha), alpha]),
                (s1, True, False, [(1. - alpha)**2, np.nan, 1.]),
                (s1, True, True, [(1. - alpha), np.nan, 1.]),
                (s1, False, False, [(1. - alpha)**2, np.nan, alpha]),
                (s1, False, True, [(1. - alpha), np.nan, alpha]),
                (s2, True, False, [np.nan, (1. - alpha)**3, np.nan, np.nan, 1., np.nan]),
                (s2, True, True, [np.nan, (1. - alpha), np.nan, np.nan, 1., np.nan]),
                (s2, False, False, [np.nan, (1. - alpha)**3, np.nan, np.nan, alpha, np.nan]),
                (s2, False, True, [np.nan, (1. - alpha), np.nan, np.nan, alpha, np.nan]),
                (s3, True, False, [(1. - alpha)**3, np.nan, (1. - alpha), 1.]),
                (s3, True, True, [(1. - alpha)**2, np.nan, (1. - alpha), 1.]),
                (s3, False, False, [(1. - alpha)**3, np.nan, (1. - alpha) * alpha, alpha * ((1. - alpha)**2 + alpha)]),
                (s3, False, True, [(1. - alpha)**2, np.nan, (1. - alpha) * alpha, alpha]),
                ]:
            expected = simple_wma(s, Series(w))
            result = mom.ewma(s, com=com, adjust=adjust, ignore_na=ignore_na)
            assert_series_equal(result, expected)
            if ignore_na is False:
                # check that ignore_na defaults to False
                result = mom.ewma(s, com=com, adjust=adjust)
                assert_series_equal(result, expected)

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

    def test_ewma_halflife_arg(self):
        A = mom.ewma(self.arr, com=13.932726172912965)
        B = mom.ewma(self.arr, halflife=10.0)
        assert_almost_equal(A, B)

        self.assertRaises(Exception, mom.ewma, self.arr, span=20, halflife=50)
        self.assertRaises(Exception, mom.ewma, self.arr, com=9.5, halflife=50)
        self.assertRaises(Exception, mom.ewma, self.arr, com=9.5, span=20, halflife=50)
        self.assertRaises(Exception, mom.ewma, self.arr)

    def test_moment_preserve_series_name(self):
        # GH 10565
        s = Series(np.arange(100), name='foo')
        s2 = mom.rolling_mean(s, 30)
        s3 = mom.rolling_sum(s, 20)
        self.assertEqual(s2.name, 'foo')
        self.assertEqual(s3.name, 'foo')

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
        s = Series(arr)

        # check min_periods
        # GH 7898
        result = func(s, 50, min_periods=2)
        self.assertTrue(np.isnan(result.values[:11]).all())
        self.assertFalse(np.isnan(result.values[11:]).any())

        for min_periods in (0, 1):
            result = func(s, 50, min_periods=min_periods)
            if func == mom.ewma:
                self.assertTrue(np.isnan(result.values[:10]).all())
                self.assertFalse(np.isnan(result.values[10:]).any())
            else:
                # ewmstd, ewmvol, ewmvar (with bias=False) require at least two values
                self.assertTrue(np.isnan(result.values[:11]).all())
                self.assertFalse(np.isnan(result.values[11:]).any())

            # check series of length 0
            result = func(Series([]), 50, min_periods=min_periods)
            assert_series_equal(result, Series([]))

            # check series of length 1
            result = func(Series([1.]), 50, min_periods=min_periods)
            if func == mom.ewma:
                assert_series_equal(result, Series([1.]))
            else:
                # ewmstd, ewmvol, ewmvar with bias=False require at least two values
                assert_series_equal(result, Series([np.NaN]))

        # pass in ints
        result2 = func(np.arange(50), span=10)
        self.assertEqual(result2.dtype, np.float_)

    def _check_ew_structures(self, func):
        series_result = func(self.series, com=10)
        tm.assertIsInstance(series_result, Series)
        frame_result = func(self.frame, com=10)
        self.assertEqual(type(frame_result), DataFrame)

# create the data only once as we are not setting it
def _create_consistency_data():

    def create_series():
       return [Series(),
               Series([np.nan]),
               Series([np.nan, np.nan]),
               Series([3.]),
               Series([np.nan, 3.]),
               Series([3., np.nan]),
               Series([1., 3.]),
               Series([2., 2.]),
               Series([3., 1.]),
               Series([5., 5., 5., 5., np.nan, np.nan, np.nan, 5., 5., np.nan, np.nan]),
               Series([np.nan, 5., 5., 5., np.nan, np.nan, np.nan, 5., 5., np.nan, np.nan]),
               Series([np.nan, np.nan, 5., 5., np.nan, np.nan, np.nan, 5., 5., np.nan, np.nan]),
               Series([np.nan, 3., np.nan, 3., 4., 5., 6., np.nan, np.nan, 7., 12., 13., 14., 15.]),
               Series([np.nan, 5., np.nan, 2., 4., 0., 9., np.nan, np.nan, 3., 12., 13., 14., 15.]),
               Series([2., 3., np.nan, 3., 4., 5., 6., np.nan, np.nan, 7., 12., 13., 14., 15.]),
               Series([2., 5., np.nan, 2., 4., 0., 9., np.nan, np.nan, 3., 12., 13., 14., 15.]),
               Series(range(10)),
               Series(range(20, 0, -2)),
              ]

    def create_dataframes():
       return [DataFrame(),
               DataFrame(columns=['a']),
               DataFrame(columns=['a', 'a']),
               DataFrame(columns=['a', 'b']),
               DataFrame(np.arange(10).reshape((5, 2))),
               DataFrame(np.arange(25).reshape((5, 5))),
               DataFrame(np.arange(25).reshape((5, 5)), columns=['a', 'b', 99, 'd', 'd']),
              ] + [DataFrame(s) for s in create_series()]

    def is_constant(x):
        values = x.values.ravel()
        return len(set(values[notnull(values)])) == 1

    def no_nans(x):
        return x.notnull().all().all()

    # data is a tuple(object, is_contant, no_nans)
    data = create_series() + create_dataframes()

    return [ (x, is_constant(x), no_nans(x)) for x in data ]
_consistency_data = _create_consistency_data()

class TestMomentsConsistency(Base):
    base_functions = [
        (lambda v: Series(v).count(), None, 'count'),
        (lambda v: Series(v).max(), None, 'max'),
        (lambda v: Series(v).min(), None, 'min'),
        (lambda v: Series(v).sum(), None, 'sum'),
        (lambda v: Series(v).mean(), None, 'mean'),
        (lambda v: Series(v).std(), 1, 'std'),
        (lambda v: Series(v).cov(Series(v)), None, 'cov'),
        (lambda v: Series(v).corr(Series(v)), None, 'corr'),
        (lambda v: Series(v).var(), 1, 'var'),
        #(lambda v: Series(v).skew(), 3, 'skew'), # restore once GH 8086 is fixed
        #(lambda v: Series(v).kurt(), 4, 'kurt'), # restore once GH 8086 is fixed
        #(lambda x, min_periods: mom.expanding_quantile(x, 0.3, min_periods=min_periods, 'quantile'),
        # lambda v: Series(v).quantile(0.3), None, 'quantile'), # restore once GH 8084 is fixed
        (lambda v: Series(v).median(), None ,'median'),
        (np.nanmax, 1, 'max'),
        (np.nanmin, 1, 'min'),
        (np.nansum, 1, 'sum'),
        ]
    if np.__version__ >= LooseVersion('1.8.0'):
        base_functions += [
            (np.nanmean, 1, 'mean'),
            (lambda v: np.nanstd(v, ddof=1), 1 ,'std'),
            (lambda v: np.nanvar(v, ddof=1), 1 ,'var'),
        ]
    if np.__version__ >= LooseVersion('1.9.0'):
        base_functions += [
            (np.nanmedian, 1, 'median'),
            ]
    no_nan_functions = [
        (np.max, None, 'max'),
        (np.min, None, 'min'),
        (np.sum, None, 'sum'),
        (np.mean, None, 'mean'),
        (lambda v: np.std(v, ddof=1), 1 ,'std'),
        (lambda v: np.var(v, ddof=1), 1 ,'var'),
        (np.median, None, 'median'),
    ]

    def _create_data(self):
        super(TestMomentsConsistency, self)._create_data()
        self.data = _consistency_data

    def setUp(self):
        self._create_data()
        warnings.simplefilter("ignore", category=FutureWarning)

    def _test_moments_consistency(self,
                                  min_periods,
                                  count, mean, mock_mean, corr,
                                  var_unbiased=None, std_unbiased=None, cov_unbiased=None,
                                  var_biased=None, std_biased=None, cov_biased=None,
                                  var_debiasing_factors=None):

        def _non_null_values(x):
            values = x.values.ravel()
            return set(values[notnull(values)].tolist())

        for (x, is_constant, no_nans) in self.data:
            assert_equal = assert_series_equal if isinstance(x, Series) else assert_frame_equal
            count_x = count(x)
            mean_x = mean(x)

            if mock_mean:
                # check that mean equals mock_mean
                expected = mock_mean(x)
                assert_equal(mean_x, expected.astype('float64'))

            # check that correlation of a series with itself is either 1 or NaN
            corr_x_x = corr(x, x)
            # self.assertTrue(_non_null_values(corr_x_x).issubset(set([1.]))) # restore once rolling_cov(x, x) is identically equal to var(x)

            if is_constant:
                exp = x.max() if isinstance(x, Series) else x.max().max()

                # check mean of constant series
                expected = x * np.nan
                expected[count_x >= max(min_periods, 1)] = exp
                assert_equal(mean_x, expected)

                # check correlation of constant series with itself is NaN
                expected[:] = np.nan
                assert_equal(corr_x_x, expected)

            if var_unbiased and var_biased and var_debiasing_factors:
                # check variance debiasing factors
                var_unbiased_x = var_unbiased(x)
                var_biased_x = var_biased(x)
                var_debiasing_factors_x = var_debiasing_factors(x)
                assert_equal(var_unbiased_x, var_biased_x * var_debiasing_factors_x)

            for (std, var, cov) in [(std_biased, var_biased, cov_biased),
                                    (std_unbiased, var_unbiased, cov_unbiased)]:

                # check that var(x), std(x), and cov(x) are all >= 0
                var_x = var(x)
                std_x = std(x)
                self.assertFalse((var_x < 0).any().any())
                self.assertFalse((std_x < 0).any().any())
                if cov:
                    cov_x_x = cov(x, x)
                    self.assertFalse((cov_x_x < 0).any().any())

                    # check that var(x) == cov(x, x)
                    assert_equal(var_x, cov_x_x)

                # check that var(x) == std(x)^2
                assert_equal(var_x, std_x * std_x)

                if var is var_biased:
                    # check that biased var(x) == mean(x^2) - mean(x)^2
                    mean_x2 = mean(x * x)
                    assert_equal(var_x, mean_x2 - (mean_x * mean_x))

                if is_constant:
                    # check that variance of constant series is identically 0
                    self.assertFalse((var_x > 0).any().any())
                    expected = x * np.nan
                    expected[count_x >= max(min_periods, 1)] = 0.
                    if var is var_unbiased:
                        expected[count_x < 2] = np.nan
                    assert_equal(var_x, expected)

                if isinstance(x, Series):
                    for (y, is_constant, no_nans) in self.data:
                        if not x.isnull().equals(y.isnull()):
                            # can only easily test two Series with similar structure
                            continue

                        # check that cor(x, y) is symmetric
                        corr_x_y = corr(x, y)
                        corr_y_x = corr(y, x)
                        assert_equal(corr_x_y, corr_y_x)

                        if cov:
                            # check that cov(x, y) is symmetric
                            cov_x_y = cov(x, y)
                            cov_y_x = cov(y, x)
                            assert_equal(cov_x_y, cov_y_x)

                            # check that cov(x, y) == (var(x+y) - var(x) - var(y)) / 2
                            var_x_plus_y = var(x + y)
                            var_y = var(y)
                            assert_equal(cov_x_y, 0.5 * (var_x_plus_y - var_x - var_y))

                            # check that corr(x, y) == cov(x, y) / (std(x) * std(y))
                            std_y = std(y)
                            assert_equal(corr_x_y, cov_x_y / (std_x * std_y))

                            if cov is cov_biased:
                                # check that biased cov(x, y) == mean(x*y) - mean(x)*mean(y)
                                mean_y = mean(y)
                                mean_x_times_y = mean(x * y)
                                assert_equal(cov_x_y, mean_x_times_y - (mean_x * mean_y))

    @slow
    def test_ewm_consistency(self):

        def _weights(s, com, adjust, ignore_na):
            if isinstance(s, DataFrame):
                if not len(s.columns):
                    return DataFrame(index=s.index, columns=s.columns)
                w = concat([ _weights(s.iloc[:, i],
                                      com=com,
                                      adjust=adjust,
                                      ignore_na=ignore_na) for i, _ in enumerate(s.columns) ],
                           axis=1)
                w.index=s.index
                w.columns=s.columns
                return w

            w = Series(np.nan, index=s.index)
            alpha = 1. / (1. + com)
            if ignore_na:
                w[s.notnull()] = _weights(s[s.notnull()], com=com, adjust=adjust, ignore_na=False)
            elif adjust:
                for i in range(len(s)):
                    if s.iat[i] == s.iat[i]:
                        w.iat[i] = pow(1. / (1. - alpha), i)
            else:
                sum_wts = 0.
                prev_i = -1
                for i in range(len(s)):
                    if s.iat[i] == s.iat[i]:
                        if prev_i == -1:
                            w.iat[i] = 1.
                        else:
                            w.iat[i] = alpha * sum_wts / pow(1. - alpha, i - prev_i)
                        sum_wts += w.iat[i]
                        prev_i = i
            return w

        def _variance_debiasing_factors(s, com, adjust, ignore_na):
            weights = _weights(s, com=com, adjust=adjust, ignore_na=ignore_na)
            cum_sum = weights.cumsum().fillna(method='ffill')
            cum_sum_sq = (weights * weights).cumsum().fillna(method='ffill')
            numerator = cum_sum * cum_sum
            denominator = numerator - cum_sum_sq
            denominator[denominator <= 0.] = np.nan
            return numerator / denominator

        def _ewma(s, com, min_periods, adjust, ignore_na):
            weights = _weights(s, com=com, adjust=adjust, ignore_na=ignore_na)
            result = s.multiply(weights).cumsum().divide(weights.cumsum()).fillna(method='ffill')
            result[mom.expanding_count(s) < (max(min_periods, 1) if min_periods else 1)] = np.nan
            return result

        com = 3.
        for min_periods in [0, 1, 2, 3, 4]:
            for adjust in [True, False]:
                for ignore_na in [False, True]:
                    # test consistency between different ewm* moments
                    self._test_moments_consistency(
                        min_periods=min_periods,
                        count=mom.expanding_count,
                        mean=lambda x: mom.ewma(x, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na),
                        mock_mean=lambda x: _ewma(x, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na),
                        corr=lambda x, y: mom.ewmcorr(x, y, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na),
                        var_unbiased=lambda x: mom.ewmvar(x, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na, bias=False),
                        std_unbiased=lambda x: mom.ewmstd(x, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na, bias=False),
                        cov_unbiased=lambda x, y: mom.ewmcov(x, y, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na, bias=False),
                        var_biased=lambda x: mom.ewmvar(x, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na, bias=True),
                        std_biased=lambda x: mom.ewmstd(x, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na, bias=True),
                        cov_biased=lambda x, y: mom.ewmcov(x, y, com=com, min_periods=min_periods, adjust=adjust, ignore_na=ignore_na, bias=True),
                        var_debiasing_factors=lambda x: _variance_debiasing_factors(x, com=com, adjust=adjust, ignore_na=ignore_na))

    @slow
    def test_expanding_consistency(self):

        # suppress warnings about empty slices, as we are deliberately testing with empty/0-length Series/DataFrames
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*(empty slice|0 for slice).*", category=RuntimeWarning)

            for min_periods in [0, 1, 2, 3, 4]:

                # test consistency between different expanding_* moments
                self._test_moments_consistency(
                    min_periods=min_periods,
                    count=mom.expanding_count,
                    mean=lambda x: mom.expanding_mean(x, min_periods=min_periods),
                    mock_mean=lambda x: mom.expanding_sum(x, min_periods=min_periods) / mom.expanding_count(x),
                    corr=lambda x, y: mom.expanding_corr(x, y, min_periods=min_periods),
                    var_unbiased=lambda x: mom.expanding_var(x, min_periods=min_periods),
                    std_unbiased=lambda x: mom.expanding_std(x, min_periods=min_periods),
                    cov_unbiased=lambda x, y: mom.expanding_cov(x, y, min_periods=min_periods),
                    var_biased=lambda x: mom.expanding_var(x, min_periods=min_periods, ddof=0),
                    std_biased=lambda x: mom.expanding_std(x, min_periods=min_periods, ddof=0),
                    cov_biased=lambda x, y: mom.expanding_cov(x, y, min_periods=min_periods, ddof=0),
                    var_debiasing_factors=lambda x: mom.expanding_count(x) / (mom.expanding_count(x) - 1.).replace(0., np.nan)
                    )

                # test consistency between expanding_xyz() and either (a) expanding_apply of Series.xyz(),
                #                                                  or (b) expanding_apply of np.nanxyz()
                for (x, is_constant, no_nans) in self.data:
                    assert_equal = assert_series_equal if isinstance(x, Series) else assert_frame_equal
                    functions = self.base_functions

                    # GH 8269
                    if no_nans:
                        functions = self.base_functions + self.no_nan_functions
                    for (f, require_min_periods, name) in functions:
                        expanding_f = getattr(mom,'expanding_{0}'.format(name))

                        if require_min_periods and (min_periods is not None) and (min_periods < require_min_periods):
                            continue

                        if expanding_f is mom.expanding_count:
                            expanding_f_result = expanding_f(x)
                            expanding_apply_f_result = mom.expanding_apply(x, func=f, min_periods=0)
                        else:
                            if expanding_f in [mom.expanding_cov, mom.expanding_corr]:
                                expanding_f_result = expanding_f(x, min_periods=min_periods, pairwise=False)
                            else:
                                expanding_f_result = expanding_f(x, min_periods=min_periods)
                            expanding_apply_f_result = mom.expanding_apply(x, func=f, min_periods=min_periods)

                        if not tm._incompat_bottleneck_version(name):
                            assert_equal(expanding_f_result, expanding_apply_f_result)

                        if (expanding_f in [mom.expanding_cov, mom.expanding_corr]) and isinstance(x, DataFrame):
                            # test pairwise=True
                            expanding_f_result = expanding_f(x, x, min_periods=min_periods, pairwise=True)
                            expected = Panel(items=x.index, major_axis=x.columns, minor_axis=x.columns)
                            for i, _ in enumerate(x.columns):
                                for j, _ in enumerate(x.columns):
                                    expected.iloc[:, i, j] = expanding_f(x.iloc[:, i], x.iloc[:, j], min_periods=min_periods)
                            assert_panel_equal(expanding_f_result, expected)

    @slow
    def test_rolling_consistency(self):

        for window in [1, 2, 3, 10, 20]:
            for min_periods in set([0, 1, 2, 3, 4, window]):
                if min_periods and (min_periods > window):
                    continue
                for center in [False, True]:

                    # test consistency between different rolling_* moments
                    self._test_moments_consistency(
                        min_periods=min_periods,
                        count=lambda x: mom.rolling_count(x, window=window, center=center),
                        mean=lambda x: mom.rolling_mean(x, window=window, min_periods=min_periods, center=center),
                        mock_mean=lambda x: mom.rolling_sum(x, window=window, min_periods=min_periods, center=center).divide(
                                            mom.rolling_count(x, window=window, center=center)),
                        corr=lambda x, y: mom.rolling_corr(x, y, window=window, min_periods=min_periods, center=center),
                        var_unbiased=lambda x: mom.rolling_var(x, window=window, min_periods=min_periods, center=center),
                        std_unbiased=lambda x: mom.rolling_std(x, window=window, min_periods=min_periods, center=center),
                        cov_unbiased=lambda x, y: mom.rolling_cov(x, y, window=window, min_periods=min_periods, center=center),
                        var_biased=lambda x: mom.rolling_var(x, window=window, min_periods=min_periods, center=center, ddof=0),
                        std_biased=lambda x: mom.rolling_std(x, window=window, min_periods=min_periods, center=center, ddof=0),
                        cov_biased=lambda x, y: mom.rolling_cov(x, y, window=window, min_periods=min_periods, center=center, ddof=0),
                        var_debiasing_factors=lambda x: mom.rolling_count(x, window=window, center=center).divide(
                                                        (mom.rolling_count(x, window=window, center=center) - 1.).replace(0., np.nan)),
                        )

                    # test consistency between rolling_xyz() and either (a) rolling_apply of Series.xyz(),
                    #                                                or (b) rolling_apply of np.nanxyz()
                    for (x, is_constant, no_nans) in self.data:

                        assert_equal = assert_series_equal if isinstance(x, Series) else assert_frame_equal
                        functions = self.base_functions

                        # GH 8269
                        if no_nans:
                            functions = self.base_functions + self.no_nan_functions
                        for (f, require_min_periods, name) in functions:
                            rolling_f = getattr(mom,'rolling_{0}'.format(name))

                            if require_min_periods and (min_periods is not None) and (min_periods < require_min_periods):
                                continue

                            if rolling_f is mom.rolling_count:
                                rolling_f_result = rolling_f(x, window=window, center=center)
                                rolling_apply_f_result = mom.rolling_apply(x, window=window, func=f,
                                                                           min_periods=0, center=center)
                            else:
                                if rolling_f in [mom.rolling_cov, mom.rolling_corr]:
                                    rolling_f_result = rolling_f(x, window=window, min_periods=min_periods, center=center, pairwise=False)
                                else:
                                    rolling_f_result = rolling_f(x, window=window, min_periods=min_periods, center=center)
                                rolling_apply_f_result = mom.rolling_apply(x, window=window, func=f,
                                                                           min_periods=min_periods, center=center)
                            if not tm._incompat_bottleneck_version(name):
                                assert_equal(rolling_f_result, rolling_apply_f_result)

                            if (rolling_f in [mom.rolling_cov, mom.rolling_corr]) and isinstance(x, DataFrame):
                                # test pairwise=True
                                rolling_f_result = rolling_f(x, x, window=window, min_periods=min_periods,
                                                             center=center, pairwise=True)
                                expected = Panel(items=x.index, major_axis=x.columns, minor_axis=x.columns)
                                for i, _ in enumerate(x.columns):
                                    for j, _ in enumerate(x.columns):
                                        expected.iloc[:, i, j] = rolling_f(x.iloc[:, i], x.iloc[:, j],
                                                                           window=window, min_periods=min_periods, center=center)
                                assert_panel_equal(rolling_f_result, expected)

    # binary moments
    def test_rolling_cov(self):
        A = self.series
        B = A + randn(len(A))

        result = mom.rolling_cov(A, B, 50, min_periods=25)
        assert_almost_equal(result[-1], np.cov(A[-50:], B[-50:])[0, 1])

    def test_rolling_cov_pairwise(self):
        self._check_pairwise_moment(mom.rolling_cov, 10, min_periods=5)

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
        self._check_pairwise_moment(mom.rolling_corr, 10, min_periods=5)

    def _check_pairwise_moment(self, func, *args, **kwargs):
        panel = func(self.frame, *args, **kwargs)

        actual = panel.ix[:, 1, 5]
        expected = func(self.frame[1], self.frame[5], *args, **kwargs)
        tm.assert_series_equal(actual, expected, check_names=False)
        self.assertEqual(actual.name, 5)

    def test_flex_binary_moment(self):
        # GH3155
        # don't blow the stack
        self.assertRaises(TypeError, mom._flex_binary_moment,5,6,None)

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
            try:
                self.assertTrue(all([np.abs(np.nan_to_num(x)) <=1 for x in res]))
            except:
                print(res)


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

    def test_ewmcov_pairwise(self):
        self._check_pairwise_moment(mom.ewmcov, span=10, min_periods=5)

    def test_ewmcorr(self):
        self._check_binary_ew(mom.ewmcorr)

    def test_ewmcorr_pairwise(self):
        self._check_pairwise_moment(mom.ewmcorr, span=10, min_periods=5)

    def _check_binary_ew(self, func):
        A = Series(randn(50), index=np.arange(50))
        B = A[2:] + randn(48)

        A[:10] = np.NaN
        B[-10:] = np.NaN

        result = func(A, B, 20, min_periods=5)
        self.assertTrue(np.isnan(result.values[:14]).all())
        self.assertFalse(np.isnan(result.values[14:]).any())

        # GH 7898
        for min_periods in (0, 1, 2):
            result = func(A, B, 20, min_periods=min_periods)
            # binary functions (ewmcov, ewmcorr) with bias=False require at least two values
            self.assertTrue(np.isnan(result.values[:11]).all())
            self.assertFalse(np.isnan(result.values[11:]).any())

            # check series of length 0
            result = func(Series([]), Series([]), 50, min_periods=min_periods)
            assert_series_equal(result, Series([]))

            # check series of length 1
            result = func(Series([1.]), Series([1.]), 50, min_periods=min_periods)
            assert_series_equal(result, Series([np.NaN]))

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

        # GH 8080
        s = Series([None, None, None])
        result = mom.expanding_apply(s, lambda x: len(x), min_periods=0)
        expected = Series([1., 2., 3.])
        assert_series_equal(result, expected)

    def test_expanding_apply_args_kwargs(self):
        def mean_w_arg(x, const):
            return np.mean(x) + const

        df = DataFrame(np.random.rand(20, 3))

        expected = mom.expanding_apply(df, np.mean) + 20.

        assert_frame_equal(mom.expanding_apply(df, mean_w_arg, args=(20,)),
                            expected)
        assert_frame_equal(mom.expanding_apply(df, mean_w_arg,
                                               kwargs={'const' : 20}),
                            expected)


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

    def test_expanding_cov_pairwise(self):
        result = mom.expanding_cov(self.frame)

        rolling_result = mom.rolling_cov(self.frame, len(self.frame),
                                         min_periods=1)

        for i in result.items:
            assert_almost_equal(result[i], rolling_result[i])

    def test_expanding_corr_pairwise(self):
        result = mom.expanding_corr(self.frame)

        rolling_result = mom.rolling_corr(self.frame, len(self.frame),
                                          min_periods=1)

        for i in result.items:
            assert_almost_equal(result[i], rolling_result[i])

    def test_expanding_cov_diff_index(self):
        # GH 7512
        s1 = Series([1, 2, 3], index=[0, 1, 2])
        s2 = Series([1, 3], index=[0, 2])
        result = mom.expanding_cov(s1, s2)
        expected = Series([None, None, 2.0])
        assert_series_equal(result, expected)

        s2a = Series([1, None, 3], index=[0, 1, 2])
        result = mom.expanding_cov(s1, s2a)
        assert_series_equal(result, expected)

        s1 = Series([7, 8, 10], index=[0, 1, 3])
        s2 = Series([7, 9, 10], index=[0, 2, 3])
        result = mom.expanding_cov(s1, s2)
        expected = Series([None, None, None, 4.5])
        assert_series_equal(result, expected)

    def test_expanding_corr_diff_index(self):
        # GH 7512
        s1 = Series([1, 2, 3], index=[0, 1, 2])
        s2 = Series([1, 3], index=[0, 2])
        result = mom.expanding_corr(s1, s2)
        expected = Series([None, None, 1.0])
        assert_series_equal(result, expected)

        s2a = Series([1, None, 3], index=[0, 1, 2])
        result = mom.expanding_corr(s1, s2a)
        assert_series_equal(result, expected)

        s1 = Series([7, 8, 10], index=[0, 1, 3])
        s2 = Series([7, 9, 10], index=[0, 2, 3])
        result = mom.expanding_corr(s1, s2)
        expected = Series([None, None, None, 1.])
        assert_series_equal(result, expected)

    def test_rolling_cov_diff_length(self):
        # GH 7512
        s1 = Series([1, 2, 3], index=[0, 1, 2])
        s2 = Series([1, 3], index=[0, 2])
        result = mom.rolling_cov(s1, s2, window=3, min_periods=2)
        expected = Series([None, None, 2.0])
        assert_series_equal(result, expected)

        s2a = Series([1, None, 3], index=[0, 1, 2])
        result = mom.rolling_cov(s1, s2a, window=3, min_periods=2)
        assert_series_equal(result, expected)

    def test_rolling_corr_diff_length(self):
        # GH 7512
        s1 = Series([1, 2, 3], index=[0, 1, 2])
        s2 = Series([1, 3], index=[0, 2])
        result = mom.rolling_corr(s1, s2, window=3, min_periods=2)
        expected = Series([None, None, 1.0])
        assert_series_equal(result, expected)

        s2a = Series([1, None, 3], index=[0, 1, 2])
        result = mom.rolling_corr(s1, s2a, window=3, min_periods=2)
        assert_series_equal(result, expected)

    def test_rolling_functions_window_non_shrinkage(self):
        # GH 7764
        s = Series(range(4))
        s_expected = Series(np.nan, index=s.index)
        df = DataFrame([[1,5], [3, 2], [3,9], [-1,0]], columns=['A','B'])
        df_expected = DataFrame(np.nan, index=df.index, columns=df.columns)
        df_expected_panel = Panel(items=df.index, major_axis=df.columns, minor_axis=df.columns)

        functions = [lambda x: mom.rolling_cov(x, x, pairwise=False, window=10, min_periods=5),
                     lambda x: mom.rolling_corr(x, x, pairwise=False, window=10, min_periods=5),
                     lambda x: mom.rolling_max(x, window=10, min_periods=5),
                     lambda x: mom.rolling_min(x, window=10, min_periods=5),
                     lambda x: mom.rolling_sum(x, window=10, min_periods=5),
                     lambda x: mom.rolling_mean(x, window=10, min_periods=5),
                     lambda x: mom.rolling_std(x, window=10, min_periods=5),
                     lambda x: mom.rolling_var(x, window=10, min_periods=5),
                     lambda x: mom.rolling_skew(x, window=10, min_periods=5),
                     lambda x: mom.rolling_kurt(x, window=10, min_periods=5),
                     lambda x: mom.rolling_quantile(x, quantile=0.5, window=10, min_periods=5),
                     lambda x: mom.rolling_median(x, window=10, min_periods=5),
                     lambda x: mom.rolling_apply(x, func=sum, window=10, min_periods=5),
                     lambda x: mom.rolling_window(x, win_type='boxcar', window=10, min_periods=5),
                    ]
        for f in functions:
            try:
                s_result = f(s)
                assert_series_equal(s_result, s_expected)

                df_result = f(df)
                assert_frame_equal(df_result, df_expected)
            except (ImportError):

                # scipy needed for rolling_window
                continue

        functions = [lambda x: mom.rolling_cov(x, x, pairwise=True, window=10, min_periods=5),
                     lambda x: mom.rolling_corr(x, x, pairwise=True, window=10, min_periods=5),
                     # rolling_corr_pairwise is depracated, so the following line should be deleted
                     # when rolling_corr_pairwise is removed.
                     lambda x: mom.rolling_corr_pairwise(x, x, window=10, min_periods=5),
                    ]
        for f in functions:
            df_result_panel = f(df)
            assert_panel_equal(df_result_panel, df_expected_panel)

    def test_moment_functions_zero_length(self):
        # GH 8056
        s = Series()
        s_expected = s
        df1 = DataFrame()
        df1_expected = df1
        df1_expected_panel = Panel(items=df1.index, major_axis=df1.columns, minor_axis=df1.columns)
        df2 = DataFrame(columns=['a'])
        df2['a'] = df2['a'].astype('float64')
        df2_expected = df2
        df2_expected_panel = Panel(items=df2.index, major_axis=df2.columns, minor_axis=df2.columns)

        functions = [lambda x: mom.expanding_count(x),
                     lambda x: mom.expanding_cov(x, x, pairwise=False, min_periods=5),
                     lambda x: mom.expanding_corr(x, x, pairwise=False, min_periods=5),
                     lambda x: mom.expanding_max(x, min_periods=5),
                     lambda x: mom.expanding_min(x, min_periods=5),
                     lambda x: mom.expanding_sum(x, min_periods=5),
                     lambda x: mom.expanding_mean(x, min_periods=5),
                     lambda x: mom.expanding_std(x, min_periods=5),
                     lambda x: mom.expanding_var(x, min_periods=5),
                     lambda x: mom.expanding_skew(x, min_periods=5),
                     lambda x: mom.expanding_kurt(x, min_periods=5),
                     lambda x: mom.expanding_quantile(x, quantile=0.5, min_periods=5),
                     lambda x: mom.expanding_median(x, min_periods=5),
                     lambda x: mom.expanding_apply(x, func=sum, min_periods=5),
                     lambda x: mom.rolling_count(x, window=10),
                     lambda x: mom.rolling_cov(x, x, pairwise=False, window=10, min_periods=5),
                     lambda x: mom.rolling_corr(x, x, pairwise=False, window=10, min_periods=5),
                     lambda x: mom.rolling_max(x, window=10, min_periods=5),
                     lambda x: mom.rolling_min(x, window=10, min_periods=5),
                     lambda x: mom.rolling_sum(x, window=10, min_periods=5),
                     lambda x: mom.rolling_mean(x, window=10, min_periods=5),
                     lambda x: mom.rolling_std(x, window=10, min_periods=5),
                     lambda x: mom.rolling_var(x, window=10, min_periods=5),
                     lambda x: mom.rolling_skew(x, window=10, min_periods=5),
                     lambda x: mom.rolling_kurt(x, window=10, min_periods=5),
                     lambda x: mom.rolling_quantile(x, quantile=0.5, window=10, min_periods=5),
                     lambda x: mom.rolling_median(x, window=10, min_periods=5),
                     lambda x: mom.rolling_apply(x, func=sum, window=10, min_periods=5),
                     lambda x: mom.rolling_window(x, win_type='boxcar', window=10, min_periods=5),
                    ]
        for f in functions:
            try:
                s_result = f(s)
                assert_series_equal(s_result, s_expected)

                df1_result = f(df1)
                assert_frame_equal(df1_result, df1_expected)

                df2_result = f(df2)
                assert_frame_equal(df2_result, df2_expected)
            except (ImportError):

                # scipy needed for rolling_window
                continue

        functions = [lambda x: mom.expanding_cov(x, x, pairwise=True, min_periods=5),
                     lambda x: mom.expanding_corr(x, x, pairwise=True, min_periods=5),
                     lambda x: mom.rolling_cov(x, x, pairwise=True, window=10, min_periods=5),
                     lambda x: mom.rolling_corr(x, x, pairwise=True, window=10, min_periods=5),
                     # rolling_corr_pairwise is depracated, so the following line should be deleted
                     # when rolling_corr_pairwise is removed.
                     lambda x: mom.rolling_corr_pairwise(x, x, window=10, min_periods=5),
                    ]
        for f in functions:
            df1_result_panel = f(df1)
            assert_panel_equal(df1_result_panel, df1_expected_panel)

            df2_result_panel = f(df2)
            assert_panel_equal(df2_result_panel, df2_expected_panel)

    def test_expanding_cov_pairwise_diff_length(self):
        # GH 7512
        df1 = DataFrame([[1,5], [3, 2], [3,9]], columns=['A','B'])
        df1a = DataFrame([[1,5], [3,9]], index=[0,2], columns=['A','B'])
        df2 = DataFrame([[5,6], [None,None], [2,1]], columns=['X','Y'])
        df2a = DataFrame([[5,6], [2,1]], index=[0,2], columns=['X','Y'])
        result1 = mom.expanding_cov(df1, df2, pairwise=True)[2]
        result2 = mom.expanding_cov(df1, df2a, pairwise=True)[2]
        result3 = mom.expanding_cov(df1a, df2, pairwise=True)[2]
        result4 = mom.expanding_cov(df1a, df2a, pairwise=True)[2]
        expected = DataFrame([[-3., -5.], [-6., -10.]], index=['A','B'], columns=['X','Y'])
        assert_frame_equal(result1, expected)
        assert_frame_equal(result2, expected)
        assert_frame_equal(result3, expected)
        assert_frame_equal(result4, expected)

    def test_expanding_corr_pairwise_diff_length(self):
        # GH 7512
        df1 = DataFrame([[1,2], [3, 2], [3,4]], columns=['A','B'])
        df1a = DataFrame([[1,2], [3,4]], index=[0,2], columns=['A','B'])
        df2 = DataFrame([[5,6], [None,None], [2,1]], columns=['X','Y'])
        df2a = DataFrame([[5,6], [2,1]], index=[0,2], columns=['X','Y'])
        result1 = mom.expanding_corr(df1, df2, pairwise=True)[2]
        result2 = mom.expanding_corr(df1, df2a, pairwise=True)[2]
        result3 = mom.expanding_corr(df1a, df2, pairwise=True)[2]
        result4 = mom.expanding_corr(df1a, df2a, pairwise=True)[2]
        expected = DataFrame([[-1.0, -1.0], [-1.0, -1.0]], index=['A','B'], columns=['X','Y'])
        assert_frame_equal(result1, expected)
        assert_frame_equal(result2, expected)
        assert_frame_equal(result3, expected)
        assert_frame_equal(result4, expected)

    def test_pairwise_stats_column_names_order(self):
        # GH 7738
        df1s = [DataFrame([[2,4],[1,2],[5,2],[8,1]], columns=[0,1]),
                DataFrame([[2,4],[1,2],[5,2],[8,1]], columns=[1,0]),
                DataFrame([[2,4],[1,2],[5,2],[8,1]], columns=[1,1]),
                DataFrame([[2,4],[1,2],[5,2],[8,1]], columns=['C','C']),
                DataFrame([[2,4],[1,2],[5,2],[8,1]], columns=[1.,0]),
                DataFrame([[2,4],[1,2],[5,2],[8,1]], columns=[0.,1]),
                DataFrame([[2,4],[1,2],[5,2],[8,1]], columns=['C',1]),
                DataFrame([[2.,4.],[1.,2.],[5.,2.],[8.,1.]], columns=[1,0.]),
                DataFrame([[2,4.],[1,2.],[5,2.],[8,1.]], columns=[0,1.]),
                DataFrame([[2,4],[1,2],[5,2],[8,1.]], columns=[1.,'X']),
               ]
        df2 = DataFrame([[None,1,1],[None,1,2],[None,3,2],[None,8,1]], columns=['Y','Z','X'])
        s = Series([1,1,3,8])

        # suppress warnings about incomparable objects, as we are deliberately testing with such column labels
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*incomparable objects.*", category=RuntimeWarning)

            # DataFrame methods (which do not call _flex_binary_moment())
            for f in [lambda x: x.cov(),
                      lambda x: x.corr(),
                     ]:
                results = [f(df) for df in df1s]
                for (df, result) in zip(df1s, results):
                    assert_index_equal(result.index, df.columns)
                    assert_index_equal(result.columns, df.columns)
                for i, result in enumerate(results):
                    if i > 0:
                        self.assert_numpy_array_equal(result, results[0])

            # DataFrame with itself, pairwise=True
            for f in [lambda x: mom.expanding_cov(x, pairwise=True),
                      lambda x: mom.expanding_corr(x, pairwise=True),
                      lambda x: mom.rolling_cov(x, window=3, pairwise=True),
                      lambda x: mom.rolling_corr(x, window=3, pairwise=True),
                      lambda x: mom.ewmcov(x, com=3, pairwise=True),
                      lambda x: mom.ewmcorr(x, com=3, pairwise=True),
                     ]:
                results = [f(df) for df in df1s]
                for (df, result) in zip(df1s, results):
                    assert_index_equal(result.items, df.index)
                    assert_index_equal(result.major_axis, df.columns)
                    assert_index_equal(result.minor_axis, df.columns)
                for i, result in enumerate(results):
                    if i > 0:
                        self.assert_numpy_array_equal(result, results[0])

            # DataFrame with itself, pairwise=False
            for f in [lambda x: mom.expanding_cov(x, pairwise=False),
                      lambda x: mom.expanding_corr(x, pairwise=False),
                      lambda x: mom.rolling_cov(x, window=3, pairwise=False),
                      lambda x: mom.rolling_corr(x, window=3, pairwise=False),
                      lambda x: mom.ewmcov(x, com=3, pairwise=False),
                      lambda x: mom.ewmcorr(x, com=3, pairwise=False),
                     ]:
                results = [f(df) for df in df1s]
                for (df, result) in zip(df1s, results):
                    assert_index_equal(result.index, df.index)
                    assert_index_equal(result.columns, df.columns)
                for i, result in enumerate(results):
                    if i > 0:
                        self.assert_numpy_array_equal(result, results[0])

            # DataFrame with another DataFrame, pairwise=True
            for f in [lambda x, y: mom.expanding_cov(x, y, pairwise=True),
                      lambda x, y: mom.expanding_corr(x, y, pairwise=True),
                      lambda x, y: mom.rolling_cov(x, y, window=3, pairwise=True),
                      lambda x, y: mom.rolling_corr(x, y, window=3, pairwise=True),
                      lambda x, y: mom.ewmcov(x, y, com=3, pairwise=True),
                      lambda x, y: mom.ewmcorr(x, y, com=3, pairwise=True),
                     ]:
                results = [f(df, df2) for df in df1s]
                for (df, result) in zip(df1s, results):
                    assert_index_equal(result.items, df.index)
                    assert_index_equal(result.major_axis, df.columns)
                    assert_index_equal(result.minor_axis, df2.columns)
                for i, result in enumerate(results):
                    if i > 0:
                        self.assert_numpy_array_equal(result, results[0])

            # DataFrame with another DataFrame, pairwise=False
            for f in [lambda x, y: mom.expanding_cov(x, y, pairwise=False),
                      lambda x, y: mom.expanding_corr(x, y, pairwise=False),
                      lambda x, y: mom.rolling_cov(x, y, window=3, pairwise=False),
                      lambda x, y: mom.rolling_corr(x, y, window=3, pairwise=False),
                      lambda x, y: mom.ewmcov(x, y, com=3, pairwise=False),
                      lambda x, y: mom.ewmcorr(x, y, com=3, pairwise=False),
                     ]:
                results = [f(df, df2) if df.columns.is_unique else None for df in df1s]
                for (df, result) in zip(df1s, results):
                    if result is not None:
                        expected_index = df.index.union(df2.index)
                        expected_columns = df.columns.union(df2.columns)
                        assert_index_equal(result.index, expected_index)
                        assert_index_equal(result.columns, expected_columns)
                    else:
                        tm.assertRaisesRegexp(ValueError, "'arg1' columns are not unique", f, df, df2)
                        tm.assertRaisesRegexp(ValueError, "'arg2' columns are not unique", f, df2, df)

            # DataFrame with a Series
            for f in [lambda x, y: mom.expanding_cov(x, y),
                      lambda x, y: mom.expanding_corr(x, y),
                      lambda x, y: mom.rolling_cov(x, y, window=3),
                      lambda x, y: mom.rolling_corr(x, y, window=3),
                      lambda x, y: mom.ewmcov(x, y, com=3),
                      lambda x, y: mom.ewmcorr(x, y, com=3),
                     ]:
                results = [f(df, s) for df in df1s] + [f(s, df) for df in df1s]
                for (df, result) in zip(df1s, results):
                    assert_index_equal(result.index, df.index)
                    assert_index_equal(result.columns, df.columns)
                for i, result in enumerate(results):
                    if i > 0:
                        self.assert_numpy_array_equal(result, results[0])

    def test_rolling_skew_edge_cases(self):

        all_nan = Series([np.NaN] * 5)

        # yields all NaN (0 variance)
        d = Series([1] * 5)
        x = mom.rolling_skew(d, window=5)
        assert_series_equal(all_nan, x)

        # yields all NaN (window too small)
        d = Series(np.random.randn(5))
        x = mom.rolling_skew(d, window=2)
        assert_series_equal(all_nan, x)

        # yields [NaN, NaN, NaN, 0.177994, 1.548824]
        d = Series([-1.50837035, -0.1297039 ,  0.19501095,
                       1.73508164,  0.41941401])
        expected = Series([np.NaN, np.NaN, np.NaN,
                              0.177994, 1.548824])
        x = mom.rolling_skew(d, window=4)
        assert_series_equal(expected, x)

    def test_rolling_kurt_edge_cases(self):

        all_nan = Series([np.NaN] * 5)

        # yields all NaN (0 variance)
        d = Series([1] * 5)
        x = mom.rolling_kurt(d, window=5)
        assert_series_equal(all_nan, x)

        # yields all NaN (window too small)
        d = Series(np.random.randn(5))
        x = mom.rolling_kurt(d, window=3)
        assert_series_equal(all_nan, x)

        # yields [NaN, NaN, NaN, 1.224307, 2.671499]
        d = Series([-1.50837035, -0.1297039 ,  0.19501095,
                    1.73508164,  0.41941401])
        expected = Series([np.NaN, np.NaN, np.NaN,
                           1.224307, 2.671499])
        x = mom.rolling_kurt(d, window=4)
        assert_series_equal(expected, x)

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
            self.assertTrue(np.isnan(result[13]))
            self.assertFalse(np.isnan(result[14]))

            arr2 = randn(20)
            result = func(arr2, min_periods=5)
            self.assertTrue(isnull(result[3]))
            self.assertTrue(notnull(result[4]))

            # min_periods=0
            result0 = func(arr, min_periods=0)
            result1 = func(arr, min_periods=1)
            assert_almost_equal(result0, result1)
        else:
            result = func(arr)
            assert_almost_equal(result[-1], static_comp(arr[:50]))

    def _check_expanding_structures(self, func):
        series_result = func(self.series)
        tm.assertIsInstance(series_result, Series)
        frame_result = func(self.frame)
        self.assertEqual(type(frame_result), DataFrame)

    def _check_expanding(self, func, static_comp, has_min_periods=True,
                         has_time_rule=True,
                         preserve_nan=True):
        self._check_expanding_ndarray(func, static_comp,
                                      has_min_periods=has_min_periods,
                                      has_time_rule=has_time_rule,
                                      preserve_nan=preserve_nan)
        self._check_expanding_structures(func)

    def test_rolling_max_gh6297(self):
        """Replicate result expected in GH #6297"""

        indices = [datetime(1975, 1, i) for i in range(1, 6)]
        # So that we can have 2 datapoints on one of the days
        indices.append(datetime(1975, 1, 3, 6, 0))
        series = Series(range(1, 7), index=indices)
        # Use floats instead of ints as values
        series = series.map(lambda x: float(x))
        # Sort chronologically
        series = series.sort_index()

        expected = Series([1.0, 2.0, 6.0, 4.0, 5.0],
                          index=[datetime(1975, 1, i, 0)
                                 for i in range(1, 6)])
        x = mom.rolling_max(series, window=1, freq='D')
        assert_series_equal(expected, x)

    def test_rolling_max_how_resample(self):

        indices = [datetime(1975, 1, i) for i in range(1, 6)]
        # So that we can have 3 datapoints on last day (4, 10, and 20)
        indices.append(datetime(1975, 1, 5, 1))
        indices.append(datetime(1975, 1, 5, 2))
        series = Series(list(range(0, 5)) + [10, 20], index=indices)
        # Use floats instead of ints as values
        series = series.map(lambda x: float(x))
        # Sort chronologically
        series = series.sort_index()

        # Default how should be max
        expected = Series([0.0, 1.0, 2.0, 3.0, 20.0],
                          index=[datetime(1975, 1, i, 0)
                                 for i in range(1, 6)])
        x = mom.rolling_max(series, window=1, freq='D')
        assert_series_equal(expected, x)

        # Now specify median (10.0)
        expected = Series([0.0, 1.0, 2.0, 3.0, 10.0],
                          index=[datetime(1975, 1, i, 0)
                                 for i in range(1, 6)])
        x = mom.rolling_max(series, window=1, freq='D', how='median')
        assert_series_equal(expected, x)

        # Now specify mean (4+10+20)/3
        v = (4.0+10.0+20.0)/3.0
        expected = Series([0.0, 1.0, 2.0, 3.0, v],
                          index=[datetime(1975, 1, i, 0)
                                 for i in range(1, 6)])
        x = mom.rolling_max(series, window=1, freq='D', how='mean')
        assert_series_equal(expected, x)


    def test_rolling_min_how_resample(self):

        indices = [datetime(1975, 1, i) for i in range(1, 6)]
        # So that we can have 3 datapoints on last day (4, 10, and 20)
        indices.append(datetime(1975, 1, 5, 1))
        indices.append(datetime(1975, 1, 5, 2))
        series = Series(list(range(0, 5)) + [10, 20], index=indices)
        # Use floats instead of ints as values
        series = series.map(lambda x: float(x))
        # Sort chronologically
        series = series.sort_index()

        # Default how should be min
        expected = Series([0.0, 1.0, 2.0, 3.0, 4.0],
                          index=[datetime(1975, 1, i, 0)
                                 for i in range(1, 6)])
        x = mom.rolling_min(series, window=1, freq='D')
        assert_series_equal(expected, x)

    def test_rolling_median_how_resample(self):

        indices = [datetime(1975, 1, i) for i in range(1, 6)]
        # So that we can have 3 datapoints on last day (4, 10, and 20)
        indices.append(datetime(1975, 1, 5, 1))
        indices.append(datetime(1975, 1, 5, 2))
        series = Series(list(range(0, 5)) + [10, 20], index=indices)
        # Use floats instead of ints as values
        series = series.map(lambda x: float(x))
        # Sort chronologically
        series = series.sort_index()

        # Default how should be median
        expected = Series([0.0, 1.0, 2.0, 3.0, 10],
                          index=[datetime(1975, 1, i, 0)
                                 for i in range(1, 6)])
        x = mom.rolling_median(series, window=1, freq='D')
        assert_series_equal(expected, x)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
