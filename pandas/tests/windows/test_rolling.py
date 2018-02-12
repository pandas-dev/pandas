import pytest
import warnings

from datetime import datetime, timedelta
from numpy.random import randn
import numpy as np

import pandas as pd
from pandas import (Series, DataFrame, bdate_range,
                    isna, notna, Timestamp, Index)
import pandas.core.window as rwindow
import pandas.tseries.offsets as offsets
from pandas.errors import UnsupportedFunctionCall
import pandas.util.testing as tm
import pandas.util._test_decorators as td
from pandas.compat import range, zip

N, K = 100, 10


def assert_equal(left, right):
    if isinstance(left, Series):
        tm.assert_series_equal(left, right)
    else:
        tm.assert_frame_equal(left, right)


class Base(object):

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


class TestRolling(Base):

    def setup_method(self, method):
        self._create_data()

    def test_doc_string(self):

        df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
        df
        df.rolling(2).sum()
        df.rolling(2, min_periods=1).sum()

    def test_constructor(self):
        # GH 12669

        for o in [self.series, self.frame]:
            c = o.rolling

            # valid
            c(window=2)
            c(window=2, min_periods=1)
            c(window=2, min_periods=1, center=True)
            c(window=2, min_periods=1, center=False)

            # GH 13383
            c(0)
            with pytest.raises(ValueError):
                c(-1)

            # not valid
            for w in [2., 'foo', np.array([2])]:
                with pytest.raises(ValueError):
                    c(window=w)
                with pytest.raises(ValueError):
                    c(window=2, min_periods=w)
                with pytest.raises(ValueError):
                    c(window=2, min_periods=1, center=w)

    @td.skip_if_no_scipy
    def test_constructor_with_win_type(self):
        # GH 13383
        for o in [self.series, self.frame]:
            c = o.rolling
            c(0, win_type='boxcar')
            with pytest.raises(ValueError):
                c(-1, win_type='boxcar')

    def test_constructor_with_timedelta_window(self):
        # GH 15440
        n = 10
        df = DataFrame({'value': np.arange(n)},
                       index=pd.date_range('2015-12-24', periods=n, freq="D"))
        expected_data = np.append([0., 1.], np.arange(3., 27., 3))
        for window in [timedelta(days=3), pd.Timedelta(days=3)]:
            result = df.rolling(window=window).sum()
            expected = DataFrame({'value': expected_data},
                                 index=pd.date_range('2015-12-24', periods=n,
                                                     freq="D"))
            tm.assert_frame_equal(result, expected)
            expected = df.rolling('3D').sum()
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        'window', [timedelta(days=3), pd.Timedelta(days=3), '3D'])
    def test_constructor_with_timedelta_window_and_minperiods(self, window):
        # GH 15305
        n = 10
        df = DataFrame({'value': np.arange(n)},
                       index=pd.date_range('2017-08-08', periods=n, freq="D"))
        expected = DataFrame(
            {'value': np.append([np.NaN, 1.], np.arange(3., 27., 3))},
            index=pd.date_range('2017-08-08', periods=n, freq="D"))
        result_roll_sum = df.rolling(window=window, min_periods=2).sum()
        result_roll_generic = df.rolling(window=window,
                                         min_periods=2).apply(sum)
        tm.assert_frame_equal(result_roll_sum, expected)
        tm.assert_frame_equal(result_roll_generic, expected)

    def test_numpy_compat(self):
        # see gh-12811
        r = rwindow.Rolling(Series([2, 4, 6]), window=2)

        msg = "numpy operations are not valid with window objects"

        for func in ('std', 'mean', 'sum', 'max', 'min', 'var'):
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(r, func), 1, 2, 3)
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(r, func), dtype=np.float64)

    def test_closed(self):
        df = DataFrame({'A': [0, 1, 2, 3, 4]})
        # closed only allowed for datetimelike
        with pytest.raises(ValueError):
            df.rolling(window=3, closed='neither')

    @pytest.mark.parametrize('roller', ['1s', 1])
    def tests_empty_df_rolling(self, roller):
        # GH 15819 Verifies that datetime and integer rolling windows can be
        # applied to empty DataFrames
        expected = DataFrame()
        result = DataFrame().rolling(roller).sum()
        tm.assert_frame_equal(result, expected)

        # Verifies that datetime and integer rolling windows can be applied to
        # empty DataFrames with datetime index
        expected = DataFrame(index=pd.DatetimeIndex([]))
        result = DataFrame(index=pd.DatetimeIndex([])).rolling(roller).sum()
        tm.assert_frame_equal(result, expected)

    def test_missing_minp_zero(self):
        # https://github.com/pandas-dev/pandas/pull/18921
        # minp=0
        x = pd.Series([np.nan])
        result = x.rolling(1, min_periods=0).sum()
        expected = pd.Series([0.0])
        tm.assert_series_equal(result, expected)

        # minp=1
        result = x.rolling(1, min_periods=1).sum()
        expected = pd.Series([np.nan])
        tm.assert_series_equal(result, expected)

    def test_missing_minp_zero_variable(self):
        # https://github.com/pandas-dev/pandas/pull/18921
        x = pd.Series([np.nan] * 4,
                      index=pd.DatetimeIndex(['2017-01-01', '2017-01-04',
                                              '2017-01-06', '2017-01-07']))
        result = x.rolling(pd.Timedelta("2d"), min_periods=0).sum()
        expected = pd.Series(0.0, index=x.index)
        tm.assert_series_equal(result, expected)

    def test_multi_index_names(self):

        # GH 16789, 16825
        cols = pd.MultiIndex.from_product([['A', 'B'], ['C', 'D', 'E']],
                                          names=['1', '2'])
        df = DataFrame(np.ones((10, 6)), columns=cols)
        result = df.rolling(3).cov()

        tm.assert_index_equal(result.columns, df.columns)
        assert result.index.names == [None, '1', '2']


class TestMoments(Base):

    def setup_method(self, method):
        self._create_data()

    def test_centered_axis_validation(self):

        # ok
        Series(np.ones(10)).rolling(window=3, center=True, axis=0).mean()

        # bad axis
        with pytest.raises(ValueError):
            Series(np.ones(10)).rolling(window=3, center=True, axis=1).mean()

        # ok ok
        DataFrame(np.ones((10, 10))).rolling(window=3, center=True,
                                             axis=0).mean()
        DataFrame(np.ones((10, 10))).rolling(window=3, center=True,
                                             axis=1).mean()

        # bad axis
        with pytest.raises(ValueError):
            (DataFrame(np.ones((10, 10)))
             .rolling(window=3, center=True, axis=2).mean())

    def test_rolling_sum(self):
        self._check_moment_func(np.nansum, name='sum',
                                zero_min_periods_equal=False)

    def test_rolling_count(self):
        counter = lambda x: np.isfinite(x).astype(float).sum()
        self._check_moment_func(counter, name='count', has_min_periods=False,
                                fill_value=0)

    def test_rolling_mean(self):
        self._check_moment_func(np.mean, name='mean')

    @td.skip_if_no_scipy
    def test_cmov_mean(self):
        # GH 8238
        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48,
                         10.63, 14.48])
        result = Series(vals).rolling(5, center=True).mean()
        expected = Series([np.nan, np.nan, 9.962, 11.27, 11.564, 12.516,
                           12.818, 12.952, np.nan, np.nan])
        tm.assert_series_equal(expected, result)

    @td.skip_if_no_scipy
    def test_cmov_window(self):
        # GH 8238
        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48,
                         10.63, 14.48])
        result = Series(vals).rolling(5, win_type='boxcar', center=True).mean()
        expected = Series([np.nan, np.nan, 9.962, 11.27, 11.564, 12.516,
                           12.818, 12.952, np.nan, np.nan])
        tm.assert_series_equal(expected, result)

    @td.skip_if_no_scipy
    def test_cmov_window_corner(self):
        # GH 8238
        # all nan
        vals = pd.Series([np.nan] * 10)
        result = vals.rolling(5, center=True, win_type='boxcar').mean()
        assert np.isnan(result).all()

        # empty
        vals = pd.Series([])
        result = vals.rolling(5, center=True, win_type='boxcar').mean()
        assert len(result) == 0

        # shorter than window
        vals = pd.Series(np.random.randn(5))
        result = vals.rolling(10, win_type='boxcar').mean()
        assert np.isnan(result).all()
        assert len(result) == 5

    @td.skip_if_no_scipy
    def test_cmov_window_frame(self):
        # Gh 8238
        vals = np.array([[12.18, 3.64], [10.18, 9.16], [13.24, 14.61],
                         [4.51, 8.11], [6.15, 11.44], [9.14, 6.21],
                         [11.31, 10.67], [2.94, 6.51], [9.42, 8.39], [12.44,
                                                                      7.34]])

        xp = np.array([[np.nan, np.nan], [np.nan, np.nan], [9.252, 9.392],
                       [8.644, 9.906], [8.87, 10.208], [6.81, 8.588],
                       [7.792, 8.644], [9.05, 7.824], [np.nan, np.nan
                                                       ], [np.nan, np.nan]])

        # DataFrame
        rs = DataFrame(vals).rolling(5, win_type='boxcar', center=True).mean()
        tm.assert_frame_equal(DataFrame(xp), rs)

        # invalid method
        with pytest.raises(AttributeError):
            (DataFrame(vals).rolling(5, win_type='boxcar', center=True)
             .std())

        # sum
        xp = np.array([[np.nan, np.nan], [np.nan, np.nan], [46.26, 46.96],
                       [43.22, 49.53], [44.35, 51.04], [34.05, 42.94],
                       [38.96, 43.22], [45.25, 39.12], [np.nan, np.nan
                                                        ], [np.nan, np.nan]])

        rs = DataFrame(vals).rolling(5, win_type='boxcar', center=True).sum()
        tm.assert_frame_equal(DataFrame(xp), rs)

    @td.skip_if_no_scipy
    def test_cmov_window_na_min_periods(self):
        # min_periods
        vals = Series(np.random.randn(10))
        vals[4] = np.nan
        vals[8] = np.nan

        xp = vals.rolling(5, min_periods=4, center=True).mean()
        rs = vals.rolling(5, win_type='boxcar', min_periods=4,
                          center=True).mean()
        tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_regular(self):
        # GH 8238
        win_types = ['triang', 'blackman', 'hamming', 'bartlett', 'bohman',
                     'blackmanharris', 'nuttall', 'barthann']

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48,
                         10.63, 14.48])
        xps = {
            'hamming': [np.nan, np.nan, 8.71384, 9.56348, 12.38009, 14.03687,
                        13.8567, 11.81473, np.nan, np.nan],
            'triang': [np.nan, np.nan, 9.28667, 10.34667, 12.00556, 13.33889,
                       13.38, 12.33667, np.nan, np.nan],
            'barthann': [np.nan, np.nan, 8.4425, 9.1925, 12.5575, 14.3675,
                         14.0825, 11.5675, np.nan, np.nan],
            'bohman': [np.nan, np.nan, 7.61599, 9.1764, 12.83559, 14.17267,
                       14.65923, 11.10401, np.nan, np.nan],
            'blackmanharris': [np.nan, np.nan, 6.97691, 9.16438, 13.05052,
                               14.02156, 15.10512, 10.74574, np.nan, np.nan],
            'nuttall': [np.nan, np.nan, 7.04618, 9.16786, 13.02671, 14.03559,
                        15.05657, 10.78514, np.nan, np.nan],
            'blackman': [np.nan, np.nan, 7.73345, 9.17869, 12.79607, 14.20036,
                         14.57726, 11.16988, np.nan, np.nan],
            'bartlett': [np.nan, np.nan, 8.4425, 9.1925, 12.5575, 14.3675,
                         14.0825, 11.5675, np.nan, np.nan]
        }

        for wt in win_types:
            xp = Series(xps[wt])
            rs = Series(vals).rolling(5, win_type=wt, center=True).mean()
            tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_regular_linear_range(self):
        # GH 8238
        win_types = ['triang', 'blackman', 'hamming', 'bartlett', 'bohman',
                     'blackmanharris', 'nuttall', 'barthann']

        vals = np.array(range(10), dtype=np.float)
        xp = vals.copy()
        xp[:2] = np.nan
        xp[-2:] = np.nan
        xp = Series(xp)

        for wt in win_types:
            rs = Series(vals).rolling(5, win_type=wt, center=True).mean()
            tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_regular_missing_data(self):
        # GH 8238
        win_types = ['triang', 'blackman', 'hamming', 'bartlett', 'bohman',
                     'blackmanharris', 'nuttall', 'barthann']

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, np.nan,
                         10.63, 14.48])
        xps = {
            'bartlett': [np.nan, np.nan, 9.70333, 10.5225, 8.4425, 9.1925,
                         12.5575, 14.3675, 15.61667, 13.655],
            'blackman': [np.nan, np.nan, 9.04582, 11.41536, 7.73345, 9.17869,
                         12.79607, 14.20036, 15.8706, 13.655],
            'barthann': [np.nan, np.nan, 9.70333, 10.5225, 8.4425, 9.1925,
                         12.5575, 14.3675, 15.61667, 13.655],
            'bohman': [np.nan, np.nan, 8.9444, 11.56327, 7.61599, 9.1764,
                       12.83559, 14.17267, 15.90976, 13.655],
            'hamming': [np.nan, np.nan, 9.59321, 10.29694, 8.71384, 9.56348,
                        12.38009, 14.20565, 15.24694, 13.69758],
            'nuttall': [np.nan, np.nan, 8.47693, 12.2821, 7.04618, 9.16786,
                        13.02671, 14.03673, 16.08759, 13.65553],
            'triang': [np.nan, np.nan, 9.33167, 9.76125, 9.28667, 10.34667,
                       12.00556, 13.82125, 14.49429, 13.765],
            'blackmanharris': [np.nan, np.nan, 8.42526, 12.36824, 6.97691,
                               9.16438, 13.05052, 14.02175, 16.1098, 13.65509]
        }

        for wt in win_types:
            xp = Series(xps[wt])
            rs = Series(vals).rolling(5, win_type=wt, min_periods=3).mean()
            tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_special(self):
        # GH 8238
        win_types = ['kaiser', 'gaussian', 'general_gaussian', 'slepian']
        kwds = [{'beta': 1.}, {'std': 1.}, {'power': 2.,
                                            'width': 2.}, {'width': 0.5}]

        vals = np.array([6.95, 15.21, 4.72, 9.12, 13.81, 13.49, 16.68, 9.48,
                         10.63, 14.48])

        xps = {
            'gaussian': [np.nan, np.nan, 8.97297, 9.76077, 12.24763, 13.89053,
                         13.65671, 12.01002, np.nan, np.nan],
            'general_gaussian': [np.nan, np.nan, 9.85011, 10.71589, 11.73161,
                                 13.08516, 12.95111, 12.74577, np.nan, np.nan],
            'slepian': [np.nan, np.nan, 9.81073, 10.89359, 11.70284, 12.88331,
                        12.96079, 12.77008, np.nan, np.nan],
            'kaiser': [np.nan, np.nan, 9.86851, 11.02969, 11.65161, 12.75129,
                       12.90702, 12.83757, np.nan, np.nan]
        }

        for wt, k in zip(win_types, kwds):
            xp = Series(xps[wt])
            rs = Series(vals).rolling(5, win_type=wt, center=True).mean(**k)
            tm.assert_series_equal(xp, rs)

    @td.skip_if_no_scipy
    def test_cmov_window_special_linear_range(self):
        # GH 8238
        win_types = ['kaiser', 'gaussian', 'general_gaussian', 'slepian']
        kwds = [{'beta': 1.}, {'std': 1.}, {'power': 2.,
                                            'width': 2.}, {'width': 0.5}]

        vals = np.array(range(10), dtype=np.float)
        xp = vals.copy()
        xp[:2] = np.nan
        xp[-2:] = np.nan
        xp = Series(xp)

        for wt, k in zip(win_types, kwds):
            rs = Series(vals).rolling(5, win_type=wt, center=True).mean(**k)
            tm.assert_series_equal(xp, rs)

    def test_rolling_median(self):
        self._check_moment_func(np.median, name='median')

    def test_rolling_min(self):
        self._check_moment_func(np.min, name='min')

        a = pd.Series([1, 2, 3, 4, 5])
        result = a.rolling(window=100, min_periods=1).min()
        expected = pd.Series(np.ones(len(a)))
        tm.assert_series_equal(result, expected)

        with pytest.raises(ValueError):
            pd.Series([1, 2, 3]).rolling(window=3, min_periods=5).min()

    def test_rolling_max(self):
        self._check_moment_func(np.max, name='max')

        a = pd.Series([1, 2, 3, 4, 5], dtype=np.float64)
        b = a.rolling(window=100, min_periods=1).max()
        tm.assert_almost_equal(a, b)

        with pytest.raises(ValueError):
            pd.Series([1, 2, 3]).rolling(window=3, min_periods=5).max()

    @pytest.mark.parametrize('q', [0.0, .1, .5, .9, 1.0])
    def test_rolling_quantile(self, q):

        def scoreatpercentile(a, per):
            values = np.sort(a, axis=0)

            idx = int(per / 1. * (values.shape[0] - 1))

            if idx == values.shape[0] - 1:
                retval = values[-1]

            else:
                qlow = float(idx) / float(values.shape[0] - 1)
                qhig = float(idx + 1) / float(values.shape[0] - 1)
                vlow = values[idx]
                vhig = values[idx + 1]
                retval = vlow + (vhig - vlow) * (per - qlow) / (qhig - qlow)

            return retval

        def quantile_func(x):
            return scoreatpercentile(x, q)

        self._check_moment_func(quantile_func, name='quantile',
                                quantile=q)

    def test_rolling_quantile_np_percentile(self):
        # #9413: Tests that rolling window's quantile default behavior
        # is analogus to Numpy's percentile
        row = 10
        col = 5
        idx = pd.date_range('20100101', periods=row, freq='B')
        df = DataFrame(np.random.rand(row * col).reshape((row, -1)), index=idx)

        df_quantile = df.quantile([0.25, 0.5, 0.75], axis=0)
        np_percentile = np.percentile(df, [25, 50, 75], axis=0)

        tm.assert_almost_equal(df_quantile.values, np.array(np_percentile))

    def test_rolling_quantile_series(self):
        # #16211: Tests that rolling window's quantile default behavior
        # is analogus to Series' quantile
        arr = np.arange(100)
        s = Series(arr)
        q1 = s.quantile(0.1)
        q2 = s.rolling(100).quantile(0.1).iloc[-1]

        tm.assert_almost_equal(q1, q2)

    def test_rolling_quantile_param(self):
        ser = Series([0.0, .1, .5, .9, 1.0])

        with pytest.raises(ValueError):
            ser.rolling(3).quantile(-0.1)

        with pytest.raises(ValueError):
            ser.rolling(3).quantile(10.0)

        with pytest.raises(TypeError):
            ser.rolling(3).quantile('foo')

    def test_rolling_apply(self):
        # suppress warnings about empty slices, as we are deliberately testing
        # with a 0-length Series
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",
                                    message=".*(empty slice|0 for slice).*",
                                    category=RuntimeWarning)

            ser = Series([])
            tm.assert_series_equal(ser,
                                   ser.rolling(10).apply(lambda x: x.mean()))

            def f(x):
                return x[np.isfinite(x)].mean()

            self._check_moment_func(np.mean, name='apply', func=f)

        # GH 8080
        s = Series([None, None, None])
        result = s.rolling(2, min_periods=0).apply(lambda x: len(x))
        expected = Series([1., 2., 2.])
        tm.assert_series_equal(result, expected)

        result = s.rolling(2, min_periods=0).apply(len)
        tm.assert_series_equal(result, expected)

    def test_rolling_apply_out_of_bounds(self):
        # #1850
        vals = pd.Series([1, 2, 3, 4])

        result = vals.rolling(10).apply(np.sum)
        assert result.isna().all()

        result = vals.rolling(10, min_periods=1).apply(np.sum)
        expected = pd.Series([1, 3, 6, 10], dtype=float)
        tm.assert_almost_equal(result, expected)

    def test_rolling_std(self):
        self._check_moment_func(lambda x: np.std(x, ddof=1),
                                name='std')
        self._check_moment_func(lambda x: np.std(x, ddof=0),
                                name='std', ddof=0)

    def test_rolling_std_1obs(self):
        vals = pd.Series([1., 2., 3., 4., 5.])

        result = vals.rolling(1, min_periods=1).std()
        expected = pd.Series([np.nan] * 5)
        tm.assert_series_equal(result, expected)

        result = vals.rolling(1, min_periods=1).std(ddof=0)
        expected = pd.Series([0.] * 5)
        tm.assert_series_equal(result, expected)

        result = (pd.Series([np.nan, np.nan, 3, 4, 5])
                    .rolling(3, min_periods=2).std())
        assert np.isnan(result[2])

    def test_rolling_std_neg_sqrt(self):
        # unit test from Bottleneck

        # Test move_nanstd for neg sqrt.

        a = pd.Series([0.0011448196318903589, 0.00028718669878572767,
                       0.00028718669878572767, 0.00028718669878572767,
                       0.00028718669878572767])
        b = a.rolling(window=3).std()
        assert np.isfinite(b[2:]).all()

        b = a.ewm(span=3).std()
        assert np.isfinite(b[2:]).all()

    def test_rolling_var(self):
        self._check_moment_func(lambda x: np.var(x, ddof=1),
                                name='var')
        self._check_moment_func(lambda x: np.var(x, ddof=0),
                                name='var', ddof=0)

    @td.skip_if_no_scipy
    def test_rolling_skew(self):
        from scipy.stats import skew
        self._check_moment_func(lambda x: skew(x, bias=False), name='skew')

    @td.skip_if_no_scipy
    def test_rolling_kurt(self):
        from scipy.stats import kurtosis
        self._check_moment_func(lambda x: kurtosis(x, bias=False),
                                name='kurt')

    def _check_moment_func(self, static_comp, name, has_min_periods=True,
                           has_center=True, has_time_rule=True,
                           fill_value=None, zero_min_periods_equal=True,
                           **kwargs):

        def get_result(obj, window, min_periods=None, center=False):
            r = obj.rolling(window=window, min_periods=min_periods,
                            center=center)
            return getattr(r, name)(**kwargs)

        series_result = get_result(self.series, window=50)
        assert isinstance(series_result, Series)
        tm.assert_almost_equal(series_result.iloc[-1],
                               static_comp(self.series[-50:]))

        frame_result = get_result(self.frame, window=50)
        assert isinstance(frame_result, DataFrame)
        tm.assert_series_equal(frame_result.iloc[-1, :],
                               self.frame.iloc[-50:, :].apply(static_comp,
                                                              axis=0),
                               check_names=False)

        # check time_rule works
        if has_time_rule:
            win = 25
            minp = 10
            series = self.series[::2].resample('B').mean()
            frame = self.frame[::2].resample('B').mean()

            if has_min_periods:
                series_result = get_result(series, window=win,
                                           min_periods=minp)
                frame_result = get_result(frame, window=win,
                                          min_periods=minp)
            else:
                series_result = get_result(series, window=win)
                frame_result = get_result(frame, window=win)

            last_date = series_result.index[-1]
            prev_date = last_date - 24 * offsets.BDay()

            trunc_series = self.series[::2].truncate(prev_date, last_date)
            trunc_frame = self.frame[::2].truncate(prev_date, last_date)

            tm.assert_almost_equal(series_result[-1],
                                   static_comp(trunc_series))

            tm.assert_series_equal(frame_result.xs(last_date),
                                   trunc_frame.apply(static_comp),
                                   check_names=False)

        # excluding NaNs correctly
        obj = Series(randn(50))
        obj[:10] = np.NaN
        obj[-10:] = np.NaN
        if has_min_periods:
            result = get_result(obj, 50, min_periods=30)
            tm.assert_almost_equal(result.iloc[-1], static_comp(obj[10:-10]))

            # min_periods is working correctly
            result = get_result(obj, 20, min_periods=15)
            assert isna(result.iloc[23])
            assert not isna(result.iloc[24])

            assert not isna(result.iloc[-6])
            assert isna(result.iloc[-5])

            obj2 = Series(randn(20))
            result = get_result(obj2, 10, min_periods=5)
            assert isna(result.iloc[3])
            assert notna(result.iloc[4])

            if zero_min_periods_equal:
                # min_periods=0 may be equivalent to min_periods=1
                result0 = get_result(obj, 20, min_periods=0)
                result1 = get_result(obj, 20, min_periods=1)
                tm.assert_almost_equal(result0, result1)
        else:
            result = get_result(obj, 50)
            tm.assert_almost_equal(result.iloc[-1], static_comp(obj[10:-10]))

        # window larger than series length (#7297)
        if has_min_periods:
            for minp in (0, len(self.series) - 1, len(self.series)):
                result = get_result(self.series, len(self.series) + 1,
                                    min_periods=minp)
                expected = get_result(self.series, len(self.series),
                                      min_periods=minp)
                nan_mask = isna(result)
                tm.assert_series_equal(nan_mask, isna(expected))

                nan_mask = ~nan_mask
                tm.assert_almost_equal(result[nan_mask],
                                       expected[nan_mask])
        else:
            result = get_result(self.series, len(self.series) + 1)
            expected = get_result(self.series, len(self.series))
            nan_mask = isna(result)
            tm.assert_series_equal(nan_mask, isna(expected))

            nan_mask = ~nan_mask
            tm.assert_almost_equal(result[nan_mask], expected[nan_mask])

        # check center=True
        if has_center:
            if has_min_periods:
                result = get_result(obj, 20, min_periods=15, center=True)
                expected = get_result(
                    pd.concat([obj, Series([np.NaN] * 9)]), 20,
                    min_periods=15)[9:].reset_index(drop=True)
            else:
                result = get_result(obj, 20, center=True)
                expected = get_result(
                    pd.concat([obj, Series([np.NaN] * 9)]),
                    20)[9:].reset_index(drop=True)

            tm.assert_series_equal(result, expected)

            # shifter index
            s = ['x%d' % x for x in range(12)]

            if has_min_periods:
                minp = 10

                series_xp = get_result(
                    self.series.reindex(list(self.series.index) + s),
                    window=25,
                    min_periods=minp).shift(-12).reindex(self.series.index)
                frame_xp = get_result(
                    self.frame.reindex(list(self.frame.index) + s),
                    window=25,
                    min_periods=minp).shift(-12).reindex(self.frame.index)

                series_rs = get_result(self.series, window=25,
                                       min_periods=minp, center=True)
                frame_rs = get_result(self.frame, window=25, min_periods=minp,
                                      center=True)

            else:
                series_xp = get_result(
                    self.series.reindex(list(self.series.index) + s),
                    window=25).shift(-12).reindex(self.series.index)
                frame_xp = get_result(
                    self.frame.reindex(list(self.frame.index) + s),
                    window=25).shift(-12).reindex(self.frame.index)

                series_rs = get_result(self.series, window=25, center=True)
                frame_rs = get_result(self.frame, window=25, center=True)

            if fill_value is not None:
                series_xp = series_xp.fillna(fill_value)
                frame_xp = frame_xp.fillna(fill_value)
            tm.assert_series_equal(series_xp, series_rs)
            tm.assert_frame_equal(frame_xp, frame_rs)


def _rolling_consistency_cases():
    for window in [1, 2, 3, 10, 20]:
        for min_periods in set([0, 1, 2, 3, 4, window]):
            if min_periods and (min_periods > window):
                continue
            for center in [False, True]:
                yield window, min_periods, center


class TestGrouperGrouping(object):

    def setup_method(self, method):
        self.series = Series(np.arange(10))
        self.frame = DataFrame({'A': [1] * 20 + [2] * 12 + [3] * 8,
                                'B': np.arange(40)})

    def test_mutated(self):

        def f():
            self.frame.groupby('A', foo=1)
        pytest.raises(TypeError, f)

        g = self.frame.groupby('A')
        assert not g.mutated
        g = self.frame.groupby('A', mutated=True)
        assert g.mutated

    def test_getitem(self):
        g = self.frame.groupby('A')
        g_mutated = self.frame.groupby('A', mutated=True)

        expected = g_mutated.B.apply(lambda x: x.rolling(2).mean())

        result = g.rolling(2).mean().B
        tm.assert_series_equal(result, expected)

        result = g.rolling(2).B.mean()
        tm.assert_series_equal(result, expected)

        result = g.B.rolling(2).mean()
        tm.assert_series_equal(result, expected)

        result = self.frame.B.groupby(self.frame.A).rolling(2).mean()
        tm.assert_series_equal(result, expected)

    def test_getitem_multiple(self):

        # GH 13174
        g = self.frame.groupby('A')
        r = g.rolling(2)
        g_mutated = self.frame.groupby('A', mutated=True)
        expected = g_mutated.B.apply(lambda x: x.rolling(2).count())

        result = r.B.count()
        tm.assert_series_equal(result, expected)

        result = r.B.count()
        tm.assert_series_equal(result, expected)

    def test_rolling(self):
        g = self.frame.groupby('A')
        r = g.rolling(window=4)

        for f in ['sum', 'mean', 'min', 'max', 'count', 'kurt', 'skew']:

            result = getattr(r, f)()
            expected = g.apply(lambda x: getattr(x.rolling(4), f)())
            tm.assert_frame_equal(result, expected)

        for f in ['std', 'var']:
            result = getattr(r, f)(ddof=1)
            expected = g.apply(lambda x: getattr(x.rolling(4), f)(ddof=1))
            tm.assert_frame_equal(result, expected)

        result = r.quantile(0.5)
        expected = g.apply(lambda x: x.rolling(4).quantile(0.5))
        tm.assert_frame_equal(result, expected)

    def test_rolling_corr_cov(self):
        g = self.frame.groupby('A')
        r = g.rolling(window=4)

        for f in ['corr', 'cov']:
            result = getattr(r, f)(self.frame)

            def func(x):
                return getattr(x.rolling(4), f)(self.frame)
            expected = g.apply(func)
            tm.assert_frame_equal(result, expected)

            result = getattr(r.B, f)(pairwise=True)

            def func(x):
                return getattr(x.B.rolling(4), f)(pairwise=True)
            expected = g.apply(func)
            tm.assert_series_equal(result, expected)

    def test_rolling_apply(self):
        g = self.frame.groupby('A')
        r = g.rolling(window=4)

        # reduction
        result = r.apply(lambda x: x.sum())
        expected = g.apply(lambda x: x.rolling(4).apply(lambda y: y.sum()))
        tm.assert_frame_equal(result, expected)


class TestRollingTS(object):

    # rolling time-series friendly
    # xref GH13327

    def setup_method(self, method):

        self.regular = DataFrame({'A': pd.date_range('20130101',
                                                     periods=5,
                                                     freq='s'),
                                  'B': range(5)}).set_index('A')

        self.ragged = DataFrame({'B': range(5)})
        self.ragged.index = [Timestamp('20130101 09:00:00'),
                             Timestamp('20130101 09:00:02'),
                             Timestamp('20130101 09:00:03'),
                             Timestamp('20130101 09:00:05'),
                             Timestamp('20130101 09:00:06')]

    def test_doc_string(self):

        df = DataFrame({'B': [0, 1, 2, np.nan, 4]},
                       index=[Timestamp('20130101 09:00:00'),
                              Timestamp('20130101 09:00:02'),
                              Timestamp('20130101 09:00:03'),
                              Timestamp('20130101 09:00:05'),
                              Timestamp('20130101 09:00:06')])
        df
        df.rolling('2s').sum()

    def test_valid(self):

        df = self.regular

        # not a valid freq
        with pytest.raises(ValueError):
            df.rolling(window='foobar')

        # not a datetimelike index
        with pytest.raises(ValueError):
            df.reset_index().rolling(window='foobar')

        # non-fixed freqs
        for freq in ['2MS', pd.offsets.MonthBegin(2)]:
            with pytest.raises(ValueError):
                df.rolling(window=freq)

        for freq in ['1D', pd.offsets.Day(2), '2ms']:
            df.rolling(window=freq)

        # non-integer min_periods
        for minp in [1.0, 'foo', np.array([1, 2, 3])]:
            with pytest.raises(ValueError):
                df.rolling(window='1D', min_periods=minp)

        # center is not implemented
        with pytest.raises(NotImplementedError):
            df.rolling(window='1D', center=True)

    def test_on(self):

        df = self.regular

        # not a valid column
        with pytest.raises(ValueError):
            df.rolling(window='2s', on='foobar')

        # column is valid
        df = df.copy()
        df['C'] = pd.date_range('20130101', periods=len(df))
        df.rolling(window='2d', on='C').sum()

        # invalid columns
        with pytest.raises(ValueError):
            df.rolling(window='2d', on='B')

        # ok even though on non-selected
        df.rolling(window='2d', on='C').B.sum()

    def test_monotonic_on(self):

        # on/index must be monotonic
        df = DataFrame({'A': pd.date_range('20130101',
                                           periods=5,
                                           freq='s'),
                        'B': range(5)})

        assert df.A.is_monotonic
        df.rolling('2s', on='A').sum()

        df = df.set_index('A')
        assert df.index.is_monotonic
        df.rolling('2s').sum()

        # non-monotonic
        df.index = reversed(df.index.tolist())
        assert not df.index.is_monotonic

        with pytest.raises(ValueError):
            df.rolling('2s').sum()

        df = df.reset_index()
        with pytest.raises(ValueError):
            df.rolling('2s', on='A').sum()

    def test_frame_on(self):

        df = DataFrame({'B': range(5),
                        'C': pd.date_range('20130101 09:00:00',
                                           periods=5,
                                           freq='3s')})

        df['A'] = [Timestamp('20130101 09:00:00'),
                   Timestamp('20130101 09:00:02'),
                   Timestamp('20130101 09:00:03'),
                   Timestamp('20130101 09:00:05'),
                   Timestamp('20130101 09:00:06')]

        # we are doing simulating using 'on'
        expected = (df.set_index('A')
                    .rolling('2s')
                    .B
                    .sum()
                    .reset_index(drop=True)
                    )

        result = (df.rolling('2s', on='A')
                    .B
                    .sum()
                  )
        tm.assert_series_equal(result, expected)

        # test as a frame
        # we should be ignoring the 'on' as an aggregation column
        # note that the expected is setting, computing, and resetting
        # so the columns need to be switched compared
        # to the actual result where they are ordered as in the
        # original
        expected = (df.set_index('A')
                      .rolling('2s')[['B']]
                      .sum()
                      .reset_index()[['B', 'A']]
                    )

        result = (df.rolling('2s', on='A')[['B']]
                    .sum()
                  )
        tm.assert_frame_equal(result, expected)

    def test_frame_on2(self):

        # using multiple aggregation columns
        df = DataFrame({'A': [0, 1, 2, 3, 4],
                        'B': [0, 1, 2, np.nan, 4],
                        'C': Index([Timestamp('20130101 09:00:00'),
                                    Timestamp('20130101 09:00:02'),
                                    Timestamp('20130101 09:00:03'),
                                    Timestamp('20130101 09:00:05'),
                                    Timestamp('20130101 09:00:06')])},
                       columns=['A', 'C', 'B'])

        expected1 = DataFrame({'A': [0., 1, 3, 3, 7],
                               'B': [0, 1, 3, np.nan, 4],
                               'C': df['C']},
                              columns=['A', 'C', 'B'])

        result = df.rolling('2s', on='C').sum()
        expected = expected1
        tm.assert_frame_equal(result, expected)

        expected = Series([0, 1, 3, np.nan, 4], name='B')
        result = df.rolling('2s', on='C').B.sum()
        tm.assert_series_equal(result, expected)

        expected = expected1[['A', 'B', 'C']]
        result = df.rolling('2s', on='C')[['A', 'B', 'C']].sum()
        tm.assert_frame_equal(result, expected)

    def test_basic_regular(self):

        df = self.regular.copy()

        df.index = pd.date_range('20130101', periods=5, freq='D')
        expected = df.rolling(window=1, min_periods=1).sum()
        result = df.rolling(window='1D').sum()
        tm.assert_frame_equal(result, expected)

        df.index = pd.date_range('20130101', periods=5, freq='2D')
        expected = df.rolling(window=1, min_periods=1).sum()
        result = df.rolling(window='2D', min_periods=1).sum()
        tm.assert_frame_equal(result, expected)

        expected = df.rolling(window=1, min_periods=1).sum()
        result = df.rolling(window='2D', min_periods=1).sum()
        tm.assert_frame_equal(result, expected)

        expected = df.rolling(window=1).sum()
        result = df.rolling(window='2D').sum()
        tm.assert_frame_equal(result, expected)

    def test_min_periods(self):

        # compare for min_periods
        df = self.regular

        # these slightly different
        expected = df.rolling(2, min_periods=1).sum()
        result = df.rolling('2s').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.rolling(2, min_periods=1).sum()
        result = df.rolling('2s', min_periods=1).sum()
        tm.assert_frame_equal(result, expected)

    def test_closed(self):

        # xref GH13965

        df = DataFrame({'A': [1] * 5},
                       index=[Timestamp('20130101 09:00:01'),
                              Timestamp('20130101 09:00:02'),
                              Timestamp('20130101 09:00:03'),
                              Timestamp('20130101 09:00:04'),
                              Timestamp('20130101 09:00:06')])

        # closed must be 'right', 'left', 'both', 'neither'
        with pytest.raises(ValueError):
            self.regular.rolling(window='2s', closed="blabla")

        expected = df.copy()
        expected["A"] = [1.0, 2, 2, 2, 1]
        result = df.rolling('2s', closed='right').sum()
        tm.assert_frame_equal(result, expected)

        # default should be 'right'
        result = df.rolling('2s').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.copy()
        expected["A"] = [1.0, 2, 3, 3, 2]
        result = df.rolling('2s', closed='both').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.copy()
        expected["A"] = [np.nan, 1.0, 2, 2, 1]
        result = df.rolling('2s', closed='left').sum()
        tm.assert_frame_equal(result, expected)

        expected = df.copy()
        expected["A"] = [np.nan, 1.0, 1, 1, np.nan]
        result = df.rolling('2s', closed='neither').sum()
        tm.assert_frame_equal(result, expected)

    def test_ragged_sum(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 3, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=2).sum()
        expected = df.copy()
        expected['B'] = [np.nan, np.nan, 3, np.nan, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 5, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s').sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 5, 7]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='4s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 6, 9]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='4s', min_periods=3).sum()
        expected = df.copy()
        expected['B'] = [np.nan, np.nan, 3, 6, 9]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).sum()
        expected = df.copy()
        expected['B'] = [0.0, 1, 3, 6, 10]
        tm.assert_frame_equal(result, expected)

    def test_ragged_mean(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).mean()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).mean()
        expected = df.copy()
        expected['B'] = [0.0, 1, 1.5, 3.0, 3.5]
        tm.assert_frame_equal(result, expected)

    def test_ragged_median(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).median()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).median()
        expected = df.copy()
        expected['B'] = [0.0, 1, 1.5, 3.0, 3.5]
        tm.assert_frame_equal(result, expected)

    def test_ragged_quantile(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).quantile(0.5)
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).quantile(0.5)
        expected = df.copy()
        expected['B'] = [0.0, 1, 1.5, 3.0, 3.5]
        tm.assert_frame_equal(result, expected)

    def test_ragged_std(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).std(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='1s', min_periods=1).std(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s', min_periods=1).std(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] + [0.5] * 4
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).std(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan, 0.707107, 1.0, 1.0, 1.290994]
        tm.assert_frame_equal(result, expected)

    def test_ragged_var(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).var(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='1s', min_periods=1).var(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='3s', min_periods=1).var(ddof=0)
        expected = df.copy()
        expected['B'] = [0.0] + [0.25] * 4
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).var(ddof=1)
        expected = df.copy()
        expected['B'] = [np.nan, 0.5, 1.0, 1.0, 1 + 2 / 3.]
        tm.assert_frame_equal(result, expected)

    def test_ragged_skew(self):

        df = self.ragged
        result = df.rolling(window='3s', min_periods=1).skew()
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).skew()
        expected = df.copy()
        expected['B'] = [np.nan] * 2 + [0.0, 0.0, 0.0]
        tm.assert_frame_equal(result, expected)

    def test_ragged_kurt(self):

        df = self.ragged
        result = df.rolling(window='3s', min_periods=1).kurt()
        expected = df.copy()
        expected['B'] = [np.nan] * 5
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).kurt()
        expected = df.copy()
        expected['B'] = [np.nan] * 4 + [-1.2]
        tm.assert_frame_equal(result, expected)

    def test_ragged_count(self):

        df = self.ragged
        result = df.rolling(window='1s', min_periods=1).count()
        expected = df.copy()
        expected['B'] = [1.0, 1, 1, 1, 1]
        tm.assert_frame_equal(result, expected)

        df = self.ragged
        result = df.rolling(window='1s').count()
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).count()
        expected = df.copy()
        expected['B'] = [1.0, 1, 2, 1, 2]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=2).count()
        expected = df.copy()
        expected['B'] = [np.nan, np.nan, 2, np.nan, 2]
        tm.assert_frame_equal(result, expected)

    def test_regular_min(self):

        df = DataFrame({'A': pd.date_range('20130101',
                                           periods=5,
                                           freq='s'),
                        'B': [0.0, 1, 2, 3, 4]}).set_index('A')
        result = df.rolling('1s').min()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        df = DataFrame({'A': pd.date_range('20130101',
                                           periods=5,
                                           freq='s'),
                        'B': [5, 4, 3, 4, 5]}).set_index('A')

        tm.assert_frame_equal(result, expected)
        result = df.rolling('2s').min()
        expected = df.copy()
        expected['B'] = [5.0, 4, 3, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling('5s').min()
        expected = df.copy()
        expected['B'] = [5.0, 4, 3, 3, 3]
        tm.assert_frame_equal(result, expected)

    def test_ragged_min(self):

        df = self.ragged

        result = df.rolling(window='1s', min_periods=1).min()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).min()
        expected = df.copy()
        expected['B'] = [0.0, 1, 1, 3, 3]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).min()
        expected = df.copy()
        expected['B'] = [0.0, 0, 0, 1, 1]
        tm.assert_frame_equal(result, expected)

    def test_perf_min(self):

        N = 10000

        dfp = DataFrame({'B': np.random.randn(N)},
                        index=pd.date_range('20130101',
                                            periods=N,
                                            freq='s'))
        expected = dfp.rolling(2, min_periods=1).min()
        result = dfp.rolling('2s').min()
        assert ((result - expected) < 0.01).all().bool()

        expected = dfp.rolling(200, min_periods=1).min()
        result = dfp.rolling('200s').min()
        assert ((result - expected) < 0.01).all().bool()

    def test_ragged_max(self):

        df = self.ragged

        result = df.rolling(window='1s', min_periods=1).max()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).max()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).max()
        expected = df.copy()
        expected['B'] = [0.0, 1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

    def test_ragged_apply(self):

        df = self.ragged

        f = lambda x: 1
        result = df.rolling(window='1s', min_periods=1).apply(f)
        expected = df.copy()
        expected['B'] = 1.
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='2s', min_periods=1).apply(f)
        expected = df.copy()
        expected['B'] = 1.
        tm.assert_frame_equal(result, expected)

        result = df.rolling(window='5s', min_periods=1).apply(f)
        expected = df.copy()
        expected['B'] = 1.
        tm.assert_frame_equal(result, expected)

    def test_all(self):

        # simple comparison of integer vs time-based windowing
        df = self.regular * 2
        er = df.rolling(window=1)
        r = df.rolling(window='1s')

        for f in ['sum', 'mean', 'count', 'median', 'std',
                  'var', 'kurt', 'skew', 'min', 'max']:

            result = getattr(r, f)()
            expected = getattr(er, f)()
            tm.assert_frame_equal(result, expected)

        result = r.quantile(0.5)
        expected = er.quantile(0.5)
        tm.assert_frame_equal(result, expected)

        result = r.apply(lambda x: 1)
        expected = er.apply(lambda x: 1)
        tm.assert_frame_equal(result, expected)

    def test_all2(self):

        # more sophisticated comparison of integer vs.
        # time-based windowing
        df = DataFrame({'B': np.arange(50)},
                       index=pd.date_range('20130101',
                                           periods=50, freq='H')
                       )
        # in-range data
        dft = df.between_time("09:00", "16:00")

        r = dft.rolling(window='5H')

        for f in ['sum', 'mean', 'count', 'median', 'std',
                  'var', 'kurt', 'skew', 'min', 'max']:

            result = getattr(r, f)()

            # we need to roll the days separately
            # to compare with a time-based roll
            # finally groupby-apply will return a multi-index
            # so we need to drop the day
            def agg_by_day(x):
                x = x.between_time("09:00", "16:00")
                return getattr(x.rolling(5, min_periods=1), f)()
            expected = df.groupby(df.index.day).apply(
                agg_by_day).reset_index(level=0, drop=True)

            tm.assert_frame_equal(result, expected)

    def test_groupby_monotonic(self):

        # GH 15130
        # we don't need to validate monotonicity when grouping

        data = [
            ['David', '1/1/2015', 100], ['David', '1/5/2015', 500],
            ['David', '5/30/2015', 50], ['David', '7/25/2015', 50],
            ['Ryan', '1/4/2014', 100], ['Ryan', '1/19/2015', 500],
            ['Ryan', '3/31/2016', 50], ['Joe', '7/1/2015', 100],
            ['Joe', '9/9/2015', 500], ['Joe', '10/15/2015', 50]]

        df = DataFrame(data=data, columns=['name', 'date', 'amount'])
        df['date'] = pd.to_datetime(df['date'])

        expected = df.set_index('date').groupby('name').apply(
            lambda x: x.rolling('180D')['amount'].sum())
        result = df.groupby('name').rolling('180D', on='date')['amount'].sum()
        tm.assert_series_equal(result, expected)

    def test_non_monotonic(self):
        # GH 13966 (similar to #15130, closed by #15175)

        dates = pd.date_range(start='2016-01-01 09:30:00',
                              periods=20, freq='s')
        df = DataFrame({'A': [1] * 20 + [2] * 12 + [3] * 8,
                        'B': np.concatenate((dates, dates)),
                        'C': np.arange(40)})

        result = df.groupby('A').rolling('4s', on='B').C.mean()
        expected = df.set_index('B').groupby('A').apply(
            lambda x: x.rolling('4s')['C'].mean())
        tm.assert_series_equal(result, expected)

        df2 = df.sort_values('B')
        result = df2.groupby('A').rolling('4s', on='B').C.mean()
        tm.assert_series_equal(result, expected)

    def test_rolling_cov_offset(self):
        # GH16058

        idx = pd.date_range('2017-01-01', periods=24, freq='1h')
        ss = Series(np.arange(len(idx)), index=idx)

        result = ss.rolling('2h').cov()
        expected = Series([np.nan] + [0.5] * (len(idx) - 1), index=idx)
        tm.assert_series_equal(result, expected)

        expected2 = ss.rolling(2, min_periods=1).cov()
        tm.assert_series_equal(result, expected2)

        result = ss.rolling('3h').cov()
        expected = Series([np.nan, 0.5] + [1.0] * (len(idx) - 2), index=idx)
        tm.assert_series_equal(result, expected)

        expected2 = ss.rolling(3, min_periods=1).cov()
        tm.assert_series_equal(result, expected2)
