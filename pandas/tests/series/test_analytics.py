# coding=utf-8
# pylint: disable-msg=E1101,W0612

from itertools import product
from distutils.version import LooseVersion

import pytest

from numpy import nan
import numpy as np
import pandas as pd

from pandas import (Series, Categorical, DataFrame, isnull, notnull,
                    bdate_range, date_range, _np_version_under1p10)
from pandas.core.index import MultiIndex
from pandas.tseries.index import Timestamp
from pandas.tseries.tdi import Timedelta
import pandas.core.config as cf

import pandas.core.nanops as nanops

from pandas.compat import lrange, range
from pandas import compat
from pandas.util.testing import (assert_series_equal, assert_almost_equal,
                                 assert_frame_equal, assert_index_equal)
import pandas.util.testing as tm

from .common import TestData


class TestSeriesAnalytics(TestData, tm.TestCase):

    def test_sum_zero(self):
        arr = np.array([])
        self.assertEqual(nanops.nansum(arr), 0)

        arr = np.empty((10, 0))
        self.assertTrue((nanops.nansum(arr, axis=1) == 0).all())

        # GH #844
        s = Series([], index=[])
        self.assertEqual(s.sum(), 0)

        df = DataFrame(np.empty((10, 0)))
        self.assertTrue((df.sum(1) == 0).all())

    def test_nansum_buglet(self):
        s = Series([1.0, np.nan], index=[0, 1])
        result = np.nansum(s)
        assert_almost_equal(result, 1)

    def test_overflow(self):
        # GH 6915
        # overflowing on the smaller int dtypes
        for dtype in ['int32', 'int64']:
            v = np.arange(5000000, dtype=dtype)
            s = Series(v)

            # no bottleneck
            result = s.sum(skipna=False)
            self.assertEqual(int(result), v.sum(dtype='int64'))
            result = s.min(skipna=False)
            self.assertEqual(int(result), 0)
            result = s.max(skipna=False)
            self.assertEqual(int(result), v[-1])

            # use bottleneck if available
            result = s.sum()
            self.assertEqual(int(result), v.sum(dtype='int64'))
            result = s.min()
            self.assertEqual(int(result), 0)
            result = s.max()
            self.assertEqual(int(result), v[-1])

        for dtype in ['float32', 'float64']:
            v = np.arange(5000000, dtype=dtype)
            s = Series(v)

            # no bottleneck
            result = s.sum(skipna=False)
            self.assertEqual(result, v.sum(dtype=dtype))
            result = s.min(skipna=False)
            self.assertTrue(np.allclose(float(result), 0.0))
            result = s.max(skipna=False)
            self.assertTrue(np.allclose(float(result), v[-1]))

            # use bottleneck if available
            result = s.sum()
            self.assertEqual(result, v.sum(dtype=dtype))
            result = s.min()
            self.assertTrue(np.allclose(float(result), 0.0))
            result = s.max()
            self.assertTrue(np.allclose(float(result), v[-1]))

    def test_sum(self):
        self._check_stat_op('sum', np.sum, check_allna=True)

    def test_sum_inf(self):
        import pandas.core.nanops as nanops

        s = Series(np.random.randn(10))
        s2 = s.copy()

        s[5:8] = np.inf
        s2[5:8] = np.nan

        self.assertTrue(np.isinf(s.sum()))

        arr = np.random.randn(100, 100).astype('f4')
        arr[:, 2] = np.inf

        with cf.option_context("mode.use_inf_as_null", True):
            assert_almost_equal(s.sum(), s2.sum())

        res = nanops.nansum(arr, axis=1)
        self.assertTrue(np.isinf(res).all())

    def test_mean(self):
        self._check_stat_op('mean', np.mean)

    def test_median(self):
        self._check_stat_op('median', np.median)

        # test with integers, test failure
        int_ts = Series(np.ones(10, dtype=int), index=lrange(10))
        self.assertAlmostEqual(np.median(int_ts), int_ts.median())

    def test_mode(self):
        # No mode should be found.
        exp = Series([], dtype=np.float64)
        tm.assert_series_equal(Series([]).mode(), exp)

        exp = Series([], dtype=np.int64)
        tm.assert_series_equal(Series([1]).mode(), exp)

        exp = Series([], dtype=np.object)
        tm.assert_series_equal(Series(['a', 'b', 'c']).mode(), exp)

        # Test numerical data types.
        exp_single = [1]
        data_single = [1] * 5 + [2] * 3

        exp_multi = [1, 3]
        data_multi = [1] * 5 + [2] * 3 + [3] * 5

        for dt in np.typecodes['AllInteger'] + np.typecodes['Float']:
            s = Series(data_single, dtype=dt)
            exp = Series(exp_single, dtype=dt)
            tm.assert_series_equal(s.mode(), exp)

            s = Series(data_multi, dtype=dt)
            exp = Series(exp_multi, dtype=dt)
            tm.assert_series_equal(s.mode(), exp)

        # Test string and object types.
        exp = ['b']
        data = ['a'] * 2 + ['b'] * 3

        s = Series(data, dtype='c')
        exp = Series(exp, dtype='c')
        tm.assert_series_equal(s.mode(), exp)

        exp = ['bar']
        data = ['foo'] * 2 + ['bar'] * 3

        for dt in [str, object]:
            s = Series(data, dtype=dt)
            exp = Series(exp, dtype=dt)
            tm.assert_series_equal(s.mode(), exp)

        # Test datetime types.
        exp = Series([], dtype="M8[ns]")
        s = Series(['2011-01-03', '2013-01-02',
                    '1900-05-03'], dtype='M8[ns]')
        tm.assert_series_equal(s.mode(), exp)

        exp = Series(['2011-01-03', '2013-01-02'], dtype='M8[ns]')
        s = Series(['2011-01-03', '2013-01-02', '1900-05-03',
                    '2011-01-03', '2013-01-02'], dtype='M8[ns]')
        tm.assert_series_equal(s.mode(), exp)

        # gh-5986: Test timedelta types.
        exp = Series([], dtype='timedelta64[ns]')
        s = Series(['1 days', '-1 days', '0 days'],
                   dtype='timedelta64[ns]')
        tm.assert_series_equal(s.mode(), exp)

        exp = Series(['2 min', '1 day'], dtype='timedelta64[ns]')
        s = Series(['1 day', '1 day', '-1 day', '-1 day 2 min',
                    '2 min', '2 min'], dtype='timedelta64[ns]')
        tm.assert_series_equal(s.mode(), exp)

        # Test mixed dtype.
        exp = Series(['foo'])
        s = Series([1, 'foo', 'foo'])
        tm.assert_series_equal(s.mode(), exp)

        # Test for uint64 overflow.
        exp = Series([2**63], dtype=np.uint64)
        s = Series([1, 2**63, 2**63], dtype=np.uint64)
        tm.assert_series_equal(s.mode(), exp)

        exp = Series([], dtype=np.uint64)
        s = Series([1, 2**63], dtype=np.uint64)
        tm.assert_series_equal(s.mode(), exp)

        # Test category dtype.
        c = Categorical([1, 2])
        exp = Categorical([], categories=[1, 2])
        exp = Series(exp, dtype='category')
        tm.assert_series_equal(Series(c).mode(), exp)

        c = Categorical([1, 'a', 'a'])
        exp = Categorical(['a'], categories=[1, 'a'])
        exp = Series(exp, dtype='category')
        tm.assert_series_equal(Series(c).mode(), exp)

        c = Categorical([1, 1, 2, 3, 3])
        exp = Categorical([1, 3], categories=[1, 2, 3])
        exp = Series(exp, dtype='category')
        tm.assert_series_equal(Series(c).mode(), exp)

    def test_prod(self):
        self._check_stat_op('prod', np.prod)

    def test_min(self):
        self._check_stat_op('min', np.min, check_objects=True)

    def test_max(self):
        self._check_stat_op('max', np.max, check_objects=True)

    def test_var_std(self):
        alt = lambda x: np.std(x, ddof=1)
        self._check_stat_op('std', alt)

        alt = lambda x: np.var(x, ddof=1)
        self._check_stat_op('var', alt)

        result = self.ts.std(ddof=4)
        expected = np.std(self.ts.values, ddof=4)
        assert_almost_equal(result, expected)

        result = self.ts.var(ddof=4)
        expected = np.var(self.ts.values, ddof=4)
        assert_almost_equal(result, expected)

        # 1 - element series with ddof=1
        s = self.ts.iloc[[0]]
        result = s.var(ddof=1)
        self.assertTrue(isnull(result))

        result = s.std(ddof=1)
        self.assertTrue(isnull(result))

    def test_sem(self):
        alt = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
        self._check_stat_op('sem', alt)

        result = self.ts.sem(ddof=4)
        expected = np.std(self.ts.values,
                          ddof=4) / np.sqrt(len(self.ts.values))
        assert_almost_equal(result, expected)

        # 1 - element series with ddof=1
        s = self.ts.iloc[[0]]
        result = s.sem(ddof=1)
        self.assertTrue(isnull(result))

    def test_skew(self):
        tm._skip_if_no_scipy()

        from scipy.stats import skew
        alt = lambda x: skew(x, bias=False)
        self._check_stat_op('skew', alt)

        # test corner cases, skew() returns NaN unless there's at least 3
        # values
        min_N = 3
        for i in range(1, min_N + 1):
            s = Series(np.ones(i))
            df = DataFrame(np.ones((i, i)))
            if i < min_N:
                self.assertTrue(np.isnan(s.skew()))
                self.assertTrue(np.isnan(df.skew()).all())
            else:
                self.assertEqual(0, s.skew())
                self.assertTrue((df.skew() == 0).all())

    def test_kurt(self):
        tm._skip_if_no_scipy()

        from scipy.stats import kurtosis
        alt = lambda x: kurtosis(x, bias=False)
        self._check_stat_op('kurt', alt)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        self.assertAlmostEqual(s.kurt(), s.kurt(level=0)['bar'])

        # test corner cases, kurt() returns NaN unless there's at least 4
        # values
        min_N = 4
        for i in range(1, min_N + 1):
            s = Series(np.ones(i))
            df = DataFrame(np.ones((i, i)))
            if i < min_N:
                self.assertTrue(np.isnan(s.kurt()))
                self.assertTrue(np.isnan(df.kurt()).all())
            else:
                self.assertEqual(0, s.kurt())
                self.assertTrue((df.kurt() == 0).all())

    def test_describe(self):
        s = Series([0, 1, 2, 3, 4], name='int_data')
        result = s.describe()
        expected = Series([5, 2, s.std(), 0, 1, 2, 3, 4],
                          name='int_data',
                          index=['count', 'mean', 'std', 'min', '25%',
                                 '50%', '75%', 'max'])
        self.assert_series_equal(result, expected)

        s = Series([True, True, False, False, False], name='bool_data')
        result = s.describe()
        expected = Series([5, 2, False, 3], name='bool_data',
                          index=['count', 'unique', 'top', 'freq'])
        self.assert_series_equal(result, expected)

        s = Series(['a', 'a', 'b', 'c', 'd'], name='str_data')
        result = s.describe()
        expected = Series([5, 4, 'a', 2], name='str_data',
                          index=['count', 'unique', 'top', 'freq'])
        self.assert_series_equal(result, expected)

    def test_argsort(self):
        self._check_accum_op('argsort', check_dtype=False)
        argsorted = self.ts.argsort()
        self.assertTrue(issubclass(argsorted.dtype.type, np.integer))

        # GH 2967 (introduced bug in 0.11-dev I think)
        s = Series([Timestamp('201301%02d' % (i + 1)) for i in range(5)])
        self.assertEqual(s.dtype, 'datetime64[ns]')
        shifted = s.shift(-1)
        self.assertEqual(shifted.dtype, 'datetime64[ns]')
        self.assertTrue(isnull(shifted[4]))

        result = s.argsort()
        expected = Series(lrange(5), dtype='int64')
        assert_series_equal(result, expected)

        result = shifted.argsort()
        expected = Series(lrange(4) + [-1], dtype='int64')
        assert_series_equal(result, expected)

    def test_argsort_stable(self):
        s = Series(np.random.randint(0, 100, size=10000))
        mindexer = s.argsort(kind='mergesort')
        qindexer = s.argsort()

        mexpected = np.argsort(s.values, kind='mergesort')
        qexpected = np.argsort(s.values, kind='quicksort')

        self.assert_series_equal(mindexer, Series(mexpected),
                                 check_dtype=False)
        self.assert_series_equal(qindexer, Series(qexpected),
                                 check_dtype=False)
        self.assertFalse(np.array_equal(qindexer, mindexer))

    def test_cumsum(self):
        self._check_accum_op('cumsum')

    def test_cumprod(self):
        self._check_accum_op('cumprod')

    def test_cummin(self):
        self.assert_numpy_array_equal(self.ts.cummin().values,
                                      np.minimum.accumulate(np.array(self.ts)))
        ts = self.ts.copy()
        ts[::2] = np.NaN
        result = ts.cummin()[1::2]
        expected = np.minimum.accumulate(ts.valid())

        self.assert_series_equal(result, expected)

    def test_cummax(self):
        self.assert_numpy_array_equal(self.ts.cummax().values,
                                      np.maximum.accumulate(np.array(self.ts)))
        ts = self.ts.copy()
        ts[::2] = np.NaN
        result = ts.cummax()[1::2]
        expected = np.maximum.accumulate(ts.valid())

        self.assert_series_equal(result, expected)

    def test_cummin_datetime64(self):
        s = pd.Series(pd.to_datetime(['NaT', '2000-1-2', 'NaT', '2000-1-1',
                                      'NaT', '2000-1-3']))

        expected = pd.Series(pd.to_datetime(['NaT', '2000-1-2', 'NaT',
                                             '2000-1-1', 'NaT', '2000-1-1']))
        result = s.cummin(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', '2000-1-2', '2000-1-1', '2000-1-1', '2000-1-1'
             ]))
        result = s.cummin(skipna=False)
        self.assert_series_equal(expected, result)

    def test_cummax_datetime64(self):
        s = pd.Series(pd.to_datetime(['NaT', '2000-1-2', 'NaT', '2000-1-1',
                                      'NaT', '2000-1-3']))

        expected = pd.Series(pd.to_datetime(['NaT', '2000-1-2', 'NaT',
                                             '2000-1-2', 'NaT', '2000-1-3']))
        result = s.cummax(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', '2000-1-2', '2000-1-2', '2000-1-2', '2000-1-3'
             ]))
        result = s.cummax(skipna=False)
        self.assert_series_equal(expected, result)

    def test_cummin_timedelta64(self):
        s = pd.Series(pd.to_timedelta(['NaT',
                                       '2 min',
                                       'NaT',
                                       '1 min',
                                       'NaT',
                                       '3 min', ]))

        expected = pd.Series(pd.to_timedelta(['NaT',
                                              '2 min',
                                              'NaT',
                                              '1 min',
                                              'NaT',
                                              '1 min', ]))
        result = s.cummin(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_timedelta(['NaT',
                                              '2 min',
                                              '2 min',
                                              '1 min',
                                              '1 min',
                                              '1 min', ]))
        result = s.cummin(skipna=False)
        self.assert_series_equal(expected, result)

    def test_cummax_timedelta64(self):
        s = pd.Series(pd.to_timedelta(['NaT',
                                       '2 min',
                                       'NaT',
                                       '1 min',
                                       'NaT',
                                       '3 min', ]))

        expected = pd.Series(pd.to_timedelta(['NaT',
                                              '2 min',
                                              'NaT',
                                              '2 min',
                                              'NaT',
                                              '3 min', ]))
        result = s.cummax(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_timedelta(['NaT',
                                              '2 min',
                                              '2 min',
                                              '2 min',
                                              '2 min',
                                              '3 min', ]))
        result = s.cummax(skipna=False)
        self.assert_series_equal(expected, result)

    def test_npdiff(self):
        pytest.skip("skipping due to Series no longer being an "
                    "ndarray")

        # no longer works as the return type of np.diff is now nd.array
        s = Series(np.arange(5))

        r = np.diff(s)
        assert_series_equal(Series([nan, 0, 0, 0, nan]), r)

    def _check_stat_op(self, name, alternate, check_objects=False,
                       check_allna=False):
        import pandas.core.nanops as nanops

        def testit():
            f = getattr(Series, name)

            # add some NaNs
            self.series[5:15] = np.NaN

            # idxmax, idxmin, min, and max are valid for dates
            if name not in ['max', 'min']:
                ds = Series(date_range('1/1/2001', periods=10))
                self.assertRaises(TypeError, f, ds)

            # skipna or no
            self.assertTrue(notnull(f(self.series)))
            self.assertTrue(isnull(f(self.series, skipna=False)))

            # check the result is correct
            nona = self.series.dropna()
            assert_almost_equal(f(nona), alternate(nona.values))
            assert_almost_equal(f(self.series), alternate(nona.values))

            allna = self.series * nan

            if check_allna:
                # xref 9422
                # bottleneck >= 1.0 give 0.0 for an allna Series sum
                try:
                    self.assertTrue(nanops._USE_BOTTLENECK)
                    import bottleneck as bn  # noqa
                    self.assertTrue(bn.__version__ >= LooseVersion('1.0'))
                    self.assertEqual(f(allna), 0.0)
                except:
                    self.assertTrue(np.isnan(f(allna)))

            # dtype=object with None, it works!
            s = Series([1, 2, 3, None, 5])
            f(s)

            # 2888
            l = [0]
            l.extend(lrange(2 ** 40, 2 ** 40 + 1000))
            s = Series(l, dtype='int64')
            assert_almost_equal(float(f(s)), float(alternate(s.values)))

            # check date range
            if check_objects:
                s = Series(bdate_range('1/1/2000', periods=10))
                res = f(s)
                exp = alternate(s)
                self.assertEqual(res, exp)

            # check on string data
            if name not in ['sum', 'min', 'max']:
                self.assertRaises(TypeError, f, Series(list('abc')))

            # Invalid axis.
            self.assertRaises(ValueError, f, self.series, axis=1)

            # Unimplemented numeric_only parameter.
            if 'numeric_only' in compat.signature(f).args:
                self.assertRaisesRegexp(NotImplementedError, name, f,
                                        self.series, numeric_only=True)

        testit()

        try:
            import bottleneck as bn  # noqa
            nanops._USE_BOTTLENECK = False
            testit()
            nanops._USE_BOTTLENECK = True
        except ImportError:
            pass

    def _check_accum_op(self, name, check_dtype=True):
        func = getattr(np, name)
        self.assert_numpy_array_equal(func(self.ts).values,
                                      func(np.array(self.ts)),
                                      check_dtype=check_dtype)

        # with missing values
        ts = self.ts.copy()
        ts[::2] = np.NaN

        result = func(ts)[1::2]
        expected = func(np.array(ts.valid()))

        self.assert_numpy_array_equal(result.values, expected,
                                      check_dtype=False)

    def test_compress(self):
        cond = [True, False, True, False, False]
        s = Series([1, -1, 5, 8, 7],
                   index=list('abcde'), name='foo')
        expected = Series(s.values.compress(cond),
                          index=list('ac'), name='foo')
        tm.assert_series_equal(s.compress(cond), expected)

    def test_numpy_compress(self):
        cond = [True, False, True, False, False]
        s = Series([1, -1, 5, 8, 7],
                   index=list('abcde'), name='foo')
        expected = Series(s.values.compress(cond),
                          index=list('ac'), name='foo')
        tm.assert_series_equal(np.compress(cond, s), expected)

        msg = "the 'axis' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.compress,
                              cond, s, axis=1)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.compress,
                              cond, s, out=s)

    def test_round(self):
        self.ts.index.name = "index_name"
        result = self.ts.round(2)
        expected = Series(np.round(self.ts.values, 2),
                          index=self.ts.index, name='ts')
        assert_series_equal(result, expected)
        self.assertEqual(result.name, self.ts.name)

    def test_numpy_round(self):
        # See gh-12600
        s = Series([1.53, 1.36, 0.06])
        out = np.round(s, decimals=0)
        expected = Series([2., 1., 0.])
        assert_series_equal(out, expected)

        msg = "the 'out' parameter is not supported"
        with tm.assertRaisesRegexp(ValueError, msg):
            np.round(s, decimals=0, out=s)

    def test_built_in_round(self):
        if not compat.PY3:
            pytest.skip(
                'build in round cannot be overriden prior to Python 3')

        s = Series([1.123, 2.123, 3.123], index=lrange(3))
        result = round(s)
        expected_rounded0 = Series([1., 2., 3.], index=lrange(3))
        self.assert_series_equal(result, expected_rounded0)

        decimals = 2
        expected_rounded = Series([1.12, 2.12, 3.12], index=lrange(3))
        result = round(s, decimals)
        self.assert_series_equal(result, expected_rounded)

    def test_prod_numpy16_bug(self):
        s = Series([1., 1., 1.], index=lrange(3))
        result = s.prod()
        self.assertNotIsInstance(result, Series)

    def test_all_any(self):
        ts = tm.makeTimeSeries()
        bool_series = ts > 0
        self.assertFalse(bool_series.all())
        self.assertTrue(bool_series.any())

        # Alternative types, with implicit 'object' dtype.
        s = Series(['abc', True])
        self.assertEqual('abc', s.any())  # 'abc' || True => 'abc'

    def test_all_any_params(self):
        # Check skipna, with implicit 'object' dtype.
        s1 = Series([np.nan, True])
        s2 = Series([np.nan, False])
        self.assertTrue(s1.all(skipna=False))  # nan && True => True
        self.assertTrue(s1.all(skipna=True))
        self.assertTrue(np.isnan(s2.any(skipna=False)))  # nan || False => nan
        self.assertFalse(s2.any(skipna=True))

        # Check level.
        s = pd.Series([False, False, True, True, False, True],
                      index=[0, 0, 1, 1, 2, 2])
        assert_series_equal(s.all(level=0), Series([False, True, False]))
        assert_series_equal(s.any(level=0), Series([False, True, True]))

        # bool_only is not implemented with level option.
        self.assertRaises(NotImplementedError, s.any, bool_only=True, level=0)
        self.assertRaises(NotImplementedError, s.all, bool_only=True, level=0)

        # bool_only is not implemented alone.
        self.assertRaises(NotImplementedError, s.any, bool_only=True)
        self.assertRaises(NotImplementedError, s.all, bool_only=True)

    def test_modulo(self):
        with np.errstate(all='ignore'):

            # GH3590, modulo as ints
            p = DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})
            result = p['first'] % p['second']
            expected = Series(p['first'].values % p['second'].values,
                              dtype='float64')
            expected.iloc[0:3] = np.nan
            assert_series_equal(result, expected)

            result = p['first'] % 0
            expected = Series(np.nan, index=p.index, name='first')
            assert_series_equal(result, expected)

            p = p.astype('float64')
            result = p['first'] % p['second']
            expected = Series(p['first'].values % p['second'].values)
            assert_series_equal(result, expected)

            p = p.astype('float64')
            result = p['first'] % p['second']
            result2 = p['second'] % p['first']
            self.assertFalse(np.array_equal(result, result2))

            # GH 9144
            s = Series([0, 1])

            result = s % 0
            expected = Series([nan, nan])
            assert_series_equal(result, expected)

            result = 0 % s
            expected = Series([nan, 0.0])
            assert_series_equal(result, expected)

    def test_ops_consistency_on_empty(self):

        # GH 7869
        # consistency on empty

        # float
        result = Series(dtype=float).sum()
        self.assertEqual(result, 0)

        result = Series(dtype=float).mean()
        self.assertTrue(isnull(result))

        result = Series(dtype=float).median()
        self.assertTrue(isnull(result))

        # timedelta64[ns]
        result = Series(dtype='m8[ns]').sum()
        self.assertEqual(result, Timedelta(0))

        result = Series(dtype='m8[ns]').mean()
        self.assertTrue(result is pd.NaT)

        result = Series(dtype='m8[ns]').median()
        self.assertTrue(result is pd.NaT)

    def test_corr(self):
        tm._skip_if_no_scipy()

        import scipy.stats as stats

        # full overlap
        self.assertAlmostEqual(self.ts.corr(self.ts), 1)

        # partial overlap
        self.assertAlmostEqual(self.ts[:15].corr(self.ts[5:]), 1)

        self.assertTrue(isnull(self.ts[:15].corr(self.ts[5:], min_periods=12)))

        ts1 = self.ts[:15].reindex(self.ts.index)
        ts2 = self.ts[5:].reindex(self.ts.index)
        self.assertTrue(isnull(ts1.corr(ts2, min_periods=12)))

        # No overlap
        self.assertTrue(np.isnan(self.ts[::2].corr(self.ts[1::2])))

        # all NA
        cp = self.ts[:10].copy()
        cp[:] = np.nan
        self.assertTrue(isnull(cp.corr(cp)))

        A = tm.makeTimeSeries()
        B = tm.makeTimeSeries()
        result = A.corr(B)
        expected, _ = stats.pearsonr(A, B)
        self.assertAlmostEqual(result, expected)

    def test_corr_rank(self):
        tm._skip_if_no_scipy()

        import scipy
        import scipy.stats as stats

        # kendall and spearman
        A = tm.makeTimeSeries()
        B = tm.makeTimeSeries()
        A[-5:] = A[:5]
        result = A.corr(B, method='kendall')
        expected = stats.kendalltau(A, B)[0]
        self.assertAlmostEqual(result, expected)

        result = A.corr(B, method='spearman')
        expected = stats.spearmanr(A, B)[0]
        self.assertAlmostEqual(result, expected)

        # these methods got rewritten in 0.8
        if scipy.__version__ < LooseVersion('0.9'):
            pytest.skip("skipping corr rank because of scipy version "
                        "{0}".format(scipy.__version__))

        # results from R
        A = Series(
            [-0.89926396, 0.94209606, -1.03289164, -0.95445587, 0.76910310, -
             0.06430576, -2.09704447, 0.40660407, -0.89926396, 0.94209606])
        B = Series(
            [-1.01270225, -0.62210117, -1.56895827, 0.59592943, -0.01680292,
             1.17258718, -1.06009347, -0.10222060, -0.89076239, 0.89372375])
        kexp = 0.4319297
        sexp = 0.5853767
        self.assertAlmostEqual(A.corr(B, method='kendall'), kexp)
        self.assertAlmostEqual(A.corr(B, method='spearman'), sexp)

    def test_cov(self):
        # full overlap
        self.assertAlmostEqual(self.ts.cov(self.ts), self.ts.std() ** 2)

        # partial overlap
        self.assertAlmostEqual(self.ts[:15].cov(self.ts[5:]),
                               self.ts[5:15].std() ** 2)

        # No overlap
        self.assertTrue(np.isnan(self.ts[::2].cov(self.ts[1::2])))

        # all NA
        cp = self.ts[:10].copy()
        cp[:] = np.nan
        self.assertTrue(isnull(cp.cov(cp)))

        # min_periods
        self.assertTrue(isnull(self.ts[:15].cov(self.ts[5:], min_periods=12)))

        ts1 = self.ts[:15].reindex(self.ts.index)
        ts2 = self.ts[5:].reindex(self.ts.index)
        self.assertTrue(isnull(ts1.cov(ts2, min_periods=12)))

    def test_count(self):
        self.assertEqual(self.ts.count(), len(self.ts))

        self.ts[::2] = np.NaN

        self.assertEqual(self.ts.count(), np.isfinite(self.ts).sum())

        mi = MultiIndex.from_arrays([list('aabbcc'), [1, 2, 2, nan, 1, 2]])
        ts = Series(np.arange(len(mi)), index=mi)

        left = ts.count(level=1)
        right = Series([2, 3, 1], index=[1, 2, nan])
        assert_series_equal(left, right)

        ts.iloc[[0, 3, 5]] = nan
        assert_series_equal(ts.count(level=1), right - 1)

    def test_dot(self):
        a = Series(np.random.randn(4), index=['p', 'q', 'r', 's'])
        b = DataFrame(np.random.randn(3, 4), index=['1', '2', '3'],
                      columns=['p', 'q', 'r', 's']).T

        result = a.dot(b)
        expected = Series(np.dot(a.values, b.values), index=['1', '2', '3'])
        assert_series_equal(result, expected)

        # Check index alignment
        b2 = b.reindex(index=reversed(b.index))
        result = a.dot(b)
        assert_series_equal(result, expected)

        # Check ndarray argument
        result = a.dot(b.values)
        self.assertTrue(np.all(result == expected.values))
        assert_almost_equal(a.dot(b['2'].values), expected['2'])

        # Check series argument
        assert_almost_equal(a.dot(b['1']), expected['1'])
        assert_almost_equal(a.dot(b2['1']), expected['1'])

        self.assertRaises(Exception, a.dot, a.values[:3])
        self.assertRaises(ValueError, a.dot, b.T)

    def test_value_counts_nunique(self):

        # basics.rst doc example
        series = Series(np.random.randn(500))
        series[20:500] = np.nan
        series[10:20] = 5000
        result = series.nunique()
        self.assertEqual(result, 11)

    def test_unique(self):

        # 714 also, dtype=float
        s = Series([1.2345] * 100)
        s[::2] = np.nan
        result = s.unique()
        self.assertEqual(len(result), 2)

        s = Series([1.2345] * 100, dtype='f4')
        s[::2] = np.nan
        result = s.unique()
        self.assertEqual(len(result), 2)

        # NAs in object arrays #714
        s = Series(['foo'] * 100, dtype='O')
        s[::2] = np.nan
        result = s.unique()
        self.assertEqual(len(result), 2)

        # decision about None
        s = Series([1, 2, 3, None, None, None], dtype=object)
        result = s.unique()
        expected = np.array([1, 2, 3, None], dtype=object)
        self.assert_numpy_array_equal(result, expected)

    def test_drop_duplicates(self):
        # check both int and object
        for s in [Series([1, 2, 3, 3]), Series(['1', '2', '3', '3'])]:
            expected = Series([False, False, False, True])
            assert_series_equal(s.duplicated(), expected)
            assert_series_equal(s.drop_duplicates(), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, False, True, False])
            assert_series_equal(s.duplicated(keep='last'), expected)
            assert_series_equal(s.drop_duplicates(keep='last'), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep='last', inplace=True)
            assert_series_equal(sc, s[~expected])

            # deprecate take_last
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.duplicated(take_last=True), expected)
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(
                    s.drop_duplicates(take_last=True), s[~expected])
            sc = s.copy()
            with tm.assert_produces_warning(FutureWarning):
                sc.drop_duplicates(take_last=True, inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, False, True, True])
            assert_series_equal(s.duplicated(keep=False), expected)
            assert_series_equal(s.drop_duplicates(keep=False), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep=False, inplace=True)
            assert_series_equal(sc, s[~expected])

        for s in [Series([1, 2, 3, 5, 3, 2, 4]),
                  Series(['1', '2', '3', '5', '3', '2', '4'])]:
            expected = Series([False, False, False, False, True, True, False])
            assert_series_equal(s.duplicated(), expected)
            assert_series_equal(s.drop_duplicates(), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, True, True, False, False, False, False])
            assert_series_equal(s.duplicated(keep='last'), expected)
            assert_series_equal(s.drop_duplicates(keep='last'), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep='last', inplace=True)
            assert_series_equal(sc, s[~expected])

            # deprecate take_last
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.duplicated(take_last=True), expected)
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(
                    s.drop_duplicates(take_last=True), s[~expected])
            sc = s.copy()
            with tm.assert_produces_warning(FutureWarning):
                sc.drop_duplicates(take_last=True, inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, True, True, False, True, True, False])
            assert_series_equal(s.duplicated(keep=False), expected)
            assert_series_equal(s.drop_duplicates(keep=False), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep=False, inplace=True)
            assert_series_equal(sc, s[~expected])

    def test_rank(self):
        tm._skip_if_no_scipy()
        from scipy.stats import rankdata

        self.ts[::2] = np.nan
        self.ts[:10][::3] = 4.

        ranks = self.ts.rank()
        oranks = self.ts.astype('O').rank()

        assert_series_equal(ranks, oranks)

        mask = np.isnan(self.ts)
        filled = self.ts.fillna(np.inf)

        # rankdata returns a ndarray
        exp = Series(rankdata(filled), index=filled.index, name='ts')
        exp[mask] = np.nan

        tm.assert_series_equal(ranks, exp)

        iseries = Series(np.arange(5).repeat(2))

        iranks = iseries.rank()
        exp = iseries.astype(float).rank()
        assert_series_equal(iranks, exp)
        iseries = Series(np.arange(5)) + 1.0
        exp = iseries / 5.0
        iranks = iseries.rank(pct=True)

        assert_series_equal(iranks, exp)

        iseries = Series(np.repeat(1, 100))
        exp = Series(np.repeat(0.505, 100))
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries[1] = np.nan
        exp = Series(np.repeat(50.0 / 99.0, 100))
        exp[1] = np.nan
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.arange(5)) + 1.0
        iseries[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.repeat(np.nan, 100))
        exp = iseries.copy()
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.arange(5)) + 1
        iseries[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        rng = date_range('1/1/1990', periods=5)
        iseries = Series(np.arange(5), rng) + 1
        iseries.iloc[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series([1e-50, 1e-100, 1e-20, 1e-2, 1e-20 + 1e-30, 1e-1])
        exp = Series([2, 1, 3, 5, 4, 6.0])
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

        # GH 5968
        iseries = Series(['3 day', '1 day 10m', '-2 day', pd.NaT],
                         dtype='m8[ns]')
        exp = Series([3, 2, 1, np.nan])
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

        values = np.array(
            [-50, -1, -1e-20, -1e-25, -1e-50, 0, 1e-40, 1e-20, 1e-10, 2, 40
             ], dtype='float64')
        random_order = np.random.permutation(len(values))
        iseries = Series(values[random_order])
        exp = Series(random_order + 1.0, dtype='float64')
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

    def test_rank_categorical(self):
        # GH issue #15420 rank incorrectly orders ordered categories

        # Test ascending/descending ranking for ordered categoricals
        exp = pd.Series([1., 2., 3., 4., 5., 6.])
        exp_desc = pd.Series([6., 5., 4., 3., 2., 1.])
        ordered = pd.Series(
            ['first', 'second', 'third', 'fourth', 'fifth', 'sixth']
        ).astype('category', ).cat.set_categories(
            ['first', 'second', 'third', 'fourth', 'fifth', 'sixth'],
            ordered=True
        )
        assert_series_equal(ordered.rank(), exp)
        assert_series_equal(ordered.rank(ascending=False), exp_desc)

        # Unordered categoricals should be ranked as objects
        unordered = pd.Series(
            ['first', 'second', 'third', 'fourth', 'fifth', 'sixth'],
        ).astype('category').cat.set_categories(
            ['first', 'second', 'third', 'fourth', 'fifth', 'sixth'],
            ordered=False
        )
        exp_unordered = pd.Series([2., 4., 6., 3., 1., 5.])
        res = unordered.rank()
        assert_series_equal(res, exp_unordered)

        # Test na_option for rank data
        na_ser = pd.Series(
            ['first', 'second', 'third', 'fourth', 'fifth', 'sixth', np.NaN]
        ).astype('category', ).cat.set_categories(
            [
                'first', 'second', 'third', 'fourth',
                'fifth', 'sixth', 'seventh'
            ],
            ordered=True
        )

        exp_top = pd.Series([2., 3., 4., 5., 6., 7., 1.])
        exp_bot = pd.Series([1., 2., 3., 4., 5., 6., 7.])
        exp_keep = pd.Series([1., 2., 3., 4., 5., 6., np.NaN])

        assert_series_equal(na_ser.rank(na_option='top'), exp_top)
        assert_series_equal(na_ser.rank(na_option='bottom'), exp_bot)
        assert_series_equal(na_ser.rank(na_option='keep'), exp_keep)

        # Test na_option for rank data with ascending False
        exp_top = pd.Series([7., 6., 5., 4., 3., 2., 1.])
        exp_bot = pd.Series([6., 5., 4., 3., 2., 1., 7.])
        exp_keep = pd.Series([6., 5., 4., 3., 2., 1., np.NaN])

        assert_series_equal(
            na_ser.rank(na_option='top', ascending=False),
            exp_top
        )
        assert_series_equal(
            na_ser.rank(na_option='bottom', ascending=False),
            exp_bot
        )
        assert_series_equal(
            na_ser.rank(na_option='keep', ascending=False),
            exp_keep
        )

        # Test with pct=True
        na_ser = pd.Series(
            ['first', 'second', 'third', 'fourth', np.NaN],
        ).astype('category').cat.set_categories(
            ['first', 'second', 'third', 'fourth'],
            ordered=True
        )
        exp_top = pd.Series([0.4, 0.6, 0.8, 1., 0.2])
        exp_bot = pd.Series([0.2, 0.4, 0.6, 0.8, 1.])
        exp_keep = pd.Series([0.25, 0.5, 0.75, 1., np.NaN])

        assert_series_equal(na_ser.rank(na_option='top', pct=True), exp_top)
        assert_series_equal(na_ser.rank(na_option='bottom', pct=True), exp_bot)
        assert_series_equal(na_ser.rank(na_option='keep', pct=True), exp_keep)

    def test_rank_signature(self):
        s = Series([0, 1])
        s.rank(method='average')
        self.assertRaises(ValueError, s.rank, 'average')

    def test_rank_inf(self):
        pytest.skip('DataFrame.rank does not currently rank '
                    'np.inf and -np.inf properly')

        values = np.array(
            [-np.inf, -50, -1, -1e-20, -1e-25, -1e-50, 0, 1e-40, 1e-20, 1e-10,
             2, 40, np.inf], dtype='float64')
        random_order = np.random.permutation(len(values))
        iseries = Series(values[random_order])
        exp = Series(random_order + 1.0, dtype='float64')
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

    def test_clip(self):
        val = self.ts.median()

        self.assertEqual(self.ts.clip_lower(val).min(), val)
        self.assertEqual(self.ts.clip_upper(val).max(), val)

        self.assertEqual(self.ts.clip(lower=val).min(), val)
        self.assertEqual(self.ts.clip(upper=val).max(), val)

        result = self.ts.clip(-0.5, 0.5)
        expected = np.clip(self.ts, -0.5, 0.5)
        assert_series_equal(result, expected)
        tm.assertIsInstance(expected, Series)

    def test_clip_types_and_nulls(self):

        sers = [Series([np.nan, 1.0, 2.0, 3.0]), Series([None, 'a', 'b', 'c']),
                Series(pd.to_datetime(
                    [np.nan, 1, 2, 3], unit='D'))]

        for s in sers:
            thresh = s[2]
            l = s.clip_lower(thresh)
            u = s.clip_upper(thresh)
            self.assertEqual(l[notnull(l)].min(), thresh)
            self.assertEqual(u[notnull(u)].max(), thresh)
            self.assertEqual(list(isnull(s)), list(isnull(l)))
            self.assertEqual(list(isnull(s)), list(isnull(u)))

    def test_clip_against_series(self):
        # GH #6966

        s = Series([1.0, 1.0, 4.0])
        threshold = Series([1.0, 2.0, 3.0])

        assert_series_equal(s.clip_lower(threshold), Series([1.0, 2.0, 4.0]))
        assert_series_equal(s.clip_upper(threshold), Series([1.0, 1.0, 3.0]))

        lower = Series([1.0, 2.0, 3.0])
        upper = Series([1.5, 2.5, 3.5])
        assert_series_equal(s.clip(lower, upper), Series([1.0, 2.0, 3.5]))
        assert_series_equal(s.clip(1.5, upper), Series([1.5, 1.5, 3.5]))

    def test_clip_with_datetimes(self):

        # GH 11838
        # naive and tz-aware datetimes

        t = Timestamp('2015-12-01 09:30:30')
        s = Series([Timestamp('2015-12-01 09:30:00'), Timestamp(
            '2015-12-01 09:31:00')])
        result = s.clip(upper=t)
        expected = Series([Timestamp('2015-12-01 09:30:00'), Timestamp(
            '2015-12-01 09:30:30')])
        assert_series_equal(result, expected)

        t = Timestamp('2015-12-01 09:30:30', tz='US/Eastern')
        s = Series([Timestamp('2015-12-01 09:30:00', tz='US/Eastern'),
                    Timestamp('2015-12-01 09:31:00', tz='US/Eastern')])
        result = s.clip(upper=t)
        expected = Series([Timestamp('2015-12-01 09:30:00', tz='US/Eastern'),
                           Timestamp('2015-12-01 09:30:30', tz='US/Eastern')])
        assert_series_equal(result, expected)

    def test_cummethods_bool(self):
        # GH 6270
        # looks like a buggy np.maximum.accumulate for numpy 1.6.1, py 3.2
        def cummin(x):
            return np.minimum.accumulate(x)

        def cummax(x):
            return np.maximum.accumulate(x)

        a = pd.Series([False, False, False, True, True, False, False])
        b = ~a
        c = pd.Series([False] * len(b))
        d = ~c
        methods = {'cumsum': np.cumsum,
                   'cumprod': np.cumprod,
                   'cummin': cummin,
                   'cummax': cummax}
        args = product((a, b, c, d), methods)
        for s, method in args:
            expected = Series(methods[method](s.values))
            result = getattr(s, method)()
            assert_series_equal(result, expected)

        e = pd.Series([False, True, nan, False])
        cse = pd.Series([0, 1, nan, 1], dtype=object)
        cpe = pd.Series([False, 0, nan, 0])
        cmin = pd.Series([False, False, nan, False])
        cmax = pd.Series([False, True, nan, True])
        expecteds = {'cumsum': cse,
                     'cumprod': cpe,
                     'cummin': cmin,
                     'cummax': cmax}

        for method in methods:
            res = getattr(e, method)()
            assert_series_equal(res, expecteds[method])

    def test_isin(self):
        s = Series(['A', 'B', 'C', 'a', 'B', 'B', 'A', 'C'])

        result = s.isin(['A', 'C'])
        expected = Series([True, False, True, False, False, False, True, True])
        assert_series_equal(result, expected)

    def test_isin_with_string_scalar(self):
        # GH4763
        s = Series(['A', 'B', 'C', 'a', 'B', 'B', 'A', 'C'])
        with tm.assertRaises(TypeError):
            s.isin('a')

        with tm.assertRaises(TypeError):
            s = Series(['aaa', 'b', 'c'])
            s.isin('aaa')

    def test_isin_with_i8(self):
        # GH 5021

        expected = Series([True, True, False, False, False])
        expected2 = Series([False, True, False, False, False])

        # datetime64[ns]
        s = Series(date_range('jan-01-2013', 'jan-05-2013'))

        result = s.isin(s[0:2])
        assert_series_equal(result, expected)

        result = s.isin(s[0:2].values)
        assert_series_equal(result, expected)

        # fails on dtype conversion in the first place
        result = s.isin(s[0:2].values.astype('datetime64[D]'))
        assert_series_equal(result, expected)

        result = s.isin([s[1]])
        assert_series_equal(result, expected2)

        result = s.isin([np.datetime64(s[1])])
        assert_series_equal(result, expected2)

        result = s.isin(set(s[0:2]))
        assert_series_equal(result, expected)

        # timedelta64[ns]
        s = Series(pd.to_timedelta(lrange(5), unit='d'))
        result = s.isin(s[0:2])
        assert_series_equal(result, expected)

    def test_timedelta64_analytics(self):
        from pandas import date_range

        # index min/max
        td = Series(date_range('2012-1-1', periods=3, freq='D')) - \
            Timestamp('20120101')

        result = td.idxmin()
        self.assertEqual(result, 0)

        result = td.idxmax()
        self.assertEqual(result, 2)

        # GH 2982
        # with NaT
        td[0] = np.nan

        result = td.idxmin()
        self.assertEqual(result, 1)

        result = td.idxmax()
        self.assertEqual(result, 2)

        # abs
        s1 = Series(date_range('20120101', periods=3))
        s2 = Series(date_range('20120102', periods=3))
        expected = Series(s2 - s1)

        # this fails as numpy returns timedelta64[us]
        # result = np.abs(s1-s2)
        # assert_frame_equal(result,expected)

        result = (s1 - s2).abs()
        assert_series_equal(result, expected)

        # max/min
        result = td.max()
        expected = Timedelta('2 days')
        self.assertEqual(result, expected)

        result = td.min()
        expected = Timedelta('1 days')
        self.assertEqual(result, expected)

    def test_idxmin(self):
        # test idxmin
        # _check_stat_op approach can not be used here because of isnull check.

        # add some NaNs
        self.series[5:15] = np.NaN

        # skipna or no
        self.assertEqual(self.series[self.series.idxmin()], self.series.min())
        self.assertTrue(isnull(self.series.idxmin(skipna=False)))

        # no NaNs
        nona = self.series.dropna()
        self.assertEqual(nona[nona.idxmin()], nona.min())
        self.assertEqual(nona.index.values.tolist().index(nona.idxmin()),
                         nona.values.argmin())

        # all NaNs
        allna = self.series * nan
        self.assertTrue(isnull(allna.idxmin()))

        # datetime64[ns]
        from pandas import date_range
        s = Series(date_range('20130102', periods=6))
        result = s.idxmin()
        self.assertEqual(result, 0)

        s[0] = np.nan
        result = s.idxmin()
        self.assertEqual(result, 1)

    def test_numpy_argmin(self):
        # argmin is aliased to idxmin
        data = np.random.randint(0, 11, size=10)
        result = np.argmin(Series(data))
        self.assertEqual(result, np.argmin(data))

        if not _np_version_under1p10:
            msg = "the 'out' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, msg, np.argmin,
                                  Series(data), out=data)

    def test_idxmax(self):
        # test idxmax
        # _check_stat_op approach can not be used here because of isnull check.

        # add some NaNs
        self.series[5:15] = np.NaN

        # skipna or no
        self.assertEqual(self.series[self.series.idxmax()], self.series.max())
        self.assertTrue(isnull(self.series.idxmax(skipna=False)))

        # no NaNs
        nona = self.series.dropna()
        self.assertEqual(nona[nona.idxmax()], nona.max())
        self.assertEqual(nona.index.values.tolist().index(nona.idxmax()),
                         nona.values.argmax())

        # all NaNs
        allna = self.series * nan
        self.assertTrue(isnull(allna.idxmax()))

        from pandas import date_range
        s = Series(date_range('20130102', periods=6))
        result = s.idxmax()
        self.assertEqual(result, 5)

        s[5] = np.nan
        result = s.idxmax()
        self.assertEqual(result, 4)

        # Float64Index
        # GH 5914
        s = pd.Series([1, 2, 3], [1.1, 2.1, 3.1])
        result = s.idxmax()
        self.assertEqual(result, 3.1)
        result = s.idxmin()
        self.assertEqual(result, 1.1)

        s = pd.Series(s.index, s.index)
        result = s.idxmax()
        self.assertEqual(result, 3.1)
        result = s.idxmin()
        self.assertEqual(result, 1.1)

    def test_numpy_argmax(self):

        # argmax is aliased to idxmax
        data = np.random.randint(0, 11, size=10)
        result = np.argmax(Series(data))
        self.assertEqual(result, np.argmax(data))

        if not _np_version_under1p10:
            msg = "the 'out' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, msg, np.argmax,
                                  Series(data), out=data)

    def test_ptp(self):
        N = 1000
        arr = np.random.randn(N)
        ser = Series(arr)
        self.assertEqual(np.ptp(ser), np.ptp(arr))

        # GH11163
        s = Series([3, 5, np.nan, -3, 10])
        self.assertEqual(s.ptp(), 13)
        self.assertTrue(pd.isnull(s.ptp(skipna=False)))

        mi = pd.MultiIndex.from_product([['a', 'b'], [1, 2, 3]])
        s = pd.Series([1, np.nan, 7, 3, 5, np.nan], index=mi)

        expected = pd.Series([6, 2], index=['a', 'b'], dtype=np.float64)
        self.assert_series_equal(s.ptp(level=0), expected)

        expected = pd.Series([np.nan, np.nan], index=['a', 'b'])
        self.assert_series_equal(s.ptp(level=0, skipna=False), expected)

        with self.assertRaises(ValueError):
            s.ptp(axis=1)

        s = pd.Series(['a', 'b', 'c', 'd', 'e'])
        with self.assertRaises(TypeError):
            s.ptp()

        with self.assertRaises(NotImplementedError):
            s.ptp(numeric_only=True)

    def test_empty_timeseries_redections_return_nat(self):
        # covers #11245
        for dtype in ('m8[ns]', 'm8[ns]', 'M8[ns]', 'M8[ns, UTC]'):
            self.assertIs(Series([], dtype=dtype).min(), pd.NaT)
            self.assertIs(Series([], dtype=dtype).max(), pd.NaT)

    def test_unique_data_ownership(self):
        # it works! #1807
        Series(Series(["a", "c", "b"]).unique()).sort_values()

    def test_repeat(self):
        s = Series(np.random.randn(3), index=['a', 'b', 'c'])

        reps = s.repeat(5)
        exp = Series(s.values.repeat(5), index=s.index.values.repeat(5))
        assert_series_equal(reps, exp)

        with tm.assert_produces_warning(FutureWarning):
            result = s.repeat(reps=5)
            assert_series_equal(result, exp)

        to_rep = [2, 3, 4]
        reps = s.repeat(to_rep)
        exp = Series(s.values.repeat(to_rep),
                     index=s.index.values.repeat(to_rep))
        assert_series_equal(reps, exp)

    def test_numpy_repeat(self):
        s = Series(np.arange(3), name='x')
        expected = Series(s.values.repeat(2), name='x',
                          index=s.index.values.repeat(2))
        assert_series_equal(np.repeat(s, 2), expected)

        msg = "the 'axis' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.repeat, s, 2, axis=0)

    def test_searchsorted(self):
        s = Series([1, 2, 3])

        idx = s.searchsorted(1, side='left')
        tm.assert_numpy_array_equal(idx, np.array([0], dtype=np.intp))

        idx = s.searchsorted(1, side='right')
        tm.assert_numpy_array_equal(idx, np.array([1], dtype=np.intp))

        with tm.assert_produces_warning(FutureWarning):
            idx = s.searchsorted(v=1, side='left')
            tm.assert_numpy_array_equal(idx, np.array([0], dtype=np.intp))

    def test_searchsorted_numeric_dtypes_scalar(self):
        s = Series([1, 2, 90, 1000, 3e9])
        r = s.searchsorted(30)
        e = 2
        self.assertEqual(r, e)

        r = s.searchsorted([30])
        e = np.array([2], dtype=np.intp)
        tm.assert_numpy_array_equal(r, e)

    def test_searchsorted_numeric_dtypes_vector(self):
        s = Series([1, 2, 90, 1000, 3e9])
        r = s.searchsorted([91, 2e6])
        e = np.array([3, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(r, e)

    def test_search_sorted_datetime64_scalar(self):
        s = Series(pd.date_range('20120101', periods=10, freq='2D'))
        v = pd.Timestamp('20120102')
        r = s.searchsorted(v)
        e = 1
        self.assertEqual(r, e)

    def test_search_sorted_datetime64_list(self):
        s = Series(pd.date_range('20120101', periods=10, freq='2D'))
        v = [pd.Timestamp('20120102'), pd.Timestamp('20120104')]
        r = s.searchsorted(v)
        e = np.array([1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(r, e)

    def test_searchsorted_sorter(self):
        # GH8490
        s = Series([3, 1, 2])
        r = s.searchsorted([0, 3], sorter=np.argsort(s))
        e = np.array([0, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(r, e)

    def test_is_unique(self):
        # GH11946
        s = Series(np.random.randint(0, 10, size=1000))
        self.assertFalse(s.is_unique)
        s = Series(np.arange(1000))
        self.assertTrue(s.is_unique)

    def test_is_monotonic(self):

        s = Series(np.random.randint(0, 10, size=1000))
        self.assertFalse(s.is_monotonic)
        s = Series(np.arange(1000))
        self.assertTrue(s.is_monotonic)
        self.assertTrue(s.is_monotonic_increasing)
        s = Series(np.arange(1000, 0, -1))
        self.assertTrue(s.is_monotonic_decreasing)

        s = Series(pd.date_range('20130101', periods=10))
        self.assertTrue(s.is_monotonic)
        self.assertTrue(s.is_monotonic_increasing)
        s = Series(list(reversed(s.tolist())))
        self.assertFalse(s.is_monotonic)
        self.assertTrue(s.is_monotonic_decreasing)

    def test_nsmallest_nlargest(self):
        # float, int, datetime64 (use i8), timedelts64 (same),
        # object that are numbers, object that are strings

        base = [3, 2, 1, 2, 5]

        s_list = [
            Series(base, dtype='int8'),
            Series(base, dtype='int16'),
            Series(base, dtype='int32'),
            Series(base, dtype='int64'),
            Series(base, dtype='float32'),
            Series(base, dtype='float64'),
            Series(base, dtype='uint8'),
            Series(base, dtype='uint16'),
            Series(base, dtype='uint32'),
            Series(base, dtype='uint64'),
            Series(base).astype('timedelta64[ns]'),
            Series(pd.to_datetime(['2003', '2002', '2001', '2002', '2005'])),
        ]

        raising = [
            Series([3., 2, 1, 2, '5'], dtype='object'),
            Series([3., 2, 1, 2, 5], dtype='object'),
            # not supported on some archs
            # Series([3., 2, 1, 2, 5], dtype='complex256'),
            Series([3., 2, 1, 2, 5], dtype='complex128'),
        ]

        for r in raising:
            dt = r.dtype
            msg = "Cannot use method 'n(larg|small)est' with dtype %s" % dt
            args = 2, len(r), 0, -1
            methods = r.nlargest, r.nsmallest
            for method, arg in product(methods, args):
                with tm.assertRaisesRegexp(TypeError, msg):
                    method(arg)

        for s in s_list:

            assert_series_equal(s.nsmallest(2), s.iloc[[2, 1]])

            assert_series_equal(s.nsmallest(2, keep='last'), s.iloc[[2, 3]])
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(
                    s.nsmallest(2, take_last=True), s.iloc[[2, 3]])

            assert_series_equal(s.nlargest(3), s.iloc[[4, 0, 1]])

            assert_series_equal(s.nlargest(3, keep='last'), s.iloc[[4, 0, 3]])
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(
                    s.nlargest(3, take_last=True), s.iloc[[4, 0, 3]])

            empty = s.iloc[0:0]
            assert_series_equal(s.nsmallest(0), empty)
            assert_series_equal(s.nsmallest(-1), empty)
            assert_series_equal(s.nlargest(0), empty)
            assert_series_equal(s.nlargest(-1), empty)

            assert_series_equal(s.nsmallest(len(s)), s.sort_values())
            assert_series_equal(s.nsmallest(len(s) + 1), s.sort_values())
            assert_series_equal(s.nlargest(len(s)), s.iloc[[4, 0, 1, 3, 2]])
            assert_series_equal(s.nlargest(len(s) + 1),
                                s.iloc[[4, 0, 1, 3, 2]])

        s = Series([3., np.nan, 1, 2, 5])
        assert_series_equal(s.nlargest(), s.iloc[[4, 0, 3, 2]])
        assert_series_equal(s.nsmallest(), s.iloc[[2, 3, 0, 4]])

        msg = 'keep must be either "first", "last"'
        with tm.assertRaisesRegexp(ValueError, msg):
            s.nsmallest(keep='invalid')
        with tm.assertRaisesRegexp(ValueError, msg):
            s.nlargest(keep='invalid')

        # GH 13412
        s = Series([1, 4, 3, 2], index=[0, 0, 1, 1])
        result = s.nlargest(3)
        expected = s.sort_values(ascending=False).head(3)
        assert_series_equal(result, expected)
        result = s.nsmallest(3)
        expected = s.sort_values().head(3)
        assert_series_equal(result, expected)

    def test_sort_index_level(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        s = Series([1, 2], mi)
        backwards = s.iloc[[1, 0]]

        res = s.sort_index(level='A')
        assert_series_equal(backwards, res)

        res = s.sort_index(level=['A', 'B'])
        assert_series_equal(backwards, res)

        res = s.sort_index(level='A', sort_remaining=False)
        assert_series_equal(s, res)

        res = s.sort_index(level=['A', 'B'], sort_remaining=False)
        assert_series_equal(s, res)

    def test_apply_categorical(self):
        values = pd.Categorical(list('ABBABCD'), categories=list('DCBA'),
                                ordered=True)
        s = pd.Series(values, name='XX', index=list('abcdefg'))
        result = s.apply(lambda x: x.lower())

        # should be categorical dtype when the number of categories are
        # the same
        values = pd.Categorical(list('abbabcd'), categories=list('dcba'),
                                ordered=True)
        exp = pd.Series(values, name='XX', index=list('abcdefg'))
        tm.assert_series_equal(result, exp)
        tm.assert_categorical_equal(result.values, exp.values)

        result = s.apply(lambda x: 'A')
        exp = pd.Series(['A'] * 7, name='XX', index=list('abcdefg'))
        tm.assert_series_equal(result, exp)
        self.assertEqual(result.dtype, np.object)

    def test_shift_int(self):
        ts = self.ts.astype(int)
        shifted = ts.shift(1)
        expected = ts.astype(float).shift(1)
        assert_series_equal(shifted, expected)

    def test_shift_categorical(self):
        # GH 9416
        s = pd.Series(['a', 'b', 'c', 'd'], dtype='category')

        assert_series_equal(s.iloc[:-1], s.shift(1).shift(-1).valid())

        sp1 = s.shift(1)
        assert_index_equal(s.index, sp1.index)
        self.assertTrue(np.all(sp1.values.codes[:1] == -1))
        self.assertTrue(np.all(s.values.codes[:-1] == sp1.values.codes[1:]))

        sn2 = s.shift(-2)
        assert_index_equal(s.index, sn2.index)
        self.assertTrue(np.all(sn2.values.codes[-2:] == -1))
        self.assertTrue(np.all(s.values.codes[2:] == sn2.values.codes[:-2]))

        assert_index_equal(s.values.categories, sp1.values.categories)
        assert_index_equal(s.values.categories, sn2.values.categories)

    def test_reshape_deprecate(self):
        x = Series(np.random.random(10), name='x')
        tm.assert_produces_warning(FutureWarning, x.reshape, x.shape)

    def test_reshape_non_2d(self):
        # see gh-4554
        with tm.assert_produces_warning(FutureWarning):
            x = Series(np.random.random(201), name='x')
            self.assertTrue(x.reshape(x.shape, ) is x)

        # see gh-2719
        with tm.assert_produces_warning(FutureWarning):
            a = Series([1, 2, 3, 4])
            result = a.reshape(2, 2)
            expected = a.values.reshape(2, 2)
            tm.assert_numpy_array_equal(result, expected)
            self.assertIsInstance(result, type(expected))

    def test_reshape_2d_return_array(self):
        x = Series(np.random.random(201), name='x')

        with tm.assert_produces_warning(FutureWarning):
            result = x.reshape((-1, 1))
            self.assertNotIsInstance(result, Series)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result2 = np.reshape(x, (-1, 1))
            self.assertNotIsInstance(result2, Series)

        with tm.assert_produces_warning(FutureWarning):
            result = x[:, None]
            expected = x.reshape((-1, 1))
            assert_almost_equal(result, expected)

    def test_reshape_bad_kwarg(self):
        a = Series([1, 2, 3, 4])

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            msg = "'foo' is an invalid keyword argument for this function"
            tm.assertRaisesRegexp(TypeError, msg, a.reshape, (2, 2), foo=2)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            msg = r"reshape\(\) got an unexpected keyword argument 'foo'"
            tm.assertRaisesRegexp(TypeError, msg, a.reshape, a.shape, foo=2)

    def test_numpy_reshape(self):
        a = Series([1, 2, 3, 4])

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = np.reshape(a, (2, 2))
            expected = a.values.reshape(2, 2)
            tm.assert_numpy_array_equal(result, expected)
            self.assertIsInstance(result, type(expected))

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = np.reshape(a, a.shape)
            tm.assert_series_equal(result, a)

    def test_unstack(self):
        from numpy import nan

        index = MultiIndex(levels=[['bar', 'foo'], ['one', 'three', 'two']],
                           labels=[[1, 1, 0, 0], [0, 1, 0, 2]])

        s = Series(np.arange(4.), index=index)
        unstacked = s.unstack()

        expected = DataFrame([[2., nan, 3.], [0., 1., nan]],
                             index=['bar', 'foo'],
                             columns=['one', 'three', 'two'])

        assert_frame_equal(unstacked, expected)

        unstacked = s.unstack(level=0)
        assert_frame_equal(unstacked, expected.T)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        exp_index = MultiIndex(levels=[['one', 'two', 'three'], [0, 1]],
                               labels=[[0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1]])
        expected = DataFrame({'bar': s.values},
                             index=exp_index).sort_index(level=0)
        unstacked = s.unstack(0)
        assert_frame_equal(unstacked, expected)

        # GH5873
        idx = pd.MultiIndex.from_arrays([[101, 102], [3.5, np.nan]])
        ts = pd.Series([1, 2], index=idx)
        left = ts.unstack()
        right = DataFrame([[nan, 1], [2, nan]], index=[101, 102],
                          columns=[nan, 3.5])
        assert_frame_equal(left, right)

        idx = pd.MultiIndex.from_arrays([['cat', 'cat', 'cat', 'dog', 'dog'
                                          ], ['a', 'a', 'b', 'a', 'b'],
                                         [1, 2, 1, 1, np.nan]])
        ts = pd.Series([1.0, 1.1, 1.2, 1.3, 1.4], index=idx)
        right = DataFrame([[1.0, 1.3], [1.1, nan], [nan, 1.4], [1.2, nan]],
                          columns=['cat', 'dog'])
        tpls = [('a', 1), ('a', 2), ('b', nan), ('b', 1)]
        right.index = pd.MultiIndex.from_tuples(tpls)
        assert_frame_equal(ts.unstack(level=0), right)

    def test_value_counts_datetime(self):
        # most dtypes are tested in test_base.py
        values = [pd.Timestamp('2011-01-01 09:00'),
                  pd.Timestamp('2011-01-01 10:00'),
                  pd.Timestamp('2011-01-01 11:00'),
                  pd.Timestamp('2011-01-01 09:00'),
                  pd.Timestamp('2011-01-01 09:00'),
                  pd.Timestamp('2011-01-01 11:00')]

        exp_idx = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 11:00',
                                    '2011-01-01 10:00'])
        exp = pd.Series([3, 2, 1], index=exp_idx, name='xxx')

        s = pd.Series(values, name='xxx')
        tm.assert_series_equal(s.value_counts(), exp)
        # check DatetimeIndex outputs the same result
        idx = pd.DatetimeIndex(values, name='xxx')
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3., 2., 1]) / 6.,
                        index=exp_idx, name='xxx')
        tm.assert_series_equal(s.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_datetime_tz(self):
        values = [pd.Timestamp('2011-01-01 09:00', tz='US/Eastern'),
                  pd.Timestamp('2011-01-01 10:00', tz='US/Eastern'),
                  pd.Timestamp('2011-01-01 11:00', tz='US/Eastern'),
                  pd.Timestamp('2011-01-01 09:00', tz='US/Eastern'),
                  pd.Timestamp('2011-01-01 09:00', tz='US/Eastern'),
                  pd.Timestamp('2011-01-01 11:00', tz='US/Eastern')]

        exp_idx = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 11:00',
                                    '2011-01-01 10:00'], tz='US/Eastern')
        exp = pd.Series([3, 2, 1], index=exp_idx, name='xxx')

        s = pd.Series(values, name='xxx')
        tm.assert_series_equal(s.value_counts(), exp)
        idx = pd.DatetimeIndex(values, name='xxx')
        tm.assert_series_equal(idx.value_counts(), exp)

        exp = pd.Series(np.array([3., 2., 1]) / 6.,
                        index=exp_idx, name='xxx')
        tm.assert_series_equal(s.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_period(self):
        values = [pd.Period('2011-01', freq='M'),
                  pd.Period('2011-02', freq='M'),
                  pd.Period('2011-03', freq='M'),
                  pd.Period('2011-01', freq='M'),
                  pd.Period('2011-01', freq='M'),
                  pd.Period('2011-03', freq='M')]

        exp_idx = pd.PeriodIndex(['2011-01', '2011-03', '2011-02'], freq='M')
        exp = pd.Series([3, 2, 1], index=exp_idx, name='xxx')

        s = pd.Series(values, name='xxx')
        tm.assert_series_equal(s.value_counts(), exp)
        # check DatetimeIndex outputs the same result
        idx = pd.PeriodIndex(values, name='xxx')
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3., 2., 1]) / 6.,
                        index=exp_idx, name='xxx')
        tm.assert_series_equal(s.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical_ordered(self):
        # most dtypes are tested in test_base.py
        values = pd.Categorical([1, 2, 3, 1, 1, 3], ordered=True)

        exp_idx = pd.CategoricalIndex([1, 3, 2], categories=[1, 2, 3],
                                      ordered=True)
        exp = pd.Series([3, 2, 1], index=exp_idx, name='xxx')

        s = pd.Series(values, name='xxx')
        tm.assert_series_equal(s.value_counts(), exp)
        # check CategoricalIndex outputs the same result
        idx = pd.CategoricalIndex(values, name='xxx')
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3., 2., 1]) / 6.,
                        index=exp_idx, name='xxx')
        tm.assert_series_equal(s.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical_not_ordered(self):
        values = pd.Categorical([1, 2, 3, 1, 1, 3], ordered=False)

        exp_idx = pd.CategoricalIndex([1, 3, 2], categories=[1, 2, 3],
                                      ordered=False)
        exp = pd.Series([3, 2, 1], index=exp_idx, name='xxx')

        s = pd.Series(values, name='xxx')
        tm.assert_series_equal(s.value_counts(), exp)
        # check CategoricalIndex outputs the same result
        idx = pd.CategoricalIndex(values, name='xxx')
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3., 2., 1]) / 6.,
                        index=exp_idx, name='xxx')
        tm.assert_series_equal(s.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)
