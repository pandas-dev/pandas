# -*- coding: utf-8 -*-

from __future__ import print_function

from datetime import timedelta, datetime
from distutils.version import LooseVersion
import sys
import nose

from numpy import nan
from numpy.random import randn
import numpy as np

from pandas.compat import lrange
from pandas import (compat, isnull, notnull, DataFrame, Series,
                    MultiIndex, date_range, Timestamp, _np_version_under1p11)
import pandas as pd
import pandas.core.common as com
import pandas.core.nanops as nanops

from pandas.util.testing import (assert_almost_equal,
                                 assert_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp)

import pandas.util.testing as tm
from pandas import _np_version_under1p9

from pandas.tests.frame.common import TestData


class TestDataFrameAnalytics(tm.TestCase, TestData):

    _multiprocess_can_split_ = True

    # ---------------------------------------------------------------------=
    # Correlation and covariance

    def test_corr_pearson(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('pearson')

    def test_corr_kendall(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('kendall')

    def test_corr_spearman(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('spearman')

    def _check_method(self, method='pearson', check_minp=False):
        if not check_minp:
            correls = self.frame.corr(method=method)
            exp = self.frame['A'].corr(self.frame['C'], method=method)
            assert_almost_equal(correls['A']['C'], exp)
        else:
            result = self.frame.corr(min_periods=len(self.frame) - 8)
            expected = self.frame.corr()
            expected.ix['A', 'B'] = expected.ix['B', 'A'] = nan
            assert_frame_equal(result, expected)

    def test_corr_non_numeric(self):
        tm._skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        # exclude non-numeric types
        result = self.mixed_frame.corr()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].corr()
        assert_frame_equal(result, expected)

    def test_corr_nooverlap(self):
        tm._skip_if_no_scipy()

        # nothing in common
        for meth in ['pearson', 'kendall', 'spearman']:
            df = DataFrame({'A': [1, 1.5, 1, np.nan, np.nan, np.nan],
                            'B': [np.nan, np.nan, np.nan, 1, 1.5, 1],
                            'C': [np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan]})
            rs = df.corr(meth)
            self.assertTrue(isnull(rs.ix['A', 'B']))
            self.assertTrue(isnull(rs.ix['B', 'A']))
            self.assertEqual(rs.ix['A', 'A'], 1)
            self.assertEqual(rs.ix['B', 'B'], 1)
            self.assertTrue(isnull(rs.ix['C', 'C']))

    def test_corr_constant(self):
        tm._skip_if_no_scipy()

        # constant --> all NA

        for meth in ['pearson', 'spearman']:
            df = DataFrame({'A': [1, 1, 1, np.nan, np.nan, np.nan],
                            'B': [np.nan, np.nan, np.nan, 1, 1, 1]})
            rs = df.corr(meth)
            self.assertTrue(isnull(rs.values).all())

    def test_corr_int(self):
        # dtypes other than float64 #1761
        df3 = DataFrame({"a": [1, 2, 3, 4], "b": [1, 2, 3, 4]})

        # it works!
        df3.cov()
        df3.corr()

    def test_corr_int_and_boolean(self):
        tm._skip_if_no_scipy()

        # when dtypes of pandas series are different
        # then ndarray will have dtype=object,
        # so it need to be properly handled
        df = DataFrame({"a": [True, False], "b": [1, 0]})

        expected = DataFrame(np.ones((2, 2)), index=[
                             'a', 'b'], columns=['a', 'b'])
        for meth in ['pearson', 'kendall', 'spearman']:
            assert_frame_equal(df.corr(meth), expected)

    def test_cov(self):
        # min_periods no NAs (corner case)
        expected = self.frame.cov()
        result = self.frame.cov(min_periods=len(self.frame))

        assert_frame_equal(expected, result)

        result = self.frame.cov(min_periods=len(self.frame) + 1)
        self.assertTrue(isnull(result.values).all())

        # with NAs
        frame = self.frame.copy()
        frame['A'][:5] = nan
        frame['B'][5:10] = nan
        result = self.frame.cov(min_periods=len(self.frame) - 8)
        expected = self.frame.cov()
        expected.ix['A', 'B'] = np.nan
        expected.ix['B', 'A'] = np.nan

        # regular
        self.frame['A'][:5] = nan
        self.frame['B'][:10] = nan
        cov = self.frame.cov()

        assert_almost_equal(cov['A']['C'],
                            self.frame['A'].cov(self.frame['C']))

        # exclude non-numeric types
        result = self.mixed_frame.cov()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].cov()
        assert_frame_equal(result, expected)

        # Single column frame
        df = DataFrame(np.linspace(0.0, 1.0, 10))
        result = df.cov()
        expected = DataFrame(np.cov(df.values.T).reshape((1, 1)),
                             index=df.columns, columns=df.columns)
        assert_frame_equal(result, expected)
        df.ix[0] = np.nan
        result = df.cov()
        expected = DataFrame(np.cov(df.values[1:].T).reshape((1, 1)),
                             index=df.columns, columns=df.columns)
        assert_frame_equal(result, expected)

    def test_corrwith(self):
        a = self.tsframe
        noise = Series(randn(len(a)), index=a.index)

        b = self.tsframe.add(noise, axis=0)

        # make sure order does not matter
        b = b.reindex(columns=b.columns[::-1], index=b.index[::-1][10:])
        del b['B']

        colcorr = a.corrwith(b, axis=0)
        assert_almost_equal(colcorr['A'], a['A'].corr(b['A']))

        rowcorr = a.corrwith(b, axis=1)
        assert_series_equal(rowcorr, a.T.corrwith(b.T, axis=0))

        dropped = a.corrwith(b, axis=0, drop=True)
        assert_almost_equal(dropped['A'], a['A'].corr(b['A']))
        self.assertNotIn('B', dropped)

        dropped = a.corrwith(b, axis=1, drop=True)
        self.assertNotIn(a.index[-1], dropped.index)

        # non time-series data
        index = ['a', 'b', 'c', 'd', 'e']
        columns = ['one', 'two', 'three', 'four']
        df1 = DataFrame(randn(5, 4), index=index, columns=columns)
        df2 = DataFrame(randn(4, 4), index=index[:4], columns=columns)
        correls = df1.corrwith(df2, axis=1)
        for row in index[:4]:
            assert_almost_equal(correls[row], df1.ix[row].corr(df2.ix[row]))

    def test_corrwith_with_objects(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame()
        cols = ['A', 'B', 'C', 'D']

        df1['obj'] = 'foo'
        df2['obj'] = 'bar'

        result = df1.corrwith(df2)
        expected = df1.ix[:, cols].corrwith(df2.ix[:, cols])
        assert_series_equal(result, expected)

        result = df1.corrwith(df2, axis=1)
        expected = df1.ix[:, cols].corrwith(df2.ix[:, cols], axis=1)
        assert_series_equal(result, expected)

    def test_corrwith_series(self):
        result = self.tsframe.corrwith(self.tsframe['A'])
        expected = self.tsframe.apply(self.tsframe['A'].corr)

        assert_series_equal(result, expected)

    def test_corrwith_matches_corrcoef(self):
        df1 = DataFrame(np.arange(10000), columns=['a'])
        df2 = DataFrame(np.arange(10000) ** 2, columns=['a'])
        c1 = df1.corrwith(df2)['a']
        c2 = np.corrcoef(df1['a'], df2['a'])[0][1]

        assert_almost_equal(c1, c2)
        self.assertTrue(c1 < 1)

    def test_bool_describe_in_mixed_frame(self):
        df = DataFrame({
            'string_data': ['a', 'b', 'c', 'd', 'e'],
            'bool_data': [True, True, False, False, False],
            'int_data': [10, 20, 30, 40, 50],
        })

        # Integer data are included in .describe() output,
        # Boolean and string data are not.
        result = df.describe()
        expected = DataFrame({'int_data': [5, 30, df.int_data.std(),
                                           10, 20, 30, 40, 50]},
                             index=['count', 'mean', 'std', 'min', '25%',
                                    '50%', '75%', 'max'])
        assert_frame_equal(result, expected)

        # Top value is a boolean value that is False
        result = df.describe(include=['bool'])

        expected = DataFrame({'bool_data': [5, 2, False, 3]},
                             index=['count', 'unique', 'top', 'freq'])
        assert_frame_equal(result, expected)

    def test_reduce_mixed_frame(self):
        # GH 6806
        df = DataFrame({
            'bool_data': [True, True, False, False, False],
            'int_data': [10, 20, 30, 40, 50],
            'string_data': ['a', 'b', 'c', 'd', 'e'],
        })
        df.reindex(columns=['bool_data', 'int_data', 'string_data'])
        test = df.sum(axis=0)
        assert_almost_equal(test.values, [2, 150, 'abcde'])
        assert_series_equal(test, df.T.sum(axis=1))

    def test_count(self):
        f = lambda s: notnull(s).sum()
        self._check_stat_op('count', f,
                            has_skipna=False,
                            has_numeric_only=True,
                            check_dtype=False,
                            check_dates=True)

        # corner case
        frame = DataFrame()
        ct1 = frame.count(1)
        tm.assertIsInstance(ct1, Series)

        ct2 = frame.count(0)
        tm.assertIsInstance(ct2, Series)

        # GH #423
        df = DataFrame(index=lrange(10))
        result = df.count(1)
        expected = Series(0, index=df.index)
        assert_series_equal(result, expected)

        df = DataFrame(columns=lrange(10))
        result = df.count(0)
        expected = Series(0, index=df.columns)
        assert_series_equal(result, expected)

        df = DataFrame()
        result = df.count()
        expected = Series(0, index=[])
        assert_series_equal(result, expected)

    def test_sum(self):
        self._check_stat_op('sum', np.sum, has_numeric_only=True)

        # mixed types (with upcasting happening)
        self._check_stat_op('sum', np.sum,
                            frame=self.mixed_float.astype('float32'),
                            has_numeric_only=True, check_dtype=False,
                            check_less_precise=True)

    def test_stat_operators_attempt_obj_array(self):
        data = {
            'a': [-0.00049987540199591344, -0.0016467257772919831,
                  0.00067695870775883013],
            'b': [-0, -0, 0.0],
            'c': [0.00031111847529610595, 0.0014902627951905339,
                  -0.00094099200035979691]
        }
        df1 = DataFrame(data, index=['foo', 'bar', 'baz'],
                        dtype='O')
        methods = ['sum', 'mean', 'prod', 'var', 'std', 'skew', 'min', 'max']

        # GH #676
        df2 = DataFrame({0: [np.nan, 2], 1: [np.nan, 3],
                         2: [np.nan, 4]}, dtype=object)

        for df in [df1, df2]:
            for meth in methods:
                self.assertEqual(df.values.dtype, np.object_)
                result = getattr(df, meth)(1)
                expected = getattr(df.astype('f8'), meth)(1)

                if not tm._incompat_bottleneck_version(meth):
                    assert_series_equal(result, expected)

    def test_mean(self):
        self._check_stat_op('mean', np.mean, check_dates=True)

    def test_product(self):
        self._check_stat_op('product', np.prod)

    def test_median(self):
        def wrapper(x):
            if isnull(x).any():
                return np.nan
            return np.median(x)

        self._check_stat_op('median', wrapper, check_dates=True)

    def test_min(self):
        self._check_stat_op('min', np.min, check_dates=True)
        self._check_stat_op('min', np.min, frame=self.intframe)

    def test_cummin(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cummin = self.tsframe.cummin()
        expected = self.tsframe.apply(Series.cummin)
        assert_frame_equal(cummin, expected)

        # axis = 1
        cummin = self.tsframe.cummin(axis=1)
        expected = self.tsframe.apply(Series.cummin, axis=1)
        assert_frame_equal(cummin, expected)

        # it works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cummin()  # noqa

        # fix issue
        cummin_xs = self.tsframe.cummin(axis=1)
        self.assertEqual(np.shape(cummin_xs), np.shape(self.tsframe))

    def test_cummax(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cummax = self.tsframe.cummax()
        expected = self.tsframe.apply(Series.cummax)
        assert_frame_equal(cummax, expected)

        # axis = 1
        cummax = self.tsframe.cummax(axis=1)
        expected = self.tsframe.apply(Series.cummax, axis=1)
        assert_frame_equal(cummax, expected)

        # it works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cummax()  # noqa

        # fix issue
        cummax_xs = self.tsframe.cummax(axis=1)
        self.assertEqual(np.shape(cummax_xs), np.shape(self.tsframe))

    def test_max(self):
        self._check_stat_op('max', np.max, check_dates=True)
        self._check_stat_op('max', np.max, frame=self.intframe)

    def test_mad(self):
        f = lambda x: np.abs(x - x.mean()).mean()
        self._check_stat_op('mad', f)

    def test_var_std(self):
        alt = lambda x: np.var(x, ddof=1)
        self._check_stat_op('var', alt)

        alt = lambda x: np.std(x, ddof=1)
        self._check_stat_op('std', alt)

        result = self.tsframe.std(ddof=4)
        expected = self.tsframe.apply(lambda x: x.std(ddof=4))
        assert_almost_equal(result, expected)

        result = self.tsframe.var(ddof=4)
        expected = self.tsframe.apply(lambda x: x.var(ddof=4))
        assert_almost_equal(result, expected)

        arr = np.repeat(np.random.random((1, 1000)), 1000, 0)
        result = nanops.nanvar(arr, axis=0)
        self.assertFalse((result < 0).any())
        if nanops._USE_BOTTLENECK:
            nanops._USE_BOTTLENECK = False
            result = nanops.nanvar(arr, axis=0)
            self.assertFalse((result < 0).any())
            nanops._USE_BOTTLENECK = True

    def test_numeric_only_flag(self):
        # GH #9201
        methods = ['sem', 'var', 'std']
        df1 = DataFrame(np.random.randn(5, 3), columns=['foo', 'bar', 'baz'])
        # set one entry to a number in str format
        df1.ix[0, 'foo'] = '100'

        df2 = DataFrame(np.random.randn(5, 3), columns=['foo', 'bar', 'baz'])
        # set one entry to a non-number str
        df2.ix[0, 'foo'] = 'a'

        for meth in methods:
            result = getattr(df1, meth)(axis=1, numeric_only=True)
            expected = getattr(df1[['bar', 'baz']], meth)(axis=1)
            assert_series_equal(expected, result)

            result = getattr(df2, meth)(axis=1, numeric_only=True)
            expected = getattr(df2[['bar', 'baz']], meth)(axis=1)
            assert_series_equal(expected, result)

            # df1 has all numbers, df2 has a letter inside
            self.assertRaises(TypeError, lambda: getattr(df1, meth)
                              (axis=1, numeric_only=False))
            self.assertRaises(TypeError, lambda: getattr(df2, meth)
                              (axis=1, numeric_only=False))

    def test_quantile(self):
        from numpy import percentile

        q = self.tsframe.quantile(0.1, axis=0)
        self.assertEqual(q['A'], percentile(self.tsframe['A'], 10))
        q = self.tsframe.quantile(0.9, axis=1)
        q = self.intframe.quantile(0.1)
        self.assertEqual(q['A'], percentile(self.intframe['A'], 10))

        # test degenerate case
        q = DataFrame({'x': [], 'y': []}).quantile(0.1, axis=0)
        assert(np.isnan(q['x']) and np.isnan(q['y']))

        # non-numeric exclusion
        df = DataFrame({'col1': ['A', 'A', 'B', 'B'], 'col2': [1, 2, 3, 4]})
        rs = df.quantile(0.5)
        xp = df.median()
        assert_series_equal(rs, xp)

        # axis
        df = DataFrame({"A": [1, 2, 3], "B": [2, 3, 4]}, index=[1, 2, 3])
        result = df.quantile(.5, axis=1)
        expected = Series([1.5, 2.5, 3.5], index=[1, 2, 3])
        assert_series_equal(result, expected)

        result = df.quantile([.5, .75], axis=1)
        expected = DataFrame({1: [1.5, 1.75], 2: [2.5, 2.75],
                              3: [3.5, 3.75]}, index=[0.5, 0.75])
        assert_frame_equal(result, expected, check_index_type=True)

        # We may want to break API in the future to change this
        # so that we exclude non-numeric along the same axis
        # See GH #7312
        df = DataFrame([[1, 2, 3],
                        ['a', 'b', 4]])
        result = df.quantile(.5, axis=1)
        expected = Series([3., 4.], index=[0, 1])
        assert_series_equal(result, expected)

    def test_quantile_axis_parameter(self):
        # GH 9543/9544

        df = DataFrame({"A": [1, 2, 3], "B": [2, 3, 4]}, index=[1, 2, 3])

        result = df.quantile(.5, axis=0)

        expected = Series([2., 3.], index=["A", "B"])
        assert_series_equal(result, expected)

        expected = df.quantile(.5, axis="index")
        assert_series_equal(result, expected)

        result = df.quantile(.5, axis=1)

        expected = Series([1.5, 2.5, 3.5], index=[1, 2, 3])
        assert_series_equal(result, expected)

        result = df.quantile(.5, axis="columns")
        assert_series_equal(result, expected)

        self.assertRaises(ValueError, df.quantile, 0.1, axis=-1)
        self.assertRaises(ValueError, df.quantile, 0.1, axis="column")

    def test_quantile_interpolation(self):
        # GH #10174
        if _np_version_under1p9:
            raise nose.SkipTest("Numpy version under 1.9")

        from numpy import percentile

        # interpolation = linear (default case)
        q = self.tsframe.quantile(0.1, axis=0, interpolation='linear')
        self.assertEqual(q['A'], percentile(self.tsframe['A'], 10))
        q = self.intframe.quantile(0.1)
        self.assertEqual(q['A'], percentile(self.intframe['A'], 10))

        # test with and without interpolation keyword
        q1 = self.intframe.quantile(0.1)
        self.assertEqual(q1['A'], np.percentile(self.intframe['A'], 10))
        assert_series_equal(q, q1)

        # interpolation method other than default linear
        df = DataFrame({"A": [1, 2, 3], "B": [2, 3, 4]}, index=[1, 2, 3])
        result = df.quantile(.5, axis=1, interpolation='nearest')
        expected = Series([1., 2., 3.], index=[1, 2, 3])
        assert_series_equal(result, expected)

        # axis
        result = df.quantile([.5, .75], axis=1, interpolation='lower')
        expected = DataFrame({1: [1., 1.], 2: [2., 2.],
                              3: [3., 3.]}, index=[0.5, 0.75])
        assert_frame_equal(result, expected)

        # test degenerate case
        df = DataFrame({'x': [], 'y': []})
        q = df.quantile(0.1, axis=0, interpolation='higher')
        assert(np.isnan(q['x']) and np.isnan(q['y']))

        # multi
        df = DataFrame([[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                       columns=['a', 'b', 'c'])
        result = df.quantile([.25, .5], interpolation='midpoint')

        # https://github.com/numpy/numpy/issues/7163
        if _np_version_under1p11:
            expected = DataFrame([[1.5, 1.5, 1.5], [2.5, 2.5, 2.5]],
                                 index=[.25, .5], columns=['a', 'b', 'c'])
        else:
            expected = DataFrame([[1.5, 1.5, 1.5], [2.0, 2.0, 2.0]],
                                 index=[.25, .5], columns=['a', 'b', 'c'])
        assert_frame_equal(result, expected)

    def test_quantile_interpolation_np_lt_1p9(self):
        # GH #10174
        if not _np_version_under1p9:
            raise nose.SkipTest("Numpy version is greater than 1.9")

        from numpy import percentile

        # interpolation = linear (default case)
        q = self.tsframe.quantile(0.1, axis=0, interpolation='linear')
        self.assertEqual(q['A'], percentile(self.tsframe['A'], 10))
        q = self.intframe.quantile(0.1)
        self.assertEqual(q['A'], percentile(self.intframe['A'], 10))

        # test with and without interpolation keyword
        q1 = self.intframe.quantile(0.1)
        self.assertEqual(q1['A'], np.percentile(self.intframe['A'], 10))
        assert_series_equal(q, q1)

        # interpolation method other than default linear
        expErrMsg = "Interpolation methods other than linear"
        df = DataFrame({"A": [1, 2, 3], "B": [2, 3, 4]}, index=[1, 2, 3])
        with assertRaisesRegexp(ValueError, expErrMsg):
            df.quantile(.5, axis=1, interpolation='nearest')

        with assertRaisesRegexp(ValueError, expErrMsg):
            df.quantile([.5, .75], axis=1, interpolation='lower')

        # test degenerate case
        df = DataFrame({'x': [], 'y': []})
        with assertRaisesRegexp(ValueError, expErrMsg):
            q = df.quantile(0.1, axis=0, interpolation='higher')

        # multi
        df = DataFrame([[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                       columns=['a', 'b', 'c'])
        with assertRaisesRegexp(ValueError, expErrMsg):
            df.quantile([.25, .5], interpolation='midpoint')

    def test_quantile_multi(self):
        df = DataFrame([[1, 1, 1], [2, 2, 2], [3, 3, 3]],
                       columns=['a', 'b', 'c'])
        result = df.quantile([.25, .5])
        expected = DataFrame([[1.5, 1.5, 1.5], [2., 2., 2.]],
                             index=[.25, .5], columns=['a', 'b', 'c'])
        assert_frame_equal(result, expected)

        # axis = 1
        result = df.quantile([.25, .5], axis=1)
        expected = DataFrame([[1.5, 1.5, 1.5], [2., 2., 2.]],
                             index=[.25, .5], columns=[0, 1, 2])

        # empty
        result = DataFrame({'x': [], 'y': []}).quantile([0.1, .9], axis=0)
        expected = DataFrame({'x': [np.nan, np.nan], 'y': [np.nan, np.nan]},
                             index=[.1, .9])
        assert_frame_equal(result, expected)

    def test_quantile_datetime(self):
        df = DataFrame({'a': pd.to_datetime(['2010', '2011']), 'b': [0, 5]})

        # exclude datetime
        result = df.quantile(.5)
        expected = Series([2.5], index=['b'])

        # datetime
        result = df.quantile(.5, numeric_only=False)
        expected = Series([Timestamp('2010-07-02 12:00:00'), 2.5],
                          index=['a', 'b'])
        assert_series_equal(result, expected)

        # datetime w/ multi
        result = df.quantile([.5], numeric_only=False)
        expected = DataFrame([[Timestamp('2010-07-02 12:00:00'), 2.5]],
                             index=[.5], columns=['a', 'b'])
        assert_frame_equal(result, expected)

        # axis = 1
        df['c'] = pd.to_datetime(['2011', '2012'])
        result = df[['a', 'c']].quantile(.5, axis=1, numeric_only=False)
        expected = Series([Timestamp('2010-07-02 12:00:00'),
                           Timestamp('2011-07-02 12:00:00')],
                          index=[0, 1])
        assert_series_equal(result, expected)

        result = df[['a', 'c']].quantile([.5], axis=1, numeric_only=False)
        expected = DataFrame([[Timestamp('2010-07-02 12:00:00'),
                               Timestamp('2011-07-02 12:00:00')]],
                             index=[0.5], columns=[0, 1])
        assert_frame_equal(result, expected)

    def test_quantile_invalid(self):
        msg = 'percentiles should all be in the interval \\[0, 1\\]'
        for invalid in [-1, 2, [0.5, -1], [0.5, 2]]:
            with tm.assertRaisesRegexp(ValueError, msg):
                self.tsframe.quantile(invalid)

    def test_cumsum(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cumsum = self.tsframe.cumsum()
        expected = self.tsframe.apply(Series.cumsum)
        assert_frame_equal(cumsum, expected)

        # axis = 1
        cumsum = self.tsframe.cumsum(axis=1)
        expected = self.tsframe.apply(Series.cumsum, axis=1)
        assert_frame_equal(cumsum, expected)

        # works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cumsum()  # noqa

        # fix issue
        cumsum_xs = self.tsframe.cumsum(axis=1)
        self.assertEqual(np.shape(cumsum_xs), np.shape(self.tsframe))

    def test_cumprod(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cumprod = self.tsframe.cumprod()
        expected = self.tsframe.apply(Series.cumprod)
        assert_frame_equal(cumprod, expected)

        # axis = 1
        cumprod = self.tsframe.cumprod(axis=1)
        expected = self.tsframe.apply(Series.cumprod, axis=1)
        assert_frame_equal(cumprod, expected)

        # fix issue
        cumprod_xs = self.tsframe.cumprod(axis=1)
        self.assertEqual(np.shape(cumprod_xs), np.shape(self.tsframe))

        # ints
        df = self.tsframe.fillna(0).astype(int)
        df.cumprod(0)
        df.cumprod(1)

        # ints32
        df = self.tsframe.fillna(0).astype(np.int32)
        df.cumprod(0)
        df.cumprod(1)

    def test_rank(self):
        tm._skip_if_no_scipy()
        from scipy.stats import rankdata

        self.frame['A'][::2] = np.nan
        self.frame['B'][::3] = np.nan
        self.frame['C'][::4] = np.nan
        self.frame['D'][::5] = np.nan

        ranks0 = self.frame.rank()
        ranks1 = self.frame.rank(1)
        mask = np.isnan(self.frame.values)

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fvals)
        exp0[mask] = np.nan

        exp1 = np.apply_along_axis(rankdata, 1, fvals)
        exp1[mask] = np.nan

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # integers
        df = DataFrame(np.random.randint(0, 5, size=40).reshape((10, 4)))

        result = df.rank()
        exp = df.astype(float).rank()
        assert_frame_equal(result, exp)

        result = df.rank(1)
        exp = df.astype(float).rank(1)
        assert_frame_equal(result, exp)

    def test_rank2(self):
        df = DataFrame([[1, 3, 2], [1, 2, 3]])
        expected = DataFrame([[1.0, 3.0, 2.0], [1, 2, 3]]) / 3.0
        result = df.rank(1, pct=True)
        assert_frame_equal(result, expected)

        df = DataFrame([[1, 3, 2], [1, 2, 3]])
        expected = df.rank(0) / 2.0
        result = df.rank(0, pct=True)
        assert_frame_equal(result, expected)

        df = DataFrame([['b', 'c', 'a'], ['a', 'c', 'b']])
        expected = DataFrame([[2.0, 3.0, 1.0], [1, 3, 2]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        expected = DataFrame([[2.0, 1.5, 1.0], [1, 1.5, 2]])
        result = df.rank(0, numeric_only=False)
        assert_frame_equal(result, expected)

        df = DataFrame([['b', np.nan, 'a'], ['a', 'c', 'b']])
        expected = DataFrame([[2.0, nan, 1.0], [1.0, 3.0, 2.0]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        expected = DataFrame([[2.0, nan, 1.0], [1.0, 1.0, 2.0]])
        result = df.rank(0, numeric_only=False)
        assert_frame_equal(result, expected)

        # f7u12, this does not work without extensive workaround
        data = [[datetime(2001, 1, 5), nan, datetime(2001, 1, 2)],
                [datetime(2000, 1, 2), datetime(2000, 1, 3),
                 datetime(2000, 1, 1)]]
        df = DataFrame(data)

        # check the rank
        expected = DataFrame([[2., nan, 1.],
                              [2., 3., 1.]])
        result = df.rank(1, numeric_only=False, ascending=True)
        assert_frame_equal(result, expected)

        expected = DataFrame([[1., nan, 2.],
                              [2., 1., 3.]])
        result = df.rank(1, numeric_only=False, ascending=False)
        assert_frame_equal(result, expected)

        # mixed-type frames
        self.mixed_frame['datetime'] = datetime.now()
        self.mixed_frame['timedelta'] = timedelta(days=1, seconds=1)

        result = self.mixed_frame.rank(1)
        expected = self.mixed_frame.rank(1, numeric_only=True)
        assert_frame_equal(result, expected)

        df = DataFrame({"a": [1e-20, -5, 1e-20 + 1e-40, 10,
                              1e60, 1e80, 1e-30]})
        exp = DataFrame({"a": [3.5, 1., 3.5, 5., 6., 7., 2.]})
        assert_frame_equal(df.rank(), exp)

    def test_rank_na_option(self):
        tm._skip_if_no_scipy()
        from scipy.stats import rankdata

        self.frame['A'][::2] = np.nan
        self.frame['B'][::3] = np.nan
        self.frame['C'][::4] = np.nan
        self.frame['D'][::5] = np.nan

        # bottom
        ranks0 = self.frame.rank(na_option='bottom')
        ranks1 = self.frame.rank(1, na_option='bottom')

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fvals)
        exp1 = np.apply_along_axis(rankdata, 1, fvals)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # top
        ranks0 = self.frame.rank(na_option='top')
        ranks1 = self.frame.rank(1, na_option='top')

        fval0 = self.frame.fillna((self.frame.min() - 1).to_dict()).values
        fval1 = self.frame.T
        fval1 = fval1.fillna((fval1.min() - 1).to_dict()).T
        fval1 = fval1.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fval0)
        exp1 = np.apply_along_axis(rankdata, 1, fval1)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # descending

        # bottom
        ranks0 = self.frame.rank(na_option='top', ascending=False)
        ranks1 = self.frame.rank(1, na_option='top', ascending=False)

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, -fvals)
        exp1 = np.apply_along_axis(rankdata, 1, -fvals)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # descending

        # top
        ranks0 = self.frame.rank(na_option='bottom', ascending=False)
        ranks1 = self.frame.rank(1, na_option='bottom', ascending=False)

        fval0 = self.frame.fillna((self.frame.min() - 1).to_dict()).values
        fval1 = self.frame.T
        fval1 = fval1.fillna((fval1.min() - 1).to_dict()).T
        fval1 = fval1.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, -fval0)
        exp1 = np.apply_along_axis(rankdata, 1, -fval1)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

    def test_rank_axis(self):
        # check if using axes' names gives the same result
        df = pd.DataFrame([[2, 1], [4, 3]])
        assert_frame_equal(df.rank(axis=0), df.rank(axis='index'))
        assert_frame_equal(df.rank(axis=1), df.rank(axis='columns'))

    def test_sem(self):
        alt = lambda x: np.std(x, ddof=1) / np.sqrt(len(x))
        self._check_stat_op('sem', alt)

        result = self.tsframe.sem(ddof=4)
        expected = self.tsframe.apply(
            lambda x: x.std(ddof=4) / np.sqrt(len(x)))
        assert_almost_equal(result, expected)

        arr = np.repeat(np.random.random((1, 1000)), 1000, 0)
        result = nanops.nansem(arr, axis=0)
        self.assertFalse((result < 0).any())
        if nanops._USE_BOTTLENECK:
            nanops._USE_BOTTLENECK = False
            result = nanops.nansem(arr, axis=0)
            self.assertFalse((result < 0).any())
            nanops._USE_BOTTLENECK = True

    def test_skew(self):
        tm._skip_if_no_scipy()
        from scipy.stats import skew

        def alt(x):
            if len(x) < 3:
                return np.nan
            return skew(x, bias=False)

        self._check_stat_op('skew', alt)

    def test_kurt(self):
        tm._skip_if_no_scipy()

        from scipy.stats import kurtosis

        def alt(x):
            if len(x) < 4:
                return np.nan
            return kurtosis(x, bias=False)

        self._check_stat_op('kurt', alt)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        df = DataFrame(np.random.randn(6, 3), index=index)

        kurt = df.kurt()
        kurt2 = df.kurt(level=0).xs('bar')
        assert_series_equal(kurt, kurt2, check_names=False)
        self.assertTrue(kurt.name is None)
        self.assertEqual(kurt2.name, 'bar')

    def _check_stat_op(self, name, alternative, frame=None, has_skipna=True,
                       has_numeric_only=False, check_dtype=True,
                       check_dates=False, check_less_precise=False):
        if frame is None:
            frame = self.frame
            # set some NAs
            frame.ix[5:10] = np.nan
            frame.ix[15:20, -2:] = np.nan

        f = getattr(frame, name)

        if check_dates:
            df = DataFrame({'b': date_range('1/1/2001', periods=2)})
            _f = getattr(df, name)
            result = _f()
            self.assertIsInstance(result, Series)

            df['a'] = lrange(len(df))
            result = getattr(df, name)()
            self.assertIsInstance(result, Series)
            self.assertTrue(len(result))

        if has_skipna:
            def skipna_wrapper(x):
                nona = x.dropna()
                if len(nona) == 0:
                    return np.nan
                return alternative(nona)

            def wrapper(x):
                return alternative(x.values)

            result0 = f(axis=0, skipna=False)
            result1 = f(axis=1, skipna=False)
            assert_series_equal(result0, frame.apply(wrapper),
                                check_dtype=check_dtype,
                                check_less_precise=check_less_precise)
            # HACK: win32
            assert_series_equal(result1, frame.apply(wrapper, axis=1),
                                check_dtype=False,
                                check_less_precise=check_less_precise)
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        result0 = f(axis=0)
        result1 = f(axis=1)
        assert_series_equal(result0, frame.apply(skipna_wrapper),
                            check_dtype=check_dtype,
                            check_less_precise=check_less_precise)
        if not tm._incompat_bottleneck_version(name):
            assert_series_equal(result1, frame.apply(skipna_wrapper, axis=1),
                                check_dtype=False,
                                check_less_precise=check_less_precise)

        # check dtypes
        if check_dtype:
            lcd_dtype = frame.values.dtype
            self.assertEqual(lcd_dtype, result0.dtype)
            self.assertEqual(lcd_dtype, result1.dtype)

        # result = f(axis=1)
        # comp = frame.apply(alternative, axis=1).reindex(result.index)
        # assert_series_equal(result, comp)

        # bad axis
        assertRaisesRegexp(ValueError, 'No axis named 2', f, axis=2)
        # make sure works on mixed-type frame
        getattr(self.mixed_frame, name)(axis=0)
        getattr(self.mixed_frame, name)(axis=1)

        if has_numeric_only:
            getattr(self.mixed_frame, name)(axis=0, numeric_only=True)
            getattr(self.mixed_frame, name)(axis=1, numeric_only=True)
            getattr(self.frame, name)(axis=0, numeric_only=False)
            getattr(self.frame, name)(axis=1, numeric_only=False)

        # all NA case
        if has_skipna:
            all_na = self.frame * np.NaN
            r0 = getattr(all_na, name)(axis=0)
            r1 = getattr(all_na, name)(axis=1)
            if not tm._incompat_bottleneck_version(name):
                self.assertTrue(np.isnan(r0).all())
                self.assertTrue(np.isnan(r1).all())

    def test_mode(self):
        df = pd.DataFrame({"A": [12, 12, 11, 12, 19, 11],
                           "B": [10, 10, 10, np.nan, 3, 4],
                           "C": [8, 8, 8, 9, 9, 9],
                           "D": np.arange(6, dtype='int64'),
                           "E": [8, 8, 1, 1, 3, 3]})
        assert_frame_equal(df[["A"]].mode(),
                           pd.DataFrame({"A": [12]}))
        expected = pd.Series([], dtype='int64', name='D').to_frame()
        assert_frame_equal(df[["D"]].mode(), expected)
        expected = pd.Series([1, 3, 8], dtype='int64', name='E').to_frame()
        assert_frame_equal(df[["E"]].mode(), expected)
        assert_frame_equal(df[["A", "B"]].mode(),
                           pd.DataFrame({"A": [12], "B": [10.]}))
        assert_frame_equal(df.mode(),
                           pd.DataFrame({"A": [12, np.nan, np.nan],
                                         "B": [10, np.nan, np.nan],
                                         "C": [8, 9, np.nan],
                                         "D": [np.nan, np.nan, np.nan],
                                         "E": [1, 3, 8]}))

        # outputs in sorted order
        df["C"] = list(reversed(df["C"]))
        com.pprint_thing(df["C"])
        com.pprint_thing(df["C"].mode())
        a, b = (df[["A", "B", "C"]].mode(),
                pd.DataFrame({"A": [12, np.nan],
                              "B": [10, np.nan],
                              "C": [8, 9]}))
        com.pprint_thing(a)
        com.pprint_thing(b)
        assert_frame_equal(a, b)
        # should work with heterogeneous types
        df = pd.DataFrame({"A": np.arange(6, dtype='int64'),
                           "B": pd.date_range('2011', periods=6),
                           "C": list('abcdef')})
        exp = pd.DataFrame({"A": pd.Series([], dtype=df["A"].dtype),
                            "B": pd.Series([], dtype=df["B"].dtype),
                            "C": pd.Series([], dtype=df["C"].dtype)})
        assert_frame_equal(df.mode(), exp)

        # and also when not empty
        df.loc[1, "A"] = 0
        df.loc[4, "B"] = df.loc[3, "B"]
        df.loc[5, "C"] = 'e'
        exp = pd.DataFrame({"A": pd.Series([0], dtype=df["A"].dtype),
                            "B": pd.Series([df.loc[3, "B"]],
                                           dtype=df["B"].dtype),
                            "C": pd.Series(['e'], dtype=df["C"].dtype)})

        assert_frame_equal(df.mode(), exp)

    def test_operators_timedelta64(self):
        from datetime import timedelta
        df = DataFrame(dict(A=date_range('2012-1-1', periods=3, freq='D'),
                            B=date_range('2012-1-2', periods=3, freq='D'),
                            C=Timestamp('20120101') -
                            timedelta(minutes=5, seconds=5)))

        diffs = DataFrame(dict(A=df['A'] - df['C'],
                               B=df['A'] - df['B']))

        # min
        result = diffs.min()
        self.assertEqual(result[0], diffs.ix[0, 'A'])
        self.assertEqual(result[1], diffs.ix[0, 'B'])

        result = diffs.min(axis=1)
        self.assertTrue((result == diffs.ix[0, 'B']).all())

        # max
        result = diffs.max()
        self.assertEqual(result[0], diffs.ix[2, 'A'])
        self.assertEqual(result[1], diffs.ix[2, 'B'])

        result = diffs.max(axis=1)
        self.assertTrue((result == diffs['A']).all())

        # abs
        result = diffs.abs()
        result2 = abs(diffs)
        expected = DataFrame(dict(A=df['A'] - df['C'],
                                  B=df['B'] - df['A']))
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        # mixed frame
        mixed = diffs.copy()
        mixed['C'] = 'foo'
        mixed['D'] = 1
        mixed['E'] = 1.
        mixed['F'] = Timestamp('20130101')

        # results in an object array
        from pandas.tseries.timedeltas import (
            _coerce_scalar_to_timedelta_type as _coerce)

        result = mixed.min()
        expected = Series([_coerce(timedelta(seconds=5 * 60 + 5)),
                           _coerce(timedelta(days=-1)),
                           'foo', 1, 1.0,
                           Timestamp('20130101')],
                          index=mixed.columns)
        assert_series_equal(result, expected)

        # excludes numeric
        result = mixed.min(axis=1)
        expected = Series([1, 1, 1.], index=[0, 1, 2])
        assert_series_equal(result, expected)

        # works when only those columns are selected
        result = mixed[['A', 'B']].min(1)
        expected = Series([timedelta(days=-1)] * 3)
        assert_series_equal(result, expected)

        result = mixed[['A', 'B']].min()
        expected = Series([timedelta(seconds=5 * 60 + 5),
                           timedelta(days=-1)], index=['A', 'B'])
        assert_series_equal(result, expected)

        # GH 3106
        df = DataFrame({'time': date_range('20130102', periods=5),
                        'time2': date_range('20130105', periods=5)})
        df['off1'] = df['time2'] - df['time']
        self.assertEqual(df['off1'].dtype, 'timedelta64[ns]')

        df['off2'] = df['time'] - df['time2']
        df._consolidate_inplace()
        self.assertTrue(df['off1'].dtype == 'timedelta64[ns]')
        self.assertTrue(df['off2'].dtype == 'timedelta64[ns]')

    def test_sum_corner(self):
        axis0 = self.empty.sum(0)
        axis1 = self.empty.sum(1)
        tm.assertIsInstance(axis0, Series)
        tm.assertIsInstance(axis1, Series)
        self.assertEqual(len(axis0), 0)
        self.assertEqual(len(axis1), 0)

    def test_sum_object(self):
        values = self.frame.values.astype(int)
        frame = DataFrame(values, index=self.frame.index,
                          columns=self.frame.columns)
        deltas = frame * timedelta(1)
        deltas.sum()

    def test_sum_bool(self):
        # ensure this works, bug report
        bools = np.isnan(self.frame)
        bools.sum(1)
        bools.sum(0)

    def test_mean_corner(self):
        # unit test when have object data
        the_mean = self.mixed_frame.mean(axis=0)
        the_sum = self.mixed_frame.sum(axis=0, numeric_only=True)
        self.assertTrue(the_sum.index.equals(the_mean.index))
        self.assertTrue(len(the_mean.index) < len(self.mixed_frame.columns))

        # xs sum mixed type, just want to know it works...
        the_mean = self.mixed_frame.mean(axis=1)
        the_sum = self.mixed_frame.sum(axis=1, numeric_only=True)
        self.assertTrue(the_sum.index.equals(the_mean.index))

        # take mean of boolean column
        self.frame['bool'] = self.frame['A'] > 0
        means = self.frame.mean(0)
        self.assertEqual(means['bool'], self.frame['bool'].values.mean())

    def test_stats_mixed_type(self):
        # don't blow up
        self.mixed_frame.std(1)
        self.mixed_frame.var(1)
        self.mixed_frame.mean(1)
        self.mixed_frame.skew(1)

    def test_median_corner(self):
        def wrapper(x):
            if isnull(x).any():
                return np.nan
            return np.median(x)

        self._check_stat_op('median', wrapper, frame=self.intframe,
                            check_dtype=False, check_dates=True)

    # Miscellanea

    def test_count_objects(self):
        dm = DataFrame(self.mixed_frame._series)
        df = DataFrame(self.mixed_frame._series)

        assert_series_equal(dm.count(), df.count())
        assert_series_equal(dm.count(1), df.count(1))

    def test_cumsum_corner(self):
        dm = DataFrame(np.arange(20).reshape(4, 5),
                       index=lrange(4), columns=lrange(5))
        # ?(wesm)
        result = dm.cumsum()  # noqa

    def test_sum_bools(self):
        df = DataFrame(index=lrange(1), columns=lrange(10))
        bools = isnull(df)
        self.assertEqual(bools.sum(axis=1)[0], 10)

    # Index of max / min

    def test_idxmin(self):
        frame = self.frame
        frame.ix[5:10] = np.nan
        frame.ix[15:20, -2:] = np.nan
        for skipna in [True, False]:
            for axis in [0, 1]:
                for df in [frame, self.intframe]:
                    result = df.idxmin(axis=axis, skipna=skipna)
                    expected = df.apply(
                        Series.idxmin, axis=axis, skipna=skipna)
                    assert_series_equal(result, expected)

        self.assertRaises(ValueError, frame.idxmin, axis=2)

    def test_idxmax(self):
        frame = self.frame
        frame.ix[5:10] = np.nan
        frame.ix[15:20, -2:] = np.nan
        for skipna in [True, False]:
            for axis in [0, 1]:
                for df in [frame, self.intframe]:
                    result = df.idxmax(axis=axis, skipna=skipna)
                    expected = df.apply(
                        Series.idxmax, axis=axis, skipna=skipna)
                    assert_series_equal(result, expected)

        self.assertRaises(ValueError, frame.idxmax, axis=2)

    # ----------------------------------------------------------------------
    # Logical reductions

    def test_any_all(self):
        self._check_bool_op('any', np.any, has_skipna=True, has_bool_only=True)
        self._check_bool_op('all', np.all, has_skipna=True, has_bool_only=True)

        df = DataFrame(randn(10, 4)) > 0
        df.any(1)
        df.all(1)
        df.any(1, bool_only=True)
        df.all(1, bool_only=True)

        # skip pathological failure cases
        # class CantNonzero(object):

        #     def __nonzero__(self):
        #         raise ValueError

        # df[4] = CantNonzero()

        # it works!
        # df.any(1)
        # df.all(1)
        # df.any(1, bool_only=True)
        # df.all(1, bool_only=True)

        # df[4][4] = np.nan
        # df.any(1)
        # df.all(1)
        # df.any(1, bool_only=True)
        # df.all(1, bool_only=True)

    def _check_bool_op(self, name, alternative, frame=None, has_skipna=True,
                       has_bool_only=False):
        if frame is None:
            frame = self.frame > 0
            # set some NAs
            frame = DataFrame(frame.values.astype(object), frame.index,
                              frame.columns)
            frame.ix[5:10] = np.nan
            frame.ix[15:20, -2:] = np.nan

        f = getattr(frame, name)

        if has_skipna:
            def skipna_wrapper(x):
                nona = x.dropna().values
                return alternative(nona)

            def wrapper(x):
                return alternative(x.values)

            result0 = f(axis=0, skipna=False)
            result1 = f(axis=1, skipna=False)
            assert_series_equal(result0, frame.apply(wrapper))
            assert_series_equal(result1, frame.apply(wrapper, axis=1),
                                check_dtype=False)  # HACK: win32
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        result0 = f(axis=0)
        result1 = f(axis=1)
        assert_series_equal(result0, frame.apply(skipna_wrapper))
        assert_series_equal(result1, frame.apply(skipna_wrapper, axis=1),
                            check_dtype=False)

        # result = f(axis=1)
        # comp = frame.apply(alternative, axis=1).reindex(result.index)
        # assert_series_equal(result, comp)

        # bad axis
        self.assertRaises(ValueError, f, axis=2)

        # make sure works on mixed-type frame
        mixed = self.mixed_frame
        mixed['_bool_'] = np.random.randn(len(mixed)) > 0
        getattr(mixed, name)(axis=0)
        getattr(mixed, name)(axis=1)

        class NonzeroFail:

            def __nonzero__(self):
                raise ValueError

        mixed['_nonzero_fail_'] = NonzeroFail()

        if has_bool_only:
            getattr(mixed, name)(axis=0, bool_only=True)
            getattr(mixed, name)(axis=1, bool_only=True)
            getattr(frame, name)(axis=0, bool_only=False)
            getattr(frame, name)(axis=1, bool_only=False)

        # all NA case
        if has_skipna:
            all_na = frame * np.NaN
            r0 = getattr(all_na, name)(axis=0)
            r1 = getattr(all_na, name)(axis=1)
            if name == 'any':
                self.assertFalse(r0.any())
                self.assertFalse(r1.any())
            else:
                self.assertTrue(r0.all())
                self.assertTrue(r1.all())

    # ----------------------------------------------------------------------
    # Top / bottom

    def test_nlargest(self):
        # GH10393
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10])})
        result = df.nlargest(5, 'a')
        expected = df.sort_values('a', ascending=False).head(5)
        assert_frame_equal(result, expected)

    def test_nlargest_multiple_columns(self):
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10]),
                           'c': np.random.permutation(10).astype('float64')})
        result = df.nlargest(5, ['a', 'b'])
        expected = df.sort_values(['a', 'b'], ascending=False).head(5)
        assert_frame_equal(result, expected)

    def test_nsmallest(self):
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10])})
        result = df.nsmallest(5, 'a')
        expected = df.sort_values('a').head(5)
        assert_frame_equal(result, expected)

    def test_nsmallest_multiple_columns(self):
        from string import ascii_lowercase
        df = pd.DataFrame({'a': np.random.permutation(10),
                           'b': list(ascii_lowercase[:10]),
                           'c': np.random.permutation(10).astype('float64')})
        result = df.nsmallest(5, ['a', 'c'])
        expected = df.sort_values(['a', 'c']).head(5)
        assert_frame_equal(result, expected)

    # ----------------------------------------------------------------------
    # Isin

    def test_isin(self):
        # GH #4211
        df = DataFrame({'vals': [1, 2, 3, 4], 'ids': ['a', 'b', 'f', 'n'],
                        'ids2': ['a', 'n', 'c', 'n']},
                       index=['foo', 'bar', 'baz', 'qux'])
        other = ['a', 'b', 'c']

        result = df.isin(other)
        expected = DataFrame([df.loc[s].isin(other) for s in df.index])
        assert_frame_equal(result, expected)

    def test_isin_empty(self):
        df = DataFrame({'A': ['a', 'b', 'c'], 'B': ['a', 'e', 'f']})
        result = df.isin([])
        expected = pd.DataFrame(False, df.index, df.columns)
        assert_frame_equal(result, expected)

    def test_isin_dict(self):
        df = DataFrame({'A': ['a', 'b', 'c'], 'B': ['a', 'e', 'f']})
        d = {'A': ['a']}

        expected = DataFrame(False, df.index, df.columns)
        expected.loc[0, 'A'] = True

        result = df.isin(d)
        assert_frame_equal(result, expected)

        # non unique columns
        df = DataFrame({'A': ['a', 'b', 'c'], 'B': ['a', 'e', 'f']})
        df.columns = ['A', 'A']
        expected = DataFrame(False, df.index, df.columns)
        expected.loc[0, 'A'] = True
        result = df.isin(d)
        assert_frame_equal(result, expected)

    def test_isin_with_string_scalar(self):
        # GH4763
        df = DataFrame({'vals': [1, 2, 3, 4], 'ids': ['a', 'b', 'f', 'n'],
                        'ids2': ['a', 'n', 'c', 'n']},
                       index=['foo', 'bar', 'baz', 'qux'])
        with tm.assertRaises(TypeError):
            df.isin('a')

        with tm.assertRaises(TypeError):
            df.isin('aaa')

    def test_isin_df(self):
        df1 = DataFrame({'A': [1, 2, 3, 4], 'B': [2, np.nan, 4, 4]})
        df2 = DataFrame({'A': [0, 2, 12, 4], 'B': [2, np.nan, 4, 5]})
        expected = DataFrame(False, df1.index, df1.columns)
        result = df1.isin(df2)
        expected['A'].loc[[1, 3]] = True
        expected['B'].loc[[0, 2]] = True
        assert_frame_equal(result, expected)

        # partial overlapping columns
        df2.columns = ['A', 'C']
        result = df1.isin(df2)
        expected['B'] = False
        assert_frame_equal(result, expected)

    def test_isin_df_dupe_values(self):
        df1 = DataFrame({'A': [1, 2, 3, 4], 'B': [2, np.nan, 4, 4]})
        # just cols duped
        df2 = DataFrame([[0, 2], [12, 4], [2, np.nan], [4, 5]],
                        columns=['B', 'B'])
        with tm.assertRaises(ValueError):
            df1.isin(df2)

        # just index duped
        df2 = DataFrame([[0, 2], [12, 4], [2, np.nan], [4, 5]],
                        columns=['A', 'B'], index=[0, 0, 1, 1])
        with tm.assertRaises(ValueError):
            df1.isin(df2)

        # cols and index:
        df2.columns = ['B', 'B']
        with tm.assertRaises(ValueError):
            df1.isin(df2)

    def test_isin_dupe_self(self):
        other = DataFrame({'A': [1, 0, 1, 0], 'B': [1, 1, 0, 0]})
        df = DataFrame([[1, 1], [1, 0], [0, 0]], columns=['A', 'A'])
        result = df.isin(other)
        expected = DataFrame(False, index=df.index, columns=df.columns)
        expected.loc[0] = True
        expected.iloc[1, 1] = True
        assert_frame_equal(result, expected)

    def test_isin_against_series(self):
        df = pd.DataFrame({'A': [1, 2, 3, 4], 'B': [2, np.nan, 4, 4]},
                          index=['a', 'b', 'c', 'd'])
        s = pd.Series([1, 3, 11, 4], index=['a', 'b', 'c', 'd'])
        expected = DataFrame(False, index=df.index, columns=df.columns)
        expected['A'].loc['a'] = True
        expected.loc['d'] = True
        result = df.isin(s)
        assert_frame_equal(result, expected)

    def test_isin_multiIndex(self):
        idx = MultiIndex.from_tuples([(0, 'a', 'foo'), (0, 'a', 'bar'),
                                      (0, 'b', 'bar'), (0, 'b', 'baz'),
                                      (2, 'a', 'foo'), (2, 'a', 'bar'),
                                      (2, 'c', 'bar'), (2, 'c', 'baz'),
                                      (1, 'b', 'foo'), (1, 'b', 'bar'),
                                      (1, 'c', 'bar'), (1, 'c', 'baz')])
        df1 = DataFrame({'A': np.ones(12),
                         'B': np.zeros(12)}, index=idx)
        df2 = DataFrame({'A': [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                         'B': [1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1]})
        # against regular index
        expected = DataFrame(False, index=df1.index, columns=df1.columns)
        result = df1.isin(df2)
        assert_frame_equal(result, expected)

        df2.index = idx
        expected = df2.values.astype(np.bool)
        expected[:, 1] = ~expected[:, 1]
        expected = DataFrame(expected, columns=['A', 'B'], index=idx)

        result = df1.isin(df2)
        assert_frame_equal(result, expected)

    # ----------------------------------------------------------------------
    # Row deduplication

    def test_drop_duplicates(self):
        df = DataFrame({'AAA': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('AAA')
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep='last')
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep=False)
        expected = df.ix[[]]
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 0)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates('AAA', take_last=True)
            expected = df.ix[[6, 7]]
            assert_frame_equal(result, expected)

        # multi column
        expected = df.ix[[0, 1, 2, 3]]
        result = df.drop_duplicates(np.array(['AAA', 'B']))
        assert_frame_equal(result, expected)
        result = df.drop_duplicates(['AAA', 'B'])
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AAA', 'B'), keep='last')
        expected = df.ix[[0, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AAA', 'B'), keep=False)
        expected = df.ix[[0]]
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(('AAA', 'B'), take_last=True)
        expected = df.ix[[0, 5, 6, 7]]
        assert_frame_equal(result, expected)

        # consider everything
        df2 = df.ix[:, ['AAA', 'B', 'C']]

        result = df2.drop_duplicates()
        # in this case only
        expected = df2.drop_duplicates(['AAA', 'B'])
        assert_frame_equal(result, expected)

        result = df2.drop_duplicates(keep='last')
        expected = df2.drop_duplicates(['AAA', 'B'], keep='last')
        assert_frame_equal(result, expected)

        result = df2.drop_duplicates(keep=False)
        expected = df2.drop_duplicates(['AAA', 'B'], keep=False)
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df2.drop_duplicates(take_last=True)
        with tm.assert_produces_warning(FutureWarning):
            expected = df2.drop_duplicates(['AAA', 'B'], take_last=True)
        assert_frame_equal(result, expected)

        # integers
        result = df.drop_duplicates('C')
        expected = df.iloc[[0, 2]]
        assert_frame_equal(result, expected)
        result = df.drop_duplicates('C', keep='last')
        expected = df.iloc[[-2, -1]]
        assert_frame_equal(result, expected)

        df['E'] = df['C'].astype('int8')
        result = df.drop_duplicates('E')
        expected = df.iloc[[0, 2]]
        assert_frame_equal(result, expected)
        result = df.drop_duplicates('E', keep='last')
        expected = df.iloc[[-2, -1]]
        assert_frame_equal(result, expected)

        # GH 11376
        df = pd.DataFrame({'x': [7, 6, 3, 3, 4, 8, 0],
                           'y': [0, 6, 5, 5, 9, 1, 2]})
        expected = df.loc[df.index != 3]
        assert_frame_equal(df.drop_duplicates(), expected)

        df = pd.DataFrame([[1, 0], [0, 2]])
        assert_frame_equal(df.drop_duplicates(), df)

        df = pd.DataFrame([[-2, 0], [0, -4]])
        assert_frame_equal(df.drop_duplicates(), df)

        x = np.iinfo(np.int64).max / 3 * 2
        df = pd.DataFrame([[-x, x], [0, x + 4]])
        assert_frame_equal(df.drop_duplicates(), df)

        df = pd.DataFrame([[-x, x], [x, x + 4]])
        assert_frame_equal(df.drop_duplicates(), df)

        # GH 11864
        df = pd.DataFrame([i] * 9 for i in range(16))
        df = df.append([[1] + [0] * 8], ignore_index=True)

        for keep in ['first', 'last', False]:
            assert_equal(df.duplicated(keep=keep).sum(), 0)

    def test_drop_duplicates_for_take_all(self):
        df = DataFrame({'AAA': ['foo', 'bar', 'baz', 'bar',
                                'foo', 'bar', 'qux', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('AAA')
        expected = df.iloc[[0, 1, 2, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep='last')
        expected = df.iloc[[2, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', keep=False)
        expected = df.iloc[[2, 6]]
        assert_frame_equal(result, expected)

        # multiple columns
        result = df.drop_duplicates(['AAA', 'B'])
        expected = df.iloc[[0, 1, 2, 3, 4, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['AAA', 'B'], keep='last')
        expected = df.iloc[[0, 1, 2, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['AAA', 'B'], keep=False)
        expected = df.iloc[[0, 1, 2, 6]]
        assert_frame_equal(result, expected)

    def test_drop_duplicates_tuple(self):
        df = DataFrame({('AA', 'AB'): ['foo', 'bar', 'foo', 'bar',
                                       'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates(('AA', 'AB'))
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AA', 'AB'), keep='last')
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AA', 'AB'), keep=False)
        expected = df.ix[[]]  # empty df
        self.assertEqual(len(result), 0)
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(('AA', 'AB'), take_last=True)
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        expected = df.ix[[0, 1, 2, 3]]
        result = df.drop_duplicates((('AA', 'AB'), 'B'))
        assert_frame_equal(result, expected)

    def test_drop_duplicates_NA(self):
        # none
        df = DataFrame({'A': [None, None, 'foo', 'bar',
                              'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('A')
        expected = df.ix[[0, 2, 3]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep='last')
        expected = df.ix[[1, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep=False)
        expected = df.ix[[]]  # empty df
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 0)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates('A', take_last=True)
        expected = df.ix[[1, 6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        result = df.drop_duplicates(['A', 'B'])
        expected = df.ix[[0, 2, 3, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['A', 'B'], keep='last')
        expected = df.ix[[1, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['A', 'B'], keep=False)
        expected = df.ix[[6]]
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(['A', 'B'], take_last=True)
        expected = df.ix[[1, 5, 6, 7]]
        assert_frame_equal(result, expected)

        # nan
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                        'D': lrange(8)})

        # single column
        result = df.drop_duplicates('C')
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep='last')
        expected = df.ix[[3, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep=False)
        expected = df.ix[[]]  # empty df
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 0)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates('C', take_last=True)
        expected = df.ix[[3, 7]]
        assert_frame_equal(result, expected)

        # multi column
        result = df.drop_duplicates(['C', 'B'])
        expected = df.ix[[0, 1, 2, 4]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['C', 'B'], keep='last')
        expected = df.ix[[1, 3, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['C', 'B'], keep=False)
        expected = df.ix[[1]]
        assert_frame_equal(result, expected)

        # deprecate take_last
        with tm.assert_produces_warning(FutureWarning):
            result = df.drop_duplicates(['C', 'B'], take_last=True)
        expected = df.ix[[1, 3, 6, 7]]
        assert_frame_equal(result, expected)

    def test_drop_duplicates_NA_for_take_all(self):
        # none
        df = DataFrame({'A': [None, None, 'foo', 'bar',
                              'foo', 'baz', 'bar', 'qux'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 2., 3, 1.]})

        # single column
        result = df.drop_duplicates('A')
        expected = df.iloc[[0, 2, 3, 5, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep='last')
        expected = df.iloc[[1, 4, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', keep=False)
        expected = df.iloc[[5, 7]]
        assert_frame_equal(result, expected)

        # nan

        # single column
        result = df.drop_duplicates('C')
        expected = df.iloc[[0, 1, 5, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep='last')
        expected = df.iloc[[3, 5, 6, 7]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', keep=False)
        expected = df.iloc[[5, 6]]
        assert_frame_equal(result, expected)

    def test_drop_duplicates_inplace(self):
        orig = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                          'B': ['one', 'one', 'two', 'two',
                                'two', 'two', 'one', 'two'],
                          'C': [1, 1, 2, 2, 2, 2, 1, 2],
                          'D': lrange(8)})

        # single column
        df = orig.copy()
        df.drop_duplicates('A', inplace=True)
        expected = orig[:2]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates('A', keep='last', inplace=True)
        expected = orig.ix[[6, 7]]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates('A', keep=False, inplace=True)
        expected = orig.ix[[]]
        result = df
        assert_frame_equal(result, expected)
        self.assertEqual(len(df), 0)

        # deprecate take_last
        df = orig.copy()
        with tm.assert_produces_warning(FutureWarning):
            df.drop_duplicates('A', take_last=True, inplace=True)
        expected = orig.ix[[6, 7]]
        result = df
        assert_frame_equal(result, expected)

        # multi column
        df = orig.copy()
        df.drop_duplicates(['A', 'B'], inplace=True)
        expected = orig.ix[[0, 1, 2, 3]]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates(['A', 'B'], keep='last', inplace=True)
        expected = orig.ix[[0, 5, 6, 7]]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates(['A', 'B'], keep=False, inplace=True)
        expected = orig.ix[[0]]
        result = df
        assert_frame_equal(result, expected)

        # deprecate take_last
        df = orig.copy()
        with tm.assert_produces_warning(FutureWarning):
            df.drop_duplicates(['A', 'B'], take_last=True, inplace=True)
        expected = orig.ix[[0, 5, 6, 7]]
        result = df
        assert_frame_equal(result, expected)

        # consider everything
        orig2 = orig.ix[:, ['A', 'B', 'C']].copy()

        df2 = orig2.copy()
        df2.drop_duplicates(inplace=True)
        # in this case only
        expected = orig2.drop_duplicates(['A', 'B'])
        result = df2
        assert_frame_equal(result, expected)

        df2 = orig2.copy()
        df2.drop_duplicates(keep='last', inplace=True)
        expected = orig2.drop_duplicates(['A', 'B'], keep='last')
        result = df2
        assert_frame_equal(result, expected)

        df2 = orig2.copy()
        df2.drop_duplicates(keep=False, inplace=True)
        expected = orig2.drop_duplicates(['A', 'B'], keep=False)
        result = df2
        assert_frame_equal(result, expected)

        # deprecate take_last
        df2 = orig2.copy()
        with tm.assert_produces_warning(FutureWarning):
            df2.drop_duplicates(take_last=True, inplace=True)
        with tm.assert_produces_warning(FutureWarning):
            expected = orig2.drop_duplicates(['A', 'B'], take_last=True)
        result = df2
        assert_frame_equal(result, expected)

    # Rounding

    def test_round(self):
        # GH 2665

        # Test that rounding an empty DataFrame does nothing
        df = DataFrame()
        assert_frame_equal(df, df.round())

        # Here's the test frame we'll be working with
        df = DataFrame(
            {'col1': [1.123, 2.123, 3.123], 'col2': [1.234, 2.234, 3.234]})

        # Default round to integer (i.e. decimals=0)
        expected_rounded = DataFrame(
            {'col1': [1., 2., 3.], 'col2': [1., 2., 3.]})
        assert_frame_equal(df.round(), expected_rounded)

        # Round with an integer
        decimals = 2
        expected_rounded = DataFrame(
            {'col1': [1.12, 2.12, 3.12], 'col2': [1.23, 2.23, 3.23]})
        assert_frame_equal(df.round(decimals), expected_rounded)

        # This should also work with np.round (since np.round dispatches to
        # df.round)
        assert_frame_equal(np.round(df, decimals), expected_rounded)

        # Round with a list
        round_list = [1, 2]
        with self.assertRaises(TypeError):
            df.round(round_list)

        # Round with a dictionary
        expected_rounded = DataFrame(
            {'col1': [1.1, 2.1, 3.1], 'col2': [1.23, 2.23, 3.23]})
        round_dict = {'col1': 1, 'col2': 2}
        assert_frame_equal(df.round(round_dict), expected_rounded)

        # Incomplete dict
        expected_partially_rounded = DataFrame(
            {'col1': [1.123, 2.123, 3.123], 'col2': [1.2, 2.2, 3.2]})
        partial_round_dict = {'col2': 1}
        assert_frame_equal(
            df.round(partial_round_dict), expected_partially_rounded)

        # Dict with unknown elements
        wrong_round_dict = {'col3': 2, 'col2': 1}
        assert_frame_equal(
            df.round(wrong_round_dict), expected_partially_rounded)

        # float input to `decimals`
        non_int_round_dict = {'col1': 1, 'col2': 0.5}
        with self.assertRaises(TypeError):
            df.round(non_int_round_dict)

        # String input
        non_int_round_dict = {'col1': 1, 'col2': 'foo'}
        with self.assertRaises(TypeError):
            df.round(non_int_round_dict)

        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        # List input
        non_int_round_dict = {'col1': 1, 'col2': [1, 2]}
        with self.assertRaises(TypeError):
            df.round(non_int_round_dict)

        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        # Non integer Series inputs
        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        non_int_round_Series = Series(non_int_round_dict)
        with self.assertRaises(TypeError):
            df.round(non_int_round_Series)

        # Negative numbers
        negative_round_dict = {'col1': -1, 'col2': -2}
        big_df = df * 100
        expected_neg_rounded = DataFrame(
            {'col1': [110., 210, 310], 'col2': [100., 200, 300]})
        assert_frame_equal(
            big_df.round(negative_round_dict), expected_neg_rounded)

        # nan in Series round
        nan_round_Series = Series({'col1': nan, 'col2': 1})

        # TODO(wesm): unused?
        expected_nan_round = DataFrame({  # noqa
            'col1': [1.123, 2.123, 3.123],
            'col2': [1.2, 2.2, 3.2]})

        if sys.version < LooseVersion('2.7'):
            # Rounding with decimal is a ValueError in Python < 2.7
            with self.assertRaises(ValueError):
                df.round(nan_round_Series)
        else:
            with self.assertRaises(TypeError):
                df.round(nan_round_Series)

        # Make sure this doesn't break existing Series.round
        assert_series_equal(df['col1'].round(1), expected_rounded['col1'])

        # named columns
        # GH 11986
        decimals = 2
        expected_rounded = DataFrame(
            {'col1': [1.12, 2.12, 3.12], 'col2': [1.23, 2.23, 3.23]})
        df.columns.name = "cols"
        expected_rounded.columns.name = "cols"
        assert_frame_equal(df.round(decimals), expected_rounded)

        # interaction of named columns & series
        assert_series_equal(df['col1'].round(decimals),
                            expected_rounded['col1'])
        assert_series_equal(df.round(decimals)['col1'],
                            expected_rounded['col1'])

    def test_round_mixed_type(self):
        # GH11885
        df = DataFrame({'col1': [1.1, 2.2, 3.3, 4.4],
                        'col2': ['1', 'a', 'c', 'f'],
                        'col3': date_range('20111111', periods=4)})
        round_0 = DataFrame({'col1': [1., 2., 3., 4.],
                             'col2': ['1', 'a', 'c', 'f'],
                             'col3': date_range('20111111', periods=4)})
        assert_frame_equal(df.round(), round_0)
        assert_frame_equal(df.round(1), df)
        assert_frame_equal(df.round({'col1': 1}), df)
        assert_frame_equal(df.round({'col1': 0}), round_0)
        assert_frame_equal(df.round({'col1': 0, 'col2': 1}), round_0)
        assert_frame_equal(df.round({'col3': 1}), df)

    def test_round_issue(self):
        # GH11611

        df = pd.DataFrame(np.random.random([3, 3]), columns=['A', 'B', 'C'],
                          index=['first', 'second', 'third'])

        dfs = pd.concat((df, df), axis=1)
        rounded = dfs.round()
        self.assertTrue(rounded.index.equals(dfs.index))

        decimals = pd.Series([1, 0, 2], index=['A', 'B', 'A'])
        self.assertRaises(ValueError, df.round, decimals)

    def test_built_in_round(self):
        if not compat.PY3:
            raise nose.SkipTest("build in round cannot be overriden "
                                "prior to Python 3")

        # GH11763
        # Here's the test frame we'll be working with
        df = DataFrame(
            {'col1': [1.123, 2.123, 3.123], 'col2': [1.234, 2.234, 3.234]})

        # Default round to integer (i.e. decimals=0)
        expected_rounded = DataFrame(
            {'col1': [1., 2., 3.], 'col2': [1., 2., 3.]})
        assert_frame_equal(round(df), expected_rounded)

    # Clip

    def test_clip(self):
        median = self.frame.median().median()

        capped = self.frame.clip_upper(median)
        self.assertFalse((capped.values > median).any())

        floored = self.frame.clip_lower(median)
        self.assertFalse((floored.values < median).any())

        double = self.frame.clip(upper=median, lower=median)
        self.assertFalse((double.values != median).any())

    def test_dataframe_clip(self):
        # GH #2747
        df = DataFrame(np.random.randn(1000, 2))

        for lb, ub in [(-1, 1), (1, -1)]:
            clipped_df = df.clip(lb, ub)

            lb, ub = min(lb, ub), max(ub, lb)
            lb_mask = df.values <= lb
            ub_mask = df.values >= ub
            mask = ~lb_mask & ~ub_mask
            self.assertTrue((clipped_df.values[lb_mask] == lb).all())
            self.assertTrue((clipped_df.values[ub_mask] == ub).all())
            self.assertTrue((clipped_df.values[mask] ==
                             df.values[mask]).all())

    def test_clip_against_series(self):
        # GH #6966

        df = DataFrame(np.random.randn(1000, 2))
        lb = Series(np.random.randn(1000))
        ub = lb + 1

        clipped_df = df.clip(lb, ub, axis=0)

        for i in range(2):
            lb_mask = df.iloc[:, i] <= lb
            ub_mask = df.iloc[:, i] >= ub
            mask = ~lb_mask & ~ub_mask

            result = clipped_df.loc[lb_mask, i]
            assert_series_equal(result, lb[lb_mask], check_names=False)
            self.assertEqual(result.name, i)

            result = clipped_df.loc[ub_mask, i]
            assert_series_equal(result, ub[ub_mask], check_names=False)
            self.assertEqual(result.name, i)

            assert_series_equal(clipped_df.loc[mask, i], df.loc[mask, i])

    def test_clip_against_frame(self):
        df = DataFrame(np.random.randn(1000, 2))
        lb = DataFrame(np.random.randn(1000, 2))
        ub = lb + 1

        clipped_df = df.clip(lb, ub)

        lb_mask = df <= lb
        ub_mask = df >= ub
        mask = ~lb_mask & ~ub_mask

        assert_frame_equal(clipped_df[lb_mask], lb[lb_mask])
        assert_frame_equal(clipped_df[ub_mask], ub[ub_mask])
        assert_frame_equal(clipped_df[mask], df[mask])

    # Matrix-like

    def test_dot(self):
        a = DataFrame(np.random.randn(3, 4), index=['a', 'b', 'c'],
                      columns=['p', 'q', 'r', 's'])
        b = DataFrame(np.random.randn(4, 2), index=['p', 'q', 'r', 's'],
                      columns=['one', 'two'])

        result = a.dot(b)
        expected = DataFrame(np.dot(a.values, b.values),
                             index=['a', 'b', 'c'],
                             columns=['one', 'two'])
        # Check alignment
        b1 = b.reindex(index=reversed(b.index))
        result = a.dot(b)
        assert_frame_equal(result, expected)

        # Check series argument
        result = a.dot(b['one'])
        assert_series_equal(result, expected['one'], check_names=False)
        self.assertTrue(result.name is None)

        result = a.dot(b1['one'])
        assert_series_equal(result, expected['one'], check_names=False)
        self.assertTrue(result.name is None)

        # can pass correct-length arrays
        row = a.ix[0].values

        result = a.dot(row)
        exp = a.dot(a.ix[0])
        assert_series_equal(result, exp)

        with assertRaisesRegexp(ValueError, 'Dot product shape mismatch'):
            a.dot(row[:-1])

        a = np.random.rand(1, 5)
        b = np.random.rand(5, 1)
        A = DataFrame(a)

        # TODO(wesm): unused
        B = DataFrame(b)  # noqa

        # it works
        result = A.dot(b)

        # unaligned
        df = DataFrame(randn(3, 4), index=[1, 2, 3], columns=lrange(4))
        df2 = DataFrame(randn(5, 3), index=lrange(5), columns=[1, 2, 3])

        assertRaisesRegexp(ValueError, 'aligned', df.dot, df2)
