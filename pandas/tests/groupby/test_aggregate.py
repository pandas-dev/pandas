# -*- coding: utf-8 -*-
from __future__ import print_function
import nose
from datetime import datetime


from pandas import date_range
from pandas.core.index import MultiIndex
from pandas.core.api import DataFrame

from pandas.core.series import Series

from pandas.util.testing import (assert_frame_equal, assert_series_equal
                                 )

from pandas.core.groupby import (SpecificationError)
from pandas.compat import (lmap, OrderedDict)
from pandas.formats.printing import pprint_thing

from pandas import compat

import pandas.core.common as com
import numpy as np

import pandas.util.testing as tm
import pandas as pd


class TestGroupByAggregate(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.ts = tm.makeTimeSeries()

        self.seriesd = tm.getSeriesData()
        self.tsd = tm.getTimeSeriesData()
        self.frame = DataFrame(self.seriesd)
        self.tsframe = DataFrame(self.tsd)

        self.df = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.random.randn(8)})

        self.df_mixed_floats = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.array(
                 np.random.randn(8), dtype='float32')})

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'], ['one', 'two',
                                                                  'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        self.mframe = DataFrame(np.random.randn(10, 3), index=index,
                                columns=['A', 'B', 'C'])

        self.three_group = DataFrame(
            {'A': ['foo', 'foo', 'foo', 'foo', 'bar', 'bar', 'bar', 'bar',
                   'foo', 'foo', 'foo'],
             'B': ['one', 'one', 'one', 'two', 'one', 'one', 'one', 'two',
                   'two', 'two', 'one'],
             'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
                   'dull', 'shiny', 'shiny', 'shiny'],
             'D': np.random.randn(11),
             'E': np.random.randn(11),
             'F': np.random.randn(11)})

    def test_agg_api(self):

        # GH 6337
        # http://stackoverflow.com/questions/21706030/pandas-groupby-agg-function-column-dtype-error
        # different api for agg when passed custom function with mixed frame

        df = DataFrame({'data1': np.random.randn(5),
                        'data2': np.random.randn(5),
                        'key1': ['a', 'a', 'b', 'b', 'a'],
                        'key2': ['one', 'two', 'one', 'two', 'one']})
        grouped = df.groupby('key1')

        def peak_to_peak(arr):
            return arr.max() - arr.min()

        expected = grouped.agg([peak_to_peak])
        expected.columns = ['data1', 'data2']
        result = grouped.agg(peak_to_peak)
        assert_frame_equal(result, expected)

    def test_agg_regression1(self):
        grouped = self.tsframe.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.agg(np.mean)
        expected = grouped.mean()
        assert_frame_equal(result, expected)

    def test_agg_datetimes_mixed(self):
        data = [[1, '2012-01-01', 1.0], [2, '2012-01-02', 2.0], [3, None, 3.0]]

        df1 = DataFrame({'key': [x[0] for x in data],
                         'date': [x[1] for x in data],
                         'value': [x[2] for x in data]})

        data = [[row[0], datetime.strptime(row[1], '%Y-%m-%d').date() if row[1]
                 else None, row[2]] for row in data]

        df2 = DataFrame({'key': [x[0] for x in data],
                         'date': [x[1] for x in data],
                         'value': [x[2] for x in data]})

        df1['weights'] = df1['value'] / df1['value'].sum()
        gb1 = df1.groupby('date').aggregate(np.sum)

        df2['weights'] = df1['value'] / df1['value'].sum()
        gb2 = df2.groupby('date').aggregate(np.sum)

        assert (len(gb1) == len(gb2))

    def test_agg_period_index(self):
        from pandas import period_range, PeriodIndex
        prng = period_range('2012-1-1', freq='M', periods=3)
        df = DataFrame(np.random.randn(3, 2), index=prng)
        rs = df.groupby(level=0).sum()
        tm.assertIsInstance(rs.index, PeriodIndex)

        # GH 3579
        index = period_range(start='1999-01', periods=5, freq='M')
        s1 = Series(np.random.rand(len(index)), index=index)
        s2 = Series(np.random.rand(len(index)), index=index)
        series = [('s1', s1), ('s2', s2)]
        df = DataFrame.from_items(series)
        grouped = df.groupby(df.index.month)
        list(grouped)

    def test_agg_dict_parameter_cast_result_dtypes(self):
        # GH 12821

        df = DataFrame(
            {'class': ['A', 'A', 'B', 'B', 'C', 'C', 'D', 'D'],
             'time': date_range('1/1/2011', periods=8, freq='H')})
        df.loc[[0, 1, 2, 5], 'time'] = None

        # test for `first` function
        exp = df.loc[[0, 3, 4, 6]].set_index('class')
        grouped = df.groupby('class')
        assert_frame_equal(grouped.first(), exp)
        assert_frame_equal(grouped.agg('first'), exp)
        assert_frame_equal(grouped.agg({'time': 'first'}), exp)
        assert_series_equal(grouped.time.first(), exp['time'])
        assert_series_equal(grouped.time.agg('first'), exp['time'])

        # test for `last` function
        exp = df.loc[[0, 3, 4, 7]].set_index('class')
        grouped = df.groupby('class')
        assert_frame_equal(grouped.last(), exp)
        assert_frame_equal(grouped.agg('last'), exp)
        assert_frame_equal(grouped.agg({'time': 'last'}), exp)
        assert_series_equal(grouped.time.last(), exp['time'])
        assert_series_equal(grouped.time.agg('last'), exp['time'])

    def test_agg_must_agg(self):
        grouped = self.df.groupby('A')['C']
        self.assertRaises(Exception, grouped.agg, lambda x: x.describe())
        self.assertRaises(Exception, grouped.agg, lambda x: x.index[:2])

    def test_agg_ser_multi_key(self):
        # TODO(wesm): unused
        ser = self.df.C  # noqa

        f = lambda x: x.sum()
        results = self.df.C.groupby([self.df.A, self.df.B]).aggregate(f)
        expected = self.df.groupby(['A', 'B']).sum()['C']
        assert_series_equal(results, expected)

    def test_agg_apply_corner(self):
        # nothing to group, all NA
        grouped = self.ts.groupby(self.ts * np.nan)
        self.assertEqual(self.ts.dtype, np.float64)

        # groupby float64 values results in Float64Index
        exp = Series([], dtype=np.float64, index=pd.Index(
            [], dtype=np.float64))
        assert_series_equal(grouped.sum(), exp)
        assert_series_equal(grouped.agg(np.sum), exp)
        assert_series_equal(grouped.apply(np.sum), exp, check_index_type=False)

        # DataFrame
        grouped = self.tsframe.groupby(self.tsframe['A'] * np.nan)
        exp_df = DataFrame(columns=self.tsframe.columns, dtype=float,
                           index=pd.Index([], dtype=np.float64))
        assert_frame_equal(grouped.sum(), exp_df, check_names=False)
        assert_frame_equal(grouped.agg(np.sum), exp_df, check_names=False)
        assert_frame_equal(grouped.apply(np.sum), exp_df.iloc[:, :0],
                           check_names=False)

    def test_agg_grouping_is_list_tuple(self):
        from pandas.core.groupby import Grouping

        df = tm.makeTimeDataFrame()

        grouped = df.groupby(lambda x: x.year)
        grouper = grouped.grouper.groupings[0].grouper
        grouped.grouper.groupings[0] = Grouping(self.ts.index, list(grouper))

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

        grouped.grouper.groupings[0] = Grouping(self.ts.index, tuple(grouper))

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    def test_aggregate_api_consistency(self):
        # GH 9052
        # make sure that the aggregates via dict
        # are consistent

        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': np.random.randn(8) + 1.0,
                        'D': np.arange(8)})

        grouped = df.groupby(['A', 'B'])
        c_mean = grouped['C'].mean()
        c_sum = grouped['C'].sum()
        d_mean = grouped['D'].mean()
        d_sum = grouped['D'].sum()

        result = grouped['D'].agg(['sum', 'mean'])
        expected = pd.concat([d_sum, d_mean],
                             axis=1)
        expected.columns = ['sum', 'mean']
        assert_frame_equal(result, expected, check_like=True)

        result = grouped.agg([np.sum, np.mean])
        expected = pd.concat([c_sum,
                              c_mean,
                              d_sum,
                              d_mean],
                             axis=1)
        expected.columns = MultiIndex.from_product([['C', 'D'],
                                                    ['sum', 'mean']])
        assert_frame_equal(result, expected, check_like=True)

        result = grouped[['D', 'C']].agg([np.sum, np.mean])
        expected = pd.concat([d_sum,
                              d_mean,
                              c_sum,
                              c_mean],
                             axis=1)
        expected.columns = MultiIndex.from_product([['D', 'C'],
                                                    ['sum', 'mean']])
        assert_frame_equal(result, expected, check_like=True)

        result = grouped.agg({'C': 'mean', 'D': 'sum'})
        expected = pd.concat([d_sum,
                              c_mean],
                             axis=1)
        assert_frame_equal(result, expected, check_like=True)

        result = grouped.agg({'C': ['mean', 'sum'],
                              'D': ['mean', 'sum']})
        expected = pd.concat([c_mean,
                              c_sum,
                              d_mean,
                              d_sum],
                             axis=1)
        expected.columns = MultiIndex.from_product([['C', 'D'],
                                                    ['mean', 'sum']])

        result = grouped[['D', 'C']].agg({'r': np.sum,
                                          'r2': np.mean})
        expected = pd.concat([d_sum,
                              c_sum,
                              d_mean,
                              c_mean],
                             axis=1)
        expected.columns = MultiIndex.from_product([['r', 'r2'],
                                                    ['D', 'C']])
        assert_frame_equal(result, expected, check_like=True)

    def test_agg_compat(self):

        # GH 12334

        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': np.random.randn(8) + 1.0,
                        'D': np.arange(8)})

        g = df.groupby(['A', 'B'])

        expected = pd.concat([g['D'].sum(),
                              g['D'].std()],
                             axis=1)
        expected.columns = MultiIndex.from_tuples([('C', 'sum'),
                                                   ('C', 'std')])
        result = g['D'].agg({'C': ['sum', 'std']})
        assert_frame_equal(result, expected, check_like=True)

        expected = pd.concat([g['D'].sum(),
                              g['D'].std()],
                             axis=1)
        expected.columns = ['C', 'D']
        result = g['D'].agg({'C': 'sum', 'D': 'std'})
        assert_frame_equal(result, expected, check_like=True)

    def test_agg_nested_dicts(self):

        # API change for disallowing these types of nested dicts
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': np.random.randn(8) + 1.0,
                        'D': np.arange(8)})

        g = df.groupby(['A', 'B'])

        def f():
            g.aggregate({'r1': {'C': ['mean', 'sum']},
                         'r2': {'D': ['mean', 'sum']}})

        self.assertRaises(SpecificationError, f)

        result = g.agg({'C': {'ra': ['mean', 'std']},
                        'D': {'rb': ['mean', 'std']}})
        expected = pd.concat([g['C'].mean(), g['C'].std(), g['D'].mean(),
                              g['D'].std()], axis=1)
        expected.columns = pd.MultiIndex.from_tuples([('ra', 'mean'), (
            'ra', 'std'), ('rb', 'mean'), ('rb', 'std')])
        assert_frame_equal(result, expected, check_like=True)

        # same name as the original column
        # GH9052
        expected = g['D'].agg({'result1': np.sum, 'result2': np.mean})
        expected = expected.rename(columns={'result1': 'D'})
        result = g['D'].agg({'D': np.sum, 'result2': np.mean})
        assert_frame_equal(result, expected, check_like=True)

    def test_agg_python_multiindex(self):
        grouped = self.mframe.groupby(['A', 'B'])

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    def test_aggregate_str_func(self):
        def _check_results(grouped):
            # single series
            result = grouped['A'].agg('std')
            expected = grouped['A'].std()
            assert_series_equal(result, expected)

            # group frame by function name
            result = grouped.aggregate('var')
            expected = grouped.var()
            assert_frame_equal(result, expected)

            # group frame by function dict
            result = grouped.agg(OrderedDict([['A', 'var'], ['B', 'std'],
                                              ['C', 'mean'], ['D', 'sem']]))
            expected = DataFrame(OrderedDict([['A', grouped['A'].var(
            )], ['B', grouped['B'].std()], ['C', grouped['C'].mean()],
                ['D', grouped['D'].sem()]]))
            assert_frame_equal(result, expected)

        by_weekday = self.tsframe.groupby(lambda x: x.weekday())
        _check_results(by_weekday)

        by_mwkday = self.tsframe.groupby([lambda x: x.month,
                                          lambda x: x.weekday()])
        _check_results(by_mwkday)

    def test_aggregate_item_by_item(self):

        df = self.df.copy()
        df['E'] = ['a'] * len(self.df)
        grouped = self.df.groupby('A')

        # API change in 0.11
        # def aggfun(ser):
        #     return len(ser + 'a')
        # result = grouped.agg(aggfun)
        # self.assertEqual(len(result.columns), 1)

        aggfun = lambda ser: ser.size
        result = grouped.agg(aggfun)
        foo = (self.df.A == 'foo').sum()
        bar = (self.df.A == 'bar').sum()
        K = len(result.columns)

        # GH5782
        # odd comparisons can result here, so cast to make easy
        exp = pd.Series(np.array([foo] * K), index=list('BCD'),
                        dtype=np.float64, name='foo')
        tm.assert_series_equal(result.xs('foo'), exp)

        exp = pd.Series(np.array([bar] * K), index=list('BCD'),
                        dtype=np.float64, name='bar')
        tm.assert_almost_equal(result.xs('bar'), exp)

        def aggfun(ser):
            return ser.size

        result = DataFrame().groupby(self.df.A).agg(aggfun)
        tm.assertIsInstance(result, DataFrame)
        self.assertEqual(len(result), 0)

    def test_agg_item_by_item_raise_typeerror(self):
        from numpy.random import randint

        df = DataFrame(randint(10, size=(20, 10)))

        def raiseException(df):
            pprint_thing('----------------------------------------')
            pprint_thing(df.to_string())
            raise TypeError

        self.assertRaises(TypeError, df.groupby(0).agg, raiseException)

    def test_series_agg_multikey(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby([lambda x: x.year, lambda x: x.month])

        result = grouped.agg(np.sum)
        expected = grouped.sum()
        assert_series_equal(result, expected)

    def test_series_agg_multi_pure_python(self):
        data = DataFrame(
            {'A': ['foo', 'foo', 'foo', 'foo', 'bar', 'bar', 'bar', 'bar',
                   'foo', 'foo', 'foo'],
             'B': ['one', 'one', 'one', 'two', 'one', 'one', 'one', 'two',
                   'two', 'two', 'one'],
             'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
                   'dull', 'shiny', 'shiny', 'shiny'],
             'D': np.random.randn(11),
             'E': np.random.randn(11),
             'F': np.random.randn(11)})

        def bad(x):
            assert (len(x.base) > 0)
            return 'foo'

        result = data.groupby(['A', 'B']).agg(bad)
        expected = data.groupby(['A', 'B']).agg(lambda x: 'foo')
        assert_frame_equal(result, expected)


def assert_fp_equal(a, b):
    assert (np.abs(a - b) < 1e-12).all()


def _check_groupby(df, result, keys, field, f=lambda x: x.sum()):
    tups = lmap(tuple, df[keys].values)
    tups = com._asarray_tuplesafe(tups)
    expected = f(df.groupby(tups)[field])
    for k, v in compat.iteritems(expected):
        assert (result[k] == v)


def test_decons():
    from pandas.core.groupby import decons_group_index, get_group_index

    def testit(label_list, shape):
        group_index = get_group_index(label_list, shape, sort=True, xnull=True)
        label_list2 = decons_group_index(group_index, shape)

        for a, b in zip(label_list, label_list2):
            assert (np.array_equal(a, b))

    shape = (4, 5, 6)
    label_list = [np.tile([0, 1, 2, 3, 0, 1, 2, 3], 100), np.tile(
        [0, 2, 4, 3, 0, 1, 2, 3], 100), np.tile(
            [5, 1, 0, 2, 3, 0, 5, 4], 100)]
    testit(label_list, shape)

    shape = (10000, 10000)
    label_list = [np.tile(np.arange(10000), 5), np.tile(np.arange(10000), 5)]
    testit(label_list, shape)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure', '-s'
                         ], exit=False)
