# -*- coding: utf-8 -*-
from __future__ import print_function
import nose

from datetime import datetime
from numpy import nan

from pandas.types.common import _ensure_platform_int
from pandas import date_range, bdate_range, Timestamp, isnull
from pandas.core.index import Index, MultiIndex, CategoricalIndex
from pandas.core.api import Categorical, DataFrame
from pandas.core.common import UnsupportedFunctionCall
from pandas.core.groupby import (SpecificationError, DataError, _nargsort,
                                 _lexsort_indexer)
from pandas.core.series import Series
from pandas.core.config import option_context
from pandas.formats.printing import pprint_thing
from pandas.util.testing import (assert_panel_equal, assert_frame_equal,
                                 assert_series_equal, assert_almost_equal,
                                 assert_index_equal, assertRaisesRegexp)
from pandas.compat import (range, long, lrange, StringIO, lmap, lzip, map, zip,
                           builtins, OrderedDict, product as cart_product)
from pandas import compat
from pandas.core.panel import Panel
from pandas.tools.merge import concat
from collections import defaultdict
from functools import partial
import pandas.core.common as com
import numpy as np

import pandas.core.nanops as nanops

import pandas.util.testing as tm
import pandas as pd


class TestGroupBy(tm.TestCase):

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

    def test_basic(self):
        def checkit(dtype):
            data = Series(np.arange(9) // 3, index=np.arange(9), dtype=dtype)

            index = np.arange(9)
            np.random.shuffle(index)
            data = data.reindex(index)

            grouped = data.groupby(lambda x: x // 3)

            for k, v in grouped:
                self.assertEqual(len(v), 3)

            agged = grouped.aggregate(np.mean)
            self.assertEqual(agged[1], 1)

            assert_series_equal(agged, grouped.agg(np.mean))  # shorthand
            assert_series_equal(agged, grouped.mean())
            assert_series_equal(grouped.agg(np.sum), grouped.sum())

            expected = grouped.apply(lambda x: x * x.sum())
            transformed = grouped.transform(lambda x: x * x.sum())
            self.assertEqual(transformed[7], 12)
            assert_series_equal(transformed, expected)

            value_grouped = data.groupby(data)
            assert_series_equal(value_grouped.aggregate(np.mean), agged,
                                check_index_type=False)

            # complex agg
            agged = grouped.aggregate([np.mean, np.std])
            agged = grouped.aggregate({'one': np.mean, 'two': np.std})

            group_constants = {0: 10, 1: 20, 2: 30}
            agged = grouped.agg(lambda x: group_constants[x.name] + x.mean())
            self.assertEqual(agged[1], 21)

            # corner cases
            self.assertRaises(Exception, grouped.aggregate, lambda x: x * 2)

        for dtype in ['int64', 'int32', 'float64', 'float32']:
            checkit(dtype)

    def test_select_bad_cols(self):
        df = DataFrame([[1, 2]], columns=['A', 'B'])
        g = df.groupby('A')
        self.assertRaises(KeyError, g.__getitem__, ['C'])  # g[['C']]

        self.assertRaises(KeyError, g.__getitem__, ['A', 'C'])  # g[['A', 'C']]
        with assertRaisesRegexp(KeyError, '^[^A]+$'):
            # A should not be referenced as a bad column...
            # will have to rethink regex if you change message!
            g[['A', 'C']]

    def test_first_last_nth(self):
        # tests for first / last / nth
        grouped = self.df.groupby('A')
        first = grouped.first()
        expected = self.df.ix[[1, 0], ['B', 'C', 'D']]
        expected.index = Index(['bar', 'foo'], name='A')
        expected = expected.sort_index()
        assert_frame_equal(first, expected)

        nth = grouped.nth(0)
        assert_frame_equal(nth, expected)

        last = grouped.last()
        expected = self.df.ix[[5, 7], ['B', 'C', 'D']]
        expected.index = Index(['bar', 'foo'], name='A')
        assert_frame_equal(last, expected)

        nth = grouped.nth(-1)
        assert_frame_equal(nth, expected)

        nth = grouped.nth(1)
        expected = self.df.ix[[2, 3], ['B', 'C', 'D']].copy()
        expected.index = Index(['foo', 'bar'], name='A')
        expected = expected.sort_index()
        assert_frame_equal(nth, expected)

        # it works!
        grouped['B'].first()
        grouped['B'].last()
        grouped['B'].nth(0)

        self.df.loc[self.df['A'] == 'foo', 'B'] = np.nan
        self.assertTrue(isnull(grouped['B'].first()['foo']))
        self.assertTrue(isnull(grouped['B'].last()['foo']))
        self.assertTrue(isnull(grouped['B'].nth(0)['foo']))

        # v0.14.0 whatsnew
        df = DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
        g = df.groupby('A')
        result = g.first()
        expected = df.iloc[[1, 2]].set_index('A')
        assert_frame_equal(result, expected)

        expected = df.iloc[[1, 2]].set_index('A')
        result = g.nth(0, dropna='any')
        assert_frame_equal(result, expected)

    def test_first_last_nth_dtypes(self):

        df = self.df_mixed_floats.copy()
        df['E'] = True
        df['F'] = 1

        # tests for first / last / nth
        grouped = df.groupby('A')
        first = grouped.first()
        expected = df.ix[[1, 0], ['B', 'C', 'D', 'E', 'F']]
        expected.index = Index(['bar', 'foo'], name='A')
        expected = expected.sort_index()
        assert_frame_equal(first, expected)

        last = grouped.last()
        expected = df.ix[[5, 7], ['B', 'C', 'D', 'E', 'F']]
        expected.index = Index(['bar', 'foo'], name='A')
        expected = expected.sort_index()
        assert_frame_equal(last, expected)

        nth = grouped.nth(1)
        expected = df.ix[[3, 2], ['B', 'C', 'D', 'E', 'F']]
        expected.index = Index(['bar', 'foo'], name='A')
        expected = expected.sort_index()
        assert_frame_equal(nth, expected)

        # GH 2763, first/last shifting dtypes
        idx = lrange(10)
        idx.append(9)
        s = Series(data=lrange(11), index=idx, name='IntCol')
        self.assertEqual(s.dtype, 'int64')
        f = s.groupby(level=0).first()
        self.assertEqual(f.dtype, 'int64')

    def test_nth(self):
        df = DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
        g = df.groupby('A')

        assert_frame_equal(g.nth(0), df.iloc[[0, 2]].set_index('A'))
        assert_frame_equal(g.nth(1), df.iloc[[1]].set_index('A'))
        assert_frame_equal(g.nth(2), df.loc[[]].set_index('A'))
        assert_frame_equal(g.nth(-1), df.iloc[[1, 2]].set_index('A'))
        assert_frame_equal(g.nth(-2), df.iloc[[0]].set_index('A'))
        assert_frame_equal(g.nth(-3), df.loc[[]].set_index('A'))
        assert_series_equal(g.B.nth(0), df.set_index('A').B.iloc[[0, 2]])
        assert_series_equal(g.B.nth(1), df.set_index('A').B.iloc[[1]])
        assert_frame_equal(g[['B']].nth(0),
                           df.ix[[0, 2], ['A', 'B']].set_index('A'))

        exp = df.set_index('A')
        assert_frame_equal(g.nth(0, dropna='any'), exp.iloc[[1, 2]])
        assert_frame_equal(g.nth(-1, dropna='any'), exp.iloc[[1, 2]])

        exp['B'] = np.nan
        assert_frame_equal(g.nth(7, dropna='any'), exp.iloc[[1, 2]])
        assert_frame_equal(g.nth(2, dropna='any'), exp.iloc[[1, 2]])

        # out of bounds, regression from 0.13.1
        # GH 6621
        df = DataFrame({'color': {0: 'green',
                                  1: 'green',
                                  2: 'red',
                                  3: 'red',
                                  4: 'red'},
                        'food': {0: 'ham',
                                 1: 'eggs',
                                 2: 'eggs',
                                 3: 'ham',
                                 4: 'pork'},
                        'two': {0: 1.5456590000000001,
                                1: -0.070345000000000005,
                                2: -2.4004539999999999,
                                3: 0.46206000000000003,
                                4: 0.52350799999999997},
                        'one': {0: 0.56573799999999996,
                                1: -0.9742360000000001,
                                2: 1.033801,
                                3: -0.78543499999999999,
                                4: 0.70422799999999997}}).set_index(['color',
                                                                     'food'])

        result = df.groupby(level=0, as_index=False).nth(2)
        expected = df.iloc[[-1]]
        assert_frame_equal(result, expected)

        result = df.groupby(level=0, as_index=False).nth(3)
        expected = df.loc[[]]
        assert_frame_equal(result, expected)

        # GH 7559
        # from the vbench
        df = DataFrame(np.random.randint(1, 10, (100, 2)), dtype='int64')
        s = df[1]
        g = df[0]
        expected = s.groupby(g).first()
        expected2 = s.groupby(g).apply(lambda x: x.iloc[0])
        assert_series_equal(expected2, expected, check_names=False)
        self.assertTrue(expected.name, 0)
        self.assertEqual(expected.name, 1)

        # validate first
        v = s[g == 1].iloc[0]
        self.assertEqual(expected.iloc[0], v)
        self.assertEqual(expected2.iloc[0], v)

        # this is NOT the same as .first (as sorted is default!)
        # as it keeps the order in the series (and not the group order)
        # related GH 7287
        expected = s.groupby(g, sort=False).first()
        result = s.groupby(g, sort=False).nth(0, dropna='all')
        assert_series_equal(result, expected)

        # doc example
        df = DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
        g = df.groupby('A')
        result = g.B.nth(0, dropna=True)
        expected = g.B.first()
        assert_series_equal(result, expected)

        # test multiple nth values
        df = DataFrame([[1, np.nan], [1, 3], [1, 4], [5, 6], [5, 7]],
                       columns=['A', 'B'])
        g = df.groupby('A')

        assert_frame_equal(g.nth(0), df.iloc[[0, 3]].set_index('A'))
        assert_frame_equal(g.nth([0]), df.iloc[[0, 3]].set_index('A'))
        assert_frame_equal(g.nth([0, 1]), df.iloc[[0, 1, 3, 4]].set_index('A'))
        assert_frame_equal(
            g.nth([0, -1]), df.iloc[[0, 2, 3, 4]].set_index('A'))
        assert_frame_equal(
            g.nth([0, 1, 2]), df.iloc[[0, 1, 2, 3, 4]].set_index('A'))
        assert_frame_equal(
            g.nth([0, 1, -1]), df.iloc[[0, 1, 2, 3, 4]].set_index('A'))
        assert_frame_equal(g.nth([2]), df.iloc[[2]].set_index('A'))
        assert_frame_equal(g.nth([3, 4]), df.loc[[]].set_index('A'))

        business_dates = pd.date_range(start='4/1/2014', end='6/30/2014',
                                       freq='B')
        df = DataFrame(1, index=business_dates, columns=['a', 'b'])
        # get the first, fourth and last two business days for each month
        key = (df.index.year, df.index.month)
        result = df.groupby(key, as_index=False).nth([0, 3, -2, -1])
        expected_dates = pd.to_datetime(
            ['2014/4/1', '2014/4/4', '2014/4/29', '2014/4/30', '2014/5/1',
             '2014/5/6', '2014/5/29', '2014/5/30', '2014/6/2', '2014/6/5',
             '2014/6/27', '2014/6/30'])
        expected = DataFrame(1, columns=['a', 'b'], index=expected_dates)
        assert_frame_equal(result, expected)

    def test_nth_multi_index(self):
        # PR 9090, related to issue 8979
        # test nth on MultiIndex, should match .first()
        grouped = self.three_group.groupby(['A', 'B'])
        result = grouped.nth(0)
        expected = grouped.first()
        assert_frame_equal(result, expected)

    def test_nth_multi_index_as_expected(self):
        # PR 9090, related to issue 8979
        # test nth on MultiIndex
        three_group = DataFrame(
            {'A': ['foo', 'foo', 'foo', 'foo', 'bar', 'bar', 'bar', 'bar',
                   'foo', 'foo', 'foo'],
             'B': ['one', 'one', 'one', 'two', 'one', 'one', 'one', 'two',
                   'two', 'two', 'one'],
             'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
                   'dull', 'shiny', 'shiny', 'shiny']})
        grouped = three_group.groupby(['A', 'B'])
        result = grouped.nth(0)
        expected = DataFrame(
            {'C': ['dull', 'dull', 'dull', 'dull']},
            index=MultiIndex.from_arrays([['bar', 'bar', 'foo', 'foo'],
                                          ['one', 'two', 'one', 'two']],
                                         names=['A', 'B']))
        assert_frame_equal(result, expected)

    def test_group_selection_cache(self):
        # GH 12839 nth, head, and tail should return same result consistently
        df = DataFrame([[1, 2], [1, 4], [5, 6]], columns=['A', 'B'])
        expected = df.iloc[[0, 2]].set_index('A')

        g = df.groupby('A')
        result1 = g.head(n=2)
        result2 = g.nth(0)
        assert_frame_equal(result1, df)
        assert_frame_equal(result2, expected)

        g = df.groupby('A')
        result1 = g.tail(n=2)
        result2 = g.nth(0)
        assert_frame_equal(result1, df)
        assert_frame_equal(result2, expected)

        g = df.groupby('A')
        result1 = g.nth(0)
        result2 = g.head(n=2)
        assert_frame_equal(result1, expected)
        assert_frame_equal(result2, df)

        g = df.groupby('A')
        result1 = g.nth(0)
        result2 = g.tail(n=2)
        assert_frame_equal(result1, expected)
        assert_frame_equal(result2, df)

    def test_grouper_index_types(self):
        # related GH5375
        # groupby misbehaving when using a Floatlike index
        df = DataFrame(np.arange(10).reshape(5, 2), columns=list('AB'))
        for index in [tm.makeFloatIndex, tm.makeStringIndex,
                      tm.makeUnicodeIndex, tm.makeIntIndex, tm.makeDateIndex,
                      tm.makePeriodIndex]:

            df.index = index(len(df))
            df.groupby(list('abcde')).apply(lambda x: x)

            df.index = list(reversed(df.index.tolist()))
            df.groupby(list('abcde')).apply(lambda x: x)

    def test_grouper_multilevel_freq(self):

        # GH 7885
        # with level and freq specified in a pd.Grouper
        from datetime import date, timedelta
        d0 = date.today() - timedelta(days=14)
        dates = date_range(d0, date.today())
        date_index = pd.MultiIndex.from_product(
            [dates, dates], names=['foo', 'bar'])
        df = pd.DataFrame(np.random.randint(0, 100, 225), index=date_index)

        # Check string level
        expected = df.reset_index().groupby([pd.Grouper(
            key='foo', freq='W'), pd.Grouper(key='bar', freq='W')]).sum()
        # reset index changes columns dtype to object
        expected.columns = pd.Index([0], dtype='int64')

        result = df.groupby([pd.Grouper(level='foo', freq='W'), pd.Grouper(
            level='bar', freq='W')]).sum()
        assert_frame_equal(result, expected)

        # Check integer level
        result = df.groupby([pd.Grouper(level=0, freq='W'), pd.Grouper(
            level=1, freq='W')]).sum()
        assert_frame_equal(result, expected)

    def test_grouper_creation_bug(self):

        # GH 8795
        df = DataFrame({'A': [0, 0, 1, 1, 2, 2], 'B': [1, 2, 3, 4, 5, 6]})
        g = df.groupby('A')
        expected = g.sum()

        g = df.groupby(pd.Grouper(key='A'))
        result = g.sum()
        assert_frame_equal(result, expected)

        result = g.apply(lambda x: x.sum())
        assert_frame_equal(result, expected)

        g = df.groupby(pd.Grouper(key='A', axis=0))
        result = g.sum()
        assert_frame_equal(result, expected)

        # GH14334
        # pd.Grouper(key=...) may be passed in a list
        df = DataFrame({'A': [0, 0, 0, 1, 1, 1],
                        'B': [1, 1, 2, 2, 3, 3],
                        'C': [1, 2, 3, 4, 5, 6]})
        # Group by single column
        expected = df.groupby('A').sum()
        g = df.groupby([pd.Grouper(key='A')])
        result = g.sum()
        assert_frame_equal(result, expected)

        # Group by two columns
        # using a combination of strings and Grouper objects
        expected = df.groupby(['A', 'B']).sum()

        # Group with two Grouper objects
        g = df.groupby([pd.Grouper(key='A'), pd.Grouper(key='B')])
        result = g.sum()
        assert_frame_equal(result, expected)

        # Group with a string and a Grouper object
        g = df.groupby(['A', pd.Grouper(key='B')])
        result = g.sum()
        assert_frame_equal(result, expected)

        # Group with a Grouper object and a string
        g = df.groupby([pd.Grouper(key='A'), 'B'])
        result = g.sum()
        assert_frame_equal(result, expected)

        # GH8866
        s = Series(np.arange(8, dtype='int64'),
                   index=pd.MultiIndex.from_product(
                       [list('ab'), range(2),
                        date_range('20130101', periods=2)],
                       names=['one', 'two', 'three']))
        result = s.groupby(pd.Grouper(level='three', freq='M')).sum()
        expected = Series([28], index=Index(
            [Timestamp('2013-01-31')], freq='M', name='three'))
        assert_series_equal(result, expected)

        # just specifying a level breaks
        result = s.groupby(pd.Grouper(level='one')).sum()
        expected = s.groupby(level='one').sum()
        assert_series_equal(result, expected)

    def test_grouper_column_and_index(self):
        # GH 14327

        # Grouping a multi-index frame by a column and an index level should
        # be equivalent to resetting the index and grouping by two columns
        idx = pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('a', 3),
                                         ('b', 1), ('b', 2), ('b', 3)])
        idx.names = ['outer', 'inner']
        df_multi = pd.DataFrame({"A": np.arange(6),
                                 'B': ['one', 'one', 'two',
                                       'two', 'one', 'one']},
                                index=idx)
        result = df_multi.groupby(['B', pd.Grouper(level='inner')]).mean()
        expected = df_multi.reset_index().groupby(['B', 'inner']).mean()
        assert_frame_equal(result, expected)

        # Test the reverse grouping order
        result = df_multi.groupby([pd.Grouper(level='inner'), 'B']).mean()
        expected = df_multi.reset_index().groupby(['inner', 'B']).mean()
        assert_frame_equal(result, expected)

        # Grouping a single-index frame by a column and the index should
        # be equivalent to resetting the index and grouping by two columns
        df_single = df_multi.reset_index('outer')
        result = df_single.groupby(['B', pd.Grouper(level='inner')]).mean()
        expected = df_single.reset_index().groupby(['B', 'inner']).mean()
        assert_frame_equal(result, expected)

        # Test the reverse grouping order
        result = df_single.groupby([pd.Grouper(level='inner'), 'B']).mean()
        expected = df_single.reset_index().groupby(['inner', 'B']).mean()
        assert_frame_equal(result, expected)

    def test_grouper_getting_correct_binner(self):

        # GH 10063
        # using a non-time-based grouper and a time-based grouper
        # and specifying levels
        df = DataFrame({'A': 1}, index=pd.MultiIndex.from_product(
            [list('ab'), date_range('20130101', periods=80)], names=['one',
                                                                     'two']))
        result = df.groupby([pd.Grouper(level='one'), pd.Grouper(
            level='two', freq='M')]).sum()
        expected = DataFrame({'A': [31, 28, 21, 31, 28, 21]},
                             index=MultiIndex.from_product(
                                 [list('ab'),
                                  date_range('20130101', freq='M', periods=3)],
                                 names=['one', 'two']))
        assert_frame_equal(result, expected)

    def test_grouper_iter(self):
        self.assertEqual(sorted(self.df.groupby('A').grouper), ['bar', 'foo'])

    def test_empty_groups(self):
        # GH # 1048
        self.assertRaises(ValueError, self.df.groupby, [])

    def test_groupby_grouper(self):
        grouped = self.df.groupby('A')

        result = self.df.groupby(grouped.grouper).mean()
        expected = grouped.mean()
        assert_frame_equal(result, expected)

    def test_groupby_duplicated_column_errormsg(self):
        # GH7511
        df = DataFrame(columns=['A', 'B', 'A', 'C'],
                       data=[range(4), range(2, 6), range(0, 8, 2)])

        self.assertRaises(ValueError, df.groupby, 'A')
        self.assertRaises(ValueError, df.groupby, ['A', 'B'])

        grouped = df.groupby('B')
        c = grouped.count()
        self.assertTrue(c.columns.nlevels == 1)
        self.assertTrue(c.columns.size == 3)

    def test_groupby_dict_mapping(self):
        # GH #679
        from pandas import Series
        s = Series({'T1': 5})
        result = s.groupby({'T1': 'T2'}).agg(sum)
        expected = s.groupby(['T2']).agg(sum)
        assert_series_equal(result, expected)

        s = Series([1., 2., 3., 4.], index=list('abcd'))
        mapping = {'a': 0, 'b': 0, 'c': 1, 'd': 1}

        result = s.groupby(mapping).mean()
        result2 = s.groupby(mapping).agg(np.mean)
        expected = s.groupby([0, 0, 1, 1]).mean()
        expected2 = s.groupby([0, 0, 1, 1]).mean()
        assert_series_equal(result, expected)
        assert_series_equal(result, result2)
        assert_series_equal(result, expected2)

    def test_groupby_grouper_f_sanity_checked(self):
        dates = date_range('01-Jan-2013', periods=12, freq='MS')
        ts = Series(np.random.randn(12), index=dates)

        # GH3035
        # index.map is used to apply grouper to the index
        # if it fails on the elements, map tries it on the entire index as
        # a sequence. That can yield invalid results that cause trouble
        # down the line.
        # the surprise comes from using key[0:6] rather then str(key)[0:6]
        # when the elements are Timestamp.
        # the result is Index[0:6], very confusing.

        self.assertRaises(AssertionError, ts.groupby, lambda key: key[0:6])

    def test_groupby_nonobject_dtype(self):
        key = self.mframe.index.labels[0]
        grouped = self.mframe.groupby(key)
        result = grouped.sum()

        expected = self.mframe.groupby(key.astype('O')).sum()
        assert_frame_equal(result, expected)

        # GH 3911, mixed frame non-conversion
        df = self.df_mixed_floats.copy()
        df['value'] = lrange(len(df))

        def max_value(group):
            return group.ix[group['value'].idxmax()]

        applied = df.groupby('A').apply(max_value)
        result = applied.get_dtype_counts().sort_values()
        expected = Series({'object': 2,
                           'float64': 2,
                           'int64': 1}).sort_values()
        assert_series_equal(result, expected)

    def test_groupby_return_type(self):

        # GH2893, return a reduced type
        df1 = DataFrame([{"val1": 1,
                          "val2": 20}, {"val1": 1,
                                        "val2": 19}, {"val1": 2,
                                                      "val2": 27}, {"val1": 2,
                                                                    "val2": 12}
                         ])

        def func(dataf):
            return dataf["val2"] - dataf["val2"].mean()

        result = df1.groupby("val1", squeeze=True).apply(func)
        tm.assertIsInstance(result, Series)

        df2 = DataFrame([{"val1": 1,
                          "val2": 20}, {"val1": 1,
                                        "val2": 19}, {"val1": 1,
                                                      "val2": 27}, {"val1": 1,
                                                                    "val2": 12}
                         ])

        def func(dataf):
            return dataf["val2"] - dataf["val2"].mean()

        result = df2.groupby("val1", squeeze=True).apply(func)
        tm.assertIsInstance(result, Series)

        # GH3596, return a consistent type (regression in 0.11 from 0.10.1)
        df = DataFrame([[1, 1], [1, 1]], columns=['X', 'Y'])
        result = df.groupby('X', squeeze=False).count()
        tm.assertIsInstance(result, DataFrame)

        # GH5592
        # inconcistent return type
        df = DataFrame(dict(A=['Tiger', 'Tiger', 'Tiger', 'Lamb', 'Lamb',
                               'Pony', 'Pony'], B=Series(
                                   np.arange(7), dtype='int64'), C=date_range(
                                       '20130101', periods=7)))

        def f(grp):
            return grp.iloc[0]

        expected = df.groupby('A').first()[['B']]
        result = df.groupby('A').apply(f)[['B']]
        assert_frame_equal(result, expected)

        def f(grp):
            if grp.name == 'Tiger':
                return None
            return grp.iloc[0]

        result = df.groupby('A').apply(f)[['B']]
        e = expected.copy()
        e.loc['Tiger'] = np.nan
        assert_frame_equal(result, e)

        def f(grp):
            if grp.name == 'Pony':
                return None
            return grp.iloc[0]

        result = df.groupby('A').apply(f)[['B']]
        e = expected.copy()
        e.loc['Pony'] = np.nan
        assert_frame_equal(result, e)

        # 5592 revisited, with datetimes
        def f(grp):
            if grp.name == 'Pony':
                return None
            return grp.iloc[0]

        result = df.groupby('A').apply(f)[['C']]
        e = df.groupby('A').first()[['C']]
        e.loc['Pony'] = pd.NaT
        assert_frame_equal(result, e)

        # scalar outputs
        def f(grp):
            if grp.name == 'Pony':
                return None
            return grp.iloc[0].loc['C']

        result = df.groupby('A').apply(f)
        e = df.groupby('A').first()['C'].copy()
        e.loc['Pony'] = np.nan
        e.name = None
        assert_series_equal(result, e)

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

    def test_get_group(self):
        wp = tm.makePanel()
        grouped = wp.groupby(lambda x: x.month, axis='major')

        gp = grouped.get_group(1)
        expected = wp.reindex(major=[x for x in wp.major_axis if x.month == 1])
        assert_panel_equal(gp, expected)

        # GH 5267
        # be datelike friendly
        df = DataFrame({'DATE': pd.to_datetime(
            ['10-Oct-2013', '10-Oct-2013', '10-Oct-2013', '11-Oct-2013',
             '11-Oct-2013', '11-Oct-2013']),
            'label': ['foo', 'foo', 'bar', 'foo', 'foo', 'bar'],
            'VAL': [1, 2, 3, 4, 5, 6]})

        g = df.groupby('DATE')
        key = list(g.groups)[0]
        result1 = g.get_group(key)
        result2 = g.get_group(Timestamp(key).to_pydatetime())
        result3 = g.get_group(str(Timestamp(key)))
        assert_frame_equal(result1, result2)
        assert_frame_equal(result1, result3)

        g = df.groupby(['DATE', 'label'])

        key = list(g.groups)[0]
        result1 = g.get_group(key)
        result2 = g.get_group((Timestamp(key[0]).to_pydatetime(), key[1]))
        result3 = g.get_group((str(Timestamp(key[0])), key[1]))
        assert_frame_equal(result1, result2)
        assert_frame_equal(result1, result3)

        # must pass a same-length tuple with multiple keys
        self.assertRaises(ValueError, lambda: g.get_group('foo'))
        self.assertRaises(ValueError, lambda: g.get_group(('foo')))
        self.assertRaises(ValueError,
                          lambda: g.get_group(('foo', 'bar', 'baz')))

    def test_get_group_empty_bins(self):
        d = pd.DataFrame([3, 1, 7, 6])
        bins = [0, 5, 10, 15]
        g = d.groupby(pd.cut(d[0], bins))

        result = g.get_group('(0, 5]')
        expected = DataFrame([3, 1], index=[0, 1])
        assert_frame_equal(result, expected)

        self.assertRaises(KeyError, lambda: g.get_group('(10, 15]'))

    def test_get_group_grouped_by_tuple(self):
        # GH 8121
        df = DataFrame([[(1, ), (1, 2), (1, ), (1, 2)]], index=['ids']).T
        gr = df.groupby('ids')
        expected = DataFrame({'ids': [(1, ), (1, )]}, index=[0, 2])
        result = gr.get_group((1, ))
        assert_frame_equal(result, expected)

        dt = pd.to_datetime(['2010-01-01', '2010-01-02', '2010-01-01',
                             '2010-01-02'])
        df = DataFrame({'ids': [(x, ) for x in dt]})
        gr = df.groupby('ids')
        result = gr.get_group(('2010-01-01', ))
        expected = DataFrame({'ids': [(dt[0], ), (dt[0], )]}, index=[0, 2])
        assert_frame_equal(result, expected)

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

    def test_grouping_error_on_multidim_input(self):
        from pandas.core.groupby import Grouping
        self.assertRaises(ValueError,
                          Grouping, self.df.index, self.df[['A', 'A']])

    def test_agg_python_multiindex(self):
        grouped = self.mframe.groupby(['A', 'B'])

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    def test_apply_describe_bug(self):
        grouped = self.mframe.groupby(level='first')
        grouped.describe()  # it works!

    def test_apply_issues(self):
        # GH 5788

        s = """2011.05.16,00:00,1.40893
2011.05.16,01:00,1.40760
2011.05.16,02:00,1.40750
2011.05.16,03:00,1.40649
2011.05.17,02:00,1.40893
2011.05.17,03:00,1.40760
2011.05.17,04:00,1.40750
2011.05.17,05:00,1.40649
2011.05.18,02:00,1.40893
2011.05.18,03:00,1.40760
2011.05.18,04:00,1.40750
2011.05.18,05:00,1.40649"""

        df = pd.read_csv(
            StringIO(s), header=None, names=['date', 'time', 'value'],
            parse_dates=[['date', 'time']])
        df = df.set_index('date_time')

        expected = df.groupby(df.index.date).idxmax()
        result = df.groupby(df.index.date).apply(lambda x: x.idxmax())
        assert_frame_equal(result, expected)

        # GH 5789
        # don't auto coerce dates
        df = pd.read_csv(
            StringIO(s), header=None, names=['date', 'time', 'value'])
        exp_idx = pd.Index(
            ['2011.05.16', '2011.05.17', '2011.05.18'
             ], dtype=object, name='date')
        expected = Series(['00:00', '02:00', '02:00'], index=exp_idx)
        result = df.groupby('date').apply(
            lambda x: x['time'][x['value'].idxmax()])
        assert_series_equal(result, expected)

    def test_time_field_bug(self):
        # Test a fix for the following error related to GH issue 11324 When
        # non-key fields in a group-by dataframe contained time-based fields
        # that were not returned by the apply function, an exception would be
        # raised.

        df = pd.DataFrame({'a': 1, 'b': [datetime.now() for nn in range(10)]})

        def func_with_no_date(batch):
            return pd.Series({'c': 2})

        def func_with_date(batch):
            return pd.Series({'c': 2, 'b': datetime(2015, 1, 1)})

        dfg_no_conversion = df.groupby(by=['a']).apply(func_with_no_date)
        dfg_no_conversion_expected = pd.DataFrame({'c': 2}, index=[1])
        dfg_no_conversion_expected.index.name = 'a'

        dfg_conversion = df.groupby(by=['a']).apply(func_with_date)
        dfg_conversion_expected = pd.DataFrame(
            {'b': datetime(2015, 1, 1),
             'c': 2}, index=[1])
        dfg_conversion_expected.index.name = 'a'

        self.assert_frame_equal(dfg_no_conversion, dfg_no_conversion_expected)
        self.assert_frame_equal(dfg_conversion, dfg_conversion_expected)

    def test_len(self):
        df = tm.makeTimeDataFrame()
        grouped = df.groupby([lambda x: x.year, lambda x: x.month,
                              lambda x: x.day])
        self.assertEqual(len(grouped), len(df))

        grouped = df.groupby([lambda x: x.year, lambda x: x.month])
        expected = len(set([(x.year, x.month) for x in df.index]))
        self.assertEqual(len(grouped), expected)

        # issue 11016
        df = pd.DataFrame(dict(a=[np.nan] * 3, b=[1, 2, 3]))
        self.assertEqual(len(df.groupby(('a'))), 0)
        self.assertEqual(len(df.groupby(('b'))), 3)
        self.assertEqual(len(df.groupby(('a', 'b'))), 3)

    def test_groups(self):
        grouped = self.df.groupby(['A'])
        groups = grouped.groups
        self.assertIs(groups, grouped.groups)  # caching works

        for k, v in compat.iteritems(grouped.groups):
            self.assertTrue((self.df.ix[v]['A'] == k).all())

        grouped = self.df.groupby(['A', 'B'])
        groups = grouped.groups
        self.assertIs(groups, grouped.groups)  # caching works
        for k, v in compat.iteritems(grouped.groups):
            self.assertTrue((self.df.ix[v]['A'] == k[0]).all())
            self.assertTrue((self.df.ix[v]['B'] == k[1]).all())

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

    def test_basic_regression(self):
        # regression
        T = [1.0 * x for x in lrange(1, 10) * 10][:1095]
        result = Series(T, lrange(0, len(T)))

        groupings = np.random.random((1100, ))
        groupings = Series(groupings, lrange(0, len(groupings))) * 10.

        grouped = result.groupby(groupings)
        grouped.mean()

    def test_transform(self):
        data = Series(np.arange(9) // 3, index=np.arange(9))

        index = np.arange(9)
        np.random.shuffle(index)
        data = data.reindex(index)

        grouped = data.groupby(lambda x: x // 3)

        transformed = grouped.transform(lambda x: x * x.sum())
        self.assertEqual(transformed[7], 12)

        # GH 8046
        # make sure that we preserve the input order

        df = DataFrame(
            np.arange(6, dtype='int64').reshape(
                3, 2), columns=["a", "b"], index=[0, 2, 1])
        key = [0, 0, 1]
        expected = df.sort_index().groupby(key).transform(
            lambda x: x - x.mean()).groupby(key).mean()
        result = df.groupby(key).transform(lambda x: x - x.mean()).groupby(
            key).mean()
        assert_frame_equal(result, expected)

        def demean(arr):
            return arr - arr.mean()

        people = DataFrame(np.random.randn(5, 5),
                           columns=['a', 'b', 'c', 'd', 'e'],
                           index=['Joe', 'Steve', 'Wes', 'Jim', 'Travis'])
        key = ['one', 'two', 'one', 'two', 'one']
        result = people.groupby(key).transform(demean).groupby(key).mean()
        expected = people.groupby(key).apply(demean).groupby(key).mean()
        assert_frame_equal(result, expected)

        # GH 8430
        df = tm.makeTimeDataFrame()
        g = df.groupby(pd.TimeGrouper('M'))
        g.transform(lambda x: x - 1)

        # GH 9700
        df = DataFrame({'a': range(5, 10), 'b': range(5)})
        result = df.groupby('a').transform(max)
        expected = DataFrame({'b': range(5)})
        tm.assert_frame_equal(result, expected)

    def test_transform_fast(self):

        df = DataFrame({'id': np.arange(100000) / 3,
                        'val': np.random.randn(100000)})

        grp = df.groupby('id')['val']

        values = np.repeat(grp.mean().values,
                           _ensure_platform_int(grp.count().values))
        expected = pd.Series(values, index=df.index, name='val')

        result = grp.transform(np.mean)
        assert_series_equal(result, expected)

        result = grp.transform('mean')
        assert_series_equal(result, expected)

        # GH 12737
        df = pd.DataFrame({'grouping': [0, 1, 1, 3], 'f': [1.1, 2.1, 3.1, 4.5],
                           'd': pd.date_range('2014-1-1', '2014-1-4'),
                           'i': [1, 2, 3, 4]},
                          columns=['grouping', 'f', 'i', 'd'])
        result = df.groupby('grouping').transform('first')

        dates = [pd.Timestamp('2014-1-1'), pd.Timestamp('2014-1-2'),
                 pd.Timestamp('2014-1-2'), pd.Timestamp('2014-1-4')]
        expected = pd.DataFrame({'f': [1.1, 2.1, 2.1, 4.5],
                                 'd': dates,
                                 'i': [1, 2, 2, 4]},
                                columns=['f', 'i', 'd'])
        assert_frame_equal(result, expected)

        # selection
        result = df.groupby('grouping')[['f', 'i']].transform('first')
        expected = expected[['f', 'i']]
        assert_frame_equal(result, expected)

        # dup columns
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=['g', 'a', 'a'])
        result = df.groupby('g').transform('first')
        expected = df.drop('g', axis=1)
        assert_frame_equal(result, expected)

    def test_transform_broadcast(self):
        grouped = self.ts.groupby(lambda x: x.month)
        result = grouped.transform(np.mean)

        self.assert_index_equal(result.index, self.ts.index)
        for _, gp in grouped:
            assert_fp_equal(result.reindex(gp.index), gp.mean())

        grouped = self.tsframe.groupby(lambda x: x.month)
        result = grouped.transform(np.mean)
        self.assert_index_equal(result.index, self.tsframe.index)
        for _, gp in grouped:
            agged = gp.mean()
            res = result.reindex(gp.index)
            for col in self.tsframe:
                assert_fp_equal(res[col], agged[col])

        # group columns
        grouped = self.tsframe.groupby({'A': 0, 'B': 0, 'C': 1, 'D': 1},
                                       axis=1)
        result = grouped.transform(np.mean)
        self.assert_index_equal(result.index, self.tsframe.index)
        self.assert_index_equal(result.columns, self.tsframe.columns)
        for _, gp in grouped:
            agged = gp.mean(1)
            res = result.reindex(columns=gp.columns)
            for idx in gp.index:
                assert_fp_equal(res.xs(idx), agged[idx])

    def test_transform_axis(self):

        # make sure that we are setting the axes
        # correctly when on axis=0 or 1
        # in the presence of a non-monotonic indexer
        # GH12713

        base = self.tsframe.iloc[0:5]
        r = len(base.index)
        c = len(base.columns)
        tso = DataFrame(np.random.randn(r, c),
                        index=base.index,
                        columns=base.columns,
                        dtype='float64')
        # monotonic
        ts = tso
        grouped = ts.groupby(lambda x: x.weekday())
        result = ts - grouped.transform('mean')
        expected = grouped.apply(lambda x: x - x.mean())
        assert_frame_equal(result, expected)

        ts = ts.T
        grouped = ts.groupby(lambda x: x.weekday(), axis=1)
        result = ts - grouped.transform('mean')
        expected = grouped.apply(lambda x: (x.T - x.mean(1)).T)
        assert_frame_equal(result, expected)

        # non-monotonic
        ts = tso.iloc[[1, 0] + list(range(2, len(base)))]
        grouped = ts.groupby(lambda x: x.weekday())
        result = ts - grouped.transform('mean')
        expected = grouped.apply(lambda x: x - x.mean())
        assert_frame_equal(result, expected)

        ts = ts.T
        grouped = ts.groupby(lambda x: x.weekday(), axis=1)
        result = ts - grouped.transform('mean')
        expected = grouped.apply(lambda x: (x.T - x.mean(1)).T)
        assert_frame_equal(result, expected)

    def test_transform_dtype(self):
        # GH 9807
        # Check transform dtype output is preserved
        df = DataFrame([[1, 3], [2, 3]])
        result = df.groupby(1).transform('mean')
        expected = DataFrame([[1.5], [1.5]])
        assert_frame_equal(result, expected)

    def test_transform_bug(self):
        # GH 5712
        # transforming on a datetime column
        df = DataFrame(dict(A=Timestamp('20130101'), B=np.arange(5)))
        result = df.groupby('A')['B'].transform(
            lambda x: x.rank(ascending=False))
        expected = Series(np.arange(5, 0, step=-1), name='B')
        assert_series_equal(result, expected)

    def test_transform_multiple(self):
        grouped = self.ts.groupby([lambda x: x.year, lambda x: x.month])

        grouped.transform(lambda x: x * 2)
        grouped.transform(np.mean)

    def test_dispatch_transform(self):
        df = self.tsframe[::5].reindex(self.tsframe.index)

        grouped = df.groupby(lambda x: x.month)

        filled = grouped.fillna(method='pad')
        fillit = lambda x: x.fillna(method='pad')
        expected = df.groupby(lambda x: x.month).transform(fillit)
        assert_frame_equal(filled, expected)

    def test_transform_select_columns(self):
        f = lambda x: x.mean()
        result = self.df.groupby('A')['C', 'D'].transform(f)

        selection = self.df[['C', 'D']]
        expected = selection.groupby(self.df['A']).transform(f)

        assert_frame_equal(result, expected)

    def test_transform_exclude_nuisance(self):

        # this also tests orderings in transform between
        # series/frame to make sure it's consistent
        expected = {}
        grouped = self.df.groupby('A')
        expected['C'] = grouped['C'].transform(np.mean)
        expected['D'] = grouped['D'].transform(np.mean)
        expected = DataFrame(expected)
        result = self.df.groupby('A').transform(np.mean)

        assert_frame_equal(result, expected)

    def test_transform_function_aliases(self):
        result = self.df.groupby('A').transform('mean')
        expected = self.df.groupby('A').transform(np.mean)
        assert_frame_equal(result, expected)

        result = self.df.groupby('A')['C'].transform('mean')
        expected = self.df.groupby('A')['C'].transform(np.mean)
        assert_series_equal(result, expected)

    def test_series_fast_transform_date(self):
        # GH 13191
        df = pd.DataFrame({'grouping': [np.nan, 1, 1, 3],
                           'd': pd.date_range('2014-1-1', '2014-1-4')})
        result = df.groupby('grouping')['d'].transform('first')
        dates = [pd.NaT, pd.Timestamp('2014-1-2'), pd.Timestamp('2014-1-2'),
                 pd.Timestamp('2014-1-4')]
        expected = pd.Series(dates, name='d')
        assert_series_equal(result, expected)

    def test_transform_length(self):
        # GH 9697
        df = pd.DataFrame({'col1': [1, 1, 2, 2], 'col2': [1, 2, 3, np.nan]})
        expected = pd.Series([3.0] * 4)

        def nsum(x):
            return np.nansum(x)

        results = [df.groupby('col1').transform(sum)['col2'],
                   df.groupby('col1')['col2'].transform(sum),
                   df.groupby('col1').transform(nsum)['col2'],
                   df.groupby('col1')['col2'].transform(nsum)]
        for result in results:
            assert_series_equal(result, expected, check_names=False)

    def test_transform_coercion(self):

        # 14457
        # when we are transforming be sure to not coerce
        # via assignment
        df = pd.DataFrame(dict(A=['a', 'a'], B=[0, 1]))
        g = df.groupby('A')

        expected = g.transform(np.mean)
        result = g.transform(lambda x: np.mean(x))
        assert_frame_equal(result, expected)

    def test_with_na(self):
        index = Index(np.arange(10))

        for dtype in ['float64', 'float32', 'int64', 'int32', 'int16', 'int8']:
            values = Series(np.ones(10), index, dtype=dtype)
            labels = Series([nan, 'foo', 'bar', 'bar', nan, nan, 'bar',
                             'bar', nan, 'foo'], index=index)

            # this SHOULD be an int
            grouped = values.groupby(labels)
            agged = grouped.agg(len)
            expected = Series([4, 2], index=['bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)

            # self.assertTrue(issubclass(agged.dtype.type, np.integer))

            # explicity return a float from my function
            def f(x):
                return float(len(x))

            agged = grouped.agg(f)
            expected = Series([4, 2], index=['bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)
            self.assertTrue(issubclass(agged.dtype.type, np.dtype(dtype).type))

    def test_groupby_transform_with_int(self):

        # GH 3740, make sure that we might upcast on item-by-item transform

        # floats
        df = DataFrame(dict(A=[1, 1, 1, 2, 2, 2], B=Series(1, dtype='float64'),
                            C=Series(
                                [1, 2, 3, 1, 2, 3], dtype='float64'), D='foo'))
        with np.errstate(all='ignore'):
            result = df.groupby('A').transform(
                lambda x: (x - x.mean()) / x.std())
        expected = DataFrame(dict(B=np.nan, C=Series(
            [-1, 0, 1, -1, 0, 1], dtype='float64')))
        assert_frame_equal(result, expected)

        # int case
        df = DataFrame(dict(A=[1, 1, 1, 2, 2, 2], B=1,
                            C=[1, 2, 3, 1, 2, 3], D='foo'))
        with np.errstate(all='ignore'):
            result = df.groupby('A').transform(
                lambda x: (x - x.mean()) / x.std())
        expected = DataFrame(dict(B=np.nan, C=[-1, 0, 1, -1, 0, 1]))
        assert_frame_equal(result, expected)

        # int that needs float conversion
        s = Series([2, 3, 4, 10, 5, -1])
        df = DataFrame(dict(A=[1, 1, 1, 2, 2, 2], B=1, C=s, D='foo'))
        with np.errstate(all='ignore'):
            result = df.groupby('A').transform(
                lambda x: (x - x.mean()) / x.std())

        s1 = s.iloc[0:3]
        s1 = (s1 - s1.mean()) / s1.std()
        s2 = s.iloc[3:6]
        s2 = (s2 - s2.mean()) / s2.std()
        expected = DataFrame(dict(B=np.nan, C=concat([s1, s2])))
        assert_frame_equal(result, expected)

        # int downcasting
        result = df.groupby('A').transform(lambda x: x * 2 / 2)
        expected = DataFrame(dict(B=1, C=[2, 3, 4, 10, 5, -1]))
        assert_frame_equal(result, expected)

    def test_indices_concatenation_order(self):

        # GH 2808

        def f1(x):
            y = x[(x.b % 2) == 1] ** 2
            if y.empty:
                multiindex = MultiIndex(levels=[[]] * 2, labels=[[]] * 2,
                                        names=['b', 'c'])
                res = DataFrame(None, columns=['a'], index=multiindex)
                return res
            else:
                y = y.set_index(['b', 'c'])
                return y

        def f2(x):
            y = x[(x.b % 2) == 1] ** 2
            if y.empty:
                return DataFrame()
            else:
                y = y.set_index(['b', 'c'])
                return y

        def f3(x):
            y = x[(x.b % 2) == 1] ** 2
            if y.empty:
                multiindex = MultiIndex(levels=[[]] * 2, labels=[[]] * 2,
                                        names=['foo', 'bar'])
                res = DataFrame(None, columns=['a', 'b'], index=multiindex)
                return res
            else:
                return y

        df = DataFrame({'a': [1, 2, 2, 2], 'b': lrange(4), 'c': lrange(5, 9)})

        df2 = DataFrame({'a': [3, 2, 2, 2], 'b': lrange(4), 'c': lrange(5, 9)})

        # correct result
        result1 = df.groupby('a').apply(f1)
        result2 = df2.groupby('a').apply(f1)
        assert_frame_equal(result1, result2)

        # should fail (not the same number of levels)
        self.assertRaises(AssertionError, df.groupby('a').apply, f2)
        self.assertRaises(AssertionError, df2.groupby('a').apply, f2)

        # should fail (incorrect shape)
        self.assertRaises(AssertionError, df.groupby('a').apply, f3)
        self.assertRaises(AssertionError, df2.groupby('a').apply, f3)

    def test_attr_wrapper(self):
        grouped = self.ts.groupby(lambda x: x.weekday())

        result = grouped.std()
        expected = grouped.agg(lambda x: np.std(x, ddof=1))
        assert_series_equal(result, expected)

        # this is pretty cool
        result = grouped.describe()
        expected = {}
        for name, gp in grouped:
            expected[name] = gp.describe()
        expected = DataFrame(expected).T
        assert_frame_equal(result.unstack(), expected)

        # get attribute
        result = grouped.dtype
        expected = grouped.agg(lambda x: x.dtype)

        # make sure raises error
        self.assertRaises(AttributeError, getattr, grouped, 'foo')

    def test_series_describe_multikey(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.describe().unstack()
        assert_series_equal(result['mean'], grouped.mean(), check_names=False)
        assert_series_equal(result['std'], grouped.std(), check_names=False)
        assert_series_equal(result['min'], grouped.min(), check_names=False)

    def test_series_describe_single(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby(lambda x: x.month)
        result = grouped.apply(lambda x: x.describe())
        expected = grouped.describe()
        assert_series_equal(result, expected)

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

    def test_series_index_name(self):
        grouped = self.df.ix[:, ['C']].groupby(self.df['A'])
        result = grouped.agg(lambda x: x.mean())
        self.assertEqual(result.index.name, 'A')

    def test_frame_describe_multikey(self):
        grouped = self.tsframe.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.describe()

        for col in self.tsframe:
            expected = grouped[col].describe()
            assert_series_equal(result[col], expected, check_names=False)

        groupedT = self.tsframe.groupby({'A': 0, 'B': 0,
                                         'C': 1, 'D': 1}, axis=1)
        result = groupedT.describe()

        for name, group in groupedT:
            assert_frame_equal(result[name], group.describe())

    def test_frame_groupby(self):
        grouped = self.tsframe.groupby(lambda x: x.weekday())

        # aggregate
        aggregated = grouped.aggregate(np.mean)
        self.assertEqual(len(aggregated), 5)
        self.assertEqual(len(aggregated.columns), 4)

        # by string
        tscopy = self.tsframe.copy()
        tscopy['weekday'] = [x.weekday() for x in tscopy.index]
        stragged = tscopy.groupby('weekday').aggregate(np.mean)
        assert_frame_equal(stragged, aggregated, check_names=False)

        # transform
        grouped = self.tsframe.head(30).groupby(lambda x: x.weekday())
        transformed = grouped.transform(lambda x: x - x.mean())
        self.assertEqual(len(transformed), 30)
        self.assertEqual(len(transformed.columns), 4)

        # transform propagate
        transformed = grouped.transform(lambda x: x.mean())
        for name, group in grouped:
            mean = group.mean()
            for idx in group.index:
                tm.assert_series_equal(transformed.xs(idx), mean,
                                       check_names=False)

        # iterate
        for weekday, group in grouped:
            self.assertEqual(group.index[0].weekday(), weekday)

        # groups / group_indices
        groups = grouped.groups
        indices = grouped.indices

        for k, v in compat.iteritems(groups):
            samething = self.tsframe.index.take(indices[k])
            self.assertTrue((samething == v).all())

    def test_grouping_is_iterable(self):
        # this code path isn't used anywhere else
        # not sure it's useful
        grouped = self.tsframe.groupby([lambda x: x.weekday(), lambda x: x.year
                                        ])

        # test it works
        for g in grouped.grouper.groupings[0]:
            pass

    def test_frame_groupby_columns(self):
        mapping = {'A': 0, 'B': 0, 'C': 1, 'D': 1}
        grouped = self.tsframe.groupby(mapping, axis=1)

        # aggregate
        aggregated = grouped.aggregate(np.mean)
        self.assertEqual(len(aggregated), len(self.tsframe))
        self.assertEqual(len(aggregated.columns), 2)

        # transform
        tf = lambda x: x - x.mean()
        groupedT = self.tsframe.T.groupby(mapping, axis=0)
        assert_frame_equal(groupedT.transform(tf).T, grouped.transform(tf))

        # iterate
        for k, v in grouped:
            self.assertEqual(len(v.columns), 2)

    def test_frame_set_name_single(self):
        grouped = self.df.groupby('A')

        result = grouped.mean()
        self.assertEqual(result.index.name, 'A')

        result = self.df.groupby('A', as_index=False).mean()
        self.assertNotEqual(result.index.name, 'A')

        result = grouped.agg(np.mean)
        self.assertEqual(result.index.name, 'A')

        result = grouped.agg({'C': np.mean, 'D': np.std})
        self.assertEqual(result.index.name, 'A')

        result = grouped['C'].mean()
        self.assertEqual(result.index.name, 'A')
        result = grouped['C'].agg(np.mean)
        self.assertEqual(result.index.name, 'A')
        result = grouped['C'].agg([np.mean, np.std])
        self.assertEqual(result.index.name, 'A')

        result = grouped['C'].agg({'foo': np.mean, 'bar': np.std})
        self.assertEqual(result.index.name, 'A')

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

    def test_multi_iter(self):
        s = Series(np.arange(6))
        k1 = np.array(['a', 'a', 'a', 'b', 'b', 'b'])
        k2 = np.array(['1', '2', '1', '2', '1', '2'])

        grouped = s.groupby([k1, k2])

        iterated = list(grouped)
        expected = [('a', '1', s[[0, 2]]), ('a', '2', s[[1]]),
                    ('b', '1', s[[4]]), ('b', '2', s[[3, 5]])]
        for i, ((one, two), three) in enumerate(iterated):
            e1, e2, e3 = expected[i]
            self.assertEqual(e1, one)
            self.assertEqual(e2, two)
            assert_series_equal(three, e3)

    def test_multi_iter_frame(self):
        k1 = np.array(['b', 'b', 'b', 'a', 'a', 'a'])
        k2 = np.array(['1', '2', '1', '2', '1', '2'])
        df = DataFrame({'v1': np.random.randn(6),
                        'v2': np.random.randn(6),
                        'k1': k1, 'k2': k2},
                       index=['one', 'two', 'three', 'four', 'five', 'six'])

        grouped = df.groupby(['k1', 'k2'])

        # things get sorted!
        iterated = list(grouped)
        idx = df.index
        expected = [('a', '1', df.ix[idx[[4]]]),
                    ('a', '2', df.ix[idx[[3, 5]]]),
                    ('b', '1', df.ix[idx[[0, 2]]]),
                    ('b', '2', df.ix[idx[[1]]])]
        for i, ((one, two), three) in enumerate(iterated):
            e1, e2, e3 = expected[i]
            self.assertEqual(e1, one)
            self.assertEqual(e2, two)
            assert_frame_equal(three, e3)

        # don't iterate through groups with no data
        df['k1'] = np.array(['b', 'b', 'b', 'a', 'a', 'a'])
        df['k2'] = np.array(['1', '1', '1', '2', '2', '2'])
        grouped = df.groupby(['k1', 'k2'])
        groups = {}
        for key, gp in grouped:
            groups[key] = gp
        self.assertEqual(len(groups), 2)

        # axis = 1
        three_levels = self.three_group.groupby(['A', 'B', 'C']).mean()
        grouped = three_levels.T.groupby(axis=1, level=(1, 2))
        for key, group in grouped:
            pass

    def test_multi_iter_panel(self):
        wp = tm.makePanel()
        grouped = wp.groupby([lambda x: x.month, lambda x: x.weekday()],
                             axis=1)

        for (month, wd), group in grouped:
            exp_axis = [x
                        for x in wp.major_axis
                        if x.month == month and x.weekday() == wd]
            expected = wp.reindex(major=exp_axis)
            assert_panel_equal(group, expected)

    def test_multi_func(self):
        col1 = self.df['A']
        col2 = self.df['B']

        grouped = self.df.groupby([col1.get, col2.get])
        agged = grouped.mean()
        expected = self.df.groupby(['A', 'B']).mean()
        assert_frame_equal(agged.ix[:, ['C', 'D']], expected.ix[:, ['C', 'D']],
                           check_names=False)  # TODO groupby get drops names

        # some "groups" with no data
        df = DataFrame({'v1': np.random.randn(6),
                        'v2': np.random.randn(6),
                        'k1': np.array(['b', 'b', 'b', 'a', 'a', 'a']),
                        'k2': np.array(['1', '1', '1', '2', '2', '2'])},
                       index=['one', 'two', 'three', 'four', 'five', 'six'])
        # only verify that it works for now
        grouped = df.groupby(['k1', 'k2'])
        grouped.agg(np.sum)

    def test_multi_key_multiple_functions(self):
        grouped = self.df.groupby(['A', 'B'])['C']

        agged = grouped.agg([np.mean, np.std])
        expected = DataFrame({'mean': grouped.agg(np.mean),
                              'std': grouped.agg(np.std)})
        assert_frame_equal(agged, expected)

    def test_frame_multi_key_function_list(self):
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

        grouped = data.groupby(['A', 'B'])
        funcs = [np.mean, np.std]
        agged = grouped.agg(funcs)
        expected = concat([grouped['D'].agg(funcs), grouped['E'].agg(funcs),
                           grouped['F'].agg(funcs)],
                          keys=['D', 'E', 'F'], axis=1)
        assert (isinstance(agged.index, MultiIndex))
        assert (isinstance(expected.index, MultiIndex))
        assert_frame_equal(agged, expected)

    def test_groupby_multiple_columns(self):
        data = self.df
        grouped = data.groupby(['A', 'B'])

        def _check_op(op):

            result1 = op(grouped)

            expected = defaultdict(dict)
            for n1, gp1 in data.groupby('A'):
                for n2, gp2 in gp1.groupby('B'):
                    expected[n1][n2] = op(gp2.ix[:, ['C', 'D']])
            expected = dict((k, DataFrame(v))
                            for k, v in compat.iteritems(expected))
            expected = Panel.fromDict(expected).swapaxes(0, 1)
            expected.major_axis.name, expected.minor_axis.name = 'A', 'B'

            # a little bit crude
            for col in ['C', 'D']:
                result_col = op(grouped[col])
                exp = expected[col]
                pivoted = result1[col].unstack()
                pivoted2 = result_col.unstack()
                assert_frame_equal(pivoted.reindex_like(exp), exp)
                assert_frame_equal(pivoted2.reindex_like(exp), exp)

        _check_op(lambda x: x.sum())
        _check_op(lambda x: x.mean())

        # test single series works the same
        result = data['C'].groupby([data['A'], data['B']]).mean()
        expected = data.groupby(['A', 'B']).mean()['C']

        assert_series_equal(result, expected)

    def test_groupby_as_index_agg(self):
        grouped = self.df.groupby('A', as_index=False)

        # single-key

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        assert_frame_equal(result, expected)

        result2 = grouped.agg(OrderedDict([['C', np.mean], ['D', np.sum]]))
        expected2 = grouped.mean()
        expected2['D'] = grouped.sum()['D']
        assert_frame_equal(result2, expected2)

        grouped = self.df.groupby('A', as_index=True)
        expected3 = grouped['C'].sum()
        expected3 = DataFrame(expected3).rename(columns={'C': 'Q'})
        result3 = grouped['C'].agg({'Q': np.sum})
        assert_frame_equal(result3, expected3)

        # multi-key

        grouped = self.df.groupby(['A', 'B'], as_index=False)

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        assert_frame_equal(result, expected)

        result2 = grouped.agg(OrderedDict([['C', np.mean], ['D', np.sum]]))
        expected2 = grouped.mean()
        expected2['D'] = grouped.sum()['D']
        assert_frame_equal(result2, expected2)

        expected3 = grouped['C'].sum()
        expected3 = DataFrame(expected3).rename(columns={'C': 'Q'})
        result3 = grouped['C'].agg({'Q': np.sum})
        assert_frame_equal(result3, expected3)

        # GH7115 & GH8112 & GH8582
        df = DataFrame(np.random.randint(0, 100, (50, 3)),
                       columns=['jim', 'joe', 'jolie'])
        ts = Series(np.random.randint(5, 10, 50), name='jim')

        gr = df.groupby(ts)
        gr.nth(0)  # invokes set_selection_from_grouper internally
        assert_frame_equal(gr.apply(sum), df.groupby(ts).apply(sum))

        for attr in ['mean', 'max', 'count', 'idxmax', 'cumsum', 'all']:
            gr = df.groupby(ts, as_index=False)
            left = getattr(gr, attr)()

            gr = df.groupby(ts.values, as_index=True)
            right = getattr(gr, attr)().reset_index(drop=True)

            assert_frame_equal(left, right)

    def test_series_groupby_nunique(self):
        from itertools import product
        from string import ascii_lowercase

        def check_nunique(df, keys):
            for sort, dropna in product((False, True), repeat=2):
                gr = df.groupby(keys, sort=sort)
                left = gr['julie'].nunique(dropna=dropna)

                gr = df.groupby(keys, sort=sort)
                right = gr['julie'].apply(Series.nunique, dropna=dropna)

                assert_series_equal(left, right)

        days = date_range('2015-08-23', periods=10)

        for n, m in product(10 ** np.arange(2, 6), (10, 100, 1000)):
            frame = DataFrame({
                'jim': np.random.choice(
                    list(ascii_lowercase), n),
                'joe': np.random.choice(days, n),
                'julie': np.random.randint(0, m, n)
            })

            check_nunique(frame, ['jim'])
            check_nunique(frame, ['jim', 'joe'])

            frame.loc[1::17, 'jim'] = None
            frame.loc[3::37, 'joe'] = None
            frame.loc[7::19, 'julie'] = None
            frame.loc[8::19, 'julie'] = None
            frame.loc[9::19, 'julie'] = None

            check_nunique(frame, ['jim'])
            check_nunique(frame, ['jim', 'joe'])

    def test_series_groupby_value_counts(self):
        from itertools import product

        def rebuild_index(df):
            arr = list(map(df.index.get_level_values, range(df.index.nlevels)))
            df.index = MultiIndex.from_arrays(arr, names=df.index.names)
            return df

        def check_value_counts(df, keys, bins):
            for isort, normalize, sort, ascending, dropna \
                    in product((False, True), repeat=5):

                kwargs = dict(normalize=normalize, sort=sort,
                              ascending=ascending, dropna=dropna, bins=bins)

                gr = df.groupby(keys, sort=isort)
                left = gr['3rd'].value_counts(**kwargs)

                gr = df.groupby(keys, sort=isort)
                right = gr['3rd'].apply(Series.value_counts, **kwargs)
                right.index.names = right.index.names[:-1] + ['3rd']

                # have to sort on index because of unstable sort on values
                left, right = map(rebuild_index, (left, right))  # xref GH9212
                assert_series_equal(left.sort_index(), right.sort_index())

        def loop(df):
            bins = None, np.arange(0, max(5, df['3rd'].max()) + 1, 2)
            keys = '1st', '2nd', ('1st', '2nd')
            for k, b in product(keys, bins):
                check_value_counts(df, k, b)

        days = date_range('2015-08-24', periods=10)

        for n, m in product((100, 1000), (5, 20)):
            frame = DataFrame({
                '1st': np.random.choice(
                    list('abcd'), n),
                '2nd': np.random.choice(days, n),
                '3rd': np.random.randint(1, m + 1, n)
            })

            loop(frame)

            frame.loc[1::11, '1st'] = nan
            frame.loc[3::17, '2nd'] = nan
            frame.loc[7::19, '3rd'] = nan
            frame.loc[8::19, '3rd'] = nan
            frame.loc[9::19, '3rd'] = nan

            loop(frame)

    def test_multiindex_passthru(self):

        # GH 7997
        # regression from 0.14.1
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        df.columns = pd.MultiIndex.from_tuples([(0, 1), (1, 1), (2, 1)])

        result = df.groupby(axis=1, level=[0, 1]).first()
        assert_frame_equal(result, df)

    def test_multiindex_negative_level(self):
        # GH 13901
        result = self.mframe.groupby(level=-1).sum()
        expected = self.mframe.groupby(level='second').sum()
        assert_frame_equal(result, expected)

        result = self.mframe.groupby(level=-2).sum()
        expected = self.mframe.groupby(level='first').sum()
        assert_frame_equal(result, expected)

        result = self.mframe.groupby(level=[-2, -1]).sum()
        expected = self.mframe
        assert_frame_equal(result, expected)

        result = self.mframe.groupby(level=[-1, 'first']).sum()
        expected = self.mframe.groupby(level=['second', 'first']).sum()
        assert_frame_equal(result, expected)

    def test_multifunc_select_col_integer_cols(self):
        df = self.df
        df.columns = np.arange(len(df.columns))

        # it works!
        df.groupby(1, as_index=False)[2].agg({'Q': np.mean})

    def test_as_index_series_return_frame(self):
        grouped = self.df.groupby('A', as_index=False)
        grouped2 = self.df.groupby(['A', 'B'], as_index=False)

        result = grouped['C'].agg(np.sum)
        expected = grouped.agg(np.sum).ix[:, ['A', 'C']]
        tm.assertIsInstance(result, DataFrame)
        assert_frame_equal(result, expected)

        result2 = grouped2['C'].agg(np.sum)
        expected2 = grouped2.agg(np.sum).ix[:, ['A', 'B', 'C']]
        tm.assertIsInstance(result2, DataFrame)
        assert_frame_equal(result2, expected2)

        result = grouped['C'].sum()
        expected = grouped.sum().ix[:, ['A', 'C']]
        tm.assertIsInstance(result, DataFrame)
        assert_frame_equal(result, expected)

        result2 = grouped2['C'].sum()
        expected2 = grouped2.sum().ix[:, ['A', 'B', 'C']]
        tm.assertIsInstance(result2, DataFrame)
        assert_frame_equal(result2, expected2)

        # corner case
        self.assertRaises(Exception, grouped['C'].__getitem__, 'D')

    def test_groupby_as_index_cython(self):
        data = self.df

        # single-key
        grouped = data.groupby('A', as_index=False)
        result = grouped.mean()
        expected = data.groupby(['A']).mean()
        expected.insert(0, 'A', expected.index)
        expected.index = np.arange(len(expected))
        assert_frame_equal(result, expected)

        # multi-key
        grouped = data.groupby(['A', 'B'], as_index=False)
        result = grouped.mean()
        expected = data.groupby(['A', 'B']).mean()

        arrays = lzip(*expected.index._tuple_index)
        expected.insert(0, 'A', arrays[0])
        expected.insert(1, 'B', arrays[1])
        expected.index = np.arange(len(expected))
        assert_frame_equal(result, expected)

    def test_groupby_as_index_series_scalar(self):
        grouped = self.df.groupby(['A', 'B'], as_index=False)

        # GH #421

        result = grouped['C'].agg(len)
        expected = grouped.agg(len).ix[:, ['A', 'B', 'C']]
        assert_frame_equal(result, expected)

    def test_groupby_as_index_corner(self):
        self.assertRaises(TypeError, self.ts.groupby, lambda x: x.weekday(),
                          as_index=False)

        self.assertRaises(ValueError, self.df.groupby, lambda x: x.lower(),
                          as_index=False, axis=1)

    def test_groupby_as_index_apply(self):
        # GH #4648 and #3417
        df = DataFrame({'item_id': ['b', 'b', 'a', 'c', 'a', 'b'],
                        'user_id': [1, 2, 1, 1, 3, 1],
                        'time': range(6)})

        g_as = df.groupby('user_id', as_index=True)
        g_not_as = df.groupby('user_id', as_index=False)

        res_as = g_as.head(2).index
        res_not_as = g_not_as.head(2).index
        exp = Index([0, 1, 2, 4])
        assert_index_equal(res_as, exp)
        assert_index_equal(res_not_as, exp)

        res_as_apply = g_as.apply(lambda x: x.head(2)).index
        res_not_as_apply = g_not_as.apply(lambda x: x.head(2)).index

        # apply doesn't maintain the original ordering
        # changed in GH5610 as the as_index=False returns a MI here
        exp_not_as_apply = MultiIndex.from_tuples([(0, 0), (0, 2), (1, 1), (
            2, 4)])
        tp = [(1, 0), (1, 2), (2, 1), (3, 4)]
        exp_as_apply = MultiIndex.from_tuples(tp, names=['user_id', None])

        assert_index_equal(res_as_apply, exp_as_apply)
        assert_index_equal(res_not_as_apply, exp_not_as_apply)

        ind = Index(list('abcde'))
        df = DataFrame([[1, 2], [2, 3], [1, 4], [1, 5], [2, 6]], index=ind)
        res = df.groupby(0, as_index=False).apply(lambda x: x).index
        assert_index_equal(res, ind)

    def test_groupby_head_tail(self):
        df = DataFrame([[1, 2], [1, 4], [5, 6]], columns=['A', 'B'])
        g_as = df.groupby('A', as_index=True)
        g_not_as = df.groupby('A', as_index=False)

        # as_index= False, much easier
        assert_frame_equal(df.loc[[0, 2]], g_not_as.head(1))
        assert_frame_equal(df.loc[[1, 2]], g_not_as.tail(1))

        empty_not_as = DataFrame(columns=df.columns,
                                 index=pd.Index([], dtype=df.index.dtype))
        empty_not_as['A'] = empty_not_as['A'].astype(df.A.dtype)
        empty_not_as['B'] = empty_not_as['B'].astype(df.B.dtype)
        assert_frame_equal(empty_not_as, g_not_as.head(0))
        assert_frame_equal(empty_not_as, g_not_as.tail(0))
        assert_frame_equal(empty_not_as, g_not_as.head(-1))
        assert_frame_equal(empty_not_as, g_not_as.tail(-1))

        assert_frame_equal(df, g_not_as.head(7))  # contains all
        assert_frame_equal(df, g_not_as.tail(7))

        # as_index=True, (used to be different)
        df_as = df

        assert_frame_equal(df_as.loc[[0, 2]], g_as.head(1))
        assert_frame_equal(df_as.loc[[1, 2]], g_as.tail(1))

        empty_as = DataFrame(index=df_as.index[:0], columns=df.columns)
        empty_as['A'] = empty_not_as['A'].astype(df.A.dtype)
        empty_as['B'] = empty_not_as['B'].astype(df.B.dtype)
        assert_frame_equal(empty_as, g_as.head(0))
        assert_frame_equal(empty_as, g_as.tail(0))
        assert_frame_equal(empty_as, g_as.head(-1))
        assert_frame_equal(empty_as, g_as.tail(-1))

        assert_frame_equal(df_as, g_as.head(7))  # contains all
        assert_frame_equal(df_as, g_as.tail(7))

        # test with selection
        assert_frame_equal(g_as[[]].head(1), df_as.loc[[0, 2], []])
        assert_frame_equal(g_as[['A']].head(1), df_as.loc[[0, 2], ['A']])
        assert_frame_equal(g_as[['B']].head(1), df_as.loc[[0, 2], ['B']])
        assert_frame_equal(g_as[['A', 'B']].head(1), df_as.loc[[0, 2]])

        assert_frame_equal(g_not_as[[]].head(1), df_as.loc[[0, 2], []])
        assert_frame_equal(g_not_as[['A']].head(1), df_as.loc[[0, 2], ['A']])
        assert_frame_equal(g_not_as[['B']].head(1), df_as.loc[[0, 2], ['B']])
        assert_frame_equal(g_not_as[['A', 'B']].head(1), df_as.loc[[0, 2]])

    def test_groupby_multiple_key(self):
        df = tm.makeTimeDataFrame()
        grouped = df.groupby([lambda x: x.year, lambda x: x.month,
                              lambda x: x.day])
        agged = grouped.sum()
        assert_almost_equal(df.values, agged.values)

        grouped = df.T.groupby([lambda x: x.year,
                                lambda x: x.month,
                                lambda x: x.day], axis=1)

        agged = grouped.agg(lambda x: x.sum())
        self.assert_index_equal(agged.index, df.columns)
        assert_almost_equal(df.T.values, agged.values)

        agged = grouped.agg(lambda x: x.sum())
        assert_almost_equal(df.T.values, agged.values)

    def test_groupby_multi_corner(self):
        # test that having an all-NA column doesn't mess you up
        df = self.df.copy()
        df['bad'] = np.nan
        agged = df.groupby(['A', 'B']).mean()

        expected = self.df.groupby(['A', 'B']).mean()
        expected['bad'] = np.nan

        assert_frame_equal(agged, expected)

    def test_omit_nuisance(self):
        grouped = self.df.groupby('A')

        result = grouped.mean()
        expected = self.df.ix[:, ['A', 'C', 'D']].groupby('A').mean()
        assert_frame_equal(result, expected)

        agged = grouped.agg(np.mean)
        exp = grouped.mean()
        assert_frame_equal(agged, exp)

        df = self.df.ix[:, ['A', 'C', 'D']]
        df['E'] = datetime.now()
        grouped = df.groupby('A')
        result = grouped.agg(np.sum)
        expected = grouped.sum()
        assert_frame_equal(result, expected)

        # won't work with axis = 1
        grouped = df.groupby({'A': 0, 'C': 0, 'D': 1, 'E': 1}, axis=1)
        result = self.assertRaises(TypeError, grouped.agg,
                                   lambda x: x.sum(0, numeric_only=False))

    def test_omit_nuisance_python_multiple(self):
        grouped = self.three_group.groupby(['A', 'B'])

        agged = grouped.agg(np.mean)
        exp = grouped.mean()
        assert_frame_equal(agged, exp)

    def test_empty_groups_corner(self):
        # handle empty groups
        df = DataFrame({'k1': np.array(['b', 'b', 'b', 'a', 'a', 'a']),
                        'k2': np.array(['1', '1', '1', '2', '2', '2']),
                        'k3': ['foo', 'bar'] * 3,
                        'v1': np.random.randn(6),
                        'v2': np.random.randn(6)})

        grouped = df.groupby(['k1', 'k2'])
        result = grouped.agg(np.mean)
        expected = grouped.mean()
        assert_frame_equal(result, expected)

        grouped = self.mframe[3:5].groupby(level=0)
        agged = grouped.apply(lambda x: x.mean())
        agged_A = grouped['A'].apply(np.mean)
        assert_series_equal(agged['A'], agged_A)
        self.assertEqual(agged.index.name, 'first')

    def test_apply_concat_preserve_names(self):
        grouped = self.three_group.groupby(['A', 'B'])

        def desc(group):
            result = group.describe()
            result.index.name = 'stat'
            return result

        def desc2(group):
            result = group.describe()
            result.index.name = 'stat'
            result = result[:len(group)]
            # weirdo
            return result

        def desc3(group):
            result = group.describe()

            # names are different
            result.index.name = 'stat_%d' % len(group)

            result = result[:len(group)]
            # weirdo
            return result

        result = grouped.apply(desc)
        self.assertEqual(result.index.names, ('A', 'B', 'stat'))

        result2 = grouped.apply(desc2)
        self.assertEqual(result2.index.names, ('A', 'B', 'stat'))

        result3 = grouped.apply(desc3)
        self.assertEqual(result3.index.names, ('A', 'B', None))

    def test_nonsense_func(self):
        df = DataFrame([0])
        self.assertRaises(Exception, df.groupby, lambda x: x + 'foo')

    def test_builtins_apply(self):  # GH8155
        df = pd.DataFrame(np.random.randint(1, 50, (1000, 2)),
                          columns=['jim', 'joe'])
        df['jolie'] = np.random.randn(1000)

        for keys in ['jim', ['jim', 'joe']]:  # single key & multi-key
            if keys == 'jim':
                continue
            for f in [max, min, sum]:
                fname = f.__name__
                result = df.groupby(keys).apply(f)
                result.shape
                ngroups = len(df.drop_duplicates(subset=keys))
                assert result.shape == (ngroups, 3), 'invalid frame shape: '\
                    '{} (expected ({}, 3))'.format(result.shape, ngroups)

                assert_frame_equal(result,  # numpy's equivalent function
                                   df.groupby(keys).apply(getattr(np, fname)))

                if f != sum:
                    expected = df.groupby(keys).agg(fname).reset_index()
                    expected.set_index(keys, inplace=True, drop=False)
                    assert_frame_equal(result, expected, check_dtype=False)

                assert_series_equal(getattr(result, fname)(),
                                    getattr(df, fname)())

    def test_cythonized_aggers(self):
        data = {'A': [0, 0, 0, 0, 1, 1, 1, 1, 1, 1., nan, nan],
                'B': ['A', 'B'] * 6,
                'C': np.random.randn(12)}
        df = DataFrame(data)
        df.loc[2:10:2, 'C'] = nan

        def _testit(name):

            op = lambda x: getattr(x, name)()

            # single column
            grouped = df.drop(['B'], axis=1).groupby('A')
            exp = {}
            for cat, group in grouped:
                exp[cat] = op(group['C'])
            exp = DataFrame({'C': exp})
            exp.index.name = 'A'
            result = op(grouped)
            assert_frame_equal(result, exp)

            # multiple columns
            grouped = df.groupby(['A', 'B'])
            expd = {}
            for (cat1, cat2), group in grouped:
                expd.setdefault(cat1, {})[cat2] = op(group['C'])
            exp = DataFrame(expd).T.stack(dropna=False)
            exp.index.names = ['A', 'B']
            exp.name = 'C'

            result = op(grouped)['C']
            if not tm._incompat_bottleneck_version(name):
                assert_series_equal(result, exp)

        _testit('count')
        _testit('sum')
        _testit('std')
        _testit('var')
        _testit('sem')
        _testit('mean')
        _testit('median')
        _testit('prod')
        _testit('min')
        _testit('max')

    def test_max_min_non_numeric(self):
        # #2700
        aa = DataFrame({'nn': [11, 11, 22, 22],
                        'ii': [1, 2, 3, 4],
                        'ss': 4 * ['mama']})

        result = aa.groupby('nn').max()
        self.assertTrue('ss' in result)

        result = aa.groupby('nn').min()
        self.assertTrue('ss' in result)

    def test_cython_agg_boolean(self):
        frame = DataFrame({'a': np.random.randint(0, 5, 50),
                           'b': np.random.randint(0, 2, 50).astype('bool')})
        result = frame.groupby('a')['b'].mean()
        expected = frame.groupby('a')['b'].agg(np.mean)

        assert_series_equal(result, expected)

    def test_cython_agg_nothing_to_agg(self):
        frame = DataFrame({'a': np.random.randint(0, 5, 50),
                           'b': ['foo', 'bar'] * 25})
        self.assertRaises(DataError, frame.groupby('a')['b'].mean)

        frame = DataFrame({'a': np.random.randint(0, 5, 50),
                           'b': ['foo', 'bar'] * 25})
        self.assertRaises(DataError, frame[['b']].groupby(frame['a']).mean)

    def test_cython_agg_nothing_to_agg_with_dates(self):
        frame = DataFrame({'a': np.random.randint(0, 5, 50),
                           'b': ['foo', 'bar'] * 25,
                           'dates': pd.date_range('now', periods=50,
                                                  freq='T')})
        with tm.assertRaisesRegexp(DataError, "No numeric types to aggregate"):
            frame.groupby('b').dates.mean()

    def test_groupby_timedelta_cython_count(self):
        df = DataFrame({'g': list('ab' * 2),
                        'delt': np.arange(4).astype('timedelta64[ns]')})
        expected = Series([
            2, 2
        ], index=pd.Index(['a', 'b'], name='g'), name='delt')
        result = df.groupby('g').delt.count()
        tm.assert_series_equal(expected, result)

    def test_cython_agg_frame_columns(self):
        # #2113
        df = DataFrame({'x': [1, 2, 3], 'y': [3, 4, 5]})

        df.groupby(level=0, axis='columns').mean()
        df.groupby(level=0, axis='columns').mean()
        df.groupby(level=0, axis='columns').mean()
        df.groupby(level=0, axis='columns').mean()

    def test_wrap_aggregated_output_multindex(self):
        df = self.mframe.T
        df['baz', 'two'] = 'peekaboo'

        keys = [np.array([0, 0, 1]), np.array([0, 0, 1])]
        agged = df.groupby(keys).agg(np.mean)
        tm.assertIsInstance(agged.columns, MultiIndex)

        def aggfun(ser):
            if ser.name == ('foo', 'one'):
                raise TypeError
            else:
                return ser.sum()

        agged2 = df.groupby(keys).aggregate(aggfun)
        self.assertEqual(len(agged2.columns) + 1, len(df.columns))

    def test_groupby_level(self):
        frame = self.mframe
        deleveled = frame.reset_index()

        result0 = frame.groupby(level=0).sum()
        result1 = frame.groupby(level=1).sum()

        expected0 = frame.groupby(deleveled['first'].values).sum()
        expected1 = frame.groupby(deleveled['second'].values).sum()

        expected0 = expected0.reindex(frame.index.levels[0])
        expected1 = expected1.reindex(frame.index.levels[1])

        self.assertEqual(result0.index.name, 'first')
        self.assertEqual(result1.index.name, 'second')

        assert_frame_equal(result0, expected0)
        assert_frame_equal(result1, expected1)
        self.assertEqual(result0.index.name, frame.index.names[0])
        self.assertEqual(result1.index.name, frame.index.names[1])

        # groupby level name
        result0 = frame.groupby(level='first').sum()
        result1 = frame.groupby(level='second').sum()
        assert_frame_equal(result0, expected0)
        assert_frame_equal(result1, expected1)

        # axis=1

        result0 = frame.T.groupby(level=0, axis=1).sum()
        result1 = frame.T.groupby(level=1, axis=1).sum()
        assert_frame_equal(result0, expected0.T)
        assert_frame_equal(result1, expected1.T)

        # raise exception for non-MultiIndex
        self.assertRaises(ValueError, self.df.groupby, level=1)

    def test_groupby_level_index_names(self):
        # GH4014 this used to raise ValueError since 'exp'>1 (in py2)
        df = DataFrame({'exp': ['A'] * 3 + ['B'] * 3,
                        'var1': lrange(6), }).set_index('exp')
        df.groupby(level='exp')
        self.assertRaises(ValueError, df.groupby, level='foo')

    def test_groupby_level_with_nas(self):
        index = MultiIndex(levels=[[1, 0], [0, 1, 2, 3]],
                           labels=[[1, 1, 1, 1, 0, 0, 0, 0], [0, 1, 2, 3, 0, 1,
                                                              2, 3]])

        # factorizing doesn't confuse things
        s = Series(np.arange(8.), index=index)
        result = s.groupby(level=0).sum()
        expected = Series([22., 6.], index=[1, 0])
        assert_series_equal(result, expected)

        index = MultiIndex(levels=[[1, 0], [0, 1, 2, 3]],
                           labels=[[1, 1, 1, 1, -1, 0, 0, 0], [0, 1, 2, 3, 0,
                                                               1, 2, 3]])

        # factorizing doesn't confuse things
        s = Series(np.arange(8.), index=index)
        result = s.groupby(level=0).sum()
        expected = Series([18., 6.], index=[1, 0])
        assert_series_equal(result, expected)

    def test_groupby_level_apply(self):
        frame = self.mframe

        result = frame.groupby(level=0).count()
        self.assertEqual(result.index.name, 'first')
        result = frame.groupby(level=1).count()
        self.assertEqual(result.index.name, 'second')

        result = frame['A'].groupby(level=0).count()
        self.assertEqual(result.index.name, 'first')

    def test_groupby_args(self):
        # PR8618 and issue 8015
        frame = self.mframe

        def j():
            frame.groupby()

        self.assertRaisesRegexp(TypeError,
                                "You have to supply one of 'by' and 'level'",
                                j)

        def k():
            frame.groupby(by=None, level=None)

        self.assertRaisesRegexp(TypeError,
                                "You have to supply one of 'by' and 'level'",
                                k)

    def test_groupby_level_mapper(self):
        frame = self.mframe
        deleveled = frame.reset_index()

        mapper0 = {'foo': 0, 'bar': 0, 'baz': 1, 'qux': 1}
        mapper1 = {'one': 0, 'two': 0, 'three': 1}

        result0 = frame.groupby(mapper0, level=0).sum()
        result1 = frame.groupby(mapper1, level=1).sum()

        mapped_level0 = np.array([mapper0.get(x) for x in deleveled['first']])
        mapped_level1 = np.array([mapper1.get(x) for x in deleveled['second']])
        expected0 = frame.groupby(mapped_level0).sum()
        expected1 = frame.groupby(mapped_level1).sum()
        expected0.index.name, expected1.index.name = 'first', 'second'

        assert_frame_equal(result0, expected0)
        assert_frame_equal(result1, expected1)

    def test_groupby_level_nonmulti(self):
        # GH 1313, GH 13901
        s = Series([1, 2, 3, 10, 4, 5, 20, 6],
                   Index([1, 2, 3, 1, 4, 5, 2, 6], name='foo'))
        expected = Series([11, 22, 3, 4, 5, 6],
                          Index(range(1, 7), name='foo'))

        result = s.groupby(level=0).sum()
        self.assert_series_equal(result, expected)
        result = s.groupby(level=[0]).sum()
        self.assert_series_equal(result, expected)
        result = s.groupby(level=-1).sum()
        self.assert_series_equal(result, expected)
        result = s.groupby(level=[-1]).sum()
        self.assert_series_equal(result, expected)

        tm.assertRaises(ValueError, s.groupby, level=1)
        tm.assertRaises(ValueError, s.groupby, level=-2)
        tm.assertRaises(ValueError, s.groupby, level=[])
        tm.assertRaises(ValueError, s.groupby, level=[0, 0])
        tm.assertRaises(ValueError, s.groupby, level=[0, 1])
        tm.assertRaises(ValueError, s.groupby, level=[1])

    def test_groupby_complex(self):
        # GH 12902
        a = Series(data=np.arange(4) * (1 + 2j), index=[0, 0, 1, 1])
        expected = Series((1 + 2j, 5 + 10j))

        result = a.groupby(level=0).sum()
        assert_series_equal(result, expected)

        result = a.sum(level=0)
        assert_series_equal(result, expected)

    def test_level_preserve_order(self):
        grouped = self.mframe.groupby(level=0)
        exp_labels = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 3], np.intp)
        assert_almost_equal(grouped.grouper.labels[0], exp_labels)

    def test_grouping_labels(self):
        grouped = self.mframe.groupby(self.mframe.index.get_level_values(0))
        exp_labels = np.array([2, 2, 2, 0, 0, 1, 1, 3, 3, 3], dtype=np.intp)
        assert_almost_equal(grouped.grouper.labels[0], exp_labels)

    def test_cython_fail_agg(self):
        dr = bdate_range('1/1/2000', periods=50)
        ts = Series(['A', 'B', 'C', 'D', 'E'] * 10, index=dr)

        grouped = ts.groupby(lambda x: x.month)
        summed = grouped.sum()
        expected = grouped.agg(np.sum)
        assert_series_equal(summed, expected)

    def test_apply_series_to_frame(self):
        def f(piece):
            with np.errstate(invalid='ignore'):
                logged = np.log(piece)
            return DataFrame({'value': piece,
                              'demeaned': piece - piece.mean(),
                              'logged': logged})

        dr = bdate_range('1/1/2000', periods=100)
        ts = Series(np.random.randn(100), index=dr)

        grouped = ts.groupby(lambda x: x.month)
        result = grouped.apply(f)

        tm.assertIsInstance(result, DataFrame)
        self.assert_index_equal(result.index, ts.index)

    def test_apply_series_yield_constant(self):
        result = self.df.groupby(['A', 'B'])['C'].apply(len)
        self.assertEqual(result.index.names[:2], ('A', 'B'))

    def test_apply_frame_yield_constant(self):
        # GH13568
        result = self.df.groupby(['A', 'B']).apply(len)
        self.assertTrue(isinstance(result, Series))
        self.assertIsNone(result.name)

        result = self.df.groupby(['A', 'B'])[['C', 'D']].apply(len)
        self.assertTrue(isinstance(result, Series))
        self.assertIsNone(result.name)

    def test_apply_frame_to_series(self):
        grouped = self.df.groupby(['A', 'B'])
        result = grouped.apply(len)
        expected = grouped.count()['C']
        self.assert_index_equal(result.index, expected.index)
        self.assert_numpy_array_equal(result.values, expected.values)

    def test_apply_frame_concat_series(self):
        def trans(group):
            return group.groupby('B')['C'].sum().sort_values()[:2]

        def trans2(group):
            grouped = group.groupby(df.reindex(group.index)['B'])
            return grouped.sum().sort_values()[:2]

        df = DataFrame({'A': np.random.randint(0, 5, 1000),
                        'B': np.random.randint(0, 5, 1000),
                        'C': np.random.randn(1000)})

        result = df.groupby('A').apply(trans)
        exp = df.groupby('A')['C'].apply(trans2)
        assert_series_equal(result, exp, check_names=False)
        self.assertEqual(result.name, 'C')

    def test_apply_transform(self):
        grouped = self.ts.groupby(lambda x: x.month)
        result = grouped.apply(lambda x: x * 2)
        expected = grouped.transform(lambda x: x * 2)
        assert_series_equal(result, expected)

    def test_apply_multikey_corner(self):
        grouped = self.tsframe.groupby([lambda x: x.year, lambda x: x.month])

        def f(group):
            return group.sort_values('A')[-5:]

        result = grouped.apply(f)
        for key, group in grouped:
            assert_frame_equal(result.ix[key], f(group))

    def test_mutate_groups(self):

        # GH3380

        mydf = DataFrame({
            'cat1': ['a'] * 8 + ['b'] * 6,
            'cat2': ['c'] * 2 + ['d'] * 2 + ['e'] * 2 + ['f'] * 2 + ['c'] * 2 +
            ['d'] * 2 + ['e'] * 2,
            'cat3': lmap(lambda x: 'g%s' % x, lrange(1, 15)),
            'val': np.random.randint(100, size=14),
        })

        def f_copy(x):
            x = x.copy()
            x['rank'] = x.val.rank(method='min')
            return x.groupby('cat2')['rank'].min()

        def f_no_copy(x):
            x['rank'] = x.val.rank(method='min')
            return x.groupby('cat2')['rank'].min()

        grpby_copy = mydf.groupby('cat1').apply(f_copy)
        grpby_no_copy = mydf.groupby('cat1').apply(f_no_copy)
        assert_series_equal(grpby_copy, grpby_no_copy)

    def test_no_mutate_but_looks_like(self):

        # GH 8467
        # first show's mutation indicator
        # second does not, but should yield the same results
        df = DataFrame({'key': [1, 1, 1, 2, 2, 2, 3, 3, 3], 'value': range(9)})

        result1 = df.groupby('key', group_keys=True).apply(lambda x: x[:].key)
        result2 = df.groupby('key', group_keys=True).apply(lambda x: x.key)
        assert_series_equal(result1, result2)

    def test_apply_chunk_view(self):
        # Low level tinkering could be unsafe, make sure not
        df = DataFrame({'key': [1, 1, 1, 2, 2, 2, 3, 3, 3],
                        'value': lrange(9)})

        # return view
        f = lambda x: x[:2]

        result = df.groupby('key', group_keys=False).apply(f)
        expected = df.take([0, 1, 3, 4, 6, 7])
        assert_frame_equal(result, expected)

    def test_apply_no_name_column_conflict(self):
        df = DataFrame({'name': [1, 1, 1, 1, 1, 1, 2, 2, 2, 2],
                        'name2': [0, 0, 0, 1, 1, 1, 0, 0, 1, 1],
                        'value': lrange(10)[::-1]})

        # it works! #2605
        grouped = df.groupby(['name', 'name2'])
        grouped.apply(lambda x: x.sort_values('value', inplace=True))

    def test_groupby_series_indexed_differently(self):
        s1 = Series([5.0, -9.0, 4.0, 100., -5., 55., 6.7],
                    index=Index(['a', 'b', 'c', 'd', 'e', 'f', 'g']))
        s2 = Series([1.0, 1.0, 4.0, 5.0, 5.0, 7.0],
                    index=Index(['a', 'b', 'd', 'f', 'g', 'h']))

        grouped = s1.groupby(s2)
        agged = grouped.mean()
        exp = s1.groupby(s2.reindex(s1.index).get).mean()
        assert_series_equal(agged, exp)

    def test_groupby_with_hier_columns(self):
        tuples = list(zip(*[['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux',
                             'qux'], ['one', 'two', 'one', 'two', 'one', 'two',
                                      'one', 'two']]))
        index = MultiIndex.from_tuples(tuples)
        columns = MultiIndex.from_tuples([('A', 'cat'), ('B', 'dog'), (
            'B', 'cat'), ('A', 'dog')])
        df = DataFrame(np.random.randn(8, 4), index=index, columns=columns)

        result = df.groupby(level=0).mean()
        self.assert_index_equal(result.columns, columns)

        result = df.groupby(level=0, axis=1).mean()
        self.assert_index_equal(result.index, df.index)

        result = df.groupby(level=0).agg(np.mean)
        self.assert_index_equal(result.columns, columns)

        result = df.groupby(level=0).apply(lambda x: x.mean())
        self.assert_index_equal(result.columns, columns)

        result = df.groupby(level=0, axis=1).agg(lambda x: x.mean(1))
        self.assert_index_equal(result.columns, Index(['A', 'B']))
        self.assert_index_equal(result.index, df.index)

        # add a nuisance column
        sorted_columns, _ = columns.sortlevel(0)
        df['A', 'foo'] = 'bar'
        result = df.groupby(level=0).mean()
        self.assert_index_equal(result.columns, df.columns[:-1])

    def test_pass_args_kwargs(self):
        from numpy import percentile

        def f(x, q=None, axis=0):
            return percentile(x, q, axis=axis)

        g = lambda x: percentile(x, 80, axis=0)

        # Series
        ts_grouped = self.ts.groupby(lambda x: x.month)
        agg_result = ts_grouped.agg(percentile, 80, axis=0)
        apply_result = ts_grouped.apply(percentile, 80, axis=0)
        trans_result = ts_grouped.transform(percentile, 80, axis=0)

        agg_expected = ts_grouped.quantile(.8)
        trans_expected = ts_grouped.transform(g)

        assert_series_equal(apply_result, agg_expected)
        assert_series_equal(agg_result, agg_expected, check_names=False)
        assert_series_equal(trans_result, trans_expected)

        agg_result = ts_grouped.agg(f, q=80)
        apply_result = ts_grouped.apply(f, q=80)
        trans_result = ts_grouped.transform(f, q=80)
        assert_series_equal(agg_result, agg_expected)
        assert_series_equal(apply_result, agg_expected)
        assert_series_equal(trans_result, trans_expected)

        # DataFrame
        df_grouped = self.tsframe.groupby(lambda x: x.month)
        agg_result = df_grouped.agg(percentile, 80, axis=0)
        apply_result = df_grouped.apply(DataFrame.quantile, .8)
        expected = df_grouped.quantile(.8)
        assert_frame_equal(apply_result, expected)
        assert_frame_equal(agg_result, expected, check_names=False)

        agg_result = df_grouped.agg(f, q=80)
        apply_result = df_grouped.apply(DataFrame.quantile, q=.8)
        assert_frame_equal(agg_result, expected, check_names=False)
        assert_frame_equal(apply_result, expected)

    def test_size(self):
        grouped = self.df.groupby(['A', 'B'])
        result = grouped.size()
        for key, group in grouped:
            self.assertEqual(result[key], len(group))

        grouped = self.df.groupby('A')
        result = grouped.size()
        for key, group in grouped:
            self.assertEqual(result[key], len(group))

        grouped = self.df.groupby('B')
        result = grouped.size()
        for key, group in grouped:
            self.assertEqual(result[key], len(group))

        df = DataFrame(np.random.choice(20, (1000, 3)), columns=list('abc'))
        for sort, key in cart_product((False, True), ('a', 'b', ['a', 'b'])):
            left = df.groupby(key, sort=sort).size()
            right = df.groupby(key, sort=sort)['c'].apply(lambda a: a.shape[0])
            assert_series_equal(left, right, check_names=False)

        # GH11699
        df = DataFrame([], columns=['A', 'B'])
        out = Series([], dtype='int64', index=Index([], name='A'))
        assert_series_equal(df.groupby('A').size(), out)

    def test_count(self):
        from string import ascii_lowercase
        n = 1 << 15
        dr = date_range('2015-08-30', periods=n // 10, freq='T')

        df = DataFrame({
            '1st': np.random.choice(
                list(ascii_lowercase), n),
            '2nd': np.random.randint(0, 5, n),
            '3rd': np.random.randn(n).round(3),
            '4th': np.random.randint(-10, 10, n),
            '5th': np.random.choice(dr, n),
            '6th': np.random.randn(n).round(3),
            '7th': np.random.randn(n).round(3),
            '8th': np.random.choice(dr, n) - np.random.choice(dr, 1),
            '9th': np.random.choice(
                list(ascii_lowercase), n)
        })

        for col in df.columns.drop(['1st', '2nd', '4th']):
            df.loc[np.random.choice(n, n // 10), col] = np.nan

        df['9th'] = df['9th'].astype('category')

        for key in '1st', '2nd', ['1st', '2nd']:
            left = df.groupby(key).count()
            right = df.groupby(key).apply(DataFrame.count).drop(key, axis=1)
            assert_frame_equal(left, right)

        # GH5610
        # count counts non-nulls
        df = pd.DataFrame([[1, 2, 'foo'], [1, nan, 'bar'], [3, nan, nan]],
                          columns=['A', 'B', 'C'])

        count_as = df.groupby('A').count()
        count_not_as = df.groupby('A', as_index=False).count()

        expected = DataFrame([[1, 2], [0, 0]], columns=['B', 'C'],
                             index=[1, 3])
        expected.index.name = 'A'
        assert_frame_equal(count_not_as, expected.reset_index())
        assert_frame_equal(count_as, expected)

        count_B = df.groupby('A')['B'].count()
        assert_series_equal(count_B, expected['B'])

    def test_count_object(self):
        df = pd.DataFrame({'a': ['a'] * 3 + ['b'] * 3, 'c': [2] * 3 + [3] * 3})
        result = df.groupby('c').a.count()
        expected = pd.Series([
            3, 3
        ], index=pd.Index([2, 3], name='c'), name='a')
        tm.assert_series_equal(result, expected)

        df = pd.DataFrame({'a': ['a', np.nan, np.nan] + ['b'] * 3,
                           'c': [2] * 3 + [3] * 3})
        result = df.groupby('c').a.count()
        expected = pd.Series([
            1, 3
        ], index=pd.Index([2, 3], name='c'), name='a')
        tm.assert_series_equal(result, expected)

    def test_count_cross_type(self):  # GH8169
        vals = np.hstack((np.random.randint(0, 5, (100, 2)), np.random.randint(
            0, 2, (100, 2))))

        df = pd.DataFrame(vals, columns=['a', 'b', 'c', 'd'])
        df[df == 2] = np.nan
        expected = df.groupby(['c', 'd']).count()

        for t in ['float32', 'object']:
            df['a'] = df['a'].astype(t)
            df['b'] = df['b'].astype(t)
            result = df.groupby(['c', 'd']).count()
            tm.assert_frame_equal(result, expected)

    def test_non_cython_api(self):

        # GH5610
        # non-cython calls should not include the grouper

        df = DataFrame(
            [[1, 2, 'foo'], [1,
                             nan,
                             'bar', ], [3, nan, 'baz']
             ], columns=['A', 'B', 'C'])
        g = df.groupby('A')
        gni = df.groupby('A', as_index=False)

        # mad
        expected = DataFrame([[0], [nan]], columns=['B'], index=[1, 3])
        expected.index.name = 'A'
        result = g.mad()
        assert_frame_equal(result, expected)

        expected = DataFrame([[0., 0.], [0, nan]], columns=['A', 'B'],
                             index=[0, 1])
        result = gni.mad()
        assert_frame_equal(result, expected)

        # describe
        expected = DataFrame(dict(B=concat(
            [df.loc[[0, 1], 'B'].describe(), df.loc[[2], 'B'].describe()],
            keys=[1, 3])))
        expected.index.names = ['A', None]
        result = g.describe()
        assert_frame_equal(result, expected)

        expected = concat(
            [df.loc[[0, 1], ['A', 'B']].describe(),
             df.loc[[2], ['A', 'B']].describe()], keys=[0, 1])
        result = gni.describe()
        assert_frame_equal(result, expected)

        # any
        expected = DataFrame([[True, True], [False, True]], columns=['B', 'C'],
                             index=[1, 3])
        expected.index.name = 'A'
        result = g.any()
        assert_frame_equal(result, expected)

        # idxmax
        expected = DataFrame([[0], [nan]], columns=['B'], index=[1, 3])
        expected.index.name = 'A'
        result = g.idxmax()
        assert_frame_equal(result, expected)

    def test_cython_api2(self):

        # this takes the fast apply path

        # cumsum (GH5614)
        df = DataFrame(
            [[1, 2, np.nan], [1, np.nan, 9], [3, 4, 9]
             ], columns=['A', 'B', 'C'])
        expected = DataFrame(
            [[2, np.nan], [np.nan, 9], [4, 9]], columns=['B', 'C'])
        result = df.groupby('A').cumsum()
        assert_frame_equal(result, expected)

        # GH 5755 - cumsum is a transformer and should ignore as_index
        result = df.groupby('A', as_index=False).cumsum()
        assert_frame_equal(result, expected)

        # GH 13994
        result = df.groupby('A').cumsum(axis=1)
        expected = df.cumsum(axis=1)
        assert_frame_equal(result, expected)
        result = df.groupby('A').cumprod(axis=1)
        expected = df.cumprod(axis=1)
        assert_frame_equal(result, expected)

    def test_grouping_ndarray(self):
        grouped = self.df.groupby(self.df['A'].values)

        result = grouped.sum()
        expected = self.df.groupby('A').sum()
        assert_frame_equal(result, expected, check_names=False
                           )  # Note: no names when grouping by value

    def test_agg_consistency(self):
        # agg with ([]) and () not consistent
        # GH 6715

        def P1(a):
            try:
                return np.percentile(a.dropna(), q=1)
            except:
                return np.nan

        import datetime as dt
        df = DataFrame({'col1': [1, 2, 3, 4],
                        'col2': [10, 25, 26, 31],
                        'date': [dt.date(2013, 2, 10), dt.date(2013, 2, 10),
                                 dt.date(2013, 2, 11), dt.date(2013, 2, 11)]})

        g = df.groupby('date')

        expected = g.agg([P1])
        expected.columns = expected.columns.levels[0]

        result = g.agg(P1)
        assert_frame_equal(result, expected)

    def test_apply_typecast_fail(self):
        df = DataFrame({'d': [1., 1., 1., 2., 2., 2.],
                        'c': np.tile(
                            ['a', 'b', 'c'], 2),
                        'v': np.arange(1., 7.)})

        def f(group):
            v = group['v']
            group['v2'] = (v - v.min()) / (v.max() - v.min())
            return group

        result = df.groupby('d').apply(f)

        expected = df.copy()
        expected['v2'] = np.tile([0., 0.5, 1], 2)

        assert_frame_equal(result, expected)

    def test_apply_multiindex_fail(self):
        index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1], [1, 2, 3, 1, 2, 3]
                                        ])
        df = DataFrame({'d': [1., 1., 1., 2., 2., 2.],
                        'c': np.tile(['a', 'b', 'c'], 2),
                        'v': np.arange(1., 7.)}, index=index)

        def f(group):
            v = group['v']
            group['v2'] = (v - v.min()) / (v.max() - v.min())
            return group

        result = df.groupby('d').apply(f)

        expected = df.copy()
        expected['v2'] = np.tile([0., 0.5, 1], 2)

        assert_frame_equal(result, expected)

    def test_apply_corner(self):
        result = self.tsframe.groupby(lambda x: x.year).apply(lambda x: x * 2)
        expected = self.tsframe * 2
        assert_frame_equal(result, expected)

    def test_apply_without_copy(self):
        # GH 5545
        # returning a non-copy in an applied function fails

        data = DataFrame({'id_field': [100, 100, 200, 300],
                          'category': ['a', 'b', 'c', 'c'],
                          'value': [1, 2, 3, 4]})

        def filt1(x):
            if x.shape[0] == 1:
                return x.copy()
            else:
                return x[x.category == 'c']

        def filt2(x):
            if x.shape[0] == 1:
                return x
            else:
                return x[x.category == 'c']

        expected = data.groupby('id_field').apply(filt1)
        result = data.groupby('id_field').apply(filt2)
        assert_frame_equal(result, expected)

    def test_apply_use_categorical_name(self):
        from pandas import qcut
        cats = qcut(self.df.C, 4)

        def get_stats(group):
            return {'min': group.min(),
                    'max': group.max(),
                    'count': group.count(),
                    'mean': group.mean()}

        result = self.df.groupby(cats).D.apply(get_stats)
        self.assertEqual(result.index.names[0], 'C')

    def test_apply_categorical_data(self):
        # GH 10138
        for ordered in [True, False]:
            dense = Categorical(list('abc'), ordered=ordered)
            # 'b' is in the categories but not in the list
            missing = Categorical(
                list('aaa'), categories=['a', 'b'], ordered=ordered)
            values = np.arange(len(dense))
            df = DataFrame({'missing': missing,
                            'dense': dense,
                            'values': values})
            grouped = df.groupby(['missing', 'dense'])

            # missing category 'b' should still exist in the output index
            idx = MultiIndex.from_product(
                [Categorical(['a', 'b'], ordered=ordered),
                 Categorical(['a', 'b', 'c'], ordered=ordered)],
                names=['missing', 'dense'])
            expected = DataFrame([0, 1, 2, np.nan, np.nan, np.nan],
                                 index=idx,
                                 columns=['values'])

            assert_frame_equal(grouped.apply(lambda x: np.mean(x)), expected)
            assert_frame_equal(grouped.mean(), expected)
            assert_frame_equal(grouped.agg(np.mean), expected)

            # but for transform we should still get back the original index
            idx = MultiIndex.from_product([['a'], ['a', 'b', 'c']],
                                          names=['missing', 'dense'])
            expected = Series(1, index=idx)
            assert_series_equal(grouped.apply(lambda x: 1), expected)

    def test_apply_corner_cases(self):
        # #535, can't use sliding iterator

        N = 1000
        labels = np.random.randint(0, 100, size=N)
        df = DataFrame({'key': labels,
                        'value1': np.random.randn(N),
                        'value2': ['foo', 'bar', 'baz', 'qux'] * (N // 4)})

        grouped = df.groupby('key')

        def f(g):
            g['value3'] = g['value1'] * 2
            return g

        result = grouped.apply(f)
        self.assertTrue('value3' in result)

    def test_transform_mixed_type(self):
        index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1], [1, 2, 3, 1, 2, 3]
                                        ])
        df = DataFrame({'d': [1., 1., 1., 2., 2., 2.],
                        'c': np.tile(['a', 'b', 'c'], 2),
                        'v': np.arange(1., 7.)}, index=index)

        def f(group):
            group['g'] = group['d'] * 2
            return group[:1]

        grouped = df.groupby('c')
        result = grouped.apply(f)

        self.assertEqual(result['d'].dtype, np.float64)

        # this is by definition a mutating operation!
        with option_context('mode.chained_assignment', None):
            for key, group in grouped:
                res = f(group)
                assert_frame_equal(res, result.ix[key])

    def test_groupby_wrong_multi_labels(self):
        from pandas import read_csv
        data = """index,foo,bar,baz,spam,data
0,foo1,bar1,baz1,spam2,20
1,foo1,bar2,baz1,spam3,30
2,foo2,bar2,baz1,spam2,40
3,foo1,bar1,baz2,spam1,50
4,foo3,bar1,baz2,spam1,60"""

        data = read_csv(StringIO(data), index_col=0)

        grouped = data.groupby(['foo', 'bar', 'baz', 'spam'])

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        assert_frame_equal(result, expected)

    def test_groupby_series_with_name(self):
        result = self.df.groupby(self.df['A']).mean()
        result2 = self.df.groupby(self.df['A'], as_index=False).mean()
        self.assertEqual(result.index.name, 'A')
        self.assertIn('A', result2)

        result = self.df.groupby([self.df['A'], self.df['B']]).mean()
        result2 = self.df.groupby([self.df['A'], self.df['B']],
                                  as_index=False).mean()
        self.assertEqual(result.index.names, ('A', 'B'))
        self.assertIn('A', result2)
        self.assertIn('B', result2)

    def test_seriesgroupby_name_attr(self):
        # GH 6265
        result = self.df.groupby('A')['C']
        self.assertEqual(result.count().name, 'C')
        self.assertEqual(result.mean().name, 'C')

        testFunc = lambda x: np.sum(x) * 2
        self.assertEqual(result.agg(testFunc).name, 'C')

    def test_consistency_name(self):
        # GH 12363

        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': np.random.randn(8) + 1.0,
                        'D': np.arange(8)})

        expected = df.groupby(['A']).B.count()
        result = df.B.groupby(df.A).count()
        assert_series_equal(result, expected)

    def test_groupby_name_propagation(self):
        # GH 6124
        def summarize(df, name=None):
            return Series({'count': 1, 'mean': 2, 'omissions': 3, }, name=name)

        def summarize_random_name(df):
            # Provide a different name for each Series.  In this case, groupby
            # should not attempt to propagate the Series name since they are
            # inconsistent.
            return Series({
                'count': 1,
                'mean': 2,
                'omissions': 3,
            }, name=df.iloc[0]['A'])

        metrics = self.df.groupby('A').apply(summarize)
        self.assertEqual(metrics.columns.name, None)
        metrics = self.df.groupby('A').apply(summarize, 'metrics')
        self.assertEqual(metrics.columns.name, 'metrics')
        metrics = self.df.groupby('A').apply(summarize_random_name)
        self.assertEqual(metrics.columns.name, None)

    def test_groupby_nonstring_columns(self):
        df = DataFrame([np.arange(10) for x in range(10)])
        grouped = df.groupby(0)
        result = grouped.mean()
        expected = df.groupby(df[0]).mean()
        assert_frame_equal(result, expected)

    def test_groupby_mixed_type_columns(self):
        # GH 13432, unorderable types in py3
        df = DataFrame([[0, 1, 2]], columns=['A', 'B', 0])
        expected = DataFrame([[1, 2]], columns=['B', 0],
                             index=Index([0], name='A'))

        result = df.groupby('A').first()
        tm.assert_frame_equal(result, expected)

        result = df.groupby('A').sum()
        tm.assert_frame_equal(result, expected)

    def test_cython_grouper_series_bug_noncontig(self):
        arr = np.empty((100, 100))
        arr.fill(np.nan)
        obj = Series(arr[:, 0], index=lrange(100))
        inds = np.tile(lrange(10), 10)

        result = obj.groupby(inds).agg(Series.median)
        self.assertTrue(result.isnull().all())

    def test_series_grouper_noncontig_index(self):
        index = Index(tm.rands_array(10, 100))

        values = Series(np.random.randn(50), index=index[::2])
        labels = np.random.randint(0, 5, 50)

        # it works!
        grouped = values.groupby(labels)

        # accessing the index elements causes segfault
        f = lambda x: len(set(map(id, x.index)))
        grouped.agg(f)

    def test_convert_objects_leave_decimal_alone(self):

        from decimal import Decimal

        s = Series(lrange(5))
        labels = np.array(['a', 'b', 'c', 'd', 'e'], dtype='O')

        def convert_fast(x):
            return Decimal(str(x.mean()))

        def convert_force_pure(x):
            # base will be length 0
            assert (len(x.base) > 0)
            return Decimal(str(x.mean()))

        grouped = s.groupby(labels)

        result = grouped.agg(convert_fast)
        self.assertEqual(result.dtype, np.object_)
        tm.assertIsInstance(result[0], Decimal)

        result = grouped.agg(convert_force_pure)
        self.assertEqual(result.dtype, np.object_)
        tm.assertIsInstance(result[0], Decimal)

    def test_fast_apply(self):
        # make sure that fast apply is correctly called
        # rather than raising any kind of error
        # otherwise the python path will be callsed
        # which slows things down
        N = 1000
        labels = np.random.randint(0, 2000, size=N)
        labels2 = np.random.randint(0, 3, size=N)
        df = DataFrame({'key': labels,
                        'key2': labels2,
                        'value1': np.random.randn(N),
                        'value2': ['foo', 'bar', 'baz', 'qux'] * (N // 4)})

        def f(g):
            return 1

        g = df.groupby(['key', 'key2'])

        grouper = g.grouper

        splitter = grouper._get_splitter(g._selected_obj, axis=g.axis)
        group_keys = grouper._get_group_keys()

        values, mutated = splitter.fast_apply(f, group_keys)
        self.assertFalse(mutated)

    def test_apply_with_mixed_dtype(self):
        # GH3480, apply with mixed dtype on axis=1 breaks in 0.11
        df = DataFrame({'foo1': ['one', 'two', 'two', 'three', 'one', 'two'],
                        'foo2': np.random.randn(6)})
        result = df.apply(lambda x: x, axis=1)
        assert_series_equal(df.get_dtype_counts(), result.get_dtype_counts())

        # GH 3610 incorrect dtype conversion with as_index=False
        df = DataFrame({"c1": [1, 2, 6, 6, 8]})
        df["c2"] = df.c1 / 2.0
        result1 = df.groupby("c2").mean().reset_index().c2
        result2 = df.groupby("c2", as_index=False).mean().c2
        assert_series_equal(result1, result2)

    def test_groupby_aggregation_mixed_dtype(self):

        # GH 6212
        expected = DataFrame({
            'v1': [5, 5, 7, np.nan, 3, 3, 4, 1],
            'v2': [55, 55, 77, np.nan, 33, 33, 44, 11]},
            index=MultiIndex.from_tuples([(1, 95), (1, 99), (2, 95), (2, 99),
                                          ('big', 'damp'),
                                          ('blue', 'dry'),
                                          ('red', 'red'), ('red', 'wet')],
                                         names=['by1', 'by2']))

        df = DataFrame({
            'v1': [1, 3, 5, 7, 8, 3, 5, np.nan, 4, 5, 7, 9],
            'v2': [11, 33, 55, 77, 88, 33, 55, np.nan, 44, 55, 77, 99],
            'by1': ["red", "blue", 1, 2, np.nan, "big", 1, 2, "red", 1, np.nan,
                    12],
            'by2': ["wet", "dry", 99, 95, np.nan, "damp", 95, 99, "red", 99,
                    np.nan, np.nan]
        })

        g = df.groupby(['by1', 'by2'])
        result = g[['v1', 'v2']].mean()
        assert_frame_equal(result, expected)

    def test_groupby_dtype_inference_empty(self):
        # GH 6733
        df = DataFrame({'x': [], 'range': np.arange(0, dtype='int64')})
        self.assertEqual(df['x'].dtype, np.float64)

        result = df.groupby('x').first()
        exp_index = Index([], name='x', dtype=np.float64)
        expected = DataFrame({'range': Series(
            [], index=exp_index, dtype='int64')})
        assert_frame_equal(result, expected, by_blocks=True)

    def test_groupby_list_infer_array_like(self):
        result = self.df.groupby(list(self.df['A'])).mean()
        expected = self.df.groupby(self.df['A']).mean()
        assert_frame_equal(result, expected, check_names=False)

        self.assertRaises(Exception, self.df.groupby, list(self.df['A'][:-1]))

        # pathological case of ambiguity
        df = DataFrame({'foo': [0, 1],
                        'bar': [3, 4],
                        'val': np.random.randn(2)})

        result = df.groupby(['foo', 'bar']).mean()
        expected = df.groupby([df['foo'], df['bar']]).mean()[['val']]

    def test_groupby_keys_same_size_as_index(self):
        # GH 11185
        freq = 's'
        index = pd.date_range(start=pd.Timestamp('2015-09-29T11:34:44-0700'),
                              periods=2, freq=freq)
        df = pd.DataFrame([['A', 10], ['B', 15]], columns=[
            'metric', 'values'
        ], index=index)
        result = df.groupby([pd.Grouper(level=0, freq=freq), 'metric']).mean()
        expected = df.set_index([df.index, 'metric'])

        assert_frame_equal(result, expected)

    def test_groupby_one_row(self):
        # GH 11741
        df1 = pd.DataFrame(np.random.randn(1, 4), columns=list('ABCD'))
        self.assertRaises(KeyError, df1.groupby, 'Z')
        df2 = pd.DataFrame(np.random.randn(2, 4), columns=list('ABCD'))
        self.assertRaises(KeyError, df2.groupby, 'Z')

    def test_groupby_nat_exclude(self):
        # GH 6992
        df = pd.DataFrame(
            {'values': np.random.randn(8),
             'dt': [np.nan, pd.Timestamp('2013-01-01'), np.nan, pd.Timestamp(
                 '2013-02-01'), np.nan, pd.Timestamp('2013-02-01'), np.nan,
                pd.Timestamp('2013-01-01')],
             'str': [np.nan, 'a', np.nan, 'a', np.nan, 'a', np.nan, 'b']})
        grouped = df.groupby('dt')

        expected = [pd.Index([1, 7]), pd.Index([3, 5])]
        keys = sorted(grouped.groups.keys())
        self.assertEqual(len(keys), 2)
        for k, e in zip(keys, expected):
            # grouped.groups keys are np.datetime64 with system tz
            # not to be affected by tz, only compare values
            tm.assert_index_equal(grouped.groups[k], e)

        # confirm obj is not filtered
        tm.assert_frame_equal(grouped.grouper.groupings[0].obj, df)
        self.assertEqual(grouped.ngroups, 2)

        expected = {
            Timestamp('2013-01-01 00:00:00'): np.array([1, 7], dtype=np.int64),
            Timestamp('2013-02-01 00:00:00'): np.array([3, 5], dtype=np.int64)
        }

        for k in grouped.indices:
            self.assert_numpy_array_equal(grouped.indices[k], expected[k])

        tm.assert_frame_equal(
            grouped.get_group(Timestamp('2013-01-01')), df.iloc[[1, 7]])
        tm.assert_frame_equal(
            grouped.get_group(Timestamp('2013-02-01')), df.iloc[[3, 5]])

        self.assertRaises(KeyError, grouped.get_group, pd.NaT)

        nan_df = DataFrame({'nan': [np.nan, np.nan, np.nan],
                            'nat': [pd.NaT, pd.NaT, pd.NaT]})
        self.assertEqual(nan_df['nan'].dtype, 'float64')
        self.assertEqual(nan_df['nat'].dtype, 'datetime64[ns]')

        for key in ['nan', 'nat']:
            grouped = nan_df.groupby(key)
            self.assertEqual(grouped.groups, {})
            self.assertEqual(grouped.ngroups, 0)
            self.assertEqual(grouped.indices, {})
            self.assertRaises(KeyError, grouped.get_group, np.nan)
            self.assertRaises(KeyError, grouped.get_group, pd.NaT)

    def test_dictify(self):
        dict(iter(self.df.groupby('A')))
        dict(iter(self.df.groupby(['A', 'B'])))
        dict(iter(self.df['C'].groupby(self.df['A'])))
        dict(iter(self.df['C'].groupby([self.df['A'], self.df['B']])))
        dict(iter(self.df.groupby('A')['C']))
        dict(iter(self.df.groupby(['A', 'B'])['C']))

    def test_sparse_friendly(self):
        sdf = self.df[['C', 'D']].to_sparse()
        panel = tm.makePanel()
        tm.add_nans(panel)

        def _check_work(gp):
            gp.mean()
            gp.agg(np.mean)
            dict(iter(gp))

        # it works!
        _check_work(sdf.groupby(lambda x: x // 2))
        _check_work(sdf['C'].groupby(lambda x: x // 2))
        _check_work(sdf.groupby(self.df['A']))

        # do this someday
        # _check_work(panel.groupby(lambda x: x.month, axis=1))

    def test_panel_groupby(self):
        self.panel = tm.makePanel()
        tm.add_nans(self.panel)
        grouped = self.panel.groupby({'ItemA': 0, 'ItemB': 0, 'ItemC': 1},
                                     axis='items')
        agged = grouped.mean()
        agged2 = grouped.agg(lambda x: x.mean('items'))

        tm.assert_panel_equal(agged, agged2)

        self.assert_index_equal(agged.items, Index([0, 1]))

        grouped = self.panel.groupby(lambda x: x.month, axis='major')
        agged = grouped.mean()

        exp = Index(sorted(list(set(self.panel.major_axis.month))))
        self.assert_index_equal(agged.major_axis, exp)

        grouped = self.panel.groupby({'A': 0, 'B': 0, 'C': 1, 'D': 1},
                                     axis='minor')
        agged = grouped.mean()
        self.assert_index_equal(agged.minor_axis, Index([0, 1]))

    def test_numpy_groupby(self):
        from pandas.core.groupby import numpy_groupby

        data = np.random.randn(100, 100)
        labels = np.random.randint(0, 10, size=100)

        df = DataFrame(data)

        result = df.groupby(labels).sum().values
        expected = numpy_groupby(data, labels)
        assert_almost_equal(result, expected)

        result = df.groupby(labels, axis=1).sum().values
        expected = numpy_groupby(data, labels, axis=1)
        assert_almost_equal(result, expected)

    def test_groupby_2d_malformed(self):
        d = DataFrame(index=lrange(2))
        d['group'] = ['g1', 'g2']
        d['zeros'] = [0, 0]
        d['ones'] = [1, 1]
        d['label'] = ['l1', 'l2']
        tmp = d.groupby(['group']).mean()
        res_values = np.array([[0, 1], [0, 1]], dtype=np.int64)
        self.assert_index_equal(tmp.columns, Index(['zeros', 'ones']))
        self.assert_numpy_array_equal(tmp.values, res_values)

    def test_int32_overflow(self):
        B = np.concatenate((np.arange(10000), np.arange(10000), np.arange(5000)
                            ))
        A = np.arange(25000)
        df = DataFrame({'A': A,
                        'B': B,
                        'C': A,
                        'D': B,
                        'E': np.random.randn(25000)})

        left = df.groupby(['A', 'B', 'C', 'D']).sum()
        right = df.groupby(['D', 'C', 'B', 'A']).sum()
        self.assertEqual(len(left), len(right))

    def test_int64_overflow(self):
        from pandas.core.groupby import _int64_overflow_possible

        B = np.concatenate((np.arange(1000), np.arange(1000), np.arange(500)))
        A = np.arange(2500)
        df = DataFrame({'A': A,
                        'B': B,
                        'C': A,
                        'D': B,
                        'E': A,
                        'F': B,
                        'G': A,
                        'H': B,
                        'values': np.random.randn(2500)})

        lg = df.groupby(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
        rg = df.groupby(['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A'])

        left = lg.sum()['values']
        right = rg.sum()['values']

        exp_index, _ = left.index.sortlevel(0)
        self.assert_index_equal(left.index, exp_index)

        exp_index, _ = right.index.sortlevel(0)
        self.assert_index_equal(right.index, exp_index)

        tups = list(map(tuple, df[['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'
                                   ]].values))
        tups = com._asarray_tuplesafe(tups)
        expected = df.groupby(tups).sum()['values']

        for k, v in compat.iteritems(expected):
            self.assertEqual(left[k], right[k[::-1]])
            self.assertEqual(left[k], v)
        self.assertEqual(len(left), len(right))

        # GH9096
        values = range(55109)
        data = pd.DataFrame.from_dict({'a': values,
                                       'b': values,
                                       'c': values,
                                       'd': values})
        grouped = data.groupby(['a', 'b', 'c', 'd'])
        self.assertEqual(len(grouped), len(values))

        arr = np.random.randint(-1 << 12, 1 << 12, (1 << 15, 5))
        i = np.random.choice(len(arr), len(arr) * 4)
        arr = np.vstack((arr, arr[i]))  # add sume duplicate rows

        i = np.random.permutation(len(arr))
        arr = arr[i]  # shuffle rows

        df = DataFrame(arr, columns=list('abcde'))
        df['jim'], df['joe'] = np.random.randn(2, len(df)) * 10
        gr = df.groupby(list('abcde'))

        # verify this is testing what it is supposed to test!
        self.assertTrue(_int64_overflow_possible(gr.grouper.shape))

        # mannually compute groupings
        jim, joe = defaultdict(list), defaultdict(list)
        for key, a, b in zip(map(tuple, arr), df['jim'], df['joe']):
            jim[key].append(a)
            joe[key].append(b)

        self.assertEqual(len(gr), len(jim))
        mi = MultiIndex.from_tuples(jim.keys(), names=list('abcde'))

        def aggr(func):
            f = lambda a: np.fromiter(map(func, a), dtype='f8')
            arr = np.vstack((f(jim.values()), f(joe.values()))).T
            res = DataFrame(arr, columns=['jim', 'joe'], index=mi)
            return res.sort_index()

        assert_frame_equal(gr.mean(), aggr(np.mean))
        assert_frame_equal(gr.median(), aggr(np.median))

    def test_groupby_sort_multi(self):
        df = DataFrame({'a': ['foo', 'bar', 'baz'],
                        'b': [3, 2, 1],
                        'c': [0, 1, 2],
                        'd': np.random.randn(3)})

        tups = lmap(tuple, df[['a', 'b', 'c']].values)
        tups = com._asarray_tuplesafe(tups)
        result = df.groupby(['a', 'b', 'c'], sort=True).sum()
        self.assert_numpy_array_equal(result.index.values, tups[[1, 2, 0]])

        tups = lmap(tuple, df[['c', 'a', 'b']].values)
        tups = com._asarray_tuplesafe(tups)
        result = df.groupby(['c', 'a', 'b'], sort=True).sum()
        self.assert_numpy_array_equal(result.index.values, tups)

        tups = lmap(tuple, df[['b', 'c', 'a']].values)
        tups = com._asarray_tuplesafe(tups)
        result = df.groupby(['b', 'c', 'a'], sort=True).sum()
        self.assert_numpy_array_equal(result.index.values, tups[[2, 1, 0]])

        df = DataFrame({'a': [0, 1, 2, 0, 1, 2],
                        'b': [0, 0, 0, 1, 1, 1],
                        'd': np.random.randn(6)})
        grouped = df.groupby(['a', 'b'])['d']
        result = grouped.sum()
        _check_groupby(df, result, ['a', 'b'], 'd')

    def test_intercept_builtin_sum(self):
        s = Series([1., 2., np.nan, 3.])
        grouped = s.groupby([0, 1, 2, 2])

        result = grouped.agg(builtins.sum)
        result2 = grouped.apply(builtins.sum)
        expected = grouped.sum()
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

    def test_column_select_via_attr(self):
        result = self.df.groupby('A').C.sum()
        expected = self.df.groupby('A')['C'].sum()
        assert_series_equal(result, expected)

        self.df['mean'] = 1.5
        result = self.df.groupby('A').mean()
        expected = self.df.groupby('A').agg(np.mean)
        assert_frame_equal(result, expected)

    def test_rank_apply(self):
        lev1 = tm.rands_array(10, 100)
        lev2 = tm.rands_array(10, 130)
        lab1 = np.random.randint(0, 100, size=500)
        lab2 = np.random.randint(0, 130, size=500)

        df = DataFrame({'value': np.random.randn(500),
                        'key1': lev1.take(lab1),
                        'key2': lev2.take(lab2)})

        result = df.groupby(['key1', 'key2']).value.rank()

        expected = []
        for key, piece in df.groupby(['key1', 'key2']):
            expected.append(piece.value.rank())
        expected = concat(expected, axis=0)
        expected = expected.reindex(result.index)
        assert_series_equal(result, expected)

        result = df.groupby(['key1', 'key2']).value.rank(pct=True)

        expected = []
        for key, piece in df.groupby(['key1', 'key2']):
            expected.append(piece.value.rank(pct=True))
        expected = concat(expected, axis=0)
        expected = expected.reindex(result.index)
        assert_series_equal(result, expected)

    def test_dont_clobber_name_column(self):
        df = DataFrame({'key': ['a', 'a', 'a', 'b', 'b', 'b'],
                        'name': ['foo', 'bar', 'baz'] * 2})

        result = df.groupby('key').apply(lambda x: x)
        assert_frame_equal(result, df)

    def test_skip_group_keys(self):
        from pandas import concat

        tsf = tm.makeTimeDataFrame()

        grouped = tsf.groupby(lambda x: x.month, group_keys=False)
        result = grouped.apply(lambda x: x.sort_values(by='A')[:3])

        pieces = []
        for key, group in grouped:
            pieces.append(group.sort_values(by='A')[:3])

        expected = concat(pieces)
        assert_frame_equal(result, expected)

        grouped = tsf['A'].groupby(lambda x: x.month, group_keys=False)
        result = grouped.apply(lambda x: x.sort_values()[:3])

        pieces = []
        for key, group in grouped:
            pieces.append(group.sort_values()[:3])

        expected = concat(pieces)
        assert_series_equal(result, expected)

    def test_no_nonsense_name(self):
        # GH #995
        s = self.frame['C'].copy()
        s.name = None

        result = s.groupby(self.frame['A']).agg(np.sum)
        self.assertIsNone(result.name)

    def test_wrap_agg_out(self):
        grouped = self.three_group.groupby(['A', 'B'])

        def func(ser):
            if ser.dtype == np.object:
                raise TypeError
            else:
                return ser.sum()

        result = grouped.aggregate(func)
        exp_grouped = self.three_group.ix[:, self.three_group.columns != 'C']
        expected = exp_grouped.groupby(['A', 'B']).aggregate(func)
        assert_frame_equal(result, expected)

    def test_multifunc_sum_bug(self):
        # GH #1065
        x = DataFrame(np.arange(9).reshape(3, 3))
        x['test'] = 0
        x['fl'] = [1.3, 1.5, 1.6]

        grouped = x.groupby('test')
        result = grouped.agg({'fl': 'sum', 2: 'size'})
        self.assertEqual(result['fl'].dtype, np.float64)

    def test_handle_dict_return_value(self):
        def f(group):
            return {'min': group.min(), 'max': group.max()}

        def g(group):
            return Series({'min': group.min(), 'max': group.max()})

        result = self.df.groupby('A')['C'].apply(f)
        expected = self.df.groupby('A')['C'].apply(g)

        tm.assertIsInstance(result, Series)
        assert_series_equal(result, expected)

    def test_getitem_list_of_columns(self):
        df = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.random.randn(8),
             'E': np.random.randn(8)})

        result = df.groupby('A')[['C', 'D']].mean()
        result2 = df.groupby('A')['C', 'D'].mean()
        result3 = df.groupby('A')[df.columns[2:4]].mean()

        expected = df.ix[:, ['A', 'C', 'D']].groupby('A').mean()

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)
        assert_frame_equal(result3, expected)

    def test_getitem_numeric_column_names(self):
        # GH #13731
        df = DataFrame({0: list('abcd') * 2,
                        2: np.random.randn(8),
                        4: np.random.randn(8),
                        6: np.random.randn(8)})
        result = df.groupby(0)[df.columns[1:3]].mean()
        result2 = df.groupby(0)[2, 4].mean()
        result3 = df.groupby(0)[[2, 4]].mean()

        expected = df.ix[:, [0, 2, 4]].groupby(0).mean()

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)
        assert_frame_equal(result3, expected)

    def test_agg_multiple_functions_maintain_order(self):
        # GH #610
        funcs = [('mean', np.mean), ('max', np.max), ('min', np.min)]
        result = self.df.groupby('A')['C'].agg(funcs)
        exp_cols = Index(['mean', 'max', 'min'])

        self.assert_index_equal(result.columns, exp_cols)

    def test_multiple_functions_tuples_and_non_tuples(self):
        # #1359

        funcs = [('foo', 'mean'), 'std']
        ex_funcs = [('foo', 'mean'), ('std', 'std')]

        result = self.df.groupby('A')['C'].agg(funcs)
        expected = self.df.groupby('A')['C'].agg(ex_funcs)
        assert_frame_equal(result, expected)

        result = self.df.groupby('A').agg(funcs)
        expected = self.df.groupby('A').agg(ex_funcs)
        assert_frame_equal(result, expected)

    def test_agg_multiple_functions_too_many_lambdas(self):
        grouped = self.df.groupby('A')
        funcs = ['mean', lambda x: x.mean(), lambda x: x.std()]

        self.assertRaises(SpecificationError, grouped.agg, funcs)

    def test_more_flexible_frame_multi_function(self):
        from pandas import concat

        grouped = self.df.groupby('A')

        exmean = grouped.agg(OrderedDict([['C', np.mean], ['D', np.mean]]))
        exstd = grouped.agg(OrderedDict([['C', np.std], ['D', np.std]]))

        expected = concat([exmean, exstd], keys=['mean', 'std'], axis=1)
        expected = expected.swaplevel(0, 1, axis=1).sortlevel(0, axis=1)

        d = OrderedDict([['C', [np.mean, np.std]], ['D', [np.mean, np.std]]])
        result = grouped.aggregate(d)

        assert_frame_equal(result, expected)

        # be careful
        result = grouped.aggregate(OrderedDict([['C', np.mean],
                                                ['D', [np.mean, np.std]]]))
        expected = grouped.aggregate(OrderedDict([['C', np.mean],
                                                  ['D', [np.mean, np.std]]]))
        assert_frame_equal(result, expected)

        def foo(x):
            return np.mean(x)

        def bar(x):
            return np.std(x, ddof=1)

        d = OrderedDict([['C', np.mean], ['D', OrderedDict(
            [['foo', np.mean], ['bar', np.std]])]])
        result = grouped.aggregate(d)

        d = OrderedDict([['C', [np.mean]], ['D', [foo, bar]]])
        expected = grouped.aggregate(d)

        assert_frame_equal(result, expected)

    def test_multi_function_flexible_mix(self):
        # GH #1268
        grouped = self.df.groupby('A')

        d = OrderedDict([['C', OrderedDict([['foo', 'mean'], [
            'bar', 'std'
        ]])], ['D', 'sum']])
        result = grouped.aggregate(d)
        d2 = OrderedDict([['C', OrderedDict([['foo', 'mean'], [
            'bar', 'std'
        ]])], ['D', ['sum']]])
        result2 = grouped.aggregate(d2)

        d3 = OrderedDict([['C', OrderedDict([['foo', 'mean'], [
            'bar', 'std'
        ]])], ['D', {'sum': 'sum'}]])
        expected = grouped.aggregate(d3)

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

    def test_agg_callables(self):
        # GH 7929
        df = DataFrame({'foo': [1, 2], 'bar': [3, 4]}).astype(np.int64)

        class fn_class(object):

            def __call__(self, x):
                return sum(x)

        equiv_callables = [sum, np.sum, lambda x: sum(x), lambda x: x.sum(),
                           partial(sum), fn_class()]

        expected = df.groupby("foo").agg(sum)
        for ecall in equiv_callables:
            result = df.groupby('foo').agg(ecall)
            assert_frame_equal(result, expected)

    def test_set_group_name(self):
        def f(group):
            assert group.name is not None
            return group

        def freduce(group):
            assert group.name is not None
            return group.sum()

        def foo(x):
            return freduce(x)

        def _check_all(grouped):
            # make sure all these work
            grouped.apply(f)
            grouped.aggregate(freduce)
            grouped.aggregate({'C': freduce, 'D': freduce})
            grouped.transform(f)

            grouped['C'].apply(f)
            grouped['C'].aggregate(freduce)
            grouped['C'].aggregate([freduce, foo])
            grouped['C'].transform(f)

        _check_all(self.df.groupby('A'))
        _check_all(self.df.groupby(['A', 'B']))

    def test_no_dummy_key_names(self):
        # GH #1291

        result = self.df.groupby(self.df['A'].values).sum()
        self.assertIsNone(result.index.name)

        result = self.df.groupby([self.df['A'].values, self.df['B'].values
                                  ]).sum()
        self.assertEqual(result.index.names, (None, None))

    def test_groupby_sort_categorical(self):
        # dataframe groupby sort was being ignored # GH 8868
        df = DataFrame([['(7.5, 10]', 10, 10],
                        ['(7.5, 10]', 8, 20],
                        ['(2.5, 5]', 5, 30],
                        ['(5, 7.5]', 6, 40],
                        ['(2.5, 5]', 4, 50],
                        ['(0, 2.5]', 1, 60],
                        ['(5, 7.5]', 7, 70]], columns=['range', 'foo', 'bar'])
        df['range'] = Categorical(df['range'], ordered=True)
        index = CategoricalIndex(['(0, 2.5]', '(2.5, 5]', '(5, 7.5]',
                                  '(7.5, 10]'], name='range', ordered=True)
        result_sort = DataFrame([[1, 60], [5, 30], [6, 40], [10, 10]],
                                columns=['foo', 'bar'], index=index)

        col = 'range'
        assert_frame_equal(result_sort, df.groupby(col, sort=True).first())
        # when categories is ordered, group is ordered by category's order
        assert_frame_equal(result_sort, df.groupby(col, sort=False).first())

        df['range'] = Categorical(df['range'], ordered=False)
        index = CategoricalIndex(['(0, 2.5]', '(2.5, 5]', '(5, 7.5]',
                                  '(7.5, 10]'], name='range')
        result_sort = DataFrame([[1, 60], [5, 30], [6, 40], [10, 10]],
                                columns=['foo', 'bar'], index=index)

        index = CategoricalIndex(['(7.5, 10]', '(2.5, 5]', '(5, 7.5]',
                                  '(0, 2.5]'],
                                 categories=['(7.5, 10]', '(2.5, 5]',
                                             '(5, 7.5]', '(0, 2.5]'],
                                 name='range')
        result_nosort = DataFrame([[10, 10], [5, 30], [6, 40], [1, 60]],
                                  index=index, columns=['foo', 'bar'])

        col = 'range'
        # this is an unordered categorical, but we allow this ####
        assert_frame_equal(result_sort, df.groupby(col, sort=True).first())
        assert_frame_equal(result_nosort, df.groupby(col, sort=False).first())

    def test_groupby_sort_categorical_datetimelike(self):
        # GH10505

        # use same data as test_groupby_sort_categorical, which category is
        # corresponding to datetime.month
        df = DataFrame({'dt': [datetime(2011, 7, 1), datetime(2011, 7, 1),
                               datetime(2011, 2, 1), datetime(2011, 5, 1),
                               datetime(2011, 2, 1), datetime(2011, 1, 1),
                               datetime(2011, 5, 1)],
                        'foo': [10, 8, 5, 6, 4, 1, 7],
                        'bar': [10, 20, 30, 40, 50, 60, 70]},
                       columns=['dt', 'foo', 'bar'])

        # ordered=True
        df['dt'] = Categorical(df['dt'], ordered=True)
        index = [datetime(2011, 1, 1), datetime(2011, 2, 1),
                 datetime(2011, 5, 1), datetime(2011, 7, 1)]
        result_sort = DataFrame(
            [[1, 60], [5, 30], [6, 40], [10, 10]], columns=['foo', 'bar'])
        result_sort.index = CategoricalIndex(index, name='dt', ordered=True)

        index = [datetime(2011, 7, 1), datetime(2011, 2, 1),
                 datetime(2011, 5, 1), datetime(2011, 1, 1)]
        result_nosort = DataFrame([[10, 10], [5, 30], [6, 40], [1, 60]],
                                  columns=['foo', 'bar'])
        result_nosort.index = CategoricalIndex(index, categories=index,
                                               name='dt', ordered=True)

        col = 'dt'
        assert_frame_equal(result_sort, df.groupby(col, sort=True).first())
        # when categories is ordered, group is ordered by category's order
        assert_frame_equal(result_sort, df.groupby(col, sort=False).first())

        # ordered = False
        df['dt'] = Categorical(df['dt'], ordered=False)
        index = [datetime(2011, 1, 1), datetime(2011, 2, 1),
                 datetime(2011, 5, 1), datetime(2011, 7, 1)]
        result_sort = DataFrame(
            [[1, 60], [5, 30], [6, 40], [10, 10]], columns=['foo', 'bar'])
        result_sort.index = CategoricalIndex(index, name='dt')

        index = [datetime(2011, 7, 1), datetime(2011, 2, 1),
                 datetime(2011, 5, 1), datetime(2011, 1, 1)]
        result_nosort = DataFrame([[10, 10], [5, 30], [6, 40], [1, 60]],
                                  columns=['foo', 'bar'])
        result_nosort.index = CategoricalIndex(index, categories=index,
                                               name='dt')

        col = 'dt'
        assert_frame_equal(result_sort, df.groupby(col, sort=True).first())
        assert_frame_equal(result_nosort, df.groupby(col, sort=False).first())

    def test_groupby_sort_multiindex_series(self):
        # series multiindex groupby sort argument was not being passed through
        # _compress_group_index
        # GH 9444
        index = MultiIndex(levels=[[1, 2], [1, 2]],
                           labels=[[0, 0, 0, 0, 1, 1], [1, 1, 0, 0, 0, 0]],
                           names=['a', 'b'])
        mseries = Series([0, 1, 2, 3, 4, 5], index=index)
        index = MultiIndex(levels=[[1, 2], [1, 2]],
                           labels=[[0, 0, 1], [1, 0, 0]], names=['a', 'b'])
        mseries_result = Series([0, 2, 4], index=index)

        result = mseries.groupby(level=['a', 'b'], sort=False).first()
        assert_series_equal(result, mseries_result)
        result = mseries.groupby(level=['a', 'b'], sort=True).first()
        assert_series_equal(result, mseries_result.sort_index())

    def test_groupby_categorical(self):
        levels = ['foo', 'bar', 'baz', 'qux']
        codes = np.random.randint(0, 4, size=100)

        cats = Categorical.from_codes(codes, levels, ordered=True)

        data = DataFrame(np.random.randn(100, 4))

        result = data.groupby(cats).mean()

        expected = data.groupby(np.asarray(cats)).mean()
        exp_idx = CategoricalIndex(levels, categories=cats.categories,
                                   ordered=True)
        expected = expected.reindex(exp_idx)

        assert_frame_equal(result, expected)

        grouped = data.groupby(cats)
        desc_result = grouped.describe()

        idx = cats.codes.argsort()
        ord_labels = np.asarray(cats).take(idx)
        ord_data = data.take(idx)

        exp_cats = Categorical(ord_labels, ordered=True,
                               categories=['foo', 'bar', 'baz', 'qux'])
        expected = ord_data.groupby(exp_cats, sort=False).describe()
        expected.index.names = [None, None]
        assert_frame_equal(desc_result, expected)

        # GH 10460
        expc = Categorical.from_codes(np.arange(4).repeat(8),
                                      levels, ordered=True)
        exp = CategoricalIndex(expc)
        self.assert_index_equal(desc_result.index.get_level_values(0), exp)
        exp = Index(['count', 'mean', 'std', 'min', '25%', '50%',
                     '75%', 'max'] * 4)
        self.assert_index_equal(desc_result.index.get_level_values(1), exp)

    def test_groupby_datetime_categorical(self):
        # GH9049: ensure backward compatibility
        levels = pd.date_range('2014-01-01', periods=4)
        codes = np.random.randint(0, 4, size=100)

        cats = Categorical.from_codes(codes, levels, ordered=True)

        data = DataFrame(np.random.randn(100, 4))
        result = data.groupby(cats).mean()

        expected = data.groupby(np.asarray(cats)).mean()
        expected = expected.reindex(levels)
        expected.index = CategoricalIndex(expected.index,
                                          categories=expected.index,
                                          ordered=True)

        assert_frame_equal(result, expected)

        grouped = data.groupby(cats)
        desc_result = grouped.describe()

        idx = cats.codes.argsort()
        ord_labels = cats.take_nd(idx)
        ord_data = data.take(idx)
        expected = ord_data.groupby(ord_labels).describe()
        expected.index.names = [None, None]
        assert_frame_equal(desc_result, expected)
        tm.assert_index_equal(desc_result.index, expected.index)
        tm.assert_index_equal(
            desc_result.index.get_level_values(0),
            expected.index.get_level_values(0))

        # GH 10460
        expc = Categorical.from_codes(
            np.arange(4).repeat(8), levels, ordered=True)
        exp = CategoricalIndex(expc)
        self.assert_index_equal(desc_result.index.get_level_values(0), exp)
        exp = Index(['count', 'mean', 'std', 'min', '25%', '50%',
                     '75%', 'max'] * 4)
        self.assert_index_equal(desc_result.index.get_level_values(1), exp)

    def test_groupby_categorical_index(self):

        levels = ['foo', 'bar', 'baz', 'qux']
        codes = np.random.randint(0, 4, size=20)
        cats = Categorical.from_codes(codes, levels, ordered=True)
        df = DataFrame(
            np.repeat(
                np.arange(20), 4).reshape(-1, 4), columns=list('abcd'))
        df['cats'] = cats

        # with a cat index
        result = df.set_index('cats').groupby(level=0).sum()
        expected = df[list('abcd')].groupby(cats.codes).sum()
        expected.index = CategoricalIndex(
            Categorical.from_codes(
                [0, 1, 2, 3], levels, ordered=True), name='cats')
        assert_frame_equal(result, expected)

        # with a cat column, should produce a cat index
        result = df.groupby('cats').sum()
        expected = df[list('abcd')].groupby(cats.codes).sum()
        expected.index = CategoricalIndex(
            Categorical.from_codes(
                [0, 1, 2, 3], levels, ordered=True), name='cats')
        assert_frame_equal(result, expected)

    def test_groupby_describe_categorical_columns(self):
        # GH 11558
        cats = pd.CategoricalIndex(['qux', 'foo', 'baz', 'bar'],
                                   categories=['foo', 'bar', 'baz', 'qux'],
                                   ordered=True)
        df = DataFrame(np.random.randn(20, 4), columns=cats)
        result = df.groupby([1, 2, 3, 4] * 5).describe()

        tm.assert_index_equal(result.columns, cats)
        tm.assert_categorical_equal(result.columns.values, cats.values)

    def test_groupby_unstack_categorical(self):
        # GH11558 (example is taken from the original issue)
        df = pd.DataFrame({'a': range(10),
                           'medium': ['A', 'B'] * 5,
                           'artist': list('XYXXY') * 2})
        df['medium'] = df['medium'].astype('category')

        gcat = df.groupby(['artist', 'medium'])['a'].count().unstack()
        result = gcat.describe()

        exp_columns = pd.CategoricalIndex(['A', 'B'], ordered=False,
                                          name='medium')
        tm.assert_index_equal(result.columns, exp_columns)
        tm.assert_categorical_equal(result.columns.values, exp_columns.values)

        result = gcat['A'] + gcat['B']
        expected = pd.Series([6, 4], index=pd.Index(['X', 'Y'], name='artist'))
        tm.assert_series_equal(result, expected)

    def test_groupby_groups_datetimeindex(self):
        # #1430
        from pandas.tseries.api import DatetimeIndex
        periods = 1000
        ind = DatetimeIndex(start='2012/1/1', freq='5min', periods=periods)
        df = DataFrame({'high': np.arange(periods),
                        'low': np.arange(periods)}, index=ind)
        grouped = df.groupby(lambda x: datetime(x.year, x.month, x.day))

        # it works!
        groups = grouped.groups
        tm.assertIsInstance(list(groups.keys())[0], datetime)

    def test_groupby_groups_datetimeindex_tz(self):
        # GH 3950
        dates = ['2011-07-19 07:00:00', '2011-07-19 08:00:00',
                 '2011-07-19 09:00:00', '2011-07-19 07:00:00',
                 '2011-07-19 08:00:00', '2011-07-19 09:00:00']
        df = DataFrame({'label': ['a', 'a', 'a', 'b', 'b', 'b'],
                        'datetime': dates,
                        'value1': np.arange(6, dtype='int64'),
                        'value2': [1, 2] * 3})
        df['datetime'] = df['datetime'].apply(
            lambda d: Timestamp(d, tz='US/Pacific'))

        exp_idx1 = pd.DatetimeIndex(['2011-07-19 07:00:00',
                                     '2011-07-19 07:00:00',
                                     '2011-07-19 08:00:00',
                                     '2011-07-19 08:00:00',
                                     '2011-07-19 09:00:00',
                                     '2011-07-19 09:00:00'],
                                    tz='US/Pacific', name='datetime')
        exp_idx2 = Index(['a', 'b'] * 3, name='label')
        exp_idx = MultiIndex.from_arrays([exp_idx1, exp_idx2])
        expected = DataFrame({'value1': [0, 3, 1, 4, 2, 5],
                              'value2': [1, 2, 2, 1, 1, 2]},
                             index=exp_idx, columns=['value1', 'value2'])

        result = df.groupby(['datetime', 'label']).sum()
        assert_frame_equal(result, expected)

        # by level
        didx = pd.DatetimeIndex(dates, tz='Asia/Tokyo')
        df = DataFrame({'value1': np.arange(6, dtype='int64'),
                        'value2': [1, 2, 3, 1, 2, 3]},
                       index=didx)

        exp_idx = pd.DatetimeIndex(['2011-07-19 07:00:00',
                                    '2011-07-19 08:00:00',
                                    '2011-07-19 09:00:00'], tz='Asia/Tokyo')
        expected = DataFrame({'value1': [3, 5, 7], 'value2': [2, 4, 6]},
                             index=exp_idx, columns=['value1', 'value2'])

        result = df.groupby(level=0).sum()
        assert_frame_equal(result, expected)

    def test_groupby_multi_timezone(self):

        # combining multiple / different timezones yields UTC

        data = """0,2000-01-28 16:47:00,America/Chicago
1,2000-01-29 16:48:00,America/Chicago
2,2000-01-30 16:49:00,America/Los_Angeles
3,2000-01-31 16:50:00,America/Chicago
4,2000-01-01 16:50:00,America/New_York"""

        df = pd.read_csv(StringIO(data), header=None,
                         names=['value', 'date', 'tz'])
        result = df.groupby('tz').date.apply(
            lambda x: pd.to_datetime(x).dt.tz_localize(x.name))

        expected = Series([Timestamp('2000-01-28 16:47:00-0600',
                                     tz='America/Chicago'),
                           Timestamp('2000-01-29 16:48:00-0600',
                                     tz='America/Chicago'),
                           Timestamp('2000-01-30 16:49:00-0800',
                                     tz='America/Los_Angeles'),
                           Timestamp('2000-01-31 16:50:00-0600',
                                     tz='America/Chicago'),
                           Timestamp('2000-01-01 16:50:00-0500',
                                     tz='America/New_York')],
                          name='date',
                          dtype=object)
        assert_series_equal(result, expected)

        tz = 'America/Chicago'
        res_values = df.groupby('tz').date.get_group(tz)
        result = pd.to_datetime(res_values).dt.tz_localize(tz)
        exp_values = Series(['2000-01-28 16:47:00', '2000-01-29 16:48:00',
                             '2000-01-31 16:50:00'],
                            index=[0, 1, 3], name='date')
        expected = pd.to_datetime(exp_values).dt.tz_localize(tz)
        assert_series_equal(result, expected)

    def test_groupby_groups_periods(self):
        dates = ['2011-07-19 07:00:00', '2011-07-19 08:00:00',
                 '2011-07-19 09:00:00', '2011-07-19 07:00:00',
                 '2011-07-19 08:00:00', '2011-07-19 09:00:00']
        df = DataFrame({'label': ['a', 'a', 'a', 'b', 'b', 'b'],
                        'period': [pd.Period(d, freq='H') for d in dates],
                        'value1': np.arange(6, dtype='int64'),
                        'value2': [1, 2] * 3})

        exp_idx1 = pd.PeriodIndex(['2011-07-19 07:00:00',
                                   '2011-07-19 07:00:00',
                                   '2011-07-19 08:00:00',
                                   '2011-07-19 08:00:00',
                                   '2011-07-19 09:00:00',
                                   '2011-07-19 09:00:00'],
                                  freq='H', name='period')
        exp_idx2 = Index(['a', 'b'] * 3, name='label')
        exp_idx = MultiIndex.from_arrays([exp_idx1, exp_idx2])
        expected = DataFrame({'value1': [0, 3, 1, 4, 2, 5],
                              'value2': [1, 2, 2, 1, 1, 2]},
                             index=exp_idx, columns=['value1', 'value2'])

        result = df.groupby(['period', 'label']).sum()
        assert_frame_equal(result, expected)

        # by level
        didx = pd.PeriodIndex(dates, freq='H')
        df = DataFrame({'value1': np.arange(6, dtype='int64'),
                        'value2': [1, 2, 3, 1, 2, 3]},
                       index=didx)

        exp_idx = pd.PeriodIndex(['2011-07-19 07:00:00',
                                  '2011-07-19 08:00:00',
                                  '2011-07-19 09:00:00'], freq='H')
        expected = DataFrame({'value1': [3, 5, 7], 'value2': [2, 4, 6]},
                             index=exp_idx, columns=['value1', 'value2'])

        result = df.groupby(level=0).sum()
        assert_frame_equal(result, expected)

    def test_groupby_reindex_inside_function(self):
        from pandas.tseries.api import DatetimeIndex

        periods = 1000
        ind = DatetimeIndex(start='2012/1/1', freq='5min', periods=periods)
        df = DataFrame({'high': np.arange(
            periods), 'low': np.arange(periods)}, index=ind)

        def agg_before(hour, func, fix=False):
            """
                Run an aggregate func on the subset of data.
            """

            def _func(data):
                d = data.select(lambda x: x.hour < 11).dropna()
                if fix:
                    data[data.index[0]]
                if len(d) == 0:
                    return None
                return func(d)

            return _func

        def afunc(data):
            d = data.select(lambda x: x.hour < 11).dropna()
            return np.max(d)

        grouped = df.groupby(lambda x: datetime(x.year, x.month, x.day))
        closure_bad = grouped.agg({'high': agg_before(11, np.max)})
        closure_good = grouped.agg({'high': agg_before(11, np.max, True)})

        assert_frame_equal(closure_bad, closure_good)

    def test_multiindex_columns_empty_level(self):
        l = [['count', 'values'], ['to filter', '']]
        midx = MultiIndex.from_tuples(l)

        df = DataFrame([[long(1), 'A']], columns=midx)

        grouped = df.groupby('to filter').groups
        self.assertEqual(grouped['A'], [0])

        grouped = df.groupby([('to filter', '')]).groups
        self.assertEqual(grouped['A'], [0])

        df = DataFrame([[long(1), 'A'], [long(2), 'B']], columns=midx)

        expected = df.groupby('to filter').groups
        result = df.groupby([('to filter', '')]).groups
        self.assertEqual(result, expected)

        df = DataFrame([[long(1), 'A'], [long(2), 'A']], columns=midx)

        expected = df.groupby('to filter').groups
        result = df.groupby([('to filter', '')]).groups
        tm.assert_dict_equal(result, expected)

    def test_cython_median(self):
        df = DataFrame(np.random.randn(1000))
        df.values[::2] = np.nan

        labels = np.random.randint(0, 50, size=1000).astype(float)
        labels[::17] = np.nan

        result = df.groupby(labels).median()
        exp = df.groupby(labels).agg(nanops.nanmedian)
        assert_frame_equal(result, exp)

        df = DataFrame(np.random.randn(1000, 5))
        rs = df.groupby(labels).agg(np.median)
        xp = df.groupby(labels).median()
        assert_frame_equal(rs, xp)

    def test_median_empty_bins(self):
        df = pd.DataFrame(np.random.randint(0, 44, 500))

        grps = range(0, 55, 5)
        bins = pd.cut(df[0], grps)

        result = df.groupby(bins).median()
        expected = df.groupby(bins).agg(lambda x: x.median())
        assert_frame_equal(result, expected)

    def test_groupby_categorical_no_compress(self):
        data = Series(np.random.randn(9))

        codes = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        cats = Categorical.from_codes(codes, [0, 1, 2], ordered=True)

        result = data.groupby(cats).mean()
        exp = data.groupby(codes).mean()

        exp.index = CategoricalIndex(exp.index, categories=cats.categories,
                                     ordered=cats.ordered)
        assert_series_equal(result, exp)

        codes = np.array([0, 0, 0, 1, 1, 1, 3, 3, 3])
        cats = Categorical.from_codes(codes, [0, 1, 2, 3], ordered=True)

        result = data.groupby(cats).mean()
        exp = data.groupby(codes).mean().reindex(cats.categories)
        exp.index = CategoricalIndex(exp.index, categories=cats.categories,
                                     ordered=cats.ordered)
        assert_series_equal(result, exp)

        cats = Categorical(["a", "a", "a", "b", "b", "b", "c", "c", "c"],
                           categories=["a", "b", "c", "d"], ordered=True)
        data = DataFrame({"a": [1, 1, 1, 2, 2, 2, 3, 4, 5], "b": cats})

        result = data.groupby("b").mean()
        result = result["a"].values
        exp = np.array([1, 2, 4, np.nan])
        self.assert_numpy_array_equal(result, exp)

    def test_groupby_non_arithmetic_agg_types(self):
        # GH9311, GH6620
        df = pd.DataFrame([{'a': 1,
                            'b': 1}, {'a': 1,
                                      'b': 2}, {'a': 2,
                                                'b': 3}, {'a': 2,
                                                          'b': 4}])

        dtypes = ['int8', 'int16', 'int32', 'int64', 'float32', 'float64']

        grp_exp = {'first': {'df': [{'a': 1,
                                     'b': 1}, {'a': 2,
                                               'b': 3}]},
                   'last': {'df': [{'a': 1,
                                    'b': 2}, {'a': 2,
                                              'b': 4}]},
                   'min': {'df': [{'a': 1,
                                   'b': 1}, {'a': 2,
                                             'b': 3}]},
                   'max': {'df': [{'a': 1,
                                   'b': 2}, {'a': 2,
                                             'b': 4}]},
                   'nth': {'df': [{'a': 1,
                                   'b': 2}, {'a': 2,
                                             'b': 4}],
                           'args': [1]},
                   'count': {'df': [{'a': 1,
                                     'b': 2}, {'a': 2,
                                               'b': 2}],
                             'out_type': 'int64'}}

        for dtype in dtypes:
            df_in = df.copy()
            df_in['b'] = df_in.b.astype(dtype)

            for method, data in compat.iteritems(grp_exp):
                if 'args' not in data:
                    data['args'] = []

                if 'out_type' in data:
                    out_type = data['out_type']
                else:
                    out_type = dtype

                exp = data['df']
                df_out = pd.DataFrame(exp)

                df_out['b'] = df_out.b.astype(out_type)
                df_out.set_index('a', inplace=True)

                grpd = df_in.groupby('a')
                t = getattr(grpd, method)(*data['args'])
                assert_frame_equal(t, df_out)

    def test_groupby_non_arithmetic_agg_intlike_precision(self):
        # GH9311, GH6620
        c = 24650000000000000

        inputs = ((Timestamp('2011-01-15 12:50:28.502376'),
                   Timestamp('2011-01-20 12:50:28.593448')), (1 + c, 2 + c))

        for i in inputs:
            df = pd.DataFrame([{'a': 1, 'b': i[0]}, {'a': 1, 'b': i[1]}])

            grp_exp = {'first': {'expected': i[0]},
                       'last': {'expected': i[1]},
                       'min': {'expected': i[0]},
                       'max': {'expected': i[1]},
                       'nth': {'expected': i[1],
                               'args': [1]},
                       'count': {'expected': 2}}

            for method, data in compat.iteritems(grp_exp):
                if 'args' not in data:
                    data['args'] = []

                grpd = df.groupby('a')
                res = getattr(grpd, method)(*data['args'])
                self.assertEqual(res.iloc[0].b, data['expected'])

    def test_groupby_first_datetime64(self):
        df = DataFrame([(1, 1351036800000000000), (2, 1351036800000000000)])
        df[1] = df[1].view('M8[ns]')

        self.assertTrue(issubclass(df[1].dtype.type, np.datetime64))

        result = df.groupby(level=0).first()
        got_dt = result[1].dtype
        self.assertTrue(issubclass(got_dt.type, np.datetime64))

        result = df[1].groupby(level=0).first()
        got_dt = result.dtype
        self.assertTrue(issubclass(got_dt.type, np.datetime64))

    def test_groupby_max_datetime64(self):
        # GH 5869
        # datetimelike dtype conversion from int
        df = DataFrame(dict(A=Timestamp('20130101'), B=np.arange(5)))
        expected = df.groupby('A')['A'].apply(lambda x: x.max())
        result = df.groupby('A')['A'].max()
        assert_series_equal(result, expected)

    def test_groupby_datetime64_32_bit(self):
        # GH 6410 / numpy 4328
        # 32-bit under 1.9-dev indexing issue

        df = DataFrame({"A": range(2), "B": [pd.Timestamp('2000-01-1')] * 2})
        result = df.groupby("A")["B"].transform(min)
        expected = Series([pd.Timestamp('2000-01-1')] * 2, name='B')
        assert_series_equal(result, expected)

    def test_groupby_categorical_unequal_len(self):
        # GH3011
        series = Series([np.nan, np.nan, 1, 1, 2, 2, 3, 3, 4, 4])
        # The raises only happens with categorical, not with series of types
        # category
        bins = pd.cut(series.dropna().values, 4)

        # len(bins) != len(series) here
        self.assertRaises(ValueError, lambda: series.groupby(bins).mean())

    def test_groupby_multiindex_missing_pair(self):
        # GH9049
        df = DataFrame({'group1': ['a', 'a', 'a', 'b'],
                        'group2': ['c', 'c', 'd', 'c'],
                        'value': [1, 1, 1, 5]})
        df = df.set_index(['group1', 'group2'])
        df_grouped = df.groupby(level=['group1', 'group2'], sort=True)

        res = df_grouped.agg('sum')
        idx = MultiIndex.from_tuples(
            [('a', 'c'), ('a', 'd'), ('b', 'c')], names=['group1', 'group2'])
        exp = DataFrame([[2], [1], [5]], index=idx, columns=['value'])

        tm.assert_frame_equal(res, exp)

    def test_groupby_multiindex_not_lexsorted(self):
        # GH 11640

        # define the lexsorted version
        lexsorted_mi = MultiIndex.from_tuples(
            [('a', ''), ('b1', 'c1'), ('b2', 'c2')], names=['b', 'c'])
        lexsorted_df = DataFrame([[1, 3, 4]], columns=lexsorted_mi)
        self.assertTrue(lexsorted_df.columns.is_lexsorted())

        # define the non-lexsorted version
        not_lexsorted_df = DataFrame(columns=['a', 'b', 'c', 'd'],
                                     data=[[1, 'b1', 'c1', 3],
                                           [1, 'b2', 'c2', 4]])
        not_lexsorted_df = not_lexsorted_df.pivot_table(
            index='a', columns=['b', 'c'], values='d')
        not_lexsorted_df = not_lexsorted_df.reset_index()
        self.assertFalse(not_lexsorted_df.columns.is_lexsorted())

        # compare the results
        tm.assert_frame_equal(lexsorted_df, not_lexsorted_df)

        expected = lexsorted_df.groupby('a').mean()
        with tm.assert_produces_warning(com.PerformanceWarning):
            result = not_lexsorted_df.groupby('a').mean()
        tm.assert_frame_equal(expected, result)

        # a transforming function should work regardless of sort
        # GH 14776
        df = DataFrame({'x': ['a', 'a', 'b', 'a'],
                        'y': [1, 1, 2, 2],
                        'z': [1, 2, 3, 4]}).set_index(['x', 'y'])
        self.assertFalse(df.index.is_lexsorted())

        for level in [0, 1, [0, 1]]:
            for sort in [False, True]:
                result = df.groupby(level=level, sort=sort).apply(
                    DataFrame.drop_duplicates)
                expected = df
                tm.assert_frame_equal(expected, result)

                result = df.sort_index().groupby(level=level, sort=sort).apply(
                    DataFrame.drop_duplicates)
                expected = df.sort_index()
                tm.assert_frame_equal(expected, result)

    def test_groupby_levels_and_columns(self):
        # GH9344, GH9049
        idx_names = ['x', 'y']
        idx = pd.MultiIndex.from_tuples(
            [(1, 1), (1, 2), (3, 4), (5, 6)], names=idx_names)
        df = pd.DataFrame(np.arange(12).reshape(-1, 3), index=idx)

        by_levels = df.groupby(level=idx_names).mean()
        # reset_index changes columns dtype to object
        by_columns = df.reset_index().groupby(idx_names).mean()

        tm.assert_frame_equal(by_levels, by_columns, check_column_type=False)

        by_columns.columns = pd.Index(by_columns.columns, dtype=np.int64)
        tm.assert_frame_equal(by_levels, by_columns)

    def test_gb_apply_list_of_unequal_len_arrays(self):

        # GH1738
        df = DataFrame({'group1': ['a', 'a', 'a', 'b', 'b', 'b', 'a', 'a', 'a',
                                   'b', 'b', 'b'],
                        'group2': ['c', 'c', 'd', 'd', 'd', 'e', 'c', 'c', 'd',
                                   'd', 'd', 'e'],
                        'weight': [1.1, 2, 3, 4, 5, 6, 2, 4, 6, 8, 1, 2],
                        'value': [7.1, 8, 9, 10, 11, 12, 8, 7, 6, 5, 4, 3]})
        df = df.set_index(['group1', 'group2'])
        df_grouped = df.groupby(level=['group1', 'group2'], sort=True)

        def noddy(value, weight):
            out = np.array(value * weight).repeat(3)
            return out

        # the kernel function returns arrays of unequal length
        # pandas sniffs the first one, sees it's an array and not
        # a list, and assumed the rest are of equal length
        # and so tries a vstack

        # don't die
        df_grouped.apply(lambda x: noddy(x.value, x.weight))

    def test_groupby_with_empty(self):
        index = pd.DatetimeIndex(())
        data = ()
        series = pd.Series(data, index)
        grouper = pd.tseries.resample.TimeGrouper('D')
        grouped = series.groupby(grouper)
        assert next(iter(grouped), None) is None

    def test_groupby_with_single_column(self):
        df = pd.DataFrame({'a': list('abssbab')})
        tm.assert_frame_equal(df.groupby('a').get_group('a'), df.iloc[[0, 5]])
        # GH 13530
        exp = pd.DataFrame([], index=pd.Index(['a', 'b', 's'], name='a'))
        tm.assert_frame_equal(df.groupby('a').count(), exp)
        tm.assert_frame_equal(df.groupby('a').sum(), exp)
        tm.assert_frame_equal(df.groupby('a').nth(1), exp)

    def test_groupby_with_small_elem(self):
        # GH 8542
        # length=2
        df = pd.DataFrame({'event': ['start', 'start'],
                           'change': [1234, 5678]},
                          index=pd.DatetimeIndex(['2014-09-10', '2013-10-10']))
        grouped = df.groupby([pd.TimeGrouper(freq='M'), 'event'])
        self.assertEqual(len(grouped.groups), 2)
        self.assertEqual(grouped.ngroups, 2)
        self.assertIn((pd.Timestamp('2014-09-30'), 'start'), grouped.groups)
        self.assertIn((pd.Timestamp('2013-10-31'), 'start'), grouped.groups)

        res = grouped.get_group((pd.Timestamp('2014-09-30'), 'start'))
        tm.assert_frame_equal(res, df.iloc[[0], :])
        res = grouped.get_group((pd.Timestamp('2013-10-31'), 'start'))
        tm.assert_frame_equal(res, df.iloc[[1], :])

        df = pd.DataFrame({'event': ['start', 'start', 'start'],
                           'change': [1234, 5678, 9123]},
                          index=pd.DatetimeIndex(['2014-09-10', '2013-10-10',
                                                  '2014-09-15']))
        grouped = df.groupby([pd.TimeGrouper(freq='M'), 'event'])
        self.assertEqual(len(grouped.groups), 2)
        self.assertEqual(grouped.ngroups, 2)
        self.assertIn((pd.Timestamp('2014-09-30'), 'start'), grouped.groups)
        self.assertIn((pd.Timestamp('2013-10-31'), 'start'), grouped.groups)

        res = grouped.get_group((pd.Timestamp('2014-09-30'), 'start'))
        tm.assert_frame_equal(res, df.iloc[[0, 2], :])
        res = grouped.get_group((pd.Timestamp('2013-10-31'), 'start'))
        tm.assert_frame_equal(res, df.iloc[[1], :])

        # length=3
        df = pd.DataFrame({'event': ['start', 'start', 'start'],
                           'change': [1234, 5678, 9123]},
                          index=pd.DatetimeIndex(['2014-09-10', '2013-10-10',
                                                  '2014-08-05']))
        grouped = df.groupby([pd.TimeGrouper(freq='M'), 'event'])
        self.assertEqual(len(grouped.groups), 3)
        self.assertEqual(grouped.ngroups, 3)
        self.assertIn((pd.Timestamp('2014-09-30'), 'start'), grouped.groups)
        self.assertIn((pd.Timestamp('2013-10-31'), 'start'), grouped.groups)
        self.assertIn((pd.Timestamp('2014-08-31'), 'start'), grouped.groups)

        res = grouped.get_group((pd.Timestamp('2014-09-30'), 'start'))
        tm.assert_frame_equal(res, df.iloc[[0], :])
        res = grouped.get_group((pd.Timestamp('2013-10-31'), 'start'))
        tm.assert_frame_equal(res, df.iloc[[1], :])
        res = grouped.get_group((pd.Timestamp('2014-08-31'), 'start'))
        tm.assert_frame_equal(res, df.iloc[[2], :])

    def test_groupby_with_timezone_selection(self):
        # GH 11616
        # Test that column selection returns output in correct timezone.
        np.random.seed(42)
        df = pd.DataFrame({
            'factor': np.random.randint(0, 3, size=60),
            'time': pd.date_range('01/01/2000 00:00', periods=60,
                                  freq='s', tz='UTC')
        })
        df1 = df.groupby('factor').max()['time']
        df2 = df.groupby('factor')['time'].max()
        tm.assert_series_equal(df1, df2)

    def test_timezone_info(self):
        # GH 11682
        # Timezone info lost when broadcasting scalar datetime to DataFrame
        tm._skip_if_no_pytz()
        import pytz

        df = pd.DataFrame({'a': [1], 'b': [datetime.now(pytz.utc)]})
        self.assertEqual(df['b'][0].tzinfo, pytz.utc)
        df = pd.DataFrame({'a': [1, 2, 3]})
        df['b'] = datetime.now(pytz.utc)
        self.assertEqual(df['b'][0].tzinfo, pytz.utc)

    def test_groupby_with_timegrouper(self):
        # GH 4161
        # TimeGrouper requires a sorted index
        # also verifies that the resultant index has the correct name
        import datetime as DT
        df_original = DataFrame({
            'Buyer': 'Carl Carl Carl Carl Joe Carl'.split(),
            'Quantity': [18, 3, 5, 1, 9, 3],
            'Date': [
                DT.datetime(2013, 9, 1, 13, 0),
                DT.datetime(2013, 9, 1, 13, 5),
                DT.datetime(2013, 10, 1, 20, 0),
                DT.datetime(2013, 10, 3, 10, 0),
                DT.datetime(2013, 12, 2, 12, 0),
                DT.datetime(2013, 9, 2, 14, 0),
            ]
        })

        # GH 6908 change target column's order
        df_reordered = df_original.sort_values(by='Quantity')

        for df in [df_original, df_reordered]:
            df = df.set_index(['Date'])

            expected = DataFrame(
                {'Quantity': np.nan},
                index=date_range('20130901 13:00:00',
                                 '20131205 13:00:00', freq='5D',
                                 name='Date', closed='left'))
            expected.iloc[[0, 6, 18], 0] = np.array(
                [24., 6., 9.], dtype='float64')

            result1 = df.resample('5D') .sum()
            assert_frame_equal(result1, expected)

            df_sorted = df.sort_index()
            result2 = df_sorted.groupby(pd.TimeGrouper(freq='5D')).sum()
            assert_frame_equal(result2, expected)

            result3 = df.groupby(pd.TimeGrouper(freq='5D')).sum()
            assert_frame_equal(result3, expected)

    def test_groupby_with_timegrouper_methods(self):
        # GH 3881
        # make sure API of timegrouper conforms

        import datetime as DT
        df_original = pd.DataFrame({
            'Branch': 'A A A A A B'.split(),
            'Buyer': 'Carl Mark Carl Joe Joe Carl'.split(),
            'Quantity': [1, 3, 5, 8, 9, 3],
            'Date': [
                DT.datetime(2013, 1, 1, 13, 0),
                DT.datetime(2013, 1, 1, 13, 5),
                DT.datetime(2013, 10, 1, 20, 0),
                DT.datetime(2013, 10, 2, 10, 0),
                DT.datetime(2013, 12, 2, 12, 0),
                DT.datetime(2013, 12, 2, 14, 0),
            ]
        })

        df_sorted = df_original.sort_values(by='Quantity', ascending=False)

        for df in [df_original, df_sorted]:
            df = df.set_index('Date', drop=False)
            g = df.groupby(pd.TimeGrouper('6M'))
            self.assertTrue(g.group_keys)
            self.assertTrue(isinstance(g.grouper, pd.core.groupby.BinGrouper))
            groups = g.groups
            self.assertTrue(isinstance(groups, dict))
            self.assertTrue(len(groups) == 3)

    def test_timegrouper_with_reg_groups(self):

        # GH 3794
        # allow combinateion of timegrouper/reg groups

        import datetime as DT

        df_original = DataFrame({
            'Branch': 'A A A A A A A B'.split(),
            'Buyer': 'Carl Mark Carl Carl Joe Joe Joe Carl'.split(),
            'Quantity': [1, 3, 5, 1, 8, 1, 9, 3],
            'Date': [
                DT.datetime(2013, 1, 1, 13, 0),
                DT.datetime(2013, 1, 1, 13, 5),
                DT.datetime(2013, 10, 1, 20, 0),
                DT.datetime(2013, 10, 2, 10, 0),
                DT.datetime(2013, 10, 1, 20, 0),
                DT.datetime(2013, 10, 2, 10, 0),
                DT.datetime(2013, 12, 2, 12, 0),
                DT.datetime(2013, 12, 2, 14, 0),
            ]
        }).set_index('Date')

        df_sorted = df_original.sort_values(by='Quantity', ascending=False)

        for df in [df_original, df_sorted]:
            expected = DataFrame({
                'Buyer': 'Carl Joe Mark'.split(),
                'Quantity': [10, 18, 3],
                'Date': [
                    DT.datetime(2013, 12, 31, 0, 0),
                    DT.datetime(2013, 12, 31, 0, 0),
                    DT.datetime(2013, 12, 31, 0, 0),
                ]
            }).set_index(['Date', 'Buyer'])

            result = df.groupby([pd.Grouper(freq='A'), 'Buyer']).sum()
            assert_frame_equal(result, expected)

            expected = DataFrame({
                'Buyer': 'Carl Mark Carl Joe'.split(),
                'Quantity': [1, 3, 9, 18],
                'Date': [
                    DT.datetime(2013, 1, 1, 0, 0),
                    DT.datetime(2013, 1, 1, 0, 0),
                    DT.datetime(2013, 7, 1, 0, 0),
                    DT.datetime(2013, 7, 1, 0, 0),
                ]
            }).set_index(['Date', 'Buyer'])
            result = df.groupby([pd.Grouper(freq='6MS'), 'Buyer']).sum()
            assert_frame_equal(result, expected)

        df_original = DataFrame({
            'Branch': 'A A A A A A A B'.split(),
            'Buyer': 'Carl Mark Carl Carl Joe Joe Joe Carl'.split(),
            'Quantity': [1, 3, 5, 1, 8, 1, 9, 3],
            'Date': [
                DT.datetime(2013, 10, 1, 13, 0),
                DT.datetime(2013, 10, 1, 13, 5),
                DT.datetime(2013, 10, 1, 20, 0),
                DT.datetime(2013, 10, 2, 10, 0),
                DT.datetime(2013, 10, 1, 20, 0),
                DT.datetime(2013, 10, 2, 10, 0),
                DT.datetime(2013, 10, 2, 12, 0),
                DT.datetime(2013, 10, 2, 14, 0),
            ]
        }).set_index('Date')

        df_sorted = df_original.sort_values(by='Quantity', ascending=False)
        for df in [df_original, df_sorted]:

            expected = DataFrame({
                'Buyer': 'Carl Joe Mark Carl Joe'.split(),
                'Quantity': [6, 8, 3, 4, 10],
                'Date': [
                    DT.datetime(2013, 10, 1, 0, 0),
                    DT.datetime(2013, 10, 1, 0, 0),
                    DT.datetime(2013, 10, 1, 0, 0),
                    DT.datetime(2013, 10, 2, 0, 0),
                    DT.datetime(2013, 10, 2, 0, 0),
                ]
            }).set_index(['Date', 'Buyer'])

            result = df.groupby([pd.Grouper(freq='1D'), 'Buyer']).sum()
            assert_frame_equal(result, expected)

            result = df.groupby([pd.Grouper(freq='1M'), 'Buyer']).sum()
            expected = DataFrame({
                'Buyer': 'Carl Joe Mark'.split(),
                'Quantity': [10, 18, 3],
                'Date': [
                    DT.datetime(2013, 10, 31, 0, 0),
                    DT.datetime(2013, 10, 31, 0, 0),
                    DT.datetime(2013, 10, 31, 0, 0),
                ]
            }).set_index(['Date', 'Buyer'])
            assert_frame_equal(result, expected)

            # passing the name
            df = df.reset_index()
            result = df.groupby([pd.Grouper(freq='1M', key='Date'), 'Buyer'
                                 ]).sum()
            assert_frame_equal(result, expected)

            with self.assertRaises(KeyError):
                df.groupby([pd.Grouper(freq='1M', key='foo'), 'Buyer']).sum()

            # passing the level
            df = df.set_index('Date')
            result = df.groupby([pd.Grouper(freq='1M', level='Date'), 'Buyer'
                                 ]).sum()
            assert_frame_equal(result, expected)
            result = df.groupby([pd.Grouper(freq='1M', level=0), 'Buyer']).sum(
            )
            assert_frame_equal(result, expected)

            with self.assertRaises(ValueError):
                df.groupby([pd.Grouper(freq='1M', level='foo'),
                            'Buyer']).sum()

            # multi names
            df = df.copy()
            df['Date'] = df.index + pd.offsets.MonthEnd(2)
            result = df.groupby([pd.Grouper(freq='1M', key='Date'), 'Buyer'
                                 ]).sum()
            expected = DataFrame({
                'Buyer': 'Carl Joe Mark'.split(),
                'Quantity': [10, 18, 3],
                'Date': [
                    DT.datetime(2013, 11, 30, 0, 0),
                    DT.datetime(2013, 11, 30, 0, 0),
                    DT.datetime(2013, 11, 30, 0, 0),
                ]
            }).set_index(['Date', 'Buyer'])
            assert_frame_equal(result, expected)

            # error as we have both a level and a name!
            with self.assertRaises(ValueError):
                df.groupby([pd.Grouper(freq='1M', key='Date',
                                       level='Date'), 'Buyer']).sum()

            # single groupers
            expected = DataFrame({'Quantity': [31],
                                  'Date': [DT.datetime(2013, 10, 31, 0, 0)
                                           ]}).set_index('Date')
            result = df.groupby(pd.Grouper(freq='1M')).sum()
            assert_frame_equal(result, expected)

            result = df.groupby([pd.Grouper(freq='1M')]).sum()
            assert_frame_equal(result, expected)

            expected = DataFrame({'Quantity': [31],
                                  'Date': [DT.datetime(2013, 11, 30, 0, 0)
                                           ]}).set_index('Date')
            result = df.groupby(pd.Grouper(freq='1M', key='Date')).sum()
            assert_frame_equal(result, expected)

            result = df.groupby([pd.Grouper(freq='1M', key='Date')]).sum()
            assert_frame_equal(result, expected)

        # GH 6764 multiple grouping with/without sort
        df = DataFrame({
            'date': pd.to_datetime([
                '20121002', '20121007', '20130130', '20130202', '20130305',
                '20121002', '20121207', '20130130', '20130202', '20130305',
                '20130202', '20130305'
            ]),
            'user_id': [1, 1, 1, 1, 1, 3, 3, 3, 5, 5, 5, 5],
            'whole_cost': [1790, 364, 280, 259, 201, 623, 90, 312, 359, 301,
                           359, 801],
            'cost1': [12, 15, 10, 24, 39, 1, 0, 90, 45, 34, 1, 12]
        }).set_index('date')

        for freq in ['D', 'M', 'A', 'Q-APR']:
            expected = df.groupby('user_id')[
                'whole_cost'].resample(
                    freq).sum().dropna().reorder_levels(
                        ['date', 'user_id']).sortlevel().astype('int64')
            expected.name = 'whole_cost'

            result1 = df.sort_index().groupby([pd.TimeGrouper(freq=freq),
                                               'user_id'])['whole_cost'].sum()
            assert_series_equal(result1, expected)

            result2 = df.groupby([pd.TimeGrouper(freq=freq), 'user_id'])[
                'whole_cost'].sum()
            assert_series_equal(result2, expected)

    def test_timegrouper_get_group(self):
        # GH 6914

        df_original = DataFrame({
            'Buyer': 'Carl Joe Joe Carl Joe Carl'.split(),
            'Quantity': [18, 3, 5, 1, 9, 3],
            'Date': [datetime(2013, 9, 1, 13, 0),
                     datetime(2013, 9, 1, 13, 5),
                     datetime(2013, 10, 1, 20, 0),
                     datetime(2013, 10, 3, 10, 0),
                     datetime(2013, 12, 2, 12, 0),
                     datetime(2013, 9, 2, 14, 0), ]
        })
        df_reordered = df_original.sort_values(by='Quantity')

        # single grouping
        expected_list = [df_original.iloc[[0, 1, 5]], df_original.iloc[[2, 3]],
                         df_original.iloc[[4]]]
        dt_list = ['2013-09-30', '2013-10-31', '2013-12-31']

        for df in [df_original, df_reordered]:
            grouped = df.groupby(pd.Grouper(freq='M', key='Date'))
            for t, expected in zip(dt_list, expected_list):
                dt = pd.Timestamp(t)
                result = grouped.get_group(dt)
                assert_frame_equal(result, expected)

        # multiple grouping
        expected_list = [df_original.iloc[[1]], df_original.iloc[[3]],
                         df_original.iloc[[4]]]
        g_list = [('Joe', '2013-09-30'), ('Carl', '2013-10-31'),
                  ('Joe', '2013-12-31')]

        for df in [df_original, df_reordered]:
            grouped = df.groupby(['Buyer', pd.Grouper(freq='M', key='Date')])
            for (b, t), expected in zip(g_list, expected_list):
                dt = pd.Timestamp(t)
                result = grouped.get_group((b, dt))
                assert_frame_equal(result, expected)

        # with index
        df_original = df_original.set_index('Date')
        df_reordered = df_original.sort_values(by='Quantity')

        expected_list = [df_original.iloc[[0, 1, 5]], df_original.iloc[[2, 3]],
                         df_original.iloc[[4]]]

        for df in [df_original, df_reordered]:
            grouped = df.groupby(pd.Grouper(freq='M'))
            for t, expected in zip(dt_list, expected_list):
                dt = pd.Timestamp(t)
                result = grouped.get_group(dt)
                assert_frame_equal(result, expected)

    def test_timegrouper_apply_return_type_series(self):
        # Using `apply` with the `TimeGrouper` should give the
        # same return type as an `apply` with a `Grouper`.
        # Issue #11742
        df = pd.DataFrame({'date': ['10/10/2000', '11/10/2000'],
                           'value': [10, 13]})
        df_dt = df.copy()
        df_dt['date'] = pd.to_datetime(df_dt['date'])

        def sumfunc_series(x):
            return pd.Series([x['value'].sum()], ('sum',))

        expected = df.groupby(pd.Grouper(key='date')).apply(sumfunc_series)
        result = (df_dt.groupby(pd.TimeGrouper(freq='M', key='date'))
                  .apply(sumfunc_series))
        assert_frame_equal(result.reset_index(drop=True),
                           expected.reset_index(drop=True))

    def test_timegrouper_apply_return_type_value(self):
        # Using `apply` with the `TimeGrouper` should give the
        # same return type as an `apply` with a `Grouper`.
        # Issue #11742
        df = pd.DataFrame({'date': ['10/10/2000', '11/10/2000'],
                           'value': [10, 13]})
        df_dt = df.copy()
        df_dt['date'] = pd.to_datetime(df_dt['date'])

        def sumfunc_value(x):
            return x.value.sum()

        expected = df.groupby(pd.Grouper(key='date')).apply(sumfunc_value)
        result = (df_dt.groupby(pd.TimeGrouper(freq='M', key='date'))
                  .apply(sumfunc_value))
        assert_series_equal(result.reset_index(drop=True),
                            expected.reset_index(drop=True))

    def test_cumcount(self):
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'])
        g = df.groupby('A')
        sg = g.A

        expected = Series([0, 1, 2, 0, 3])

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())

    def test_cumcount_empty(self):
        ge = DataFrame().groupby(level=0)
        se = Series().groupby(level=0)

        # edge case, as this is usually considered float
        e = Series(dtype='int64')

        assert_series_equal(e, ge.cumcount())
        assert_series_equal(e, se.cumcount())

    def test_cumcount_dupe_index(self):
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'],
                       index=[0] * 5)
        g = df.groupby('A')
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=[0] * 5)

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())

    def test_cumcount_mi(self):
        mi = MultiIndex.from_tuples([[0, 1], [1, 2], [2, 2], [2, 2], [1, 0]])
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'],
                       index=mi)
        g = df.groupby('A')
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=mi)

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())

    def test_cumcount_groupby_not_col(self):
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'],
                       index=[0] * 5)
        g = df.groupby([0, 0, 0, 1, 0])
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=[0] * 5)

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())

    def test_filter_series(self):
        s = pd.Series([1, 3, 20, 5, 22, 24, 7])
        expected_odd = pd.Series([1, 3, 5, 7], index=[0, 1, 3, 6])
        expected_even = pd.Series([20, 22, 24], index=[2, 4, 5])
        grouper = s.apply(lambda x: x % 2)
        grouped = s.groupby(grouper)
        assert_series_equal(
            grouped.filter(lambda x: x.mean() < 10), expected_odd)
        assert_series_equal(
            grouped.filter(lambda x: x.mean() > 10), expected_even)
        # Test dropna=False.
        assert_series_equal(
            grouped.filter(lambda x: x.mean() < 10, dropna=False),
            expected_odd.reindex(s.index))
        assert_series_equal(
            grouped.filter(lambda x: x.mean() > 10, dropna=False),
            expected_even.reindex(s.index))

    def test_filter_single_column_df(self):
        df = pd.DataFrame([1, 3, 20, 5, 22, 24, 7])
        expected_odd = pd.DataFrame([1, 3, 5, 7], index=[0, 1, 3, 6])
        expected_even = pd.DataFrame([20, 22, 24], index=[2, 4, 5])
        grouper = df[0].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        assert_frame_equal(
            grouped.filter(lambda x: x.mean() < 10), expected_odd)
        assert_frame_equal(
            grouped.filter(lambda x: x.mean() > 10), expected_even)
        # Test dropna=False.
        assert_frame_equal(
            grouped.filter(lambda x: x.mean() < 10, dropna=False),
            expected_odd.reindex(df.index))
        assert_frame_equal(
            grouped.filter(lambda x: x.mean() > 10, dropna=False),
            expected_even.reindex(df.index))

    def test_filter_multi_column_df(self):
        df = pd.DataFrame({'A': [1, 12, 12, 1], 'B': [1, 1, 1, 1]})
        grouper = df['A'].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        expected = pd.DataFrame({'A': [12, 12], 'B': [1, 1]}, index=[1, 2])
        assert_frame_equal(
            grouped.filter(lambda x: x['A'].sum() - x['B'].sum() > 10),
            expected)

    def test_filter_mixed_df(self):
        df = pd.DataFrame({'A': [1, 12, 12, 1], 'B': 'a b c d'.split()})
        grouper = df['A'].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        expected = pd.DataFrame({'A': [12, 12], 'B': ['b', 'c']}, index=[1, 2])
        assert_frame_equal(
            grouped.filter(lambda x: x['A'].sum() > 10), expected)

    def test_filter_out_all_groups(self):
        s = pd.Series([1, 3, 20, 5, 22, 24, 7])
        grouper = s.apply(lambda x: x % 2)
        grouped = s.groupby(grouper)
        assert_series_equal(grouped.filter(lambda x: x.mean() > 1000), s[[]])
        df = pd.DataFrame({'A': [1, 12, 12, 1], 'B': 'a b c d'.split()})
        grouper = df['A'].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        assert_frame_equal(
            grouped.filter(lambda x: x['A'].sum() > 1000), df.ix[[]])

    def test_filter_out_no_groups(self):
        s = pd.Series([1, 3, 20, 5, 22, 24, 7])
        grouper = s.apply(lambda x: x % 2)
        grouped = s.groupby(grouper)
        filtered = grouped.filter(lambda x: x.mean() > 0)
        assert_series_equal(filtered, s)
        df = pd.DataFrame({'A': [1, 12, 12, 1], 'B': 'a b c d'.split()})
        grouper = df['A'].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        filtered = grouped.filter(lambda x: x['A'].mean() > 0)
        assert_frame_equal(filtered, df)

    def test_filter_out_all_groups_in_df(self):
        # GH12768
        df = pd.DataFrame({'a': [1, 1, 2], 'b': [1, 2, 0]})
        res = df.groupby('a')
        res = res.filter(lambda x: x['b'].sum() > 5, dropna=False)
        expected = pd.DataFrame({'a': [nan] * 3, 'b': [nan] * 3})
        assert_frame_equal(expected, res)

        df = pd.DataFrame({'a': [1, 1, 2], 'b': [1, 2, 0]})
        res = df.groupby('a')
        res = res.filter(lambda x: x['b'].sum() > 5, dropna=True)
        expected = pd.DataFrame({'a': [], 'b': []}, dtype="int64")
        assert_frame_equal(expected, res)

    def test_filter_condition_raises(self):
        def raise_if_sum_is_zero(x):
            if x.sum() == 0:
                raise ValueError
            else:
                return x.sum() > 0

        s = pd.Series([-1, 0, 1, 2])
        grouper = s.apply(lambda x: x % 2)
        grouped = s.groupby(grouper)
        self.assertRaises(TypeError,
                          lambda: grouped.filter(raise_if_sum_is_zero))

    def test_filter_with_axis_in_groupby(self):
        # issue 11041
        index = pd.MultiIndex.from_product([range(10), [0, 1]])
        data = pd.DataFrame(
            np.arange(100).reshape(-1, 20), columns=index, dtype='int64')
        result = data.groupby(level=0,
                              axis=1).filter(lambda x: x.iloc[0, 0] > 10)
        expected = data.iloc[:, 12:20]
        assert_frame_equal(result, expected)

    def test_filter_bad_shapes(self):
        df = DataFrame({'A': np.arange(8),
                        'B': list('aabbbbcc'),
                        'C': np.arange(8)})
        s = df['B']
        g_df = df.groupby('B')
        g_s = s.groupby(s)

        f = lambda x: x
        self.assertRaises(TypeError, lambda: g_df.filter(f))
        self.assertRaises(TypeError, lambda: g_s.filter(f))

        f = lambda x: x == 1
        self.assertRaises(TypeError, lambda: g_df.filter(f))
        self.assertRaises(TypeError, lambda: g_s.filter(f))

        f = lambda x: np.outer(x, x)
        self.assertRaises(TypeError, lambda: g_df.filter(f))
        self.assertRaises(TypeError, lambda: g_s.filter(f))

    def test_filter_nan_is_false(self):
        df = DataFrame({'A': np.arange(8),
                        'B': list('aabbbbcc'),
                        'C': np.arange(8)})
        s = df['B']
        g_df = df.groupby(df['B'])
        g_s = s.groupby(s)

        f = lambda x: np.nan
        assert_frame_equal(g_df.filter(f), df.loc[[]])
        assert_series_equal(g_s.filter(f), s[[]])

    def test_filter_against_workaround(self):
        np.random.seed(0)
        # Series of ints
        s = Series(np.random.randint(0, 100, 1000))
        grouper = s.apply(lambda x: np.round(x, -1))
        grouped = s.groupby(grouper)
        f = lambda x: x.mean() > 10
        old_way = s[grouped.transform(f).astype('bool')]
        new_way = grouped.filter(f)
        assert_series_equal(new_way.sort_values(), old_way.sort_values())

        # Series of floats
        s = 100 * Series(np.random.random(1000))
        grouper = s.apply(lambda x: np.round(x, -1))
        grouped = s.groupby(grouper)
        f = lambda x: x.mean() > 10
        old_way = s[grouped.transform(f).astype('bool')]
        new_way = grouped.filter(f)
        assert_series_equal(new_way.sort_values(), old_way.sort_values())

        # Set up DataFrame of ints, floats, strings.
        from string import ascii_lowercase
        letters = np.array(list(ascii_lowercase))
        N = 1000
        random_letters = letters.take(np.random.randint(0, 26, N))
        df = DataFrame({'ints': Series(np.random.randint(0, 100, N)),
                        'floats': N / 10 * Series(np.random.random(N)),
                        'letters': Series(random_letters)})

        # Group by ints; filter on floats.
        grouped = df.groupby('ints')
        old_way = df[grouped.floats.
                     transform(lambda x: x.mean() > N / 20).astype('bool')]
        new_way = grouped.filter(lambda x: x['floats'].mean() > N / 20)
        assert_frame_equal(new_way, old_way)

        # Group by floats (rounded); filter on strings.
        grouper = df.floats.apply(lambda x: np.round(x, -1))
        grouped = df.groupby(grouper)
        old_way = df[grouped.letters.
                     transform(lambda x: len(x) < N / 10).astype('bool')]
        new_way = grouped.filter(lambda x: len(x.letters) < N / 10)
        assert_frame_equal(new_way, old_way)

        # Group by strings; filter on ints.
        grouped = df.groupby('letters')
        old_way = df[grouped.ints.
                     transform(lambda x: x.mean() > N / 20).astype('bool')]
        new_way = grouped.filter(lambda x: x['ints'].mean() > N / 20)
        assert_frame_equal(new_way, old_way)

    def test_filter_using_len(self):
        # BUG GH4447
        df = DataFrame({'A': np.arange(8),
                        'B': list('aabbbbcc'),
                        'C': np.arange(8)})
        grouped = df.groupby('B')
        actual = grouped.filter(lambda x: len(x) > 2)
        expected = DataFrame(
            {'A': np.arange(2, 6),
             'B': list('bbbb'),
             'C': np.arange(2, 6)}, index=np.arange(2, 6))
        assert_frame_equal(actual, expected)

        actual = grouped.filter(lambda x: len(x) > 4)
        expected = df.ix[[]]
        assert_frame_equal(actual, expected)

        # Series have always worked properly, but we'll test anyway.
        s = df['B']
        grouped = s.groupby(s)
        actual = grouped.filter(lambda x: len(x) > 2)
        expected = Series(4 * ['b'], index=np.arange(2, 6), name='B')
        assert_series_equal(actual, expected)

        actual = grouped.filter(lambda x: len(x) > 4)
        expected = s[[]]
        assert_series_equal(actual, expected)

    def test_filter_maintains_ordering(self):
        # Simple case: index is sequential. #4621
        df = DataFrame({'pid': [1, 1, 1, 2, 2, 3, 3, 3],
                        'tag': [23, 45, 62, 24, 45, 34, 25, 62]})
        s = df['pid']
        grouped = df.groupby('tag')
        actual = grouped.filter(lambda x: len(x) > 1)
        expected = df.iloc[[1, 2, 4, 7]]
        assert_frame_equal(actual, expected)

        grouped = s.groupby(df['tag'])
        actual = grouped.filter(lambda x: len(x) > 1)
        expected = s.iloc[[1, 2, 4, 7]]
        assert_series_equal(actual, expected)

        # Now index is sequentially decreasing.
        df.index = np.arange(len(df) - 1, -1, -1)
        s = df['pid']
        grouped = df.groupby('tag')
        actual = grouped.filter(lambda x: len(x) > 1)
        expected = df.iloc[[1, 2, 4, 7]]
        assert_frame_equal(actual, expected)

        grouped = s.groupby(df['tag'])
        actual = grouped.filter(lambda x: len(x) > 1)
        expected = s.iloc[[1, 2, 4, 7]]
        assert_series_equal(actual, expected)

        # Index is shuffled.
        SHUFFLED = [4, 6, 7, 2, 1, 0, 5, 3]
        df.index = df.index[SHUFFLED]
        s = df['pid']
        grouped = df.groupby('tag')
        actual = grouped.filter(lambda x: len(x) > 1)
        expected = df.iloc[[1, 2, 4, 7]]
        assert_frame_equal(actual, expected)

        grouped = s.groupby(df['tag'])
        actual = grouped.filter(lambda x: len(x) > 1)
        expected = s.iloc[[1, 2, 4, 7]]
        assert_series_equal(actual, expected)

    def test_filter_multiple_timestamp(self):
        # GH 10114
        df = DataFrame({'A': np.arange(5, dtype='int64'),
                        'B': ['foo', 'bar', 'foo', 'bar', 'bar'],
                        'C': Timestamp('20130101')})

        grouped = df.groupby(['B', 'C'])

        result = grouped['A'].filter(lambda x: True)
        assert_series_equal(df['A'], result)

        result = grouped['A'].transform(len)
        expected = Series([2, 3, 2, 3, 3], name='A')
        assert_series_equal(result, expected)

        result = grouped.filter(lambda x: True)
        assert_frame_equal(df, result)

        result = grouped.transform('sum')
        expected = DataFrame({'A': [2, 8, 2, 8, 8]})
        assert_frame_equal(result, expected)

        result = grouped.transform(len)
        expected = DataFrame({'A': [2, 3, 2, 3, 3]})
        assert_frame_equal(result, expected)

    def test_filter_and_transform_with_non_unique_int_index(self):
        # GH4620
        index = [1, 1, 1, 2, 1, 1, 0, 1]
        df = DataFrame({'pid': [1, 1, 1, 2, 2, 3, 3, 3],
                        'tag': [23, 45, 62, 24, 45, 34, 25, 62]}, index=index)
        grouped_df = df.groupby('tag')
        ser = df['pid']
        grouped_ser = ser.groupby(df['tag'])
        expected_indexes = [1, 2, 4, 7]

        # Filter DataFrame
        actual = grouped_df.filter(lambda x: len(x) > 1)
        expected = df.iloc[expected_indexes]
        assert_frame_equal(actual, expected)

        actual = grouped_df.filter(lambda x: len(x) > 1, dropna=False)
        expected = df.copy()
        expected.iloc[[0, 3, 5, 6]] = np.nan
        assert_frame_equal(actual, expected)

        # Filter Series
        actual = grouped_ser.filter(lambda x: len(x) > 1)
        expected = ser.take(expected_indexes)
        assert_series_equal(actual, expected)

        actual = grouped_ser.filter(lambda x: len(x) > 1, dropna=False)
        NA = np.nan
        expected = Series([NA, 1, 1, NA, 2, NA, NA, 3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index, name='pid')
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_multiple_non_unique_int_index(self):
        # GH4620
        index = [1, 1, 1, 2, 0, 0, 0, 1]
        df = DataFrame({'pid': [1, 1, 1, 2, 2, 3, 3, 3],
                        'tag': [23, 45, 62, 24, 45, 34, 25, 62]}, index=index)
        grouped_df = df.groupby('tag')
        ser = df['pid']
        grouped_ser = ser.groupby(df['tag'])
        expected_indexes = [1, 2, 4, 7]

        # Filter DataFrame
        actual = grouped_df.filter(lambda x: len(x) > 1)
        expected = df.iloc[expected_indexes]
        assert_frame_equal(actual, expected)

        actual = grouped_df.filter(lambda x: len(x) > 1, dropna=False)
        expected = df.copy()
        expected.iloc[[0, 3, 5, 6]] = np.nan
        assert_frame_equal(actual, expected)

        # Filter Series
        actual = grouped_ser.filter(lambda x: len(x) > 1)
        expected = ser.take(expected_indexes)
        assert_series_equal(actual, expected)

        actual = grouped_ser.filter(lambda x: len(x) > 1, dropna=False)
        NA = np.nan
        expected = Series([NA, 1, 1, NA, 2, NA, NA, 3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index, name='pid')
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_non_unique_float_index(self):
        # GH4620
        index = np.array([1, 1, 1, 2, 1, 1, 0, 1], dtype=float)
        df = DataFrame({'pid': [1, 1, 1, 2, 2, 3, 3, 3],
                        'tag': [23, 45, 62, 24, 45, 34, 25, 62]}, index=index)
        grouped_df = df.groupby('tag')
        ser = df['pid']
        grouped_ser = ser.groupby(df['tag'])
        expected_indexes = [1, 2, 4, 7]

        # Filter DataFrame
        actual = grouped_df.filter(lambda x: len(x) > 1)
        expected = df.iloc[expected_indexes]
        assert_frame_equal(actual, expected)

        actual = grouped_df.filter(lambda x: len(x) > 1, dropna=False)
        expected = df.copy()
        expected.iloc[[0, 3, 5, 6]] = np.nan
        assert_frame_equal(actual, expected)

        # Filter Series
        actual = grouped_ser.filter(lambda x: len(x) > 1)
        expected = ser.take(expected_indexes)
        assert_series_equal(actual, expected)

        actual = grouped_ser.filter(lambda x: len(x) > 1, dropna=False)
        NA = np.nan
        expected = Series([NA, 1, 1, NA, 2, NA, NA, 3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index, name='pid')
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_non_unique_timestamp_index(self):
        # GH4620
        t0 = Timestamp('2013-09-30 00:05:00')
        t1 = Timestamp('2013-10-30 00:05:00')
        t2 = Timestamp('2013-11-30 00:05:00')
        index = [t1, t1, t1, t2, t1, t1, t0, t1]
        df = DataFrame({'pid': [1, 1, 1, 2, 2, 3, 3, 3],
                        'tag': [23, 45, 62, 24, 45, 34, 25, 62]}, index=index)
        grouped_df = df.groupby('tag')
        ser = df['pid']
        grouped_ser = ser.groupby(df['tag'])
        expected_indexes = [1, 2, 4, 7]

        # Filter DataFrame
        actual = grouped_df.filter(lambda x: len(x) > 1)
        expected = df.iloc[expected_indexes]
        assert_frame_equal(actual, expected)

        actual = grouped_df.filter(lambda x: len(x) > 1, dropna=False)
        expected = df.copy()
        expected.iloc[[0, 3, 5, 6]] = np.nan
        assert_frame_equal(actual, expected)

        # Filter Series
        actual = grouped_ser.filter(lambda x: len(x) > 1)
        expected = ser.take(expected_indexes)
        assert_series_equal(actual, expected)

        actual = grouped_ser.filter(lambda x: len(x) > 1, dropna=False)
        NA = np.nan
        expected = Series([NA, 1, 1, NA, 2, NA, NA, 3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index, name='pid')
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_non_unique_string_index(self):
        # GH4620
        index = list('bbbcbbab')
        df = DataFrame({'pid': [1, 1, 1, 2, 2, 3, 3, 3],
                        'tag': [23, 45, 62, 24, 45, 34, 25, 62]}, index=index)
        grouped_df = df.groupby('tag')
        ser = df['pid']
        grouped_ser = ser.groupby(df['tag'])
        expected_indexes = [1, 2, 4, 7]

        # Filter DataFrame
        actual = grouped_df.filter(lambda x: len(x) > 1)
        expected = df.iloc[expected_indexes]
        assert_frame_equal(actual, expected)

        actual = grouped_df.filter(lambda x: len(x) > 1, dropna=False)
        expected = df.copy()
        expected.iloc[[0, 3, 5, 6]] = np.nan
        assert_frame_equal(actual, expected)

        # Filter Series
        actual = grouped_ser.filter(lambda x: len(x) > 1)
        expected = ser.take(expected_indexes)
        assert_series_equal(actual, expected)

        actual = grouped_ser.filter(lambda x: len(x) > 1, dropna=False)
        NA = np.nan
        expected = Series([NA, 1, 1, NA, 2, NA, NA, 3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index, name='pid')
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_has_access_to_grouped_cols(self):
        df = DataFrame([[1, 2], [1, 3], [5, 6]], columns=['A', 'B'])
        g = df.groupby('A')
        # previously didn't have access to col A #????
        filt = g.filter(lambda x: x['A'].sum() == 2)
        assert_frame_equal(filt, df.iloc[[0, 1]])

    def test_filter_enforces_scalarness(self):
        df = pd.DataFrame([
            ['best', 'a', 'x'],
            ['worst', 'b', 'y'],
            ['best', 'c', 'x'],
            ['best', 'd', 'y'],
            ['worst', 'd', 'y'],
            ['worst', 'd', 'y'],
            ['best', 'd', 'z'],
        ], columns=['a', 'b', 'c'])
        with tm.assertRaisesRegexp(TypeError, 'filter function returned a.*'):
            df.groupby('c').filter(lambda g: g['a'] == 'best')

    def test_filter_non_bool_raises(self):
        df = pd.DataFrame([
            ['best', 'a', 1],
            ['worst', 'b', 1],
            ['best', 'c', 1],
            ['best', 'd', 1],
            ['worst', 'd', 1],
            ['worst', 'd', 1],
            ['best', 'd', 1],
        ], columns=['a', 'b', 'c'])
        with tm.assertRaisesRegexp(TypeError, 'filter function returned a.*'):
            df.groupby('a').filter(lambda g: g.c.mean())

    def test_fill_constistency(self):

        # GH9221
        # pass thru keyword arguments to the generated wrapper
        # are set if the passed kw is None (only)
        df = DataFrame(index=pd.MultiIndex.from_product(
            [['value1', 'value2'], date_range('2014-01-01', '2014-01-06')]),
            columns=Index(
            ['1', '2'], name='id'))
        df['1'] = [np.nan, 1, np.nan, np.nan, 11, np.nan, np.nan, 2, np.nan,
                   np.nan, 22, np.nan]
        df['2'] = [np.nan, 3, np.nan, np.nan, 33, np.nan, np.nan, 4, np.nan,
                   np.nan, 44, np.nan]

        expected = df.groupby(level=0, axis=0).fillna(method='ffill')
        result = df.T.groupby(level=0, axis=1).fillna(method='ffill').T
        assert_frame_equal(result, expected)

    def test_index_label_overlaps_location(self):
        # checking we don't have any label/location confusion in the
        # the wake of GH5375
        df = DataFrame(list('ABCDE'), index=[2, 0, 2, 1, 1])
        g = df.groupby(list('ababb'))
        actual = g.filter(lambda x: len(x) > 2)
        expected = df.iloc[[1, 3, 4]]
        assert_frame_equal(actual, expected)

        ser = df[0]
        g = ser.groupby(list('ababb'))
        actual = g.filter(lambda x: len(x) > 2)
        expected = ser.take([1, 3, 4])
        assert_series_equal(actual, expected)

        # ... and again, with a generic Index of floats
        df.index = df.index.astype(float)
        g = df.groupby(list('ababb'))
        actual = g.filter(lambda x: len(x) > 2)
        expected = df.iloc[[1, 3, 4]]
        assert_frame_equal(actual, expected)

        ser = df[0]
        g = ser.groupby(list('ababb'))
        actual = g.filter(lambda x: len(x) > 2)
        expected = ser.take([1, 3, 4])
        assert_series_equal(actual, expected)

    def test_groupby_selection_with_methods(self):
        # some methods which require DatetimeIndex
        rng = pd.date_range('2014', periods=len(self.df))
        self.df.index = rng

        g = self.df.groupby(['A'])[['C']]
        g_exp = self.df[['C']].groupby(self.df['A'])
        # TODO check groupby with > 1 col ?

        # methods which are called as .foo()
        methods = ['count',
                   'corr',
                   'cummax',
                   'cummin',
                   'cumprod',
                   'describe',
                   'rank',
                   'quantile',
                   'diff',
                   'shift',
                   'all',
                   'any',
                   'idxmin',
                   'idxmax',
                   'ffill',
                   'bfill',
                   'pct_change',
                   'tshift']

        for m in methods:
            res = getattr(g, m)()
            exp = getattr(g_exp, m)()
            assert_frame_equal(res, exp)  # should always be frames!

        # methods which aren't just .foo()
        assert_frame_equal(g.fillna(0), g_exp.fillna(0))
        assert_frame_equal(g.dtypes, g_exp.dtypes)
        assert_frame_equal(g.apply(lambda x: x.sum()),
                           g_exp.apply(lambda x: x.sum()))

        assert_frame_equal(g.resample('D').mean(), g_exp.resample('D').mean())
        assert_frame_equal(g.resample('D').ohlc(),
                           g_exp.resample('D').ohlc())

        assert_frame_equal(g.filter(lambda x: len(x) == 3),
                           g_exp.filter(lambda x: len(x) == 3))

    def test_groupby_whitelist(self):
        from string import ascii_lowercase
        letters = np.array(list(ascii_lowercase))
        N = 10
        random_letters = letters.take(np.random.randint(0, 26, N))
        df = DataFrame({'floats': N / 10 * Series(np.random.random(N)),
                        'letters': Series(random_letters)})
        s = df.floats

        df_whitelist = frozenset([
            'last',
            'first',
            'mean',
            'sum',
            'min',
            'max',
            'head',
            'tail',
            'cumsum',
            'cumprod',
            'cummin',
            'cummax',
            'cumcount',
            'resample',
            'describe',
            'rank',
            'quantile',
            'fillna',
            'mad',
            'any',
            'all',
            'take',
            'idxmax',
            'idxmin',
            'shift',
            'tshift',
            'ffill',
            'bfill',
            'pct_change',
            'skew',
            'plot',
            'boxplot',
            'hist',
            'median',
            'dtypes',
            'corrwith',
            'corr',
            'cov',
            'diff',
        ])
        s_whitelist = frozenset([
            'last',
            'first',
            'mean',
            'sum',
            'min',
            'max',
            'head',
            'tail',
            'cumsum',
            'cumprod',
            'cummin',
            'cummax',
            'cumcount',
            'resample',
            'describe',
            'rank',
            'quantile',
            'fillna',
            'mad',
            'any',
            'all',
            'take',
            'idxmax',
            'idxmin',
            'shift',
            'tshift',
            'ffill',
            'bfill',
            'pct_change',
            'skew',
            'plot',
            'hist',
            'median',
            'dtype',
            'corr',
            'cov',
            'diff',
            'unique',
            # 'nlargest', 'nsmallest',
        ])

        for obj, whitelist in zip((df, s), (df_whitelist, s_whitelist)):
            gb = obj.groupby(df.letters)
            self.assertEqual(whitelist, gb._apply_whitelist)
            for m in whitelist:
                getattr(type(gb), m)

    AGG_FUNCTIONS = ['sum', 'prod', 'min', 'max', 'median', 'mean', 'skew',
                     'mad', 'std', 'var', 'sem']
    AGG_FUNCTIONS_WITH_SKIPNA = ['skew', 'mad']

    def test_groupby_whitelist_deprecations(self):
        from string import ascii_lowercase
        letters = np.array(list(ascii_lowercase))
        N = 10
        random_letters = letters.take(np.random.randint(0, 26, N))
        df = DataFrame({'floats': N / 10 * Series(np.random.random(N)),
                        'letters': Series(random_letters)})

        # 10711 deprecated
        with tm.assert_produces_warning(FutureWarning):
            df.groupby('letters').irow(0)
        with tm.assert_produces_warning(FutureWarning):
            df.groupby('letters').floats.irow(0)

    def test_regression_whitelist_methods(self):

        # GH6944
        # explicity test the whitelest methods
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'], ['one', 'two',
                                                                  'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        raw_frame = DataFrame(np.random.randn(10, 3), index=index,
                              columns=Index(['A', 'B', 'C'], name='exp'))
        raw_frame.ix[1, [1, 2]] = np.nan
        raw_frame.ix[7, [0, 1]] = np.nan

        for op, level, axis, skipna in cart_product(self.AGG_FUNCTIONS,
                                                    lrange(2), lrange(2),
                                                    [True, False]):

            if axis == 0:
                frame = raw_frame
            else:
                frame = raw_frame.T

            if op in self.AGG_FUNCTIONS_WITH_SKIPNA:
                grouped = frame.groupby(level=level, axis=axis)
                result = getattr(grouped, op)(skipna=skipna)
                expected = getattr(frame, op)(level=level, axis=axis,
                                              skipna=skipna)
                assert_frame_equal(result, expected)
            else:
                grouped = frame.groupby(level=level, axis=axis)
                result = getattr(grouped, op)()
                expected = getattr(frame, op)(level=level, axis=axis)
                assert_frame_equal(result, expected)

    def test_groupby_blacklist(self):
        from string import ascii_lowercase
        letters = np.array(list(ascii_lowercase))
        N = 10
        random_letters = letters.take(np.random.randint(0, 26, N))
        df = DataFrame({'floats': N / 10 * Series(np.random.random(N)),
                        'letters': Series(random_letters)})
        s = df.floats

        blacklist = [
            'eval', 'query', 'abs', 'where',
            'mask', 'align', 'groupby', 'clip', 'astype',
            'at', 'combine', 'consolidate', 'convert_objects',
        ]
        to_methods = [method for method in dir(df) if method.startswith('to_')]

        blacklist.extend(to_methods)

        # e.g., to_csv
        defined_but_not_allowed = ("(?:^Cannot.+{0!r}.+{1!r}.+try using the "
                                   "'apply' method$)")

        # e.g., query, eval
        not_defined = "(?:^{1!r} object has no attribute {0!r}$)"
        fmt = defined_but_not_allowed + '|' + not_defined
        for bl in blacklist:
            for obj in (df, s):
                gb = obj.groupby(df.letters)
                msg = fmt.format(bl, type(gb).__name__)
                with tm.assertRaisesRegexp(AttributeError, msg):
                    getattr(gb, bl)

    def test_tab_completion(self):
        grp = self.mframe.groupby(level='second')
        results = set([v for v in dir(grp) if not v.startswith('_')])
        expected = set(
            ['A', 'B', 'C', 'agg', 'aggregate', 'apply', 'boxplot', 'filter',
             'first', 'get_group', 'groups', 'hist', 'indices', 'last', 'max',
             'mean', 'median', 'min', 'name', 'ngroups', 'nth', 'ohlc', 'plot',
             'prod', 'size', 'std', 'sum', 'transform', 'var', 'sem', 'count',
             'head', 'irow', 'describe', 'cummax', 'quantile', 'rank',
             'cumprod', 'tail', 'resample', 'cummin', 'fillna', 'cumsum',
             'cumcount', 'all', 'shift', 'skew', 'bfill', 'ffill', 'take',
             'tshift', 'pct_change', 'any', 'mad', 'corr', 'corrwith', 'cov',
             'dtypes', 'ndim', 'diff', 'idxmax', 'idxmin',
             'ffill', 'bfill', 'pad', 'backfill', 'rolling', 'expanding'])
        self.assertEqual(results, expected)

    def test_lexsort_indexer(self):
        keys = [[nan] * 5 + list(range(100)) + [nan] * 5]
        # orders=True, na_position='last'
        result = _lexsort_indexer(keys, orders=True, na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=True, na_position='first'
        result = _lexsort_indexer(keys, orders=True, na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=False, na_position='last'
        result = _lexsort_indexer(keys, orders=False, na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

        # orders=False, na_position='first'
        result = _lexsort_indexer(keys, orders=False, na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp, dtype=np.intp))

    def test_nargsort(self):
        # np.argsort(items) places NaNs last
        items = [nan] * 5 + list(range(100)) + [nan] * 5
        # np.argsort(items2) may not place NaNs first
        items2 = np.array(items, dtype='O')

        try:
            # GH 2785; due to a regression in NumPy1.6.2
            np.argsort(np.array([[1, 2], [1, 3], [1, 2]], dtype='i'))
            np.argsort(items2, kind='mergesort')
        except TypeError:
            raise nose.SkipTest('requested sort not available for type')

        # mergesort is the most difficult to get right because we want it to be
        # stable.

        # According to numpy/core/tests/test_multiarray, """The number of
        # sorted items must be greater than ~50 to check the actual algorithm
        # because quick and merge sort fall over to insertion sort for small
        # arrays."""

        # mergesort, ascending=True, na_position='last'
        result = _nargsort(items, kind='mergesort', ascending=True,
                           na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='first'
        result = _nargsort(items, kind='mergesort', ascending=True,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='last'
        result = _nargsort(items, kind='mergesort', ascending=False,
                           na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='first'
        result = _nargsort(items, kind='mergesort', ascending=False,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='last'
        result = _nargsort(items2, kind='mergesort', ascending=True,
                           na_position='last')
        exp = list(range(5, 105)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=True, na_position='first'
        result = _nargsort(items2, kind='mergesort', ascending=True,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(5, 105))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='last'
        result = _nargsort(items2, kind='mergesort', ascending=False,
                           na_position='last')
        exp = list(range(104, 4, -1)) + list(range(5)) + list(range(105, 110))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

        # mergesort, ascending=False, na_position='first'
        result = _nargsort(items2, kind='mergesort', ascending=False,
                           na_position='first')
        exp = list(range(5)) + list(range(105, 110)) + list(range(104, 4, -1))
        tm.assert_numpy_array_equal(result, np.array(exp), check_dtype=False)

    def test_datetime_count(self):
        df = DataFrame({'a': [1, 2, 3] * 2,
                        'dates': pd.date_range('now', periods=6, freq='T')})
        result = df.groupby('a').dates.count()
        expected = Series([
            2, 2, 2
        ], index=Index([1, 2, 3], name='a'), name='dates')
        tm.assert_series_equal(result, expected)

    def test_lower_int_prec_count(self):
        df = DataFrame({'a': np.array(
            [0, 1, 2, 100], np.int8),
            'b': np.array(
            [1, 2, 3, 6], np.uint32),
            'c': np.array(
            [4, 5, 6, 8], np.int16),
            'grp': list('ab' * 2)})
        result = df.groupby('grp').count()
        expected = DataFrame({'a': [2, 2],
                              'b': [2, 2],
                              'c': [2, 2]}, index=pd.Index(list('ab'),
                                                           name='grp'))
        tm.assert_frame_equal(result, expected)

    def test_count_uses_size_on_exception(self):
        class RaisingObjectException(Exception):
            pass

        class RaisingObject(object):

            def __init__(self, msg='I will raise inside Cython'):
                super(RaisingObject, self).__init__()
                self.msg = msg

            def __eq__(self, other):
                # gets called in Cython to check that raising calls the method
                raise RaisingObjectException(self.msg)

        df = DataFrame({'a': [RaisingObject() for _ in range(4)],
                        'grp': list('ab' * 2)})
        result = df.groupby('grp').count()
        expected = DataFrame({'a': [2, 2]}, index=pd.Index(
            list('ab'), name='grp'))
        tm.assert_frame_equal(result, expected)

    def test__cython_agg_general(self):
        ops = [('mean', np.mean),
               ('median', np.median),
               ('var', np.var),
               ('add', np.sum),
               ('prod', np.prod),
               ('min', np.min),
               ('max', np.max),
               ('first', lambda x: x.iloc[0]),
               ('last', lambda x: x.iloc[-1]), ]
        df = DataFrame(np.random.randn(1000))
        labels = np.random.randint(0, 50, size=1000).astype(float)

        for op, targop in ops:
            result = df.groupby(labels)._cython_agg_general(op)
            expected = df.groupby(labels).agg(targop)
            try:
                tm.assert_frame_equal(result, expected)
            except BaseException as exc:
                exc.args += ('operation: %s' % op, )
                raise

    def test_cython_agg_empty_buckets(self):
        ops = [('mean', np.mean),
               ('median', lambda x: np.median(x) if len(x) > 0 else np.nan),
               ('var', lambda x: np.var(x, ddof=1)),
               ('add', lambda x: np.sum(x) if len(x) > 0 else np.nan),
               ('prod', np.prod),
               ('min', np.min),
               ('max', np.max), ]

        df = pd.DataFrame([11, 12, 13])
        grps = range(0, 55, 5)

        for op, targop in ops:
            result = df.groupby(pd.cut(df[0], grps))._cython_agg_general(op)
            expected = df.groupby(pd.cut(df[0], grps)).agg(lambda x: targop(x))
            try:
                tm.assert_frame_equal(result, expected)
            except BaseException as exc:
                exc.args += ('operation: %s' % op,)
                raise

    def test_cython_group_transform_algos(self):
        # GH 4095
        dtypes = [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint32,
                  np.uint64, np.float32, np.float64]

        ops = [(pd.algos.group_cumprod_float64, np.cumproduct, [np.float64]),
               (pd.algos.group_cumsum, np.cumsum, dtypes)]

        for pd_op, np_op, dtypes in ops:
            for dtype in dtypes:
                data = np.array([[1], [2], [3], [4]], dtype=dtype)
                ans = np.zeros_like(data)
                accum = np.array([[0]], dtype=dtype)
                labels = np.array([0, 0, 0, 0], dtype=np.int64)
                pd_op(ans, data, labels, accum)
                self.assert_numpy_array_equal(np_op(data), ans[:, 0],
                                              check_dtype=False)

        # with nans
        labels = np.array([0, 0, 0, 0, 0], dtype=np.int64)

        data = np.array([[1], [2], [3], [np.nan], [4]], dtype='float64')
        accum = np.array([[0.0]])
        actual = np.zeros_like(data)
        actual.fill(np.nan)
        pd.algos.group_cumprod_float64(actual, data, labels, accum)
        expected = np.array([1, 2, 6, np.nan, 24], dtype='float64')
        self.assert_numpy_array_equal(actual[:, 0], expected)

        accum = np.array([[0.0]])
        actual = np.zeros_like(data)
        actual.fill(np.nan)
        pd.algos.group_cumsum(actual, data, labels, accum)
        expected = np.array([1, 3, 6, np.nan, 10], dtype='float64')
        self.assert_numpy_array_equal(actual[:, 0], expected)

        # timedelta
        data = np.array([np.timedelta64(1, 'ns')] * 5, dtype='m8[ns]')[:, None]
        accum = np.array([[0]], dtype='int64')
        actual = np.zeros_like(data, dtype='int64')
        pd.algos.group_cumsum(actual, data.view('int64'), labels, accum)
        expected = np.array([np.timedelta64(1, 'ns'), np.timedelta64(
            2, 'ns'), np.timedelta64(3, 'ns'), np.timedelta64(4, 'ns'),
            np.timedelta64(5, 'ns')])
        self.assert_numpy_array_equal(actual[:, 0].view('m8[ns]'), expected)

    def test_cython_transform(self):
        # GH 4095
        ops = [(('cumprod',
                 ()), lambda x: x.cumprod()), (('cumsum', ()),
                                               lambda x: x.cumsum()),
               (('shift', (-1, )),
                lambda x: x.shift(-1)), (('shift',
                                          (1, )), lambda x: x.shift())]

        s = Series(np.random.randn(1000))
        s_missing = s.copy()
        s_missing.iloc[2:10] = np.nan
        labels = np.random.randint(0, 50, size=1000).astype(float)

        # series
        for (op, args), targop in ops:
            for data in [s, s_missing]:
                # print(data.head())
                expected = data.groupby(labels).transform(targop)

                tm.assert_series_equal(expected,
                                       data.groupby(labels).transform(op,
                                                                      *args))
                tm.assert_series_equal(expected, getattr(
                    data.groupby(labels), op)(*args))

        strings = list('qwertyuiopasdfghjklz')
        strings_missing = strings[:]
        strings_missing[5] = np.nan
        df = DataFrame({'float': s,
                        'float_missing': s_missing,
                        'int': [1, 1, 1, 1, 2] * 200,
                        'datetime': pd.date_range('1990-1-1', periods=1000),
                        'timedelta': pd.timedelta_range(1, freq='s',
                                                        periods=1000),
                        'string': strings * 50,
                        'string_missing': strings_missing * 50})
        df['cat'] = df['string'].astype('category')

        df2 = df.copy()
        df2.index = pd.MultiIndex.from_product([range(100), range(10)])

        # DataFrame - Single and MultiIndex,
        # group by values, index level, columns
        for df in [df, df2]:
            for gb_target in [dict(by=labels), dict(level=0), dict(by='string')
                              ]:  # dict(by='string_missing')]:
                # dict(by=['int','string'])]:

                gb = df.groupby(**gb_target)
                # whitelisted methods set the selection before applying
                # bit a of hack to make sure the cythonized shift
                # is equivalent to pre 0.17.1 behavior
                if op == 'shift':
                    gb._set_group_selection()

                for (op, args), targop in ops:
                    if op != 'shift' and 'int' not in gb_target:
                        # numeric apply fastpath promotes dtype so have
                        # to apply seperately and concat
                        i = gb[['int']].apply(targop)
                        f = gb[['float', 'float_missing']].apply(targop)
                        expected = pd.concat([f, i], axis=1)
                    else:
                        expected = gb.apply(targop)

                    expected = expected.sort_index(axis=1)
                    tm.assert_frame_equal(expected,
                                          gb.transform(op, *args).sort_index(
                                              axis=1))
                    tm.assert_frame_equal(expected, getattr(gb, op)(*args))
                    # individual columns
                    for c in df:
                        if c not in ['float', 'int', 'float_missing'
                                     ] and op != 'shift':
                            self.assertRaises(DataError, gb[c].transform, op)
                            self.assertRaises(DataError, getattr(gb[c], op))
                        else:
                            expected = gb[c].apply(targop)
                            expected.name = c
                            tm.assert_series_equal(expected,
                                                   gb[c].transform(op, *args))
                            tm.assert_series_equal(expected,
                                                   getattr(gb[c], op)(*args))

    def test_groupby_cumprod(self):
        # GH 4095
        df = pd.DataFrame({'key': ['b'] * 10, 'value': 2})

        actual = df.groupby('key')['value'].cumprod()
        expected = df.groupby('key')['value'].apply(lambda x: x.cumprod())
        expected.name = 'value'
        tm.assert_series_equal(actual, expected)

        df = pd.DataFrame({'key': ['b'] * 100, 'value': 2})
        actual = df.groupby('key')['value'].cumprod()
        # if overflows, groupby product casts to float
        # while numpy passes back invalid values
        df['value'] = df['value'].astype(float)
        expected = df.groupby('key')['value'].apply(lambda x: x.cumprod())
        expected.name = 'value'
        tm.assert_series_equal(actual, expected)

    def test_ops_general(self):
        ops = [('mean', np.mean),
               ('median', np.median),
               ('std', np.std),
               ('var', np.var),
               ('sum', np.sum),
               ('prod', np.prod),
               ('min', np.min),
               ('max', np.max),
               ('first', lambda x: x.iloc[0]),
               ('last', lambda x: x.iloc[-1]),
               ('count', np.size), ]
        try:
            from scipy.stats import sem
        except ImportError:
            pass
        else:
            ops.append(('sem', sem))
        df = DataFrame(np.random.randn(1000))
        labels = np.random.randint(0, 50, size=1000).astype(float)

        for op, targop in ops:
            result = getattr(df.groupby(labels), op)().astype(float)
            expected = df.groupby(labels).agg(targop)
            try:
                tm.assert_frame_equal(result, expected)
            except BaseException as exc:
                exc.args += ('operation: %s' % op, )
                raise

    def test_max_nan_bug(self):
        raw = """,Date,app,File
2013-04-23,2013-04-23 00:00:00,,log080001.log
2013-05-06,2013-05-06 00:00:00,,log.log
2013-05-07,2013-05-07 00:00:00,OE,xlsx"""

        df = pd.read_csv(StringIO(raw), parse_dates=[0])
        gb = df.groupby('Date')
        r = gb[['File']].max()
        e = gb['File'].max().to_frame()
        tm.assert_frame_equal(r, e)
        self.assertFalse(r['File'].isnull().any())

    def test_nlargest(self):
        a = Series([1, 3, 5, 7, 2, 9, 0, 4, 6, 10])
        b = Series(list('a' * 5 + 'b' * 5))
        gb = a.groupby(b)
        r = gb.nlargest(3)
        e = Series([
            7, 5, 3, 10, 9, 6
        ], index=MultiIndex.from_arrays([list('aaabbb'), [3, 2, 1, 9, 5, 8]]))
        tm.assert_series_equal(r, e)

        a = Series([1, 1, 3, 2, 0, 3, 3, 2, 1, 0])
        gb = a.groupby(b)
        e = Series([
            3, 2, 1, 3, 3, 2
        ], index=MultiIndex.from_arrays([list('aaabbb'), [2, 3, 1, 6, 5, 7]]))
        assert_series_equal(gb.nlargest(3, keep='last'), e)
        with tm.assert_produces_warning(FutureWarning):
            assert_series_equal(gb.nlargest(3, take_last=True), e)

    def test_nsmallest(self):
        a = Series([1, 3, 5, 7, 2, 9, 0, 4, 6, 10])
        b = Series(list('a' * 5 + 'b' * 5))
        gb = a.groupby(b)
        r = gb.nsmallest(3)
        e = Series([
            1, 2, 3, 0, 4, 6
        ], index=MultiIndex.from_arrays([list('aaabbb'), [0, 4, 1, 6, 7, 8]]))
        tm.assert_series_equal(r, e)

        a = Series([1, 1, 3, 2, 0, 3, 3, 2, 1, 0])
        gb = a.groupby(b)
        e = Series([
            0, 1, 1, 0, 1, 2
        ], index=MultiIndex.from_arrays([list('aaabbb'), [4, 1, 0, 9, 8, 7]]))
        assert_series_equal(gb.nsmallest(3, keep='last'), e)
        with tm.assert_produces_warning(FutureWarning):
            assert_series_equal(gb.nsmallest(3, take_last=True), e)

    def test_transform_doesnt_clobber_ints(self):
        # GH 7972
        n = 6
        x = np.arange(n)
        df = DataFrame({'a': x // 2, 'b': 2.0 * x, 'c': 3.0 * x})
        df2 = DataFrame({'a': x // 2 * 1.0, 'b': 2.0 * x, 'c': 3.0 * x})

        gb = df.groupby('a')
        result = gb.transform('mean')

        gb2 = df2.groupby('a')
        expected = gb2.transform('mean')
        tm.assert_frame_equal(result, expected)

    def test_groupby_categorical_two_columns(self):

        # https://github.com/pandas-dev/pandas/issues/8138
        d = {'cat':
             pd.Categorical(["a", "b", "a", "b"], categories=["a", "b", "c"],
                            ordered=True),
             'ints': [1, 1, 2, 2],
             'val': [10, 20, 30, 40]}
        test = pd.DataFrame(d)

        # Grouping on a single column
        groups_single_key = test.groupby("cat")
        res = groups_single_key.agg('mean')

        exp_index = pd.CategoricalIndex(["a", "b", "c"], name="cat",
                                        ordered=True)
        exp = DataFrame({"ints": [1.5, 1.5, np.nan], "val": [20, 30, np.nan]},
                        index=exp_index)
        tm.assert_frame_equal(res, exp)

        # Grouping on two columns
        groups_double_key = test.groupby(["cat", "ints"])
        res = groups_double_key.agg('mean')
        exp = DataFrame({"val": [10, 30, 20, 40, np.nan, np.nan],
                         "cat": pd.Categorical(["a", "a", "b", "b", "c", "c"],
                                               ordered=True),
                         "ints": [1, 2, 1, 2, 1, 2]}).set_index(["cat", "ints"
                                                                 ])
        tm.assert_frame_equal(res, exp)

        # GH 10132
        for key in [('a', 1), ('b', 2), ('b', 1), ('a', 2)]:
            c, i = key
            result = groups_double_key.get_group(key)
            expected = test[(test.cat == c) & (test.ints == i)]
            assert_frame_equal(result, expected)

        d = {'C1': [3, 3, 4, 5], 'C2': [1, 2, 3, 4], 'C3': [10, 100, 200, 34]}
        test = pd.DataFrame(d)
        values = pd.cut(test['C1'], [1, 2, 3, 6])
        values.name = "cat"
        groups_double_key = test.groupby([values, 'C2'])

        res = groups_double_key.agg('mean')
        nan = np.nan
        idx = MultiIndex.from_product(
            [Categorical(["(1, 2]", "(2, 3]", "(3, 6]"], ordered=True),
             [1, 2, 3, 4]],
            names=["cat", "C2"])
        exp = DataFrame({"C1": [nan, nan, nan, nan, 3, 3,
                                nan, nan, nan, nan, 4, 5],
                         "C3": [nan, nan, nan, nan, 10, 100,
                                nan, nan, nan, nan, 200, 34]}, index=idx)
        tm.assert_frame_equal(res, exp)

    def test_groupby_multi_categorical_as_index(self):
        # GH13204
        df = DataFrame({'cat': Categorical([1, 2, 2], [1, 2, 3]),
                        'A': [10, 11, 11],
                        'B': [101, 102, 103]})
        result = df.groupby(['cat', 'A'], as_index=False).sum()
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10, 11, 10, 11, 10, 11],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])
        tm.assert_frame_equal(result, expected)

        # function grouper
        f = lambda r: df.loc[r, 'A']
        result = df.groupby(['cat', f], as_index=False).sum()
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10.0, nan, nan, 22.0, nan, nan],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])
        tm.assert_frame_equal(result, expected)

        # another not in-axis grouper (conflicting names in index)
        s = Series(['a', 'b', 'b'], name='cat')
        result = df.groupby(['cat', s], as_index=False).sum()
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10.0, nan, nan, 22.0, nan, nan],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])
        tm.assert_frame_equal(result, expected)

        # is original index dropped?
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10, 11, 10, 11, 10, 11],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])

        for name in [None, 'X', 'B', 'cat']:
            df.index = Index(list("abc"), name=name)
            result = df.groupby(['cat', 'A'], as_index=False).sum()
            tm.assert_frame_equal(result, expected, check_index_type=True)

    def test_groupby_preserve_categorical_dtype(self):
        # GH13743, GH13854
        df = DataFrame({'A': [1, 2, 1, 1, 2],
                        'B': [10, 16, 22, 28, 34],
                        'C1': Categorical(list("abaab"),
                                          categories=list("bac"),
                                          ordered=False),
                        'C2': Categorical(list("abaab"),
                                          categories=list("bac"),
                                          ordered=True)})
        # single grouper
        exp_full = DataFrame({'A': [2.0, 1.0, np.nan],
                              'B': [25.0, 20.0, np.nan],
                              'C1': Categorical(list("bac"),
                                                categories=list("bac"),
                                                ordered=False),
                              'C2': Categorical(list("bac"),
                                                categories=list("bac"),
                                                ordered=True)})
        for col in ['C1', 'C2']:
            result1 = df.groupby(by=col, as_index=False).mean()
            result2 = df.groupby(by=col, as_index=True).mean().reset_index()
            expected = exp_full.reindex(columns=result1.columns)
            tm.assert_frame_equal(result1, expected)
            tm.assert_frame_equal(result2, expected)

        # multiple grouper
        exp_full = DataFrame({'A': [1, 1, 1, 2, 2, 2],
                              'B': [np.nan, 20.0, np.nan, 25.0, np.nan,
                                    np.nan],
                              'C1': Categorical(list("bacbac"),
                                                categories=list("bac"),
                                                ordered=False),
                              'C2': Categorical(list("bacbac"),
                                                categories=list("bac"),
                                                ordered=True)})
        for cols in [['A', 'C1'], ['A', 'C2']]:
            result1 = df.groupby(by=cols, as_index=False).mean()
            result2 = df.groupby(by=cols, as_index=True).mean().reset_index()
            expected = exp_full.reindex(columns=result1.columns)
            tm.assert_frame_equal(result1, expected)
            tm.assert_frame_equal(result2, expected)

    def test_groupby_apply_all_none(self):
        # Tests to make sure no errors if apply function returns all None
        # values. Issue 9684.
        test_df = DataFrame({'groups': [0, 0, 1, 1],
                             'random_vars': [8, 7, 4, 5]})

        def test_func(x):
            pass

        result = test_df.groupby('groups').apply(test_func)
        expected = DataFrame()
        tm.assert_frame_equal(result, expected)

    def test_groupby_apply_none_first(self):
        # GH 12824. Tests if apply returns None first.
        test_df1 = DataFrame({'groups': [1, 1, 1, 2], 'vars': [0, 1, 2, 3]})
        test_df2 = DataFrame({'groups': [1, 2, 2, 2], 'vars': [0, 1, 2, 3]})

        def test_func(x):
            if x.shape[0] < 2:
                return None
            return x.iloc[[0, -1]]

        result1 = test_df1.groupby('groups').apply(test_func)
        result2 = test_df2.groupby('groups').apply(test_func)
        index1 = MultiIndex.from_arrays([[1, 1], [0, 2]],
                                        names=['groups', None])
        index2 = MultiIndex.from_arrays([[2, 2], [1, 3]],
                                        names=['groups', None])
        expected1 = DataFrame({'groups': [1, 1], 'vars': [0, 2]},
                              index=index1)
        expected2 = DataFrame({'groups': [2, 2], 'vars': [1, 3]},
                              index=index2)
        tm.assert_frame_equal(result1, expected1)
        tm.assert_frame_equal(result2, expected2)

    def test_first_last_max_min_on_time_data(self):
        # GH 10295
        # Verify that NaT is not in the result of max, min, first and last on
        # Dataframe with datetime or timedelta values.
        from datetime import timedelta as td
        df_test = DataFrame(
            {'dt': [nan, '2015-07-24 10:10', '2015-07-25 11:11',
                    '2015-07-23 12:12', nan],
             'td': [nan, td(days=1), td(days=2), td(days=3), nan]})
        df_test.dt = pd.to_datetime(df_test.dt)
        df_test['group'] = 'A'
        df_ref = df_test[df_test.dt.notnull()]

        grouped_test = df_test.groupby('group')
        grouped_ref = df_ref.groupby('group')

        assert_frame_equal(grouped_ref.max(), grouped_test.max())
        assert_frame_equal(grouped_ref.min(), grouped_test.min())
        assert_frame_equal(grouped_ref.first(), grouped_test.first())
        assert_frame_equal(grouped_ref.last(), grouped_test.last())

    def test_groupby_preserves_sort(self):
        # Test to ensure that groupby always preserves sort order of original
        # object. Issue #8588 and #9651

        df = DataFrame(
            {'int_groups': [3, 1, 0, 1, 0, 3, 3, 3],
             'string_groups': ['z', 'a', 'z', 'a', 'a', 'g', 'g', 'g'],
             'ints': [8, 7, 4, 5, 2, 9, 1, 1],
             'floats': [2.3, 5.3, 6.2, -2.4, 2.2, 1.1, 1.1, 5],
             'strings': ['z', 'd', 'a', 'e', 'word', 'word2', '42', '47']})

        # Try sorting on different types and with different group types
        for sort_column in ['ints', 'floats', 'strings', ['ints', 'floats'],
                            ['ints', 'strings']]:
            for group_column in ['int_groups', 'string_groups',
                                 ['int_groups', 'string_groups']]:

                df = df.sort_values(by=sort_column)

                g = df.groupby(group_column)

                def test_sort(x):
                    assert_frame_equal(x, x.sort_values(by=sort_column))

                g.apply(test_sort)

    def test_nunique_with_object(self):
        # GH 11077
        data = pd.DataFrame(
            [[100, 1, 'Alice'],
             [200, 2, 'Bob'],
             [300, 3, 'Charlie'],
             [-400, 4, 'Dan'],
             [500, 5, 'Edith']],
            columns=['amount', 'id', 'name']
        )

        result = data.groupby(['id', 'amount'])['name'].nunique()
        index = MultiIndex.from_arrays([data.id, data.amount])
        expected = pd.Series([1] * 5, name='name', index=index)
        tm.assert_series_equal(result, expected)

    def test_nunique_with_empty_series(self):
        # GH 12553
        data = pd.Series(name='name')
        result = data.groupby(level=0).nunique()
        expected = pd.Series(name='name', dtype='int64')
        tm.assert_series_equal(result, expected)

    def test_transform_with_non_scalar_group(self):
        # GH 10165
        cols = pd.MultiIndex.from_tuples([
            ('syn', 'A'), ('mis', 'A'), ('non', 'A'),
            ('syn', 'C'), ('mis', 'C'), ('non', 'C'),
            ('syn', 'T'), ('mis', 'T'), ('non', 'T'),
            ('syn', 'G'), ('mis', 'G'), ('non', 'G')])
        df = pd.DataFrame(np.random.randint(1, 10, (4, 12)),
                          columns=cols,
                          index=['A', 'C', 'G', 'T'])
        self.assertRaisesRegexp(ValueError, 'transform must return a scalar '
                                'value for each group.*', df.groupby
                                (axis=1, level=1).transform,
                                lambda z: z.div(z.sum(axis=1), axis=0))

    def test_numpy_compat(self):
        # see gh-12811
        df = pd.DataFrame({'A': [1, 2, 1], 'B': [1, 2, 3]})
        g = df.groupby('A')

        msg = "numpy operations are not valid with groupby"

        for func in ('mean', 'var', 'std', 'cumprod', 'cumsum'):
            tm.assertRaisesRegexp(UnsupportedFunctionCall, msg,
                                  getattr(g, func), 1, 2, 3)
            tm.assertRaisesRegexp(UnsupportedFunctionCall, msg,
                                  getattr(g, func), foo=1)

    def test_grouping_string_repr(self):
        # GH 13394
        mi = MultiIndex.from_arrays([list("AAB"), list("aba")])
        df = DataFrame([[1, 2, 3]], columns=mi)
        gr = df.groupby(df[('A', 'a')])

        result = gr.grouper.groupings[0].__repr__()
        expected = "Grouping(('A', 'a'))"
        tm.assert_equal(result, expected)

    def test_group_shift_with_null_key(self):
        # This test is designed to replicate the segfault in issue #13813.
        n_rows = 1200

        # Generate a moderately large dataframe with occasional missing
        # values in column `B`, and then group by [`A`, `B`]. This should
        # force `-1` in `labels` array of `g.grouper.group_info` exactly
        # at those places, where the group-by key is partilly missing.
        df = DataFrame([(i % 12, i % 3 if i % 3 else np.nan, i)
                        for i in range(n_rows)], dtype=float,
                       columns=["A", "B", "Z"], index=None)
        g = df.groupby(["A", "B"])

        expected = DataFrame([(i + 12 if i % 3 and i < n_rows - 12
                              else np.nan)
                             for i in range(n_rows)], dtype=float,
                             columns=["Z"], index=None)
        result = g.shift(-1)

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
