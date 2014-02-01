from __future__ import print_function
import nose

from numpy.testing.decorators import slow

from datetime import datetime
from numpy import nan

from pandas import date_range,bdate_range, Timestamp
from pandas.core.index import Index, MultiIndex, Int64Index
from pandas.core.common import rands
from pandas.core.api import Categorical, DataFrame
from pandas.core.groupby import SpecificationError, DataError
from pandas.core.series import Series
from pandas.util.testing import (assert_panel_equal, assert_frame_equal,
                                 assert_series_equal, assert_almost_equal,
                                 assert_index_equal)
from pandas.compat import(
    range, long, lrange, StringIO, lmap, lzip, map, zip, builtins, OrderedDict
)
from pandas import compat
from pandas.core.panel import Panel
from pandas.tools.merge import concat
from collections import defaultdict
import pandas.core.common as com
import numpy as np

import pandas.core.nanops as nanops

import pandas.util.testing as tm
import pandas as pd

def commonSetUp(self):
    self.dateRange = bdate_range('1/1/2005', periods=250)
    self.stringIndex = Index([rands(8).upper() for x in range(250)])

    self.groupId = Series([x[0] for x in self.stringIndex],
                          index=self.stringIndex)
    self.groupDict = dict((k, v) for k, v in compat.iteritems(self.groupId))

    self.columnIndex = Index(['A', 'B', 'C', 'D', 'E'])

    randMat = np.random.randn(250, 5)
    self.stringMatrix = DataFrame(randMat, columns=self.columnIndex,
                                  index=self.stringIndex)

    self.timeMatrix = DataFrame(randMat, columns=self.columnIndex,
                                index=self.dateRange)


class TestGroupBy(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.ts = tm.makeTimeSeries()

        self.seriesd = tm.getSeriesData()
        self.tsd = tm.getTimeSeriesData()
        self.frame = DataFrame(self.seriesd)
        self.tsframe = DataFrame(self.tsd)

        self.df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                                   'foo', 'bar', 'foo', 'foo'],
                             'B': ['one', 'one', 'two', 'three',
                                   'two', 'two', 'one', 'three'],
                             'C': np.random.randn(8),
                             'D': np.random.randn(8)})

        self.df_mixed_floats = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                                                'foo', 'bar', 'foo', 'foo'],
                                          'B': ['one', 'one', 'two', 'three',
                                                'two', 'two', 'one', 'three'],
                                          'C': np.random.randn(8),
                                          'D': np.array(np.random.randn(8),dtype='float32')})

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        self.mframe = DataFrame(np.random.randn(10, 3), index=index,
                                columns=['A', 'B', 'C'])

        self.three_group = DataFrame({'A': ['foo', 'foo', 'foo', 'foo',
                                            'bar', 'bar', 'bar', 'bar',
                                            'foo', 'foo', 'foo'],
                                      'B': ['one', 'one', 'one', 'two',
                                            'one', 'one', 'one', 'two',
                                            'two', 'two', 'one'],
                                      'C': ['dull', 'dull', 'shiny', 'dull',
                                            'dull', 'shiny', 'shiny', 'dull',
                                            'shiny', 'shiny', 'shiny'],
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
            assert_series_equal(grouped.agg(np.sum),grouped.sum())

            transformed = grouped.transform(lambda x: x * x.sum())
            self.assertEqual(transformed[7], 12)

            value_grouped = data.groupby(data)
            assert_series_equal(value_grouped.aggregate(np.mean), agged)

            # complex agg
            agged = grouped.aggregate([np.mean, np.std])
            agged = grouped.aggregate({'one': np.mean,
                                       'two': np.std})

            group_constants = {
                0: 10,
                1: 20,
                2: 30
                }
            agged = grouped.agg(lambda x: group_constants[x.name] + x.mean())
            self.assertEqual(agged[1], 21)

            # corner cases
            self.assertRaises(Exception, grouped.aggregate, lambda x: x * 2)


        for dtype in ['int64','int32','float64','float32']:
            checkit(dtype)

    def test_first_last_nth(self):
        # tests for first / last / nth
        grouped = self.df.groupby('A')
        first = grouped.first()
        expected = self.df.ix[[1, 0], ['B', 'C', 'D']]
        expected.index = ['bar', 'foo']
        assert_frame_equal(first, expected, check_names=False)

        last = grouped.last()
        expected = self.df.ix[[5, 7], ['B', 'C', 'D']]
        expected.index = ['bar', 'foo']
        assert_frame_equal(last, expected, check_names=False)

        nth = grouped.nth(1)
        expected = self.df.ix[[3, 2], ['B', 'C', 'D']]
        expected.index = ['bar', 'foo']
        assert_frame_equal(nth, expected, check_names=False)

        # it works!
        grouped['B'].first()
        grouped['B'].last()
        grouped['B'].nth(0)

        self.df['B'][self.df['A'] == 'foo'] = np.nan
        self.assert_(com.isnull(grouped['B'].first()['foo']))
        self.assert_(com.isnull(grouped['B'].last()['foo']))
        self.assert_(com.isnull(grouped['B'].nth(0)['foo']))

    def test_first_last_nth_dtypes(self):

        df = self.df_mixed_floats.copy()
        df['E'] = True
        df['F'] = 1

        # tests for first / last / nth
        grouped = df.groupby('A')
        first = grouped.first()
        expected = df.ix[[1, 0], ['B', 'C', 'D', 'E', 'F']]
        expected.index = ['bar', 'foo']
        assert_frame_equal(first, expected, check_names=False)

        last = grouped.last()
        expected = df.ix[[5, 7], ['B', 'C', 'D', 'E', 'F']]
        expected.index = ['bar', 'foo']
        assert_frame_equal(last, expected, check_names=False)

        nth = grouped.nth(1)
        expected = df.ix[[3, 2], ['B', 'C', 'D', 'E', 'F']]
        expected.index = ['bar', 'foo']
        assert_frame_equal(nth, expected, check_names=False)

        # GH 2763, first/last shifting dtypes
        idx = lrange(10)
        idx.append(9)
        s = Series(data=lrange(11), index=idx, name='IntCol')
        self.assert_(s.dtype == 'int64')
        f = s.groupby(level=0).first()
        self.assert_(f.dtype == 'int64')

    def test_grouper_index_types(self):
        # related GH5375
        # groupby misbehaving when using a Floatlike index
        df = DataFrame(np.arange(10).reshape(5,2),columns=list('AB'))
        for index in [ tm.makeFloatIndex, tm.makeStringIndex,
                       tm.makeUnicodeIndex, tm.makeIntIndex,
                       tm.makeDateIndex, tm.makePeriodIndex ]:

            df.index = index(len(df))
            df.groupby(list('abcde')).apply(lambda x: x)

            df.index = list(reversed(df.index.tolist()))
            df.groupby(list('abcde')).apply(lambda x: x)

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

    def test_groupby_bounds_check(self):
        import pandas as pd
        # groupby_X is code-generated, so if one variant
        # does, the rest probably do to
        a = np.array([1,2],dtype='object')
        b = np.array([1,2,3],dtype='object')
        self.assertRaises(AssertionError, pd.algos.groupby_object,a, b)

    def test_groupby_grouper_f_sanity_checked(self):
        import pandas as pd
        dates = date_range('01-Jan-2013', periods=12, freq='MS')
        ts = pd.TimeSeries(np.random.randn(12), index=dates)

        # GH3035
        # index.map is used to apply grouper to the index
        # if it fails on the elements, map tries it on the entire index as
        # a sequence. That can yield invalid results that cause trouble
        # down the line.
        # the surprise comes from using key[0:6] rather then str(key)[0:6]
        # when the elements are Timestamp.
        # the result is Index[0:6], very confusing.

        self.assertRaises(AssertionError, ts.groupby,lambda key: key[0:6])

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
        result = applied.get_dtype_counts()
        result.sort()
        expected = Series({ 'object' : 2, 'float64' : 2, 'int64' : 1 })
        expected.sort()
        assert_series_equal(result,expected)

    def test_groupby_return_type(self):

        # GH2893, return a reduced type
        df1 = DataFrame([{"val1": 1, "val2" : 20}, {"val1":1, "val2": 19},
                         {"val1":2, "val2": 27}, {"val1":2, "val2": 12}])

        def func(dataf):
            return dataf["val2"]  - dataf["val2"].mean()

        result = df1.groupby("val1", squeeze=True).apply(func)
        tm.assert_isinstance(result,Series)

        df2 = DataFrame([{"val1": 1, "val2" : 20}, {"val1":1, "val2": 19},
                         {"val1":1, "val2": 27}, {"val1":1, "val2": 12}])
        def func(dataf):
            return dataf["val2"]  - dataf["val2"].mean()

        result = df2.groupby("val1", squeeze=True).apply(func)
        tm.assert_isinstance(result,Series)

        # GH3596, return a consistent type (regression in 0.11 from 0.10.1)
        df = DataFrame([[1,1],[1,1]],columns=['X','Y'])
        result = df.groupby('X',squeeze=False).count()
        tm.assert_isinstance(result,DataFrame)

        # GH5592
        # inconcistent return type
        df = DataFrame(dict(A = [ 'Tiger', 'Tiger', 'Tiger', 'Lamb', 'Lamb', 'Pony', 'Pony' ],
                            B = Series(np.arange(7),dtype='int64'),
                            C = date_range('20130101',periods=7)))

        def f(grp):
            return grp.iloc[0]
        expected = df.groupby('A').first()[['B']]
        result = df.groupby('A').apply(f)[['B']]
        assert_frame_equal(result,expected)

        def f(grp):
            if grp.name == 'Tiger':
                return None
            return grp.iloc[0]
        result = df.groupby('A').apply(f)[['B']]
        e = expected.copy()
        e.loc['Tiger'] = np.nan
        assert_frame_equal(result,e)

        def f(grp):
            if grp.name == 'Pony':
                return None
            return grp.iloc[0]
        result = df.groupby('A').apply(f)[['B']]
        e = expected.copy()
        e.loc['Pony'] = np.nan
        assert_frame_equal(result,e)

        # 5592 revisited, with datetimes
        def f(grp):
            if grp.name == 'Pony':
                return None
            return grp.iloc[0]
        result = df.groupby('A').apply(f)[['C']]
        e = df.groupby('A').first()[['C']]
        e.loc['Pony'] = np.nan
        assert_frame_equal(result,e)

        # scalar outputs
        def f(grp):
            if grp.name == 'Pony':
                return None
            return grp.iloc[0].loc['C']
        result = df.groupby('A').apply(f)
        e = df.groupby('A').first()['C'].copy()
        e.loc['Pony'] = np.nan
        e.name = None
        assert_series_equal(result,e)

    def test_agg_regression1(self):
        grouped = self.tsframe.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.agg(np.mean)
        expected = grouped.mean()
        assert_frame_equal(result, expected)

    def test_agg_datetimes_mixed(self):
        data = [[1, '2012-01-01', 1.0],
                [2, '2012-01-02', 2.0],
                [3, None, 3.0]]

        df1 = DataFrame({'key': [x[0] for x in data],
                         'date': [x[1] for x in data],
                         'value': [x[2] for x in data]})

        data = [[row[0], datetime.strptime(row[1], '%Y-%m-%d').date()
                if row[1] else None, row[2]] for row in data]

        df2 = DataFrame({'key': [x[0] for x in data],
                         'date': [x[1] for x in data],
                         'value': [x[2] for x in data]})

        df1['weights'] = df1['value'] / df1['value'].sum()
        gb1 = df1.groupby('date').aggregate(np.sum)

        df2['weights'] = df1['value'] / df1['value'].sum()
        gb2 = df2.groupby('date').aggregate(np.sum)

        assert(len(gb1) == len(gb2))

    def test_agg_period_index(self):
        from pandas import period_range, PeriodIndex
        prng = period_range('2012-1-1', freq='M', periods=3)
        df = DataFrame(np.random.randn(3, 2), index=prng)
        rs = df.groupby(level=0).sum()
        tm.assert_isinstance(rs.index, PeriodIndex)

        # GH 3579
        index = period_range(start='1999-01', periods=5, freq='M')
        s1 = Series(np.random.rand(len(index)), index=index)
        s2 = Series(np.random.rand(len(index)), index=index)
        series = [('s1', s1), ('s2',s2)]
        df = DataFrame.from_items(series)
        grouped = df.groupby(df.index.month)
        list(grouped)

    def test_agg_must_agg(self):
        grouped = self.df.groupby('A')['C']
        self.assertRaises(Exception, grouped.agg, lambda x: x.describe())
        self.assertRaises(Exception, grouped.agg, lambda x: x.index[:2])

    def test_agg_ser_multi_key(self):
        ser = self.df.C
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

    def test_agg_apply_corner(self):
        # nothing to group, all NA
        grouped = self.ts.groupby(self.ts * np.nan)

        assert_series_equal(grouped.sum(), Series([]))
        assert_series_equal(grouped.agg(np.sum), Series([]))
        assert_series_equal(grouped.apply(np.sum), Series([]))

        # DataFrame
        grouped = self.tsframe.groupby(self.tsframe['A'] * np.nan)
        exp_df = DataFrame(columns=self.tsframe.columns, dtype=float)
        assert_frame_equal(grouped.sum(), exp_df, check_names=False)
        assert_frame_equal(grouped.agg(np.sum), exp_df, check_names=False)
        assert_frame_equal(grouped.apply(np.sum), DataFrame({}, dtype=float))

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

    def test_agg_python_multiindex(self):
        grouped = self.mframe.groupby(['A', 'B'])

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    def test_apply_describe_bug(self):
        grouped = self.mframe.groupby(level='first')
        result = grouped.describe()  # it works!

    def test_apply_issues(self):
        # GH 5788

        s="""2011.05.16,00:00,1.40893
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

        df = pd.read_csv(StringIO(s), header=None, names=['date', 'time', 'value'], parse_dates=[['date', 'time']])
        df = df.set_index('date_time')

        expected = df.groupby(df.index.date).idxmax()
        result = df.groupby(df.index.date).apply(lambda x: x.idxmax())
        assert_frame_equal(result,expected)

        # GH 5789
        # don't auto coerce dates
        df = pd.read_csv(StringIO(s), header=None, names=['date', 'time', 'value'])
        expected = Series(['00:00','02:00','02:00'],index=['2011.05.16','2011.05.17','2011.05.18'])
        result = df.groupby('date').apply(lambda x: x['time'][x['value'].idxmax()])
        assert_series_equal(result,expected)

    def test_len(self):
        df = tm.makeTimeDataFrame()
        grouped = df.groupby([lambda x: x.year,
                              lambda x: x.month,
                              lambda x: x.day])
        self.assertEquals(len(grouped), len(df))

        grouped = df.groupby([lambda x: x.year,
                              lambda x: x.month])
        expected = len(set([(x.year, x.month) for x in df.index]))
        self.assertEquals(len(grouped), expected)

    def test_groups(self):
        grouped = self.df.groupby(['A'])
        groups = grouped.groups
        self.assert_(groups is grouped.groups)  # caching works

        for k, v in compat.iteritems(grouped.groups):
            self.assert_((self.df.ix[v]['A'] == k).all())

        grouped = self.df.groupby(['A', 'B'])
        groups = grouped.groups
        self.assert_(groups is grouped.groups)  # caching works
        for k, v in compat.iteritems(grouped.groups):
            self.assert_((self.df.ix[v]['A'] == k[0]).all())
            self.assert_((self.df.ix[v]['B'] == k[1]).all())

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
            result = grouped.agg(
                OrderedDict([['A', 'var'], ['B', 'std'], ['C', 'mean']]))
            expected = DataFrame(OrderedDict([['A', grouped['A'].var()],
                                              ['B', grouped['B'].std()],
                                              ['C', grouped['C'].mean()]]))
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
        assert_almost_equal(result.xs('foo'), np.array([foo] * K).astype('float64'))
        assert_almost_equal(result.xs('bar'), np.array([bar] * K).astype('float64'))

        def aggfun(ser):
            return ser.size
        result = DataFrame().groupby(self.df.A).agg(aggfun)
        tm.assert_isinstance(result, DataFrame)
        self.assertEqual(len(result), 0)

    def test_agg_item_by_item_raise_typeerror(self):
        from numpy.random import randint

        df = DataFrame(randint(10, size=(20, 10)))

        def raiseException(df):
            print('----------------------------------------')
            print(df.to_string())
            raise TypeError

        self.assertRaises(TypeError, df.groupby(0).agg,
                          raiseException)

    def test_basic_regression(self):
        # regression
        T = [1.0 * x for x in lrange(1, 10) * 10][:1095]
        result = Series(T, lrange(0, len(T)))

        groupings = np.random.random((1100,))
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

    def test_transform_broadcast(self):
        grouped = self.ts.groupby(lambda x: x.month)
        result = grouped.transform(np.mean)

        self.assert_(result.index.equals(self.ts.index))
        for _, gp in grouped:
            assert_fp_equal(result.reindex(gp.index), gp.mean())

        grouped = self.tsframe.groupby(lambda x: x.month)
        result = grouped.transform(np.mean)
        self.assert_(result.index.equals(self.tsframe.index))
        for _, gp in grouped:
            agged = gp.mean()
            res = result.reindex(gp.index)
            for col in self.tsframe:
                assert_fp_equal(res[col], agged[col])

        # group columns
        grouped = self.tsframe.groupby({'A': 0, 'B': 0, 'C': 1, 'D': 1},
                                       axis=1)
        result = grouped.transform(np.mean)
        self.assert_(result.index.equals(self.tsframe.index))
        self.assert_(result.columns.equals(self.tsframe.columns))
        for _, gp in grouped:
            agged = gp.mean(1)
            res = result.reindex(columns=gp.columns)
            for idx in gp.index:
                assert_fp_equal(res.xs(idx), agged[idx])

    def test_transform_bug(self):
        # GH 5712
        # transforming on a datetime column
        df = DataFrame(dict(A = Timestamp('20130101'), B = np.arange(5)))
        result = df.groupby('A')['B'].transform(lambda x: x.rank(ascending=False))
        expected = Series(np.arange(5,0,step=-1),name='B')
        assert_series_equal(result,expected)

    def test_transform_multiple(self):
        grouped = self.ts.groupby([lambda x: x.year, lambda x: x.month])

        transformed = grouped.transform(lambda x: x * 2)
        broadcasted = grouped.transform(np.mean)

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

    def test_with_na(self):
        index = Index(np.arange(10))

        for dtype in ['float64','float32','int64','int32','int16','int8']:
            values = Series(np.ones(10), index, dtype=dtype)
            labels = Series([nan, 'foo', 'bar', 'bar', nan, nan, 'bar',
                             'bar', nan, 'foo'], index=index)


            # this SHOULD be an int
            grouped = values.groupby(labels)
            agged = grouped.agg(len)
            expected = Series([4, 2], index=['bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)
            #self.assert_(issubclass(agged.dtype.type, np.integer))

            # explicity return a float from my function
            def f(x):
                return float(len(x))

            agged = grouped.agg(f)
            expected = Series([4, 2], index=['bar', 'foo'])

            assert_series_equal(agged, expected, check_dtype=False)
            self.assert_(issubclass(agged.dtype.type, np.dtype(dtype).type))

    def test_groupby_transform_with_int(self):

        # GH 3740, make sure that we might upcast on item-by-item transform

        # floats
        df = DataFrame(dict(A = [1,1,1,2,2,2], B = Series(1,dtype='float64'), C = Series([1,2,3,1,2,3],dtype='float64'), D = 'foo'))
        result = df.groupby('A').transform(lambda x: (x-x.mean())/x.std())
        expected = DataFrame(dict(B = np.nan, C = Series([-1,0,1,-1,0,1],dtype='float64')))
        assert_frame_equal(result,expected)

        # int case
        df = DataFrame(dict(A = [1,1,1,2,2,2], B = 1, C = [1,2,3,1,2,3], D = 'foo'))
        result = df.groupby('A').transform(lambda x: (x-x.mean())/x.std())
        expected = DataFrame(dict(B = np.nan, C = [-1,0,1,-1,0,1]))
        assert_frame_equal(result,expected)

        # int that needs float conversion
        s = Series([2,3,4,10,5,-1])
        df = DataFrame(dict(A = [1,1,1,2,2,2], B = 1, C = s, D = 'foo'))
        result = df.groupby('A').transform(lambda x: (x-x.mean())/x.std())

        s1 = s.iloc[0:3]
        s1 = (s1-s1.mean())/s1.std()
        s2 = s.iloc[3:6]
        s2 = (s2-s2.mean())/s2.std()
        expected = DataFrame(dict(B = np.nan, C = concat([s1,s2])))
        assert_frame_equal(result,expected)

        # int downcasting
        result = df.groupby('A').transform(lambda x: x*2/2)
        expected = DataFrame(dict(B = 1, C = [2,3,4,10,5,-1]))
        assert_frame_equal(result,expected)

    def test_indices_concatenation_order(self):

        # GH 2808

        def f1(x):
            y = x[(x.b % 2) == 1]**2
            if y.empty:
                multiindex = MultiIndex(
                    levels = [[]]*2,
                    labels = [[]]*2,
                    names = ['b', 'c']
                    )
                res = DataFrame(None,
                                   columns=['a'],
                                   index=multiindex)
                return res
            else:
                y = y.set_index(['b','c'])
                return y

        def f2(x):
            y = x[(x.b % 2) == 1]**2
            if y.empty:
                return DataFrame()
            else:
                y = y.set_index(['b','c'])
                return y

        def f3(x):
            y = x[(x.b % 2) == 1]**2
            if y.empty:
                multiindex = MultiIndex(
                    levels = [[]]*2,
                    labels = [[]]*2,
                    names = ['foo', 'bar']
                    )
                res = DataFrame(None,
                                columns=['a','b'],
                                index=multiindex)
                return res
            else:
                return y

        df = DataFrame({'a':[1,2,2,2],
                        'b':lrange(4),
                        'c':lrange(5,9)})

        df2 = DataFrame({'a':[3,2,2,2],
                         'b':lrange(4),
                         'c':lrange(5,9)})


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
        assert_series_equal(result['mean'], grouped.mean())
        assert_series_equal(result['std'], grouped.std())
        assert_series_equal(result['min'], grouped.min())

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
        data = DataFrame({'A': ['foo', 'foo', 'foo', 'foo',
                                'bar', 'bar', 'bar', 'bar',
                                'foo', 'foo', 'foo'],
                          'B': ['one', 'one', 'one', 'two',
                                'one', 'one', 'one', 'two',
                                'two', 'two', 'one'],
                          'C': ['dull', 'dull', 'shiny', 'dull',
                                'dull', 'shiny', 'shiny', 'dull',
                                'shiny', 'shiny', 'shiny'],
                          'D': np.random.randn(11),
                          'E': np.random.randn(11),
                          'F': np.random.randn(11)})

        def bad(x):
            assert(len(x.base) > 0)
            return 'foo'

        result = data.groupby(['A', 'B']).agg(bad)
        expected = data.groupby(['A', 'B']).agg(lambda x: 'foo')
        assert_frame_equal(result, expected)

    def test_series_index_name(self):
        grouped = self.df.ix[:, ['C']].groupby(self.df['A'])
        result = grouped.agg(lambda x: x.mean())
        self.assertEqual(result.index.name, 'A')

    def test_frame_describe_multikey(self):
        grouped = self.tsframe.groupby([lambda x: x.year,
                                        lambda x: x.month])
        result = grouped.describe()

        for col in self.tsframe:
            expected = grouped[col].describe()
            assert_series_equal(result[col], expected)

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
        transformed = grouped.transform(lambda x: x - x.mean())
        self.assertEqual(len(transformed), 30)
        self.assertEqual(len(transformed.columns), 4)

        # transform propagate
        transformed = grouped.transform(lambda x: x.mean())
        for name, group in grouped:
            mean = group.mean()
            for idx in group.index:
                assert_almost_equal(transformed.xs(idx), mean)

        # iterate
        for weekday, group in grouped:
            self.assert_(group.index[0].weekday() == weekday)

        # groups / group_indices
        groups = grouped.groups
        indices = grouped.indices

        for k, v in compat.iteritems(groups):
            samething = self.tsframe.index.take(indices[k])
            self.assertTrue((samething == v).all())

    def test_grouping_is_iterable(self):
        # this code path isn't used anywhere else
        # not sure it's useful
        grouped = self.tsframe.groupby([lambda x: x.weekday(),
                                        lambda x: x.year])

        # test it works
        for g in grouped.grouper.groupings[0]:
            pass

    def test_frame_groupby_columns(self):
        mapping = {
            'A': 0, 'B': 0, 'C': 1, 'D': 1
        }
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
        self.assert_(result.index.name == 'A')

        result = self.df.groupby('A', as_index=False).mean()
        self.assert_(result.index.name != 'A')

        result = grouped.agg(np.mean)
        self.assert_(result.index.name == 'A')

        result = grouped.agg({'C': np.mean, 'D': np.std})
        self.assert_(result.index.name == 'A')

        result = grouped['C'].mean()
        self.assert_(result.index.name == 'A')
        result = grouped['C'].agg(np.mean)
        self.assert_(result.index.name == 'A')
        result = grouped['C'].agg([np.mean, np.std])
        self.assert_(result.index.name == 'A')

        result = grouped['C'].agg({'foo': np.mean, 'bar': np.std})
        self.assert_(result.index.name == 'A')

    def test_multi_iter(self):
        s = Series(np.arange(6))
        k1 = np.array(['a', 'a', 'a', 'b', 'b', 'b'])
        k2 = np.array(['1', '2', '1', '2', '1', '2'])

        grouped = s.groupby([k1, k2])

        iterated = list(grouped)
        expected = [('a', '1', s[[0, 2]]),
                    ('a', '2', s[[1]]),
                    ('b', '1', s[[4]]),
                    ('b', '2', s[[3, 5]])]
        for i, ((one, two), three) in enumerate(iterated):
            e1, e2, e3 = expected[i]
            self.assert_(e1 == one)
            self.assert_(e2 == two)
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
            self.assert_(e1 == one)
            self.assert_(e2 == two)
            assert_frame_equal(three, e3)

        # don't iterate through groups with no data
        df['k1'] = np.array(['b', 'b', 'b', 'a', 'a', 'a'])
        df['k2'] = np.array(['1', '1', '1', '2', '2', '2'])
        grouped = df.groupby(['k1', 'k2'])
        groups = {}
        for key, gp in grouped:
            groups[key] = gp
        self.assertEquals(len(groups), 2)

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
            exp_axis = [x for x in wp.major_axis
                        if x.month == month and x.weekday() == wd]
            expected = wp.reindex(major=exp_axis)
            assert_panel_equal(group, expected)

    def test_multi_func(self):
        col1 = self.df['A']
        col2 = self.df['B']

        grouped = self.df.groupby([col1.get, col2.get])
        agged = grouped.mean()
        expected = self.df.groupby(['A', 'B']).mean()
        assert_frame_equal(agged.ix[:, ['C', 'D']],
                           expected.ix[:, ['C', 'D']],
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
        data = DataFrame({'A': ['foo', 'foo', 'foo', 'foo',
                                'bar', 'bar', 'bar', 'bar',
                                'foo', 'foo', 'foo'],
                          'B': ['one', 'one', 'one', 'two',
                                'one', 'one', 'one', 'two',
                                'two', 'two', 'one'],
                          'C': ['dull', 'dull', 'shiny', 'dull',
                                'dull', 'shiny', 'shiny', 'dull',
                                'shiny', 'shiny', 'shiny'],
                          'D': np.random.randn(11),
                          'E': np.random.randn(11),
                          'F': np.random.randn(11)})

        grouped = data.groupby(['A', 'B'])
        funcs = [np.mean, np.std]
        agged = grouped.agg(funcs)
        expected = concat([grouped['D'].agg(funcs), grouped['E'].agg(funcs),
                           grouped['F'].agg(funcs)],
                          keys=['D', 'E', 'F'], axis=1)
        assert(isinstance(agged.index, MultiIndex))
        assert(isinstance(expected.index, MultiIndex))
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
            expected = dict((k, DataFrame(v)) for k, v in compat.iteritems(expected))
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

    def test_multifunc_select_col_integer_cols(self):
        df = self.df
        df.columns = np.arange(len(df.columns))

        # it works!
        result = df.groupby(1, as_index=False)[2].agg({'Q': np.mean})

    def test_as_index_series_return_frame(self):
        grouped = self.df.groupby('A', as_index=False)
        grouped2 = self.df.groupby(['A', 'B'], as_index=False)

        result = grouped['C'].agg(np.sum)
        expected = grouped.agg(np.sum).ix[:, ['A', 'C']]
        tm.assert_isinstance(result, DataFrame)
        assert_frame_equal(result, expected)

        result2 = grouped2['C'].agg(np.sum)
        expected2 = grouped2.agg(np.sum).ix[:, ['A', 'B', 'C']]
        tm.assert_isinstance(result2, DataFrame)
        assert_frame_equal(result2, expected2)

        result = grouped['C'].sum()
        expected = grouped.sum().ix[:, ['A', 'C']]
        tm.assert_isinstance(result, DataFrame)
        assert_frame_equal(result, expected)

        result2 = grouped2['C'].sum()
        expected2 = grouped2.sum().ix[:, ['A', 'B', 'C']]
        tm.assert_isinstance(result2, DataFrame)
        assert_frame_equal(result2, expected2)

        # corner case
        self.assertRaises(Exception, grouped['C'].__getitem__,
                          'D')

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
        self.assertRaises(TypeError, self.ts.groupby,
                          lambda x: x.weekday(), as_index=False)

        self.assertRaises(ValueError, self.df.groupby,
                          lambda x: x.lower(), as_index=False, axis=1)

    def test_groupby_as_index_apply(self):
        # GH #4648 and #3417
        df = DataFrame({'item_id': ['b', 'b', 'a', 'c', 'a', 'b'],
                        'user_id': [1,2,1,1,3,1],
                        'time': range(6)})

        g_as = df.groupby('user_id', as_index=True)
        g_not_as = df.groupby('user_id', as_index=False)

        res_as = g_as.head(2).index
        exp_as = MultiIndex.from_tuples([(1, 0), (2, 1), (1, 2), (3, 4)])
        assert_index_equal(res_as, exp_as)

        res_not_as = g_not_as.head(2).index
        exp_not_as = Index([0, 1, 2, 4])
        assert_index_equal(res_not_as, exp_not_as)

        res_as_apply = g_as.apply(lambda x: x.head(2)).index
        res_not_as_apply = g_not_as.apply(lambda x: x.head(2)).index

        # apply doesn't maintain the original ordering
        exp_not_as_apply = Index([0, 2, 1, 4])
        exp_as_apply = MultiIndex.from_tuples([(1, 0), (1, 2), (2, 1), (3, 4)])

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

        empty_not_as = DataFrame(columns=df.columns)
        assert_frame_equal(empty_not_as, g_not_as.head(0))
        assert_frame_equal(empty_not_as, g_not_as.tail(0))
        assert_frame_equal(empty_not_as, g_not_as.head(-1))
        assert_frame_equal(empty_not_as, g_not_as.tail(-1))

        assert_frame_equal(df, g_not_as.head(7)) # contains all
        assert_frame_equal(df, g_not_as.tail(7))

        # as_index=True, yuck
        # prepend the A column as an index, in a roundabout way
        df_as = df.copy()
        df_as.index = df.set_index('A', append=True,
                                        drop=False).index.swaplevel(0, 1)

        assert_frame_equal(df_as.loc[[0, 2]], g_as.head(1))
        assert_frame_equal(df_as.loc[[1, 2]], g_as.tail(1))

        empty_as = DataFrame(index=df_as.index[:0], columns=df.columns)
        assert_frame_equal(empty_as, g_as.head(0))
        assert_frame_equal(empty_as, g_as.tail(0))
        assert_frame_equal(empty_as, g_as.head(-1))
        assert_frame_equal(empty_as, g_as.tail(-1))

        assert_frame_equal(df_as, g_as.head(7)) # contains all
        assert_frame_equal(df_as, g_as.tail(7))

    def test_groupby_multiple_key(self):
        df = tm.makeTimeDataFrame()
        grouped = df.groupby([lambda x: x.year,
                              lambda x: x.month,
                              lambda x: x.day])
        agged = grouped.sum()
        assert_almost_equal(df.values, agged.values)

        grouped = df.T.groupby([lambda x: x.year,
                                lambda x: x.month,
                                lambda x: x.day], axis=1)

        agged = grouped.agg(lambda x: x.sum(1))
        self.assert_(agged.index.equals(df.columns))
        assert_almost_equal(df.T.values, agged.values)

        agged = grouped.agg(lambda x: x.sum(1))
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
                                   lambda x: x.sum(1, numeric_only=False))

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
        self.assertEquals(agged.index.name, 'first')

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
        self.assertEquals(result.index.names, ('A', 'B', 'stat'))

        result2 = grouped.apply(desc2)
        self.assertEquals(result2.index.names, ('A', 'B', 'stat'))

        result3 = grouped.apply(desc3)
        self.assertEquals(result3.index.names, ('A', 'B', None))

    def test_nonsense_func(self):
        df = DataFrame([0])
        self.assertRaises(Exception, df.groupby, lambda x: x + 'foo')

    def test_cythonized_aggers(self):
        data = {'A': [0, 0, 0, 0, 1, 1, 1, 1, 1, 1., nan, nan],
                'B': ['A', 'B'] * 6,
                'C': np.random.randn(12)}
        df = DataFrame(data)
        df['C'][2:10:2] = nan

        def _testit(op):
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
            result = op(grouped)['C']
            assert_series_equal(result, exp)

        _testit(lambda x: x.sum())
        _testit(lambda x: x.mean())
        _testit(lambda x: x.prod())
        _testit(lambda x: x.min())
        _testit(lambda x: x.max())

    def test_max_min_non_numeric(self):
        # #2700
        aa = DataFrame({'nn':[11,11,22,22],'ii':[1,2,3,4],'ss':4*['mama']})

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

    def test_cython_agg_frame_columns(self):
        # #2113
        df = DataFrame({'x': [1, 2, 3], 'y': [3, 4, 5]})

        result = df.groupby(level=0, axis='columns').mean()
        result = df.groupby(level=0, axis='columns').mean()
        result = df.groupby(level=0, axis='columns').mean()
        _ = df.groupby(level=0, axis='columns').mean()

    def test_wrap_aggregated_output_multindex(self):
        df = self.mframe.T
        df['baz', 'two'] = 'peekaboo'

        keys = [np.array([0, 0, 1]), np.array([0, 0, 1])]
        agged = df.groupby(keys).agg(np.mean)
        tm.assert_isinstance(agged.columns, MultiIndex)

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

        self.assert_(result0.index.name == 'first')
        self.assert_(result1.index.name == 'second')

        assert_frame_equal(result0, expected0)
        assert_frame_equal(result1, expected1)
        self.assertEquals(result0.index.name, frame.index.names[0])
        self.assertEquals(result1.index.name, frame.index.names[1])

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
        ## GH4014 this used to raise ValueError since 'exp'>1 (in py2)
        df = DataFrame({'exp' : ['A']*3 + ['B']*3, 'var1' : lrange(6),}).set_index('exp')
        df.groupby(level='exp')
        self.assertRaises(ValueError, df.groupby, level='foo')

    def test_groupby_level_with_nas(self):
        index = MultiIndex(levels=[[1, 0], [0, 1, 2, 3]],
                           labels=[[1, 1, 1, 1, 0, 0, 0, 0],
                                   [0, 1, 2, 3, 0, 1, 2, 3]])

        # factorizing doesn't confuse things
        s = Series(np.arange(8.), index=index)
        result = s.groupby(level=0).sum()
        expected = Series([22., 6.], index=[1, 0])
        assert_series_equal(result, expected)

        index = MultiIndex(levels=[[1, 0], [0, 1, 2, 3]],
                           labels=[[1, 1, 1, 1, -1, 0, 0, 0],
                                   [0, 1, 2, 3, 0, 1, 2, 3]])

        # factorizing doesn't confuse things
        s = Series(np.arange(8.), index=index)
        result = s.groupby(level=0).sum()
        expected = Series([18., 6.], index=[1, 0])
        assert_series_equal(result, expected)

    def test_groupby_level_apply(self):
        frame = self.mframe

        result = frame.groupby(level=0).count()
        self.assert_(result.index.name == 'first')
        result = frame.groupby(level=1).count()
        self.assert_(result.index.name == 'second')

        result = frame['A'].groupby(level=0).count()
        self.assert_(result.index.name == 'first')

    def test_groupby_level_mapper(self):
        frame = self.mframe
        deleveled = frame.reset_index()

        mapper0 = {'foo': 0, 'bar': 0,
                   'baz': 1, 'qux': 1}
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

    def test_groupby_level_0_nonmulti(self):
        # #1313
        a = Series([1, 2, 3, 10, 4, 5, 20, 6], Index([1, 2, 3, 1,
                   4, 5, 2, 6], name='foo'))

        result = a.groupby(level=0).sum()
        self.assertEquals(result.index.name, a.index.name)

    def test_level_preserve_order(self):
        grouped = self.mframe.groupby(level=0)
        exp_labels = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 3])
        assert_almost_equal(grouped.grouper.labels[0], exp_labels)

    def test_grouping_labels(self):
        grouped = self.mframe.groupby(self.mframe.index.get_level_values(0))
        exp_labels = np.array([2, 2, 2, 0, 0, 1, 1, 3, 3, 3])
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
            return DataFrame({'value': piece,
                              'demeaned': piece - piece.mean(),
                              'logged': np.log(piece)})

        dr = bdate_range('1/1/2000', periods=100)
        ts = Series(np.random.randn(100), index=dr)

        grouped = ts.groupby(lambda x: x.month)
        result = grouped.apply(f)

        tm.assert_isinstance(result, DataFrame)
        self.assert_(result.index.equals(ts.index))

    def test_apply_series_yield_constant(self):
        result = self.df.groupby(['A', 'B'])['C'].apply(len)
        self.assertEquals(result.index.names[:2], ('A', 'B'))

    def test_apply_frame_to_series(self):
        grouped = self.df.groupby(['A', 'B'])
        result = grouped.apply(len)
        expected = grouped.count()['C']
        self.assert_(result.index.equals(expected.index))
        self.assert_(np.array_equal(result.values, expected.values))

    def test_apply_frame_concat_series(self):
        def trans(group):
            return group.groupby('B')['C'].sum().order()[:2]

        def trans2(group):
            grouped = group.groupby(df.reindex(group.index)['B'])
            return grouped.sum().order()[:2]

        df = DataFrame({'A': np.random.randint(0, 5, 1000),
                        'B': np.random.randint(0, 5, 1000),
                        'C': np.random.randn(1000)})

        result = df.groupby('A').apply(trans)
        exp = df.groupby('A')['C'].apply(trans2)
        assert_series_equal(result, exp)

    def test_apply_transform(self):
        grouped = self.ts.groupby(lambda x: x.month)
        result = grouped.apply(lambda x: x * 2)
        expected = grouped.transform(lambda x: x * 2)
        assert_series_equal(result, expected)

    def test_apply_multikey_corner(self):
        grouped = self.tsframe.groupby([lambda x: x.year,
                                        lambda x: x.month])

        def f(group):
            return group.sort('A')[-5:]

        result = grouped.apply(f)
        for key, group in grouped:
            assert_frame_equal(result.ix[key], f(group))

    def test_mutate_groups(self):

        # GH3380

        mydf = DataFrame({
                'cat1' : ['a'] * 8 + ['b'] * 6,
                'cat2' : ['c'] * 2 + ['d'] * 2 + ['e'] * 2 + ['f'] * 2 + ['c'] * 2 + ['d'] * 2 + ['e'] * 2,
                'cat3' : lmap(lambda x: 'g%s' % x, lrange(1,15)),
                'val' : np.random.randint(100, size=14),
                })

        def f_copy(x):
            x = x.copy()
            x['rank'] = x.val.rank(method='min')
            return x.groupby('cat2')['rank'].min()

        def f_no_copy(x):
            x['rank'] = x.val.rank(method='min')
            return x.groupby('cat2')['rank'].min()

        grpby_copy    = mydf.groupby('cat1').apply(f_copy)
        grpby_no_copy = mydf.groupby('cat1').apply(f_no_copy)
        assert_series_equal(grpby_copy,grpby_no_copy)

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
        grouped.apply(lambda x: x.sort('value'))

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
        tuples = list(zip(*[['bar', 'bar', 'baz', 'baz',
                        'foo', 'foo', 'qux', 'qux'],
                       ['one', 'two', 'one', 'two',
                        'one', 'two', 'one', 'two']]))
        index = MultiIndex.from_tuples(tuples)
        columns = MultiIndex.from_tuples([('A', 'cat'), ('B', 'dog'),
                                          ('B', 'cat'), ('A', 'dog')])
        df = DataFrame(np.random.randn(8, 4), index=index,
                       columns=columns)

        result = df.groupby(level=0).mean()
        self.assert_(result.columns.equals(columns))

        result = df.groupby(level=0, axis=1).mean()
        self.assert_(result.index.equals(df.index))

        result = df.groupby(level=0).agg(np.mean)
        self.assert_(result.columns.equals(columns))

        result = df.groupby(level=0).apply(lambda x: x.mean())
        self.assert_(result.columns.equals(columns))

        result = df.groupby(level=0, axis=1).agg(lambda x: x.mean(1))
        self.assert_(result.columns.equals(Index(['A', 'B'])))
        self.assert_(result.index.equals(df.index))

        # add a nuisance column
        sorted_columns, _ = columns.sortlevel(0)
        df['A', 'foo'] = 'bar'
        result = df.groupby(level=0).mean()
        self.assert_(result.columns.equals(df.columns[:-1]))

    def test_pass_args_kwargs(self):
        from pandas.compat.scipy import scoreatpercentile

        def f(x, q=None):
            return scoreatpercentile(x, q)
        g = lambda x: scoreatpercentile(x, 80)

        # Series
        ts_grouped = self.ts.groupby(lambda x: x.month)
        agg_result = ts_grouped.agg(scoreatpercentile, 80)
        apply_result = ts_grouped.apply(scoreatpercentile, 80)
        trans_result = ts_grouped.transform(scoreatpercentile, 80)

        agg_expected = ts_grouped.quantile(.8)
        trans_expected = ts_grouped.transform(g)

        assert_series_equal(apply_result, agg_expected)
        assert_series_equal(agg_result, agg_expected)
        assert_series_equal(trans_result, trans_expected)

        agg_result = ts_grouped.agg(f, q=80)
        apply_result = ts_grouped.apply(f, q=80)
        trans_result = ts_grouped.transform(f, q=80)
        assert_series_equal(agg_result, agg_expected)
        assert_series_equal(apply_result, agg_expected)
        assert_series_equal(trans_result, trans_expected)

        # DataFrame
        df_grouped = self.tsframe.groupby(lambda x: x.month)
        agg_result = df_grouped.agg(scoreatpercentile, 80)
        apply_result = df_grouped.apply(DataFrame.quantile, .8)
        expected = df_grouped.quantile(.8)
        assert_frame_equal(apply_result, expected)
        assert_frame_equal(agg_result, expected)

        agg_result = df_grouped.agg(f, q=80)
        apply_result = df_grouped.apply(DataFrame.quantile, q=.8)
        assert_frame_equal(agg_result, expected)
        assert_frame_equal(apply_result, expected)

    # def test_cython_na_bug(self):
    #     values = np.random.randn(10)
    #     shape = (5, 5)
    #     label_list = [np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2], dtype=np.int32),
    # np.array([1, 2, 3, 4, 0, 1, 2, 3, 3, 4], dtype=np.int32)]

    #     lib.group_aggregate(values, label_list, shape)

    def test_size(self):
        grouped = self.df.groupby(['A', 'B'])
        result = grouped.size()
        for key, group in grouped:
            self.assertEquals(result[key], len(group))

        grouped = self.df.groupby('A')
        result = grouped.size()
        for key, group in grouped:
            self.assertEquals(result[key], len(group))

        grouped = self.df.groupby('B')
        result = grouped.size()
        for key, group in grouped:
            self.assertEquals(result[key], len(group))

    def test_grouping_ndarray(self):
        grouped = self.df.groupby(self.df['A'].values)

        result = grouped.sum()
        expected = self.df.groupby('A').sum()
        assert_frame_equal(result, expected, check_names=False)  # Note: no names when grouping by value

    def test_apply_typecast_fail(self):
        df = DataFrame({'d': [1., 1., 1., 2., 2., 2.],
                        'c': np.tile(['a', 'b', 'c'], 2),
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
        index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1],
                                        [1, 2, 3, 1, 2, 3]])
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

        data = DataFrame({'id_field' : [100, 100, 200, 300], 'category' : ['a','b','c','c'], 'value' : [1,2,3,4]})

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
        assert_frame_equal(result,expected)

    def test_apply_use_categorical_name(self):
        from pandas import qcut
        cats = qcut(self.df.C, 4)

        def get_stats(group):
            return {'min': group.min(), 'max': group.max(),
                    'count': group.count(), 'mean': group.mean()}

        result = self.df.groupby(cats).D.apply(get_stats)
        self.assertEquals(result.index.names[0], 'C')

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
        index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1],
                                        [1, 2, 3, 1, 2, 3]])
        df = DataFrame({'d': [1., 1., 1., 2., 2., 2.],
                        'c': np.tile(['a', 'b', 'c'], 2),
                        'v': np.arange(1., 7.)}, index=index)

        def f(group):
            group['g'] = group['d'] * 2
            return group[:1]

        grouped = df.groupby('c')
        result = grouped.apply(f)

        self.assert_(result['d'].dtype == np.float64)

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
        self.assertEquals(result.index.name, 'A')
        self.assert_('A' in result2)

        result = self.df.groupby([self.df['A'], self.df['B']]).mean()
        result2 = self.df.groupby([self.df['A'], self.df['B']],
                                  as_index=False).mean()
        self.assertEquals(result.index.names, ('A', 'B'))
        self.assert_('A' in result2)
        self.assert_('B' in result2)

    def test_groupby_nonstring_columns(self):
        df = DataFrame([np.arange(10) for x in range(10)])
        grouped = df.groupby(0)
        result = grouped.mean()
        expected = df.groupby(df[0]).mean()
        del expected[0]
        assert_frame_equal(result, expected)

    def test_cython_grouper_series_bug_noncontig(self):
        arr = np.empty((100, 100))
        arr.fill(np.nan)
        obj = Series(arr[:, 0], index=lrange(100))
        inds = np.tile(lrange(10), 10)

        result = obj.groupby(inds).agg(Series.median)
        self.assert_(result.isnull().all())

    def test_series_grouper_noncontig_index(self):
        index = Index([tm.rands(10) for _ in range(100)])

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
            assert(len(x.base) > 0)
            return Decimal(str(x.mean()))

        grouped = s.groupby(labels)

        result = grouped.agg(convert_fast)
        self.assert_(result.dtype == np.object_)
        tm.assert_isinstance(result[0], Decimal)

        result = grouped.agg(convert_force_pure)
        self.assert_(result.dtype == np.object_)
        tm.assert_isinstance(result[0], Decimal)

    def test_apply_with_mixed_dtype(self):
        # GH3480, apply with mixed dtype on axis=1 breaks in 0.11
        df = DataFrame({'foo1' : ['one', 'two', 'two', 'three', 'one', 'two'],
                        'foo2' : np.random.randn(6)})
        result = df.apply(lambda x: x, axis=1)
        assert_series_equal(df.get_dtype_counts(), result.get_dtype_counts())


        # GH 3610 incorrect dtype conversion with as_index=False
        df = DataFrame({"c1" : [1,2,6,6,8]})
        df["c2"] = df.c1/2.0
        result1 = df.groupby("c2").mean().reset_index().c2
        result2 = df.groupby("c2", as_index=False).mean().c2
        assert_series_equal(result1,result2)

    def test_groupby_aggregation_mixed_dtype(self):

        # GH 6212
        expected = DataFrame({
            'v1': [5,5,7,np.nan,3,3,4,1],
            'v2': [55,55,77,np.nan,33,33,44,11]},
            index=MultiIndex.from_tuples([(1,95),(1,99),(2,95),(2,99),('big','damp'),
                                          ('blue','dry'),('red','red'),('red','wet')],
                                         names=['by1','by2']))

        df = DataFrame({
            'v1': [1,3,5,7,8,3,5,np.nan,4,5,7,9],
            'v2': [11,33,55,77,88,33,55,np.nan,44,55,77,99],
            'by1': ["red", "blue", 1, 2, np.nan, "big", 1, 2, "red", 1, np.nan, 12],
            'by2': ["wet", "dry", 99, 95, np.nan, "damp", 95, 99, "red", 99, np.nan,
                    np.nan]
            })

        g = df.groupby(['by1','by2'])
        result = g[['v1','v2']].mean()
        assert_frame_equal(result,expected)

    def test_groupby_list_infer_array_like(self):
        result = self.df.groupby(list(self.df['A'])).mean()
        expected = self.df.groupby(self.df['A']).mean()
        assert_frame_equal(result, expected, check_names=False)

        self.assertRaises(Exception, self.df.groupby, list(self.df['A'][:-1]))

        # pathological case of ambiguity
        df = DataFrame({'foo': [0, 1], 'bar': [3, 4],
                        'val': np.random.randn(2)})

        result = df.groupby(['foo', 'bar']).mean()
        expected = df.groupby([df['foo'], df['bar']]).mean()[['val']]

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

        self.assert_(np.array_equal(agged.items, [0, 1]))

        grouped = self.panel.groupby(lambda x: x.month, axis='major')
        agged = grouped.mean()

        self.assert_(np.array_equal(agged.major_axis, [1, 2]))

        grouped = self.panel.groupby({'A': 0, 'B': 0, 'C': 1, 'D': 1},
                                     axis='minor')
        agged = grouped.mean()
        self.assert_(np.array_equal(agged.minor_axis, [0, 1]))

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
        res_values = np.array([[0., 1.], [0., 1.]])
        self.assert_(np.array_equal(tmp.columns, ['zeros', 'ones']))
        self.assert_(np.array_equal(tmp.values, res_values))

    def test_int32_overflow(self):
        B = np.concatenate((np.arange(10000), np.arange(10000),
                            np.arange(5000)))
        A = np.arange(25000)
        df = DataFrame({'A': A, 'B': B,
                        'C': A, 'D': B,
                        'E': np.random.randn(25000)})

        left = df.groupby(['A', 'B', 'C', 'D']).sum()
        right = df.groupby(['D', 'C', 'B', 'A']).sum()
        self.assert_(len(left) == len(right))

    def test_int64_overflow(self):
        B = np.concatenate((np.arange(1000), np.arange(1000),
                            np.arange(500)))
        A = np.arange(2500)
        df = DataFrame({'A': A, 'B': B,
                        'C': A, 'D': B,
                        'E': A, 'F': B,
                        'G': A, 'H': B,
                        'values': np.random.randn(2500)})

        lg = df.groupby(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
        rg = df.groupby(['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A'])

        left = lg.sum()['values']
        right = rg.sum()['values']

        exp_index, _ = left.index.sortlevel(0)
        self.assert_(left.index.equals(exp_index))

        exp_index, _ = right.index.sortlevel(0)
        self.assert_(right.index.equals(exp_index))

        tups = list(map(tuple, df[['A', 'B', 'C', 'D',
                              'E', 'F', 'G', 'H']].values))
        tups = com._asarray_tuplesafe(tups)
        expected = df.groupby(tups).sum()['values']

        for k, v in compat.iteritems(expected):
            self.assert_(left[k] == right[k[::-1]] == v)
        self.assert_(len(left) == len(right))

    def test_groupby_sort_multi(self):
        df = DataFrame({'a': ['foo', 'bar', 'baz'],
                        'b': [3, 2, 1],
                        'c': [0, 1, 2],
                        'd': np.random.randn(3)})

        tups = lmap(tuple, df[['a', 'b', 'c']].values)
        tups = com._asarray_tuplesafe(tups)
        result = df.groupby(['a', 'b', 'c'], sort=True).sum()
        self.assert_(np.array_equal(result.index.values,
                                    tups[[1, 2, 0]]))

        tups = lmap(tuple, df[['c', 'a', 'b']].values)
        tups = com._asarray_tuplesafe(tups)
        result = df.groupby(['c', 'a', 'b'], sort=True).sum()
        self.assert_(np.array_equal(result.index.values, tups))

        tups = lmap(tuple, df[['b', 'c', 'a']].values)
        tups = com._asarray_tuplesafe(tups)
        result = df.groupby(['b', 'c', 'a'], sort=True).sum()
        self.assert_(np.array_equal(result.index.values,
                                    tups[[2, 1, 0]]))

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
        lev1 = np.array([rands(10) for _ in range(100)], dtype=object)
        lev2 = np.array([rands(10) for _ in range(130)], dtype=object)
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

    def test_dont_clobber_name_column(self):
        df = DataFrame({'key': ['a', 'a', 'a', 'b', 'b', 'b'],
                        'name': ['foo', 'bar', 'baz'] * 2})

        result = df.groupby('key').apply(lambda x: x)
        assert_frame_equal(result, df)

    def test_skip_group_keys(self):
        from pandas import concat

        tsf = tm.makeTimeDataFrame()

        grouped = tsf.groupby(lambda x: x.month, group_keys=False)
        result = grouped.apply(lambda x: x.sort_index(by='A')[:3])

        pieces = []
        for key, group in grouped:
            pieces.append(group.sort_index(by='A')[:3])

        expected = concat(pieces)
        assert_frame_equal(result, expected)

        grouped = tsf['A'].groupby(lambda x: x.month, group_keys=False)
        result = grouped.apply(lambda x: x.order()[:3])

        pieces = []
        for key, group in grouped:
            pieces.append(group.order()[:3])

        expected = concat(pieces)
        assert_series_equal(result, expected)

    def test_no_nonsense_name(self):
        # GH #995
        s = self.frame['C'].copy()
        s.name = None

        result = s.groupby(self.frame['A']).agg(np.sum)
        self.assert_(result.name is None)

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
        self.assert_(result['fl'].dtype == np.float64)

    def test_handle_dict_return_value(self):
        def f(group):
            return {'min': group.min(), 'max': group.max()}

        def g(group):
            return Series({'min': group.min(), 'max': group.max()})

        result = self.df.groupby('A')['C'].apply(f)
        expected = self.df.groupby('A')['C'].apply(g)

        tm.assert_isinstance(result, Series)
        assert_series_equal(result, expected)

    def test_getitem_list_of_columns(self):
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'three',
                              'two', 'two', 'one', 'three'],
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

    def test_agg_multiple_functions_maintain_order(self):
        # GH #610
        funcs = [('mean', np.mean), ('max', np.max), ('min', np.min)]
        result = self.df.groupby('A')['C'].agg(funcs)
        exp_cols = ['mean', 'max', 'min']

        self.assert_(np.array_equal(result.columns, exp_cols))

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
        d = OrderedDict([['C', np.mean],
                         ['D', OrderedDict([['foo', np.mean],
                                            ['bar', np.std]])]])
        result = grouped.aggregate(d)

        d = OrderedDict([['C', [np.mean]], ['D', [foo, bar]]])
        expected = grouped.aggregate(d)

        assert_frame_equal(result, expected)

    def test_multi_function_flexible_mix(self):
        # GH #1268
        grouped = self.df.groupby('A')

        d = OrderedDict([['C', OrderedDict([['foo', 'mean'],
                                            [
                                            'bar', 'std']])],
                         ['D', 'sum']])
        result = grouped.aggregate(d)
        d2 = OrderedDict([['C', OrderedDict([['foo', 'mean'],
                                             [
                                             'bar', 'std']])],
                          ['D', ['sum']]])
        result2 = grouped.aggregate(d2)

        d3 = OrderedDict([['C', OrderedDict([['foo', 'mean'],
                                             [
                                             'bar', 'std']])],
                          ['D', {'sum': 'sum'}]])
        expected = grouped.aggregate(d3)

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

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
        self.assert_(result.index.name is None)

        result = self.df.groupby([self.df['A'].values,
                                  self.df['B'].values]).sum()
        self.assert_(result.index.names == (None, None))

    def test_groupby_categorical(self):
        levels = ['foo', 'bar', 'baz', 'qux']
        labels = np.random.randint(0, 4, size=100)

        cats = Categorical(labels, levels, name='myfactor')

        data = DataFrame(np.random.randn(100, 4))

        result = data.groupby(cats).mean()

        expected = data.groupby(np.asarray(cats)).mean()
        expected = expected.reindex(levels)
        expected.index.name = 'myfactor'

        assert_frame_equal(result, expected)
        self.assert_(result.index.name == cats.name)

        grouped = data.groupby(cats)
        desc_result = grouped.describe()

        idx = cats.labels.argsort()
        ord_labels = np.asarray(cats).take(idx)
        ord_data = data.take(idx)
        expected = ord_data.groupby(ord_labels, sort=False).describe()
        expected.index.names = ['myfactor', None]
        assert_frame_equal(desc_result, expected)

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
        tm.assert_isinstance(list(groups.keys())[0], datetime)

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
        self.assert_(np.array_equal(grouped['A'], [0]))

        grouped = df.groupby([('to filter', '')]).groups
        self.assert_(np.array_equal(grouped['A'], [0]))

        df = DataFrame([[long(1), 'A'], [long(2), 'B']], columns=midx)

        expected = df.groupby('to filter').groups
        result = df.groupby([('to filter', '')]).groups
        self.assertEquals(result, expected)

        df = DataFrame([[long(1), 'A'], [long(2), 'A']], columns=midx)

        expected = df.groupby('to filter').groups
        result = df.groupby([('to filter', '')]).groups
        self.assertEquals(result, expected)

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

    def test_groupby_categorical_no_compress(self):
        data = Series(np.random.randn(9))

        labels = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        cats = Categorical(labels, [0, 1, 2])

        result = data.groupby(cats).mean()
        exp = data.groupby(labels).mean()
        assert_series_equal(result, exp)

        labels = np.array([0, 0, 0, 1, 1, 1, 3, 3, 3])
        cats = Categorical(labels, [0, 1, 2, 3])

        result = data.groupby(cats).mean()
        exp = data.groupby(labels).mean().reindex(cats.levels)
        assert_series_equal(result, exp)

    def test_groupby_first_datetime64(self):
        df = DataFrame([(1, 1351036800000000000), (2, 1351036800000000000)])
        df[1] = df[1].view('M8[ns]')

        self.assert_(issubclass(df[1].dtype.type, np.datetime64))

        result = df.groupby(level=0).first()
        got_dt = result[1].dtype
        self.assert_(issubclass(got_dt.type, np.datetime64))

        result = df[1].groupby(level=0).first()
        got_dt = result.dtype
        self.assert_(issubclass(got_dt.type, np.datetime64))

    def test_groupby_max_datetime64(self):
        # GH 5869
        # datetimelike dtype conversion from int
        df = DataFrame(dict(A = Timestamp('20130101'), B = np.arange(5)))
        expected = df.groupby('A')['A'].apply(lambda x: x.max())
        result = df.groupby('A')['A'].max()
        assert_series_equal(result,expected)

    def test_groupby_categorical_unequal_len(self):
        import pandas as pd
        #GH3011
        series = Series([np.nan, np.nan, 1, 1, 2, 2, 3, 3, 4, 4])
        bins =  pd.cut(series.dropna(), 4)

        # len(bins) != len(series) here
        self.assertRaises(AssertionError,lambda : series.groupby(bins).mean())

    def test_gb_apply_list_of_unequal_len_arrays(self):

        # GH1738
        df = DataFrame({'group1': ['a','a','a','b','b','b','a','a','a','b','b','b'],
                               'group2': ['c','c','d','d','d','e','c','c','d','d','d','e'],
                               'weight': [1.1,2,3,4,5,6,2,4,6,8,1,2],
                               'value': [7.1,8,9,10,11,12,8,7,6,5,4,3]
        })
        df = df.set_index(['group1', 'group2'])
        df_grouped = df.groupby(level=['group1','group2'], sort=True)

        def noddy(value, weight):
            out = np.array( value * weight ).repeat(3)
            return out

        # the kernel function returns arrays of unequal length
        # pandas sniffs the first one, sees it's an array and not
        # a list, and assumed the rest are of equal length
        # and so tries a vstack

        # don't die
        no_toes = df_grouped.apply(lambda x: noddy(x.value, x.weight ))

    def test_groupby_with_empty(self):
        import pandas as pd
        index = pd.DatetimeIndex(())
        data = ()
        series = pd.Series(data, index)
        grouper = pd.tseries.resample.TimeGrouper('D')
        grouped = series.groupby(grouper)
        assert next(iter(grouped), None) is None

    def test_cumcount(self):
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'])
        g = df.groupby('A')
        sg = g.A

        expected = Series([0, 1, 2, 0, 3])

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())

    def test_cumcount_empty(self):
        ge = DataFrame().groupby()
        se = Series().groupby()

        e = Series(dtype='int64')  # edge case, as this is usually considered float

        assert_series_equal(e, ge.cumcount())
        assert_series_equal(e, se.cumcount())

    def test_cumcount_dupe_index(self):
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'], index=[0] * 5)
        g = df.groupby('A')
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=[0] * 5)

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())

    def test_cumcount_mi(self):
        mi = MultiIndex.from_tuples([[0, 1], [1, 2], [2, 2], [2, 2], [1, 0]])
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'], index=mi)
        g = df.groupby('A')
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=mi)

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())

    def test_cumcount_groupby_not_col(self):
        df = DataFrame([['a'], ['a'], ['a'], ['b'], ['a']], columns=['A'], index=[0] * 5)
        g = df.groupby([0, 0, 0, 1, 0])
        sg = g.A

        expected = Series([0, 1, 2, 0, 3], index=[0] * 5)

        assert_series_equal(expected, g.cumcount())
        assert_series_equal(expected, sg.cumcount())


    def test_filter_series(self):
        import pandas as pd
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
        import pandas as pd
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
        import pandas as pd
        df = pd.DataFrame({'A': [1, 12, 12, 1], 'B': [1, 1, 1, 1]})
        grouper = df['A'].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        expected = pd.DataFrame({'A': [12, 12], 'B': [1, 1]}, index=[1, 2])
        assert_frame_equal(
            grouped.filter(lambda x: x['A'].sum() - x['B'].sum() > 10), expected)

    def test_filter_mixed_df(self):
        import pandas as pd
        df = pd.DataFrame({'A': [1, 12, 12, 1], 'B': 'a b c d'.split()})
        grouper = df['A'].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        expected = pd.DataFrame({'A': [12, 12], 'B': ['b', 'c']},
                                index=[1, 2])
        assert_frame_equal(
            grouped.filter(lambda x: x['A'].sum() > 10), expected)

    def test_filter_out_all_groups(self):
        import pandas as pd
        s = pd.Series([1, 3, 20, 5, 22, 24, 7])
        grouper = s.apply(lambda x: x % 2)
        grouped = s.groupby(grouper)
        assert_series_equal(
            grouped.filter(lambda x: x.mean() > 1000), s[[]])
        df = pd.DataFrame({'A': [1, 12, 12, 1], 'B': 'a b c d'.split()})
        grouper = df['A'].apply(lambda x: x % 2)
        grouped = df.groupby(grouper)
        assert_frame_equal(
            grouped.filter(lambda x: x['A'].sum() > 1000), df.ix[[]])

    def test_filter_out_no_groups(self):
        import pandas as pd
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

    def test_filter_condition_raises(self):
        import pandas as pd
        def raise_if_sum_is_zero(x):
            if x.sum() == 0:
                raise ValueError
            else:
                return x.sum() > 0
        s = pd.Series([-1,0,1,2])
        grouper = s.apply(lambda x: x % 2)
        grouped = s.groupby(grouper)
        self.assertRaises(TypeError,
                          lambda: grouped.filter(raise_if_sum_is_zero))

    def test_filter_bad_shapes(self):
        df = DataFrame({'A': np.arange(8), 'B': list('aabbbbcc'), 'C': np.arange(8)})
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
        df = DataFrame({'A': np.arange(8), 'B': list('aabbbbcc'), 'C': np.arange(8)})
        s = df['B']
        g_df = df.groupby(df['B'])
        g_s = s.groupby(s)

        f = lambda x: np.nan
        assert_frame_equal(g_df.filter(f), df.loc[[]])
        assert_series_equal(g_s.filter(f), s[[]])

    def test_filter_against_workaround(self):
        np.random.seed(0)
        # Series of ints
        s = Series(np.random.randint(0,100,1000))
        grouper = s.apply(lambda x: np.round(x, -1))
        grouped = s.groupby(grouper)
        f = lambda x: x.mean() > 10
        old_way = s[grouped.transform(f).astype('bool')]
        new_way = grouped.filter(f)
        assert_series_equal(new_way.order(), old_way.order())

        # Series of floats
        s = 100*Series(np.random.random(1000))
        grouper = s.apply(lambda x: np.round(x, -1))
        grouped = s.groupby(grouper)
        f = lambda x: x.mean() > 10
        old_way = s[grouped.transform(f).astype('bool')]
        new_way = grouped.filter(f)
        assert_series_equal(new_way.order(), old_way.order())

        # Set up DataFrame of ints, floats, strings.
        from string import ascii_lowercase
        letters = np.array(list(ascii_lowercase))
        N = 1000
        random_letters = letters.take(np.random.randint(0, 26, N))
        df = DataFrame({'ints': Series(np.random.randint(0, 100, N)),
                        'floats': N/10*Series(np.random.random(N)),
                        'letters': Series(random_letters)})

        # Group by ints; filter on floats.
        grouped = df.groupby('ints')
        old_way = df[grouped.floats.\
            transform(lambda x: x.mean() > N/20).astype('bool')]
        new_way = grouped.filter(lambda x: x['floats'].mean() > N/20)
        assert_frame_equal(new_way, old_way)

        # Group by floats (rounded); filter on strings.
        grouper = df.floats.apply(lambda x: np.round(x, -1))
        grouped = df.groupby(grouper)
        old_way = df[grouped.letters.\
            transform(lambda x: len(x) < N/10).astype('bool')]
        new_way = grouped.filter(
            lambda x: len(x.letters) < N/10)
        assert_frame_equal(new_way, old_way)

        # Group by strings; filter on ints.
        grouped = df.groupby('letters')
        old_way = df[grouped.ints.\
            transform(lambda x: x.mean() > N/20).astype('bool')]
        new_way = grouped.filter(lambda x: x['ints'].mean() > N/20)
        assert_frame_equal(new_way, old_way)

    def test_filter_using_len(self):
        # BUG GH4447
        df = DataFrame({'A': np.arange(8), 'B': list('aabbbbcc'), 'C': np.arange(8)})
        grouped = df.groupby('B')
        actual = grouped.filter(lambda x: len(x) > 2)
        expected = DataFrame({'A': np.arange(2, 6), 'B': list('bbbb'), 'C': np.arange(2, 6)}, index=np.arange(2, 6))
        assert_frame_equal(actual, expected)

        actual = grouped.filter(lambda x: len(x) > 4)
        expected = df.ix[[]]
        assert_frame_equal(actual, expected)

        # Series have always worked properly, but we'll test anyway.
        s = df['B']
        grouped = s.groupby(s)
        actual = grouped.filter(lambda x: len(x) > 2)
        expected = Series(4*['b'], index=np.arange(2, 6))
        assert_series_equal(actual, expected)

        actual = grouped.filter(lambda x: len(x) > 4)
        expected = s[[]]
        assert_series_equal(actual, expected)

    def test_filter_maintains_ordering(self):
        # Simple case: index is sequential. #4621
        df = DataFrame({'pid' : [1,1,1,2,2,3,3,3],
                        'tag' : [23,45,62,24,45,34,25,62]})
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

    def test_filter_and_transform_with_non_unique_int_index(self):
        # GH4620
        index = [1, 1, 1, 2, 1, 1, 0, 1]
        df = DataFrame({'pid' : [1,1,1,2,2,3,3,3],
                       'tag' : [23,45,62,24,45,34,25,62]}, index=index)
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
        expected = Series([NA,1,1,NA,2,NA,NA,3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index)
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_multiple_non_unique_int_index(self):
        # GH4620
        index = [1, 1, 1, 2, 0, 0, 0, 1]
        df = DataFrame({'pid' : [1,1,1,2,2,3,3,3],
                       'tag' : [23,45,62,24,45,34,25,62]}, index=index)
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
        expected = Series([NA,1,1,NA,2,NA,NA,3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index)
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_non_unique_float_index(self):
        # GH4620
        index = np.array([1, 1, 1, 2, 1, 1, 0, 1], dtype=float)
        df = DataFrame({'pid' : [1,1,1,2,2,3,3,3],
                       'tag' : [23,45,62,24,45,34,25,62]}, index=index)
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
        expected = Series([NA,1,1,NA,2,NA,NA,3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index)
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_non_unique_float_index(self):
        # GH4620
        index = np.array([1, 1, 1, 2, 0, 0, 0, 1], dtype=float)
        df = DataFrame({'pid' : [1,1,1,2,2,3,3,3],
                       'tag' : [23,45,62,24,45,34,25,62]}, index=index)
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
        expected = Series([NA,1,1,NA,2,NA,NA,3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index)
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
        df = DataFrame({'pid' : [1,1,1,2,2,3,3,3],
                       'tag' : [23,45,62,24,45,34,25,62]}, index=index)
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
        expected = Series([NA,1,1,NA,2,NA,NA,3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index)
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

    def test_filter_and_transform_with_non_unique_string_index(self):
        # GH4620
        index = list('bbbcbbab')
        df = DataFrame({'pid' : [1,1,1,2,2,3,3,3],
                       'tag' : [23,45,62,24,45,34,25,62]}, index=index)
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
        expected = Series([NA,1,1,NA,2,NA,NA,3], index, name='pid')
        # ^ made manually because this can get confusing!
        assert_series_equal(actual, expected)

        # Transform Series
        actual = grouped_ser.transform(len)
        expected = Series([1, 2, 2, 1, 2, 1, 1, 2], index)
        assert_series_equal(actual, expected)

        # Transform (a column from) DataFrameGroupBy
        actual = grouped_df.pid.transform(len)
        assert_series_equal(actual, expected)

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

    def test_groupby_whitelist(self):
        from string import ascii_lowercase
        letters = np.array(list(ascii_lowercase))
        N = 10
        random_letters = letters.take(np.random.randint(0, 26, N))
        df = DataFrame({'floats': N / 10 * Series(np.random.random(N)),
                        'letters': Series(random_letters)})
        s = df.floats

        df_whitelist = frozenset([
            'last', 'first',
            'mean', 'sum', 'min', 'max',
            'head', 'tail',
            'cumsum', 'cumprod', 'cummin', 'cummax', 'cumcount',
            'resample',
            'describe',
            'rank', 'quantile', 'count',
            'fillna',
            'mad',
            'any', 'all',
            'irow', 'take',
            'idxmax', 'idxmin',
            'shift', 'tshift',
            'ffill', 'bfill',
            'pct_change', 'skew',
            'plot', 'boxplot', 'hist',
            'median', 'dtypes',
            'corrwith', 'corr', 'cov',
            'diff',
        ])
        s_whitelist = frozenset([
            'last', 'first',
            'mean', 'sum', 'min', 'max',
            'head', 'tail',
            'cumsum', 'cumprod', 'cummin', 'cummax', 'cumcount',
            'resample',
            'describe',
            'rank', 'quantile', 'count',
            'fillna',
            'mad',
            'any', 'all',
            'irow', 'take',
            'idxmax', 'idxmin',
            'shift', 'tshift',
            'ffill', 'bfill',
            'pct_change', 'skew',
            'plot', 'hist',
            'median', 'dtype',
            'corr', 'cov',
            'value_counts',
            'diff',
            'unique', 'nunique',
        ])

        for obj, whitelist in zip((df, s),
                                  (df_whitelist, s_whitelist)):
            gb = obj.groupby(df.letters)
            self.assertEqual(whitelist, gb._apply_whitelist)
            for m in whitelist:
                getattr(gb, m)

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

    def test_series_groupby_plotting_nominally_works(self):
        try:
            import matplotlib as mpl
            mpl.use('Agg')
        except ImportError:
            raise nose.SkipTest("matplotlib not installed")
        n = 10
        weight = Series(np.random.normal(166, 20, size=n))
        height = Series(np.random.normal(60, 10, size=n))
        gender = tm.choice(['male', 'female'], size=n)

        weight.groupby(gender).plot()
        tm.close()
        height.groupby(gender).hist()
        tm.close()

    @slow
    def test_frame_groupby_plot_boxplot(self):
        try:
            import matplotlib.pyplot as plt
            import matplotlib as mpl
            mpl.use('Agg')
        except ImportError:
            raise nose.SkipTest("matplotlib not installed")
        tm.close()

        n = 10
        weight = Series(np.random.normal(166, 20, size=n))
        height = Series(np.random.normal(60, 10, size=n))
        gender = tm.choice(['male', 'female'], size=n)
        df = DataFrame({'height': height, 'weight': weight, 'gender': gender})
        gb = df.groupby('gender')

        res = gb.plot()
        self.assertEqual(len(plt.get_fignums()), 2)
        self.assertEqual(len(res), 2)
        tm.close()

        res = gb.boxplot()
        self.assertEqual(len(plt.get_fignums()), 1)
        self.assertEqual(len(res), 2)
        tm.close()

        with tm.assertRaisesRegexp(TypeError, '.*str.+float'):
            gb.hist()

    @slow
    def test_frame_groupby_hist(self):
        try:
            import matplotlib.pyplot as plt
            import matplotlib as mpl
            mpl.use('Agg')
        except ImportError:
            raise nose.SkipTest("matplotlib not installed")
        tm.close()

        n = 10
        weight = Series(np.random.normal(166, 20, size=n))
        height = Series(np.random.normal(60, 10, size=n))
        gender_int = tm.choice([0, 1], size=n)
        df_int = DataFrame({'height': height, 'weight': weight,
                            'gender': gender_int})
        gb = df_int.groupby('gender')
        axes = gb.hist()
        self.assertEqual(len(axes), 2)
        self.assertEqual(len(plt.get_fignums()), 2)
        tm.close()

    def test_tab_completion(self):
        grp = self.mframe.groupby(level='second')
        results = set([v for v in dir(grp) if not v.startswith('_')])
        expected = set(['A','B','C',
            'agg','aggregate','apply','boxplot','filter','first','get_group',
            'groups','hist','indices','last','max','mean','median',
            'min','name','ngroups','nth','ohlc','plot', 'prod',
            'size','std','sum','transform','var', 'count', 'head', 'describe',
            'cummax', 'quantile', 'rank', 'cumprod', 'tail',
            'resample', 'cummin', 'fillna', 'cumsum', 'cumcount',
            'all', 'shift', 'skew', 'bfill', 'irow', 'ffill',
            'take', 'tshift', 'pct_change', 'any', 'mad', 'corr', 'corrwith',
            'cov', 'dtypes', 'diff', 'idxmax', 'idxmin'
        ])
        self.assertEqual(results, expected)

def assert_fp_equal(a, b):
    assert (np.abs(a - b) < 1e-12).all()


def _check_groupby(df, result, keys, field, f=lambda x: x.sum()):
    tups = lmap(tuple, df[keys].values)
    tups = com._asarray_tuplesafe(tups)
    expected = f(df.groupby(tups)[field])
    for k, v in compat.iteritems(expected):
        assert(result[k] == v)


def test_decons():
    from pandas.core.groupby import decons_group_index, get_group_index

    def testit(label_list, shape):
        group_index = get_group_index(label_list, shape)
        label_list2 = decons_group_index(group_index, shape)

        for a, b in zip(label_list, label_list2):
            assert(np.array_equal(a, b))

    shape = (4, 5, 6)
    label_list = [np.tile([0, 1, 2, 3, 0, 1, 2, 3], 100),
                  np.tile([0, 2, 4, 3, 0, 1, 2, 3], 100),
                  np.tile([5, 1, 0, 2, 3, 0, 5, 4], 100)]
    testit(label_list, shape)

    shape = (10000, 10000)
    label_list = [np.tile(np.arange(10000), 5),
                  np.tile(np.arange(10000), 5)]
    testit(label_list, shape)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure',
                         '-s'], exit=False)
