import nose
import unittest

from datetime import datetime
from numpy import nan

from pandas.core.daterange import DateRange
from pandas.core.index import Index, MultiIndex
from pandas.core.common import rands
from pandas.core.frame import DataFrame
from pandas.core.groupby import GroupByError
from pandas.core.series import Series
from pandas.util.testing import (assert_panel_equal, assert_frame_equal,
                                 assert_series_equal, assert_almost_equal)
from pandas.core.panel import Panel
from pandas.tools.merge import concat
from collections import defaultdict
import pandas.core.datetools as dt
import numpy as np

import pandas.util.testing as tm

def commonSetUp(self):
    self.dateRange = DateRange('1/1/2005', periods=250, offset=dt.bday)
    self.stringIndex = Index([rands(8).upper() for x in xrange(250)])

    self.groupId = Series([x[0] for x in self.stringIndex],
                              index=self.stringIndex)
    self.groupDict = dict((k, v) for k, v in self.groupId.iteritems())

    self.columnIndex = Index(['A', 'B', 'C', 'D', 'E'])

    randMat = np.random.randn(250, 5)
    self.stringMatrix = DataFrame(randMat, columns=self.columnIndex,
                                  index=self.stringIndex)

    self.timeMatrix = DataFrame(randMat, columns=self.columnIndex,
                                index=self.dateRange)

class TestGroupBy(unittest.TestCase):

    def setUp(self):
        self.ts = tm.makeTimeSeries()

        self.seriesd = tm.getSeriesData()
        self.tsd = tm.getTimeSeriesData()
        self.frame = DataFrame(self.seriesd)
        self.tsframe = DataFrame(self.tsd)

        self.df = DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                                    'foo', 'bar', 'foo', 'foo'],
                             'B' : ['one', 'one', 'two', 'three',
                                    'two', 'two', 'one', 'three'],
                             'C' : np.random.randn(8),
                             'D' : np.random.randn(8)})

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        self.mframe = DataFrame(np.random.randn(10, 3), index=index,
                                columns=['A', 'B', 'C'])

        self.three_group = DataFrame({'A' : ['foo', 'foo', 'foo', 'foo',
                                             'bar', 'bar', 'bar', 'bar',
                                             'foo', 'foo', 'foo'],
                                      'B' : ['one', 'one', 'one', 'two',
                                             'one', 'one', 'one', 'two',
                                             'two', 'two', 'one'],
                                      'C' : ['dull', 'dull', 'shiny', 'dull',
                                             'dull', 'shiny', 'shiny', 'dull',
                                             'shiny', 'shiny', 'shiny'],
                                      'D' : np.random.randn(11),
                                      'E' : np.random.randn(11),
                                      'F' : np.random.randn(11)})

    def test_basic(self):
        data = Series(np.arange(9) // 3, index=np.arange(9))

        index = np.arange(9)
        np.random.shuffle(index)
        data = data.reindex(index)

        grouped = data.groupby(lambda x: x // 3)

        for k, v in grouped:
            self.assertEqual(len(v), 3)

        agged = grouped.aggregate(np.mean)
        self.assertEqual(agged[1], 1)

        assert_series_equal(agged, grouped.agg(np.mean)) # shorthand
        assert_series_equal(agged, grouped.mean())

        # Cython only returning floating point for now...
        assert_series_equal(grouped.agg(np.sum).astype(float),
                            grouped.sum())

        transformed = grouped.transform(lambda x: x * x.sum())
        self.assertEqual(transformed[7], 12)

        value_grouped = data.groupby(data)
        assert_series_equal(value_grouped.aggregate(np.mean), agged)

        # complex agg
        agged = grouped.aggregate([np.mean, np.std])
        agged = grouped.aggregate({'one' : np.mean,
                                   'two' : np.std})

        group_constants = {
            0 : 10,
            1 : 20,
            2 : 30
        }
        agged = grouped.agg(lambda x: group_constants[x.name] + x.mean())
        self.assertEqual(agged[1], 21)

        # corner cases
        self.assertRaises(Exception, grouped.aggregate, lambda x: x * 2)

    def test_groupby_dict_mapping(self):
        # GH #679
        from pandas import Series
        s = Series({'T1': 5})
        result = s.groupby({'T1': 'T2'}).agg(sum)
        expected = s.groupby(['T2']).agg(sum)
        assert_series_equal(result, expected)

        s = Series([1., 2., 3., 4.], index=list('abcd'))
        mapping = {'a' : 0, 'b' : 0, 'c' : 1, 'd' : 1}

        result = s.groupby(mapping).mean()
        result2 = s.groupby(mapping).agg(np.mean)
        expected = s.groupby([0, 0, 1, 1]).mean()
        expected2 = s.groupby([0, 0, 1, 1]).mean()
        assert_series_equal(result, expected)
        assert_series_equal(result, result2)
        assert_series_equal(result, expected2)

    def test_groupby_nonobject_dtype(self):
        key = self.mframe.index.labels[0]
        grouped = self.mframe.groupby(key)
        result = grouped.sum()

        expected = self.mframe.groupby(key.astype('O')).sum()
        assert_frame_equal(result, expected)

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

        df1['weights'] = df1['value']/df1['value'].sum()
        gb1 = df1.groupby('date').aggregate(np.sum)

        df2['weights'] = df1['value']/df1['value'].sum()
        gb2 = df2.groupby('date').aggregate(np.sum)

        assert(len(gb1) == len(gb2))

    def test_agg_must_agg(self):
        grouped = self.df.groupby('A')['C']
        self.assertRaises(Exception, grouped.agg, lambda x: x.describe())
        self.assertRaises(Exception, grouped.agg, lambda x: x.index[:2])

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
        assert_frame_equal(grouped.sum(),
                           DataFrame(columns=self.tsframe.columns))
        assert_frame_equal(grouped.agg(np.sum), DataFrame({}))
        assert_frame_equal(grouped.apply(np.sum), DataFrame({}))

    def test_agg_python_multiindex(self):
        grouped = self.mframe.groupby(['A', 'B'])

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        tm.assert_frame_equal(result, expected)

    def test_apply_describe_bug(self):
        grouped = self.mframe.groupby(level='first')
        result = grouped.describe() # it works!

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
        self.assert_(groups is grouped.groups) # caching works

        for k, v in grouped.groups.iteritems():
            self.assert_((self.df.ix[v]['A'] == k).all())

        grouped = self.df.groupby(['A', 'B'])
        groups = grouped.groups
        self.assert_(groups is grouped.groups) # caching works
        for k, v in grouped.groups.iteritems():
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
            result = grouped.agg({'A' : 'var', 'B' : 'std', 'C' : 'mean'})
            expected = DataFrame({'A' : grouped['A'].var(),
                                  'B' : grouped['B'].std(),
                                  'C' : grouped['C'].mean()})
            assert_frame_equal(result, expected)

        by_weekday = self.tsframe.groupby(lambda x: x.weekday())
        _check_results(by_weekday)

        by_mwkday = self.tsframe.groupby([lambda x: x.month,
                                          lambda x: x.weekday()])
        _check_results(by_mwkday)

    def test_basic_regression(self):
        # regression
        T = [1.0*x for x in range(1,10) *10][:1095]
        result = Series(T, range(0, len(T)))

        groupings = np.random.random((1100,))
        groupings = Series(groupings, range(0, len(groupings))) * 10.

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
            self.assert_((result.reindex(gp.index) == gp.mean()).all())

        grouped = self.tsframe.groupby(lambda x: x.month)
        result = grouped.transform(np.mean)
        self.assert_(result.index.equals(self.tsframe.index))
        for _, gp in grouped:
            agged = gp.mean()
            res = result.reindex(gp.index)
            for col in self.tsframe:
                self.assert_((res[col] == agged[col]).all())

        # group columns
        grouped = self.tsframe.groupby({'A' : 0, 'B' : 0, 'C' : 1, 'D' : 1},
                                       axis=1)
        result = grouped.transform(np.mean)
        self.assert_(result.index.equals(self.tsframe.index))
        self.assert_(result.columns.equals(self.tsframe.columns))
        for _, gp in grouped:
            agged = gp.mean(1)
            res = result.reindex(columns=gp.columns)
            for idx in gp.index:
                self.assert_((res.xs(idx) == agged[idx]).all())

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

    def test_with_na(self):
        index = Index(np.arange(10))
        values = Series(np.ones(10), index)
        labels = Series([nan, 'foo', 'bar', 'bar', nan, nan, 'bar',
                         'bar', nan, 'foo'], index=index)

        grouped = values.groupby(labels)
        agged = grouped.agg(len)
        expected = Series([4, 2], index=['bar', 'foo'])

        assert_series_equal(agged, expected, check_dtype=False)
        self.assert_(issubclass(agged.dtype.type, np.integer))

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
        assert_frame_equal(result, expected)

        # get attribute
        result = grouped.dtype
        expected = grouped.agg(lambda x: x.dtype)

        # make sure raises error
        self.assertRaises(AttributeError, getattr, grouped, 'foo')

    def test_series_describe_multikey(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby([lambda x: x.year, lambda x: x.month])
        result = grouped.describe()
        assert_series_equal(result['mean'], grouped.mean())
        assert_series_equal(result['std'], grouped.std())
        assert_series_equal(result['min'], grouped.min())

    def test_series_describe_single(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby(lambda x: x.month)
        result = grouped.apply(lambda x: x.describe())
        expected = grouped.describe()
        assert_frame_equal(result, expected)

    def test_series_agg_multikey(self):
        ts = tm.makeTimeSeries()
        grouped = ts.groupby([lambda x: x.year, lambda x: x.month])

        result = grouped.agg(np.sum)
        expected = grouped.sum()
        assert_series_equal(result, expected)

    def test_series_agg_multi_pure_python(self):
        data = DataFrame({'A' : ['foo', 'foo', 'foo', 'foo',
                                 'bar', 'bar', 'bar', 'bar',
                                 'foo', 'foo', 'foo'],
                          'B' : ['one', 'one', 'one', 'two',
                                 'one', 'one', 'one', 'two',
                                 'two', 'two', 'one'],
                          'C' : ['dull', 'dull', 'shiny', 'dull',
                                 'dull', 'shiny', 'shiny', 'dull',
                                 'shiny', 'shiny', 'shiny'],
                          'D' : np.random.randn(11),
                          'E' : np.random.randn(11),
                          'F' : np.random.randn(11)})

        def bad(x):
            assert(len(x.base) == len(x))
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
            assert_frame_equal(result[col].unstack(), expected)

        groupedT = self.tsframe.groupby({'A' : 0, 'B' : 0,
                                         'C' : 1, 'D' : 1}, axis=1)
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
        assert_frame_equal(stragged, aggregated)

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
        groups = grouped.primary.groups
        indices = grouped.primary.indices

        for k, v in groups.iteritems():
            samething = self.tsframe.index.take(indices[k])
            self.assert_(np.array_equal(v, samething))

    def test_frame_groupby_columns(self):
        mapping = {
            'A' : 0, 'B' : 0, 'C' : 1, 'D' : 1
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

        result = grouped.agg({'C' : np.mean, 'D' : np.std})
        self.assert_(result.index.name == 'A')

        result = grouped['C'].mean()
        self.assert_(result.index.name == 'A')
        result = grouped['C'].agg(np.mean)
        self.assert_(result.index.name == 'A')
        result = grouped['C'].agg([np.mean, np.std])
        self.assert_(result.index.name == 'A')

        result = grouped['C'].agg({'foo' : np.mean, 'bar' : np.std})
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
        df = DataFrame({'v1' : np.random.randn(6),
                        'v2' : np.random.randn(6),
                        'k1' : k1, 'k2' : k2},
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
                           expected.ix[:, ['C', 'D']])

        # some "groups" with no data
        df = DataFrame({'v1' : np.random.randn(6),
                        'v2' : np.random.randn(6),
                        'k1' : np.array(['b', 'b', 'b', 'a', 'a', 'a']),
                        'k2' : np.array(['1', '1', '1', '2', '2', '2'])},
                       index=['one', 'two', 'three', 'four', 'five', 'six'])
        # only verify that it works for now
        grouped = df.groupby(['k1', 'k2'])
        grouped.agg(np.sum)

    def test_multi_key_multiple_functions(self):
        grouped = self.df.groupby(['A', 'B'])['C']

        agged = grouped.agg([np.mean, np.std])
        expected = DataFrame({'mean' : grouped.agg(np.mean),
                              'std' : grouped.agg(np.std)})
        assert_frame_equal(agged, expected)

    def test_frame_multi_key_function_list(self):
        data = DataFrame({'A' : ['foo', 'foo', 'foo', 'foo',
                                 'bar', 'bar', 'bar', 'bar',
                                 'foo', 'foo', 'foo'],
                          'B' : ['one', 'one', 'one', 'two',
                                 'one', 'one', 'one', 'two',
                                 'two', 'two', 'one'],
                          'C' : ['dull', 'dull', 'shiny', 'dull',
                                 'dull', 'shiny', 'shiny', 'dull',
                                 'shiny', 'shiny', 'shiny'],
                          'D' : np.random.randn(11),
                          'E' : np.random.randn(11),
                          'F' : np.random.randn(11)})

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
            expected = dict((k, DataFrame(v)) for k, v in expected.iteritems())
            expected = Panel.fromDict(expected).swapaxes(0, 1)

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

        result2 = grouped.agg({'C' : np.mean, 'D' : np.sum})
        expected2 = grouped.mean()
        expected2['D'] = grouped.sum()['D']
        assert_frame_equal(result2, expected2)

        # multi-key

        grouped = self.df.groupby(['A', 'B'], as_index=False)

        result = grouped.agg(np.mean)
        expected = grouped.mean()
        assert_frame_equal(result, expected)

        result2 = grouped.agg({'C' : np.mean, 'D' : np.sum})
        expected2 = grouped.mean()
        expected2['D'] = grouped.sum()['D']
        assert_frame_equal(result2, expected2)

    def test_as_index_series_return_frame(self):
        grouped = self.df.groupby('A', as_index=False)
        grouped2 = self.df.groupby(['A', 'B'], as_index=False)

        result = grouped['C'].agg(np.sum)
        expected = grouped.agg(np.sum).ix[:, ['A', 'C']]
        self.assert_(isinstance(result, DataFrame))
        assert_frame_equal(result, expected)

        result2 = grouped2['C'].agg(np.sum)
        expected2 = grouped2.agg(np.sum).ix[:, ['A', 'B', 'C']]
        self.assert_(isinstance(result2, DataFrame))
        assert_frame_equal(result2, expected2)

        result = grouped['C'].sum()
        expected = grouped.sum().ix[:, ['A', 'C']]
        self.assert_(isinstance(result, DataFrame))
        assert_frame_equal(result, expected)

        result2 = grouped2['C'].sum()
        expected2 = grouped2.sum().ix[:, ['A', 'B', 'C']]
        self.assert_(isinstance(result2, DataFrame))
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

        arrays = zip(*expected.index.get_tuple_index())
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
        grouped = df.groupby({'A' : 0, 'C' : 0, 'D' : 1, 'E' : 1}, axis=1)
        result = self.assertRaises(TypeError, grouped.agg, np.sum)

    def test_omit_nuisance_python_multiple(self):
        grouped = self.three_group.groupby(['A', 'B'])

        agged = grouped.agg(np.mean)
        exp = grouped.mean()
        assert_frame_equal(agged, exp)

    def test_empty_groups_corner(self):
        # handle empty groups
        df = DataFrame({'k1' : np.array(['b', 'b', 'b', 'a', 'a', 'a']),
                        'k2' : np.array(['1', '1', '1', '2', '2', '2']),
                        'k3' : ['foo', 'bar'] * 3,
                        'v1' : np.random.randn(6),
                        'v2' : np.random.randn(6)})

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
        self.assertEquals(result.index.names, ['A', 'B', 'stat'])

        result2 = grouped.apply(desc2)
        self.assertEquals(result2.index.names, ['A', 'B', 'stat'])

        result3 = grouped.apply(desc3)
        self.assertEquals(result3.index.names, ['A', 'B', None])

    def test_nonsense_func(self):
        df = DataFrame([0])
        self.assertRaises(Exception, df.groupby, lambda x: x + 'foo')

    def test_cythonized_aggers(self):
        data = {'A' : [0, 0, 0, 0, 1, 1, 1, 1, 1, 1., nan, nan],
                'B' : ['A', 'B'] * 6,
                'C' : np.random.randn(12)}
        df = DataFrame(data)
        df['C'][2:10:2] = nan

        def _testit(op):
            # single column
            grouped = df.drop(['B'], axis=1).groupby('A')
            exp = {}
            for cat, group in grouped:
                exp[cat] = op(group['C'])
            exp = DataFrame({'C' : exp})
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

    def test_cython_agg_boolean(self):
        frame = DataFrame({'a': np.random.randint(0, 5, 50),
                           'b': np.random.randint(0, 2, 50).astype('bool')})
        result = frame.groupby('a')['b'].mean()
        expected = frame.groupby('a')['b'].agg(np.mean)

        assert_series_equal(result, expected)

    def test_cython_agg_nothing_to_agg(self):
        frame = DataFrame({'a': np.random.randint(0, 5, 50),
                           'b': ['foo', 'bar'] * 25})
        self.assertRaises(GroupByError, frame.groupby('a')['b'].mean)

        frame = DataFrame({'a': np.random.randint(0, 5, 50),
                           'b': ['foo', 'bar'] * 25})
        self.assertRaises(GroupByError, frame[['b']].groupby(frame['a']).mean)

    def test_wrap_aggregated_output_multindex(self):
        df = self.mframe.T
        df['baz', 'two'] = 'peekaboo'

        keys = [np.array([0, 0, 1]), np.array([0, 0, 1])]
        agged = df.groupby(keys).agg(np.mean)
        self.assert_(isinstance(agged.columns, MultiIndex))

    def test_grouping_attrs(self):
        deleveled = self.mframe.reset_index()
        grouped = deleveled.groupby(['first', 'second'])

        for i, ping in enumerate(grouped.groupings):
            the_counts = self.mframe.groupby(level=i).count()['A']
            other_counts = Series(ping.counts, ping.group_index)
            assert_almost_equal(the_counts,
                                other_counts.reindex(the_counts.index))

        # compute counts when group by level
        grouped = self.mframe.groupby(level=0)
        ping = grouped.groupings[0]
        the_counts = grouped.size()
        other_counts = Series(ping.counts, ping.group_index)
        assert_almost_equal(the_counts,
                            other_counts.reindex(the_counts.index))

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
        self.assertRaises(ValueError, self.df.groupby, level=0)

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

        mapper0 = {'foo' : 0, 'bar' : 0,
                   'baz' : 1, 'qux' : 1}
        mapper1 = {'one' : 0, 'two' : 0, 'three' : 1}

        result0 = frame.groupby(mapper0, level=0).sum()
        result1 = frame.groupby(mapper1, level=1).sum()

        mapped_level0 = np.array([mapper0.get(x) for x in deleveled['first']])
        mapped_level1 = np.array([mapper1.get(x) for x in deleveled['second']])
        expected0 = frame.groupby(mapped_level0).sum()
        expected1 = frame.groupby(mapped_level1).sum()

        assert_frame_equal(result0, expected0)
        assert_frame_equal(result1, expected1)

    def test_level_preserve_order(self):
        grouped = self.mframe.groupby(level=0)
        exp_labels = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 3])
        assert_almost_equal(grouped.groupings[0].labels, exp_labels)

    def test_grouping_labels(self):
        grouped = self.mframe.groupby(self.mframe.index.get_level_values(0))
        exp_labels = np.array([2, 2, 2, 0, 0, 1, 1, 3, 3, 3])
        assert_almost_equal(grouped.groupings[0].labels, exp_labels)

    def test_cython_fail_agg(self):
        dr = DateRange('1/1/2000', periods=50)
        ts = Series(['A', 'B', 'C', 'D', 'E'] * 10, index=dr)

        grouped = ts.groupby(lambda x: x.month)
        summed = grouped.sum()
        expected = grouped.agg(np.sum)
        assert_series_equal(summed, expected)

    def test_apply_series_to_frame(self):
        def f(piece):
            return DataFrame({'value' : piece,
                              'demeaned' : piece - piece.mean(),
                              'logged' : np.log(piece)})

        dr = DateRange('1/1/2000', periods=100)
        ts = Series(np.random.randn(100), index=dr)

        grouped = ts.groupby(lambda x: x.month)
        result = grouped.apply(f)

        self.assert_(isinstance(result, DataFrame))
        self.assert_(result.index.equals(ts.index))

    def test_apply_frame_to_series(self):
        grouped = self.df.groupby(['A', 'B'])
        result = grouped.apply(len)
        expected = grouped.count()['C']
        self.assert_(result.index.equals(expected.index))
        self.assert_(np.array_equal(result.values, expected.values))

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

    def test_groupby_series_indexed_differently(self):
        s1 = Series([5.0,-9.0,4.0,100.,-5.,55.,6.7],
                    index=Index(['a','b','c','d','e','f','g']))
        s2 = Series([1.0,1.0,4.0,5.0,5.0,7.0],
                    index=Index(['a','b','d','f','g','h']))

        grouped = s1.groupby(s2)
        agged = grouped.mean()
        exp = s1.groupby(s2.reindex(s1.index).get).mean()
        assert_series_equal(agged, exp)

    def test_groupby_with_hier_columns(self):
        tuples = zip(*[['bar', 'bar', 'baz', 'baz',
                        'foo', 'foo', 'qux', 'qux'],
                       ['one', 'two', 'one', 'two',
                        'one', 'two', 'one', 'two']])
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
        self.assert_(result.columns.equals(sorted_columns))

    def test_pass_args_kwargs(self):
        from scipy.stats import scoreatpercentile

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
    #                   np.array([1, 2, 3, 4, 0, 1, 2, 3, 3, 4], dtype=np.int32)]

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
        assert_frame_equal(result, expected)

    def test_apply_typecast_fail(self):
        df = DataFrame({'d' : [1.,1.,1.,2.,2.,2.],
                        'c' : np.tile(['a','b','c'], 2),
                        'v' : np.arange(1., 7.)})

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
        df = DataFrame({'d' : [1.,1.,1.,2.,2.,2.],
                        'c' : np.tile(['a','b','c'], 2),
                        'v' : np.arange(1., 7.)}, index=index)

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

    def test_transform_mixed_type(self):
        index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1],
                                        [1, 2, 3, 1, 2, 3]])
        df = DataFrame({'d' : [1.,1.,1.,2.,2.,2.],
                        'c' : np.tile(['a','b','c'], 2),
                        'v' : np.arange(1., 7.)}, index=index)

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
        from cStringIO import StringIO
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
        self.assertEquals(result.index.names, ['A', 'B'])
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
        obj = Series(arr[:, 0], index=range(100))
        inds = np.tile(range(10), 10)

        result = obj.groupby(inds).agg(Series.median)
        self.assert_(result.isnull().all())

    def test_convert_objects_leave_decimal_alone(self):
        from decimal import Decimal

        s = Series(range(5))
        labels = np.array(['a', 'b', 'c', 'd', 'e'], dtype='O')

        def convert_fast(x):
            return Decimal(str(x.mean()))

        def convert_force_pure(x):
            # base will be length 0
            assert(len(x.base) == len(x))
            return Decimal(str(x.mean()))

        grouped = s.groupby(labels)

        result = grouped.agg(convert_fast)
        self.assert_(result.dtype == np.object_)
        self.assert_(isinstance(result[0], Decimal))

        result = grouped.agg(convert_force_pure)
        self.assert_(result.dtype == np.object_)
        self.assert_(isinstance(result[0], Decimal))

    def test_groupby_list_infer_array_like(self):
        result = self.df.groupby(list(self.df['A'])).mean()
        expected = self.df.groupby(self.df['A']).mean()
        assert_frame_equal(result, expected)

        self.assertRaises(Exception, self.df.groupby, list(self.df['A'][:-1]))

        # pathological case of ambiguity
        df = DataFrame({'foo' : [0, 1], 'bar' : [3, 4],
                        'val' : np.random.randn(2)})

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
        grouped = self.panel.groupby({'ItemA' : 0, 'ItemB' : 0, 'ItemC' : 1},
                                     axis='items')
        agged = grouped.agg(np.mean)
        self.assert_(np.array_equal(agged.items, [0, 1]))

        grouped = self.panel.groupby(lambda x: x.month, axis='major')
        agged = grouped.agg(np.mean)

        self.assert_(np.array_equal(agged.major_axis, [1, 2]))

        grouped = self.panel.groupby({'A' : 0, 'B' : 0, 'C' : 1, 'D' : 1},
                                     axis='minor')
        agged = grouped.agg(np.mean)
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
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
