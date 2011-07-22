import unittest

from numpy import nan

from pandas.core.daterange import DateRange
from pandas.core.index import Index
from pandas.core.common import rands, groupby
from pandas.core.frame import DataFrame
from pandas.core.series import Series
from pandas.util.testing import (assert_panel_equal, assert_frame_equal,
                                 assert_series_equal, assert_almost_equal)
from pandas.core.panel import WidePanel
from collections import defaultdict
import pandas.core.datetools as dt
import numpy as np

import pandas.util.testing as tm

# unittest.TestCase

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


class GroupByTestCase(unittest.TestCase):
    setUp = commonSetUp

    def test_python_grouper(self):
        groupFunc = self.groupDict.get
        groups = groupby(self.stringIndex, groupFunc)
        setDict = dict((k, set(v)) for k, v in groups.iteritems())
        for idx in self.stringIndex:
            key = groupFunc(idx)
            groupSet = setDict[key]
            assert(idx in groupSet)

class TestSeriesGroupBy(unittest.TestCase):

    def test_groupby(self):
        data = Series(np.arange(9) / 3, index=np.arange(9))

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
        assert_series_equal(grouped.agg(np.sum), grouped.sum())

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
        agged = grouped.agg(lambda x: group_constants[x.groupName] + x.mean())
        self.assertEqual(agged[1], 21)

        # corner cases
        self.assertRaises(Exception, grouped._aggregate_named,
                          lambda x: x * 2)

    def test_groupby_transform(self):
        data = Series(np.arange(9) / 3, index=np.arange(9))

        index = np.arange(9)
        np.random.shuffle(index)
        data = data.reindex(index)

        grouped = data.groupby(lambda x: x // 3)

        transformed = grouped.transform(lambda x: x * x.sum())
        self.assertEqual(transformed[7], 12)

        transformed = grouped.transform(np.mean)
        for name, group in grouped:
            mean = group.mean()
            for idx in group.index:
                self.assertEqual(transformed[idx], mean)

    def test_groupby_with_na(self):
        index = Index(np.arange(10))
        values = Series(np.ones(10), index)
        labels = Series([nan, 'foo', 'bar', 'bar', nan, nan, 'bar',
                         'bar', nan, 'foo'], index=index)

        grouped = values.groupby(labels)
        agged = grouped.agg(len)
        expected = Series([4, 2], index=['bar', 'foo'])
        assert_series_equal(agged, expected)

class TestDataFrameGroupBy(unittest.TestCase):

    def setUp(self):
        self.seriesd = tm.getSeriesData()
        self.tsd = tm.getTimeSeriesData()
        self.frame = DataFrame(self.seriesd)
        self.tsframe = DataFrame(self.tsd)

    def test_groupby(self):
        grouped = self.tsframe.groupby(lambda x: x.weekday())

        # aggregate
        aggregated = grouped.aggregate(np.mean)
        self.assertEqual(len(aggregated), 5)
        self.assertEqual(len(aggregated.columns), 4)

        # by string
        tscopy = self.tsframe.copy()
        tscopy['weekday'] = [x.weekday() for x in tscopy.index]
        stragged = tscopy.groupby('weekday').aggregate(np.mean)
        del stragged['weekday']
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

    def test_groupby_columns(self):
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

        # tgroupby
        grouping = {
            'A' : 0,
            'B' : 1,
            'C' : 0,
            'D' : 1
        }

        grouped = self.frame.tgroupby(grouping.get, np.mean)
        self.assertEqual(len(grouped), len(self.frame.index))
        self.assertEqual(len(grouped.columns), 2)

    def test_groupby_multiple_columns(self):
        data = DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                                 'foo', 'bar', 'foo', 'foo'],
                          'B' : ['one', 'one', 'two', 'three',
                                 'two', 'two', 'one', 'three'],
                          'C' : np.random.randn(8),
                          'D' : np.random.randn(8)})

        result1 = data.groupby('A', 'B').sum()
        result1_col = data.groupby('A', 'B')['C'].sum()
        expected = defaultdict(dict)
        for n1, gp1 in data.groupby('A'):
            for n2, gp2 in gp1.groupby('B'):
                expected[n1][n2] = gp2.sum()
        expected = dict((k, DataFrame(v)) for k, v in expected.iteritems())
        expected = WidePanel.fromDict(expected)
        assert_panel_equal(result1, expected)
        assert_panel_equal(result1['C'], expected['C'])

        # result2 = data.groupby('B', 'A').sum()
        # assert_panel_equal(result2, expected2)

class TestPanelGroupBy(unittest.TestCase):

    def setUp(self):
        self.panel = tm.makeWidePanel()
        tm.add_nans(self.panel)

    def test_groupby(self):
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


class TestAggregate(unittest.TestCase):
    setUp = commonSetUp

class TestTransform(unittest.TestCase):
    setUp = commonSetUp

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
