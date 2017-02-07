# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import datetime

import numpy as np
import pandas as pd

from pandas import Index, Series
from pandas.core.index import MultiIndex, RangeIndex

from pandas.compat import lrange, range, zip
from pandas.util.testing import assert_series_equal, assert_frame_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesAlterAxes(TestData, tm.TestCase):

    def test_setindex(self):
        # wrong type
        series = self.series.copy()
        self.assertRaises(TypeError, setattr, series, 'index', None)

        # wrong length
        series = self.series.copy()
        self.assertRaises(Exception, setattr, series, 'index',
                          np.arange(len(series) - 1))

        # works
        series = self.series.copy()
        series.index = np.arange(len(series))
        tm.assertIsInstance(series.index, Index)

    def test_rename(self):
        renamer = lambda x: x.strftime('%Y%m%d')
        renamed = self.ts.rename(renamer)
        self.assertEqual(renamed.index[0], renamer(self.ts.index[0]))

        # dict
        rename_dict = dict(zip(self.ts.index, renamed.index))
        renamed2 = self.ts.rename(rename_dict)
        assert_series_equal(renamed, renamed2)

        # partial dict
        s = Series(np.arange(4), index=['a', 'b', 'c', 'd'], dtype='int64')
        renamed = s.rename({'b': 'foo', 'd': 'bar'})
        self.assert_index_equal(renamed.index, Index(['a', 'foo', 'c', 'bar']))

        # index with name
        renamer = Series(np.arange(4),
                         index=Index(['a', 'b', 'c', 'd'], name='name'),
                         dtype='int64')
        renamed = renamer.rename({})
        self.assertEqual(renamed.index.name, renamer.index.name)

    def test_rename_by_series(self):
        s = Series(range(5), name='foo')
        renamer = Series({1: 10, 2: 20})
        result = s.rename(renamer)
        expected = Series(range(5), index=[0, 10, 20, 3, 4], name='foo')
        tm.assert_series_equal(result, expected)

    def test_rename_set_name(self):
        s = Series(range(4), index=list('abcd'))
        for name in ['foo', 123, 123., datetime(2001, 11, 11), ('foo',)]:
            result = s.rename(name)
            self.assertEqual(result.name, name)
            self.assert_numpy_array_equal(result.index.values, s.index.values)
            self.assertTrue(s.name is None)

    def test_rename_set_name_inplace(self):
        s = Series(range(3), index=list('abc'))
        for name in ['foo', 123, 123., datetime(2001, 11, 11), ('foo',)]:
            s.rename(name, inplace=True)
            self.assertEqual(s.name, name)

            exp = np.array(['a', 'b', 'c'], dtype=np.object_)
            self.assert_numpy_array_equal(s.index.values, exp)

    def test_set_name_attribute(self):
        s = Series([1, 2, 3])
        s2 = Series([1, 2, 3], name='bar')
        for name in [7, 7., 'name', datetime(2001, 1, 1), (1,), u"\u05D0"]:
            s.name = name
            self.assertEqual(s.name, name)
            s2.name = name
            self.assertEqual(s2.name, name)

    def test_set_name(self):
        s = Series([1, 2, 3])
        s2 = s._set_name('foo')
        self.assertEqual(s2.name, 'foo')
        self.assertTrue(s.name is None)
        self.assertTrue(s is not s2)

    def test_rename_inplace(self):
        renamer = lambda x: x.strftime('%Y%m%d')
        expected = renamer(self.ts.index[0])

        self.ts.rename(renamer, inplace=True)
        self.assertEqual(self.ts.index[0], expected)

    def test_set_index_makes_timeseries(self):
        idx = tm.makeDateIndex(10)

        s = Series(lrange(10))
        s.index = idx

        with tm.assert_produces_warning(FutureWarning):
            self.assertTrue(s.is_time_series)
        self.assertTrue(s.index.is_all_dates)

    def test_reset_index(self):
        df = tm.makeDataFrame()[:5]
        ser = df.stack()
        ser.index.names = ['hash', 'category']

        ser.name = 'value'
        df = ser.reset_index()
        self.assertIn('value', df)

        df = ser.reset_index(name='value2')
        self.assertIn('value2', df)

        # check inplace
        s = ser.reset_index(drop=True)
        s2 = ser
        s2.reset_index(drop=True, inplace=True)
        assert_series_equal(s, s2)

        # level
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        rs = s.reset_index(level=1)
        self.assertEqual(len(rs.columns), 2)

        rs = s.reset_index(level=[0, 2], drop=True)
        self.assert_index_equal(rs.index, Index(index.get_level_values(1)))
        tm.assertIsInstance(rs, Series)

    def test_reset_index_range(self):
        # GH 12071
        s = pd.Series(range(2), name='A', dtype='int64')
        series_result = s.reset_index()
        tm.assertIsInstance(series_result.index, RangeIndex)
        series_expected = pd.DataFrame([[0, 0], [1, 1]],
                                       columns=['index', 'A'],
                                       index=RangeIndex(stop=2))
        assert_frame_equal(series_result, series_expected)

    def test_reorder_levels(self):
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]],
                           names=['L0', 'L1', 'L2'])
        s = Series(np.arange(6), index=index)

        # no change, position
        result = s.reorder_levels([0, 1, 2])
        assert_series_equal(s, result)

        # no change, labels
        result = s.reorder_levels(['L0', 'L1', 'L2'])
        assert_series_equal(s, result)

        # rotate, position
        result = s.reorder_levels([1, 2, 0])
        e_idx = MultiIndex(levels=[['one', 'two', 'three'], [0, 1], ['bar']],
                           labels=[[0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L1', 'L2', 'L0'])
        expected = Series(np.arange(6), index=e_idx)
        assert_series_equal(result, expected)

        result = s.reorder_levels([0, 0, 0])
        e_idx = MultiIndex(levels=[['bar'], ['bar'], ['bar']],
                           labels=[[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L0', 'L0', 'L0'])
        expected = Series(range(6), index=e_idx)
        assert_series_equal(result, expected)

        result = s.reorder_levels(['L0', 'L0', 'L0'])
        assert_series_equal(result, expected)
