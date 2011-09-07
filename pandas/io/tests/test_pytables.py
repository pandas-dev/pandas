import nose
import unittest
import os
import sys

import numpy as np

from pandas import (Series, DataFrame, Panel, LongPanel, DateRange)
from pandas.io.pytables import HDFStore
import pandas.util.testing as tm

try:
    import tables
except ImportError:
    raise nose.SkipTest('no pytables')

class TesttHDFStore(unittest.TestCase):
    path = '__test__.h5'
    scratchpath = '__scratch__.h5'

    def setUp(self):
        self.store = HDFStore(self.path)

    def tearDown(self):
        self.store.close()
        os.remove(self.path)

    def test_len(self):
        self.store['a'] = tm.makeTimeSeries()
        self.store['b'] = tm.makeStringSeries()
        self.store['c'] = tm.makeDataFrame()
        self.store['d'] = tm.makePanel()
        self.assertEquals(len(self.store), 4)

    def test_repr(self):
        repr(self.store)
        self.store['a'] = tm.makeTimeSeries()
        self.store['b'] = tm.makeStringSeries()
        self.store['c'] = tm.makeDataFrame()
        self.store['d'] = tm.makePanel()
        repr(self.store)

    def test_reopen_handle(self):
        self.store['a'] = tm.makeTimeSeries()
        self.store.open('w', warn=False)
        self.assert_(self.store.handle.isopen)
        self.assertEquals(len(self.store), 0)

    def test_flush(self):
        self.store['a'] = tm.makeTimeSeries()
        self.store.flush()

    def test_get(self):
        self.store['a'] = tm.makeTimeSeries()
        left = self.store.get('a')
        right = self.store['a']
        tm.assert_series_equal(left, right)

        self.assertRaises(AttributeError, self.store.get, 'b')

    def test_put(self):
        ts = tm.makeTimeSeries()
        df = tm.makeTimeDataFrame()
        self.store['a'] = ts
        self.store['b'] = df[:10]
        self.store.put('c', df[:10], table=True)

        # not OK, not a table
        self.assertRaises(ValueError, self.store.put, 'b', df[10:], append=True)

        # node does not currently exist, test _is_table_type returns False in
        # this case
        self.assertRaises(ValueError, self.store.put, 'f', df[10:], append=True)

        # OK
        self.store.put('c', df[10:], append=True)

        # overwrite table
        self.store.put('c', df[:10], table=True, append=False)
        tm.assert_frame_equal(df[:10], self.store['c'])

    def test_put_compression(self):
        df = tm.makeTimeDataFrame()
        self.store.put('c', df, table=True, compression='blosc')
        tm.assert_frame_equal(self.store['c'], df)

        self.store.put('c', df, table=True, compression='zlib')
        tm.assert_frame_equal(self.store['c'], df)

        # can't compress if table=False
        self.assertRaises(ValueError, self.store.put, 'b', df,
                          table=False, compression='blosc')

    def test_put_integer(self):
        # non-date, non-string index
        df = DataFrame(np.random.randn(50, 100))
        self._check_roundtrip(df, tm.assert_frame_equal)

    def test_append(self):
        df = tm.makeTimeDataFrame()
        self.store.put('c', df[:10], table=True)
        self.store.append('c', df[10:])
        tm.assert_frame_equal(self.store['c'], df)

    def test_remove(self):
        ts = tm.makeTimeSeries()
        df = tm.makeDataFrame()
        self.store['a'] = ts
        self.store['b'] = df
        self.store.remove('a')
        self.assertEquals(len(self.store), 1)
        tm.assert_frame_equal(df, self.store['b'])

        self.store.remove('b')
        self.assertEquals(len(self.store), 0)

    def test_remove_crit(self):
        wp = tm.makePanel()
        self.store.put('wp', wp, table=True)
        date = wp.major_axis[len(wp.major_axis) // 2]

        crit1 = {
            'field' : 'index',
            'op' : '>',
            'value' : date
        }
        crit2 = {
            'field' : 'column',
            'value' : ['A', 'D']
        }
        self.store.remove('wp', where=[crit1])
        self.store.remove('wp', where=[crit2])
        result = self.store['wp']
        expected = wp.truncate(after=date).reindex(minor=['B', 'C'])
        tm.assert_panel_equal(result, expected)

    def test_series(self):
        s = tm.makeStringSeries()
        self._check_roundtrip(s, tm.assert_series_equal)

        ts = tm.makeTimeSeries()
        self._check_roundtrip(ts, tm.assert_series_equal)

    def test_timeseries_preepoch(self):
        if sys.version_info[0] == 2 and sys.version_info[1] < 7:
            raise nose.SkipTest

        dr = DateRange('1/1/1940', '1/1/1960')
        ts = Series(np.random.randn(len(dr)), index=dr)
        self._check_roundtrip(ts, tm.assert_series_equal)

    def test_frame(self):
        df = tm.makeDataFrame()

        # put in some random NAs
        df.values[0, 0] = np.nan
        df.values[5, 3] = np.nan

        self._check_roundtrip_table(df, tm.assert_frame_equal)
        self._check_roundtrip(df, tm.assert_frame_equal)

        self._check_roundtrip_table(df, tm.assert_frame_equal,
                                    compression=True)
        self._check_roundtrip(df, tm.assert_frame_equal,
                                    compression=True)

        tdf = tm.makeTimeDataFrame()
        self._check_roundtrip(tdf, tm.assert_frame_equal)
        self._check_roundtrip(tdf, tm.assert_frame_equal,
                              compression=True)

        # not consolidated
        df['foo'] = np.random.randn(len(df))
        self.store['df'] = df
        recons = self.store['df']
        self.assert_(recons._data.is_consolidated())

    def test_store_mixed(self):
        def _make_one():
            df = tm.makeDataFrame()
            df['obj1'] = 'foo'
            df['obj2'] = 'bar'
            df['bool1'] = df['A'] > 0
            df['bool2'] = df['B'] > 0
            df['int1'] = 1
            df['int2'] = 2
            return df.consolidate()

        df1 = _make_one()
        df2 = _make_one()

        self._check_roundtrip(df1, tm.assert_frame_equal)
        self._check_roundtrip(df2, tm.assert_frame_equal)

        self.store['obj'] = df1
        tm.assert_frame_equal(self.store['obj'], df1)
        self.store['obj'] = df2
        tm.assert_frame_equal(self.store['obj'], df2)

        # storing in Table not yet supported
        self.assertRaises(Exception, self.store.put, 'foo',
                          df1, table=True)

        # check that can store Series of all of these types
        self._check_roundtrip(df1['obj1'], tm.assert_series_equal)
        self._check_roundtrip(df1['bool1'], tm.assert_series_equal)
        self._check_roundtrip(df1['int1'], tm.assert_series_equal)

        # try with compression
        self._check_roundtrip(df1['obj1'], tm.assert_series_equal,
                              compression=True)
        self._check_roundtrip(df1['bool1'], tm.assert_series_equal,
                              compression=True)
        self._check_roundtrip(df1['int1'], tm.assert_series_equal,
                              compression=True)
        self._check_roundtrip(df1, tm.assert_frame_equal,
                              compression=True)

    def test_wide(self):
        wp = tm.makePanel()
        self._check_roundtrip(wp, tm.assert_panel_equal)

    def test_wide_table(self):
        wp = tm.makePanel()
        self._check_roundtrip_table(wp, tm.assert_panel_equal)

    def test_long(self):
        def _check(left, right):
            tm.assert_panel_equal(left.to_wide(),
                                  right.to_wide())

        wp = tm.makePanel()
        self._check_roundtrip(wp.to_long(), _check)

    def test_longpanel(self):
        pass

    def test_overwrite_node(self):
        self.store['a'] = tm.makeTimeDataFrame()
        ts = tm.makeTimeSeries()
        self.store['a'] = ts

        tm.assert_series_equal(self.store['a'], ts)

    def test_panel_select(self):
        wp = tm.makePanel()
        self.store.put('wp', wp, table=True)
        date = wp.major_axis[len(wp.major_axis) // 2]

        crit1 = {
            'field' : 'index',
            'op' : '>=',
            'value' : date
        }
        crit2 = {
            'field' : 'column',
            'value' : ['A', 'D']
        }

        result = self.store.select('wp', [crit1, crit2])
        expected = wp.truncate(before=date).reindex(minor=['A', 'D'])
        tm.assert_panel_equal(result, expected)

    def test_frame_select(self):
        df = tm.makeTimeDataFrame()
        self.store.put('frame', df, table=True)
        date = df.index[len(df) // 2]

        crit1 = {
            'field' : 'index',
            'op' : '>=',
            'value' : date
        }
        crit2 = {
            'field' : 'column',
            'value' : ['A', 'D']
        }
        crit3 = {
            'field' : 'column',
            'value' : 'A'
        }

        result = self.store.select('frame', [crit1, crit2])
        expected = df.ix[date:, ['A', 'D']]
        tm.assert_frame_equal(result, expected)

        result = self.store.select('frame', [crit3])
        expected = df.ix[:, ['A']]
        tm.assert_frame_equal(result, expected)

        # can't select if not written as table
        self.store['frame'] = df
        self.assertRaises(Exception, self.store.select,
                          'frame', [crit1, crit2])

    def test_select_filter_corner(self):
        df = DataFrame(np.random.randn(50, 100))
        df.index = ['%.3d' % c for c in df.index]
        df.columns = ['%.3d' % c for c in df.columns]
        self.store.put('frame', df, table=True)

        crit = {
            'field' : 'column',
            'value' : df.columns[:75]
        }
        result = self.store.select('frame', [crit])
        tm.assert_frame_equal(result, df.ix[:, df.columns[:75]])

    def _check_roundtrip(self, obj, comparator, compression=False):
        options = {}
        if compression:
            options['complib'] = 'blosc'

        store = HDFStore(self.scratchpath, 'w', **options)
        try:
            store['obj'] = obj
            retrieved = store['obj']
            comparator(retrieved, obj)
        finally:
            store.close()
            os.remove(self.scratchpath)

    def _check_roundtrip_table(self, obj, comparator, compression=False):
        options = {}
        if compression:
            options['complib'] = 'blosc'

        store = HDFStore(self.scratchpath, 'w', **options)
        try:
            store.put('obj', obj, table=True)
            retrieved = store['obj']
            sorted_obj = _test_sort(obj)
            comparator(retrieved, sorted_obj)
        finally:
            store.close()
            os.remove(self.scratchpath)

    def test_legacy_read(self):
        pth = curpath()
        store = HDFStore(os.path.join(pth, 'legacy.h5'), 'r')
        store['a']
        store['b']
        store['c']
        store['d']
        store.close()

def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

def _test_sort(obj):
    if isinstance(obj, DataFrame):
        return obj.reindex(sorted(obj.index))
    elif isinstance(obj, Panel):
        return obj.reindex(major=sorted(obj.major_axis))
    else:
        raise ValueError('type not supported here')

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

