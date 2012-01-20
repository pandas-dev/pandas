import nose
import unittest
import os
import sys

import numpy as np

from pandas import Series, DataFrame, Panel, DateRange, MultiIndex
from pandas.io.pytables import HDFStore
import pandas.util.testing as tm

try:
    import tables
except ImportError:
    raise nose.SkipTest('no pytables')

from distutils.version import LooseVersion

_default_compressor = LooseVersion(tables.__version__) >= '2.2' \
                      and 'blosc' or 'zlib'

class TesttHDFStore(unittest.TestCase):
    path = '__test__.h5'
    scratchpath = '__scratch__.h5'

    def setUp(self):
        self.store = HDFStore(self.path)

    def tearDown(self):
        self.store.close()
        os.remove(self.path)

    def test_len_keys(self):
        self.store['a'] = tm.makeTimeSeries()
        self.store['b'] = tm.makeStringSeries()
        self.store['c'] = tm.makeDataFrame()
        self.store['d'] = tm.makePanel()
        self.assertEquals(len(self.store), 4)
        self.assert_(set(self.store.keys()) == set(['a', 'b', 'c', 'd']))

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

        self.store.put('c', df, table=True, compression='zlib')
        tm.assert_frame_equal(self.store['c'], df)

        # can't compress if table=False
        self.assertRaises(ValueError, self.store.put, 'b', df,
                          table=False, compression='zlib')

    def test_put_compression_blosc(self):
        tm.skip_if_no_package('tables', '2.2', app='blosc support')
        df = tm.makeTimeDataFrame()

        # can't compress if table=False
        self.assertRaises(ValueError, self.store.put, 'b', df,
                          table=False, compression='blosc')

        self.store.put('c', df, table=True, compression='blosc')
        tm.assert_frame_equal(self.store['c'], df)

    def test_put_integer(self):
        # non-date, non-string index
        df = DataFrame(np.random.randn(50, 100))
        self._check_roundtrip(df, tm.assert_frame_equal)

    def test_append(self):
        df = tm.makeTimeDataFrame()
        self.store.put('c', df[:10], table=True)
        self.store.append('c', df[10:])
        tm.assert_frame_equal(self.store['c'], df)

    def test_append_diff_item_order(self):
        wp = tm.makePanel()
        wp1 = wp.ix[:, :10, :]
        wp2 = wp.ix[['ItemC', 'ItemB', 'ItemA'], 10:, :]

        self.store.put('panel', wp1, table=True)
        self.assertRaises(Exception, self.store.put, 'panel', wp2,
                          append=True)

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

    def test_remove_where_not_exist(self):
        crit1 = {
            'field' : 'index',
            'op' : '>',
            'value' : 'foo'
        }
        self.store.remove('a', where=[crit1])

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

    def test_float_index(self):
        # GH #454
        index = np.random.randn(10)
        s = Series(np.random.randn(10), index=index)
        self._check_roundtrip(s, tm.assert_series_equal)

    def test_tuple_index(self):
        # GH #492
        col = np.arange(10)
        idx = [(0.,1.), (2., 3.), (4., 5.)]
        data = np.random.randn(30).reshape((3, 10))
        DF = DataFrame(data, index=idx, columns=col)
        self._check_roundtrip(DF, tm.assert_frame_equal)

    def test_timeseries_preepoch(self):
        if sys.version_info[0] == 2 and sys.version_info[1] < 7:
            raise nose.SkipTest

        dr = DateRange('1/1/1940', '1/1/1960')
        ts = Series(np.random.randn(len(dr)), index=dr)
        try:
            self._check_roundtrip(ts, tm.assert_series_equal)
        except OverflowError:
            raise nose.SkipTest('known failer on some windows platforms')

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

        # empty
        self.assertRaises(ValueError, self._check_roundtrip, df[:0],
                          tm.assert_frame_equal)

    def test_can_serialize_dates(self):
        rng = [x.date() for x in DateRange('1/1/2000', '1/30/2000')]
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)
        self._check_roundtrip(frame, tm.assert_frame_equal)

    def test_store_hierarchical(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['foo', 'bar'])
        frame = DataFrame(np.random.randn(10, 3), index=index,
                          columns=['A', 'B', 'C'])

        self._check_roundtrip(frame, tm.assert_frame_equal)
        self._check_roundtrip(frame.T, tm.assert_frame_equal)
        self._check_roundtrip(frame['A'], tm.assert_series_equal)

        # check that the names are stored
        try:
            store = HDFStore(self.scratchpath)
            store['frame'] = frame
            recons = store['frame']
            assert(recons.index.names == ['foo', 'bar'])
        finally:
            store.close()
            os.remove(self.scratchpath)

    def test_store_index_name(self):
        df = tm.makeDataFrame()
        df.index.name = 'foo'
        try:
            store = HDFStore(self.scratchpath)
            store['frame'] = df
            recons = store['frame']
            assert(recons.index.name == 'foo')
        finally:
            store.close()
            os.remove(self.scratchpath)

    def test_store_series_name(self):
        df = tm.makeDataFrame()
        series = df['A']

        try:
            store = HDFStore(self.scratchpath)
            store['series'] = series
            recons = store['series']
            assert(recons.name == 'A')
        finally:
            store.close()
            os.remove(self.scratchpath)

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

    def test_wide_table_dups(self):
        wp = tm.makePanel()
        try:
            store = HDFStore(self.scratchpath)
            store._quiet = True
            store.put('panel', wp, table=True)
            store.put('panel', wp, table=True, append=True)
            recons = store['panel']
            tm.assert_panel_equal(recons, wp)
        finally:
            store.close()
            os.remove(self.scratchpath)

    def test_long(self):
        def _check(left, right):
            tm.assert_panel_equal(left.to_panel(), right.to_panel())

        wp = tm.makePanel()
        self._check_roundtrip(wp.to_frame(), _check)

        # empty
        self.assertRaises(ValueError, self._check_roundtrip, wp.to_frame()[:0],
                          _check)

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
            options['complib'] = _default_compressor

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
            options['complib'] = _default_compressor

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

