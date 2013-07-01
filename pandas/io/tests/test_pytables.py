import nose
import unittest
import os
import sys
import warnings

import datetime
import numpy as np

import pandas
from pandas import (Series, DataFrame, Panel, MultiIndex, bdate_range,
                    date_range, Index)
from pandas.io.pytables import (HDFStore, get_store, Term, read_hdf,
                                IncompatibilityWarning, PerformanceWarning,
                                AttributeConflictWarning)
import pandas.util.testing as tm
from pandas.tests.test_series import assert_series_equal
from pandas.tests.test_frame import assert_frame_equal
from pandas import concat, Timestamp
from pandas.util import py3compat

from numpy.testing.decorators import slow

try:
    import tables
except ImportError:
    raise nose.SkipTest('no pytables')

from distutils.version import LooseVersion

_default_compressor = LooseVersion(tables.__version__) >= '2.2' \
    and 'blosc' or 'zlib'

_multiprocess_can_split_ = False

# contextmanager to ensure the file cleanup
def safe_remove(path):
    if path is not None:
        import os
        try:
            os.remove(path)
        except:
            pass

def safe_close(store):
    try:
        if store is not None:
            store.close()
    except:
        pass

from contextlib import contextmanager

@contextmanager
def ensure_clean(path, mode='a', complevel=None, complib=None,
              fletcher32=False):
    store = HDFStore(path, mode=mode, complevel=complevel,
                     complib=complib, fletcher32=False)
    try:
        yield store
    finally:
        safe_close(store)
        if mode == 'w' or mode == 'a':
            safe_remove(path)

# set these parameters so we don't have file sharing
tables.parameters.MAX_NUMEXPR_THREADS = 1
tables.parameters.MAX_BLOSC_THREADS   = 1
tables.parameters.MAX_THREADS   = 1

def _maybe_remove(store, key):
    """For tests using tables, try removing the table to be sure there is
    no content from previous tests using the same table name."""
    try:
        store.remove(key)
    except:
        pass


class TestHDFStore(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        self.path = '__%s__.h5' % tm.rands(10)

    def tearDown(self):
        pass

    def test_factory_fun(self):
        try:
            with get_store(self.path) as tbl:
                raise ValueError('blah')
        except ValueError:
            pass
        finally:
            safe_remove(self.path)

        try:
            with get_store(self.path) as tbl:
                tbl['a'] = tm.makeDataFrame()

            with get_store(self.path) as tbl:
                self.assertEquals(len(tbl), 1)
                self.assertEquals(type(tbl['a']), DataFrame)
        finally:
            safe_remove(self.path)

    def test_conv_read_write(self):

        try:

            def roundtrip(key, obj,**kwargs):
                obj.to_hdf(self.path, key,**kwargs)
                return read_hdf(self.path, key)

            o = tm.makeTimeSeries()
            assert_series_equal(o, roundtrip('series',o))

            o = tm.makeStringSeries()
            assert_series_equal(o, roundtrip('string_series',o))

            o = tm.makeDataFrame()
            assert_frame_equal(o, roundtrip('frame',o))

            o = tm.makePanel()
            tm.assert_panel_equal(o, roundtrip('panel',o))

            # table
            df = DataFrame(dict(A=range(5), B=range(5)))
            df.to_hdf(self.path,'table',append=True)
            result = read_hdf(self.path, 'table', where = ['index>2'])
            assert_frame_equal(df[df.index>2],result)

        finally:
            safe_remove(self.path)

    def test_keys(self):

        with ensure_clean(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store['b'] = tm.makeStringSeries()
            store['c'] = tm.makeDataFrame()
            store['d'] = tm.makePanel()
            store['foo/bar'] = tm.makePanel()
            self.assertEquals(len(store), 5)
            self.assert_(set(
                    store.keys()) == set(['/a', '/b', '/c', '/d', '/foo/bar']))

    def test_repr(self):

        with ensure_clean(self.path) as store:
            repr(store)
            store['a'] = tm.makeTimeSeries()
            store['b'] = tm.makeStringSeries()
            store['c'] = tm.makeDataFrame()
            store['d'] = tm.makePanel()
            store['foo/bar'] = tm.makePanel()
            store.append('e', tm.makePanel())

            df = tm.makeDataFrame()
            df['obj1'] = 'foo'
            df['obj2'] = 'bar'
            df['bool1'] = df['A'] > 0
            df['bool2'] = df['B'] > 0
            df['bool3'] = True
            df['int1'] = 1
            df['int2'] = 2
            df['timestamp1'] = Timestamp('20010102')
            df['timestamp2'] = Timestamp('20010103')
            df['datetime1']  = datetime.datetime(2001,1,2,0,0)
            df['datetime2']  = datetime.datetime(2001,1,3,0,0)
            df.ix[3:6,['obj1']] = np.nan
            df = df.consolidate().convert_objects()
            store['df'] = df

            # make a random group in hdf space
            store._handle.createGroup(store._handle.root,'bah')

            repr(store)
            str(store)

    def test_contains(self):

        with ensure_clean(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store['b'] = tm.makeDataFrame()
            store['foo/bar'] = tm.makeDataFrame()
            self.assert_('a' in store)
            self.assert_('b' in store)
            self.assert_('c' not in store)
            self.assert_('foo/bar' in store)
            self.assert_('/foo/bar' in store)
            self.assert_('/foo/b' not in store)
            self.assert_('bar' not in store)

            # GH 2694
            warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)
            store['node())'] = tm.makeDataFrame()
            self.assert_('node())' in store)
            warnings.filterwarnings('always', category=tables.NaturalNameWarning)

    def test_versioning(self):

        with ensure_clean(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store['b'] = tm.makeDataFrame()
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df1')
            store.append('df1', df[:10])
            store.append('df1', df[10:])
            self.assert_(store.root.a._v_attrs.pandas_version == '0.10.1')
            self.assert_(store.root.b._v_attrs.pandas_version == '0.10.1')
            self.assert_(store.root.df1._v_attrs.pandas_version == '0.10.1')

            # write a file and wipe its versioning
            _maybe_remove(store, 'df2')
            store.append('df2', df)

            # this is an error because its table_type is appendable, but no version
            # info
            store.get_node('df2')._v_attrs.pandas_version = None
            self.assertRaises(Exception, store.select, 'df2')

    def test_reopen_handle(self):

        with ensure_clean(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store.open('w', warn=False)
            self.assert_(store._handle.isopen)
            self.assertEquals(len(store), 0)

    def test_flush(self):

        with ensure_clean(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store.flush()

    def test_get(self):

        with ensure_clean(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            left = store.get('a')
            right = store['a']
            tm.assert_series_equal(left, right)

            left = store.get('/a')
            right = store['/a']
            tm.assert_series_equal(left, right)

            self.assertRaises(KeyError, store.get, 'b')

    def test_getattr(self):

        with ensure_clean(self.path) as store:

            s = tm.makeTimeSeries()
            store['a'] = s

            # test attribute access
            result = store.a
            tm.assert_series_equal(result, s)
            result = getattr(store,'a')
            tm.assert_series_equal(result, s)

            df = tm.makeTimeDataFrame()
            store['df'] = df
            result = store.df
            tm.assert_frame_equal(result, df)

            # errors
            self.assertRaises(AttributeError, getattr, store, 'd')

            for x in ['mode','path','handle','complib']:
                self.assertRaises(AttributeError, getattr, store, x)

            # not stores
            for x in ['mode','path','handle','complib']:
                getattr(store,"_%s" % x)

    def test_put(self):

        with ensure_clean(self.path) as store:

            ts = tm.makeTimeSeries()
            df = tm.makeTimeDataFrame()
            store['a'] = ts
            store['b'] = df[:10]
            store['foo/bar/bah'] = df[:10]
            store['foo'] = df[:10]
            store['/foo'] = df[:10]
            store.put('c', df[:10], table=True)

            # not OK, not a table
            self.assertRaises(
                ValueError, store.put, 'b', df[10:], append=True)

            # node does not currently exist, test _is_table_type returns False in
            # this case
            # _maybe_remove(store, 'f')
            # self.assertRaises(ValueError, store.put, 'f', df[10:], append=True)

            # can't put to a table (use append instead)
            self.assertRaises(ValueError, store.put, 'c', df[10:], append=True)

            # overwrite table
            store.put('c', df[:10], table=True, append=False)
            tm.assert_frame_equal(df[:10], store['c'])

    def test_put_string_index(self):

        with ensure_clean(self.path) as store:

            index = Index(
                ["I am a very long string index: %s" % i for i in range(20)])
            s = Series(np.arange(20), index=index)
            df = DataFrame({'A': s, 'B': s})

            store['a'] = s
            tm.assert_series_equal(store['a'], s)

            store['b'] = df
            tm.assert_frame_equal(store['b'], df)

            # mixed length
            index = Index(['abcdefghijklmnopqrstuvwxyz1234567890'] + ["I am a very long string index: %s" % i for i in range(20)])
            s = Series(np.arange(21), index=index)
            df = DataFrame({'A': s, 'B': s})
            store['a'] = s
            tm.assert_series_equal(store['a'], s)

            store['b'] = df
            tm.assert_frame_equal(store['b'], df)

    def test_put_compression(self):

        with ensure_clean(self.path) as store:
            df = tm.makeTimeDataFrame()

            store.put('c', df, table=True, complib='zlib')
            tm.assert_frame_equal(store['c'], df)

            # can't compress if table=False
            self.assertRaises(ValueError, store.put, 'b', df,
                              table=False, complib='zlib')

    def test_put_compression_blosc(self):
        tm.skip_if_no_package('tables', '2.2', app='blosc support')
        df = tm.makeTimeDataFrame()

        with ensure_clean(self.path) as store:

            # can't compress if table=False
            self.assertRaises(ValueError, store.put, 'b', df,
                              table=False, complib='blosc')

            store.put('c', df, table=True, complib='blosc')
            tm.assert_frame_equal(store['c'], df)

    def test_put_integer(self):
        # non-date, non-string index
        df = DataFrame(np.random.randn(50, 100))
        self._check_roundtrip(df, tm.assert_frame_equal)

    def test_put_mixed_type(self):
        df = tm.makeTimeDataFrame()
        df['obj1'] = 'foo'
        df['obj2'] = 'bar'
        df['bool1'] = df['A'] > 0
        df['bool2'] = df['B'] > 0
        df['bool3'] = True
        df['int1'] = 1
        df['int2'] = 2
        df['timestamp1'] = Timestamp('20010102')
        df['timestamp2'] = Timestamp('20010103')
        df['datetime1'] = datetime.datetime(2001, 1, 2, 0, 0)
        df['datetime2'] = datetime.datetime(2001, 1, 3, 0, 0)
        df.ix[3:6, ['obj1']] = np.nan
        df = df.consolidate().convert_objects()

        with ensure_clean(self.path) as store:
            _maybe_remove(store, 'df')
            warnings.filterwarnings('ignore', category=PerformanceWarning)
            store.put('df',df)
            expected = store.get('df')
            tm.assert_frame_equal(expected,df)
            warnings.filterwarnings('always', category=PerformanceWarning)

    def test_append(self):

        with ensure_clean(self.path) as store:
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df1')
            store.append('df1', df[:10])
            store.append('df1', df[10:])
            tm.assert_frame_equal(store['df1'], df)

            _maybe_remove(store, 'df2')
            store.put('df2', df[:10], table=True)
            store.append('df2', df[10:])
            tm.assert_frame_equal(store['df2'], df)

            _maybe_remove(store, 'df3')
            store.append('/df3', df[:10])
            store.append('/df3', df[10:])
            tm.assert_frame_equal(store['df3'], df)

            # this is allowed by almost always don't want to do it
            warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)
            _maybe_remove(store, '/df3 foo')
            store.append('/df3 foo', df[:10])
            store.append('/df3 foo', df[10:])
            tm.assert_frame_equal(store['df3 foo'], df)
            warnings.filterwarnings('always', category=tables.NaturalNameWarning)

            # panel
            wp = tm.makePanel()
            _maybe_remove(store, 'wp1')
            store.append('wp1', wp.ix[:, :10, :])
            store.append('wp1', wp.ix[:, 10:, :])
            tm.assert_panel_equal(store['wp1'], wp)

            # ndim
            p4d = tm.makePanel4D()
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :])
            store.append('p4d', p4d.ix[:, :, 10:, :])
            tm.assert_panel4d_equal(store['p4d'], p4d)

            # test using axis labels
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=[
                    'items', 'major_axis', 'minor_axis'])
            store.append('p4d', p4d.ix[:, :, 10:, :], axes=[
                    'items', 'major_axis', 'minor_axis'])
            tm.assert_panel4d_equal(store['p4d'], p4d)

            # test using differnt number of items on each axis
            p4d2 = p4d.copy()
            p4d2['l4'] = p4d['l1']
            p4d2['l5'] = p4d['l1']
            _maybe_remove(store, 'p4d2')
            store.append(
                'p4d2', p4d2, axes=['items', 'major_axis', 'minor_axis'])
            tm.assert_panel4d_equal(store['p4d2'], p4d2)

            # test using differt order of items on the non-index axes
            _maybe_remove(store, 'wp1')
            wp_append1 = wp.ix[:, :10, :]
            store.append('wp1', wp_append1)
            wp_append2 = wp.ix[:, 10:, :].reindex(items=wp.items[::-1])
            store.append('wp1', wp_append2)
            tm.assert_panel_equal(store['wp1'], wp)

            # dtype issues - mizxed type in a single object column
            df = DataFrame(data=[[1, 2], [0, 1], [1, 2], [0, 0]])
            df['mixed_column'] = 'testing'
            df.ix[2, 'mixed_column'] = np.nan
            _maybe_remove(store, 'df')
            store.append('df', df)
            tm.assert_frame_equal(store['df'], df)

            # uints - test storage of uints
            uint_data = DataFrame({'u08' : Series(np.random.random_integers(0, high=255, size=5), dtype=np.uint8),
                                   'u16' : Series(np.random.random_integers(0, high=65535, size=5), dtype=np.uint16),
                                   'u32' : Series(np.random.random_integers(0, high=2**30, size=5), dtype=np.uint32),
                                   'u64' : Series([2**58, 2**59, 2**60, 2**61, 2**62], dtype=np.uint64)},
                                  index=np.arange(5))
            _maybe_remove(store, 'uints')
            store.append('uints', uint_data)
            tm.assert_frame_equal(store['uints'], uint_data)

            # uints - test storage of uints in indexable columns
            _maybe_remove(store, 'uints')
            store.append('uints', uint_data, data_columns=['u08','u16','u32']) # 64-bit indices not yet supported
            tm.assert_frame_equal(store['uints'], uint_data)

    def test_encoding(self):

        if sys.byteorder != 'little':
            raise nose.SkipTest('system byteorder is not little, skipping test_encoding!')

        with ensure_clean(self.path) as store:
            df = DataFrame(dict(A='foo',B='bar'),index=range(5))
            df.loc[2,'A'] = np.nan
            df.loc[3,'B'] = np.nan
            _maybe_remove(store, 'df')
            store.append('df', df, encoding='ascii')
            tm.assert_frame_equal(store['df'], df)

            expected = df.reindex(columns=['A'])
            result = store.select('df',Term('columns=A',encoding='ascii'))
            tm.assert_frame_equal(result,expected)

    def test_append_some_nans(self):

        with ensure_clean(self.path) as store:
            df = DataFrame({'A' : Series(np.random.randn(20)).astype('int32'),
                            'A1' : np.random.randn(20),
                            'A2' : np.random.randn(20),
                            'B' : 'foo', 'C' : 'bar', 'D' : Timestamp("20010101"), 'E' : datetime.datetime(2001,1,2,0,0) },
                           index=np.arange(20))
            # some nans
            _maybe_remove(store, 'df1')
            df.ix[0:15,['A1','B','D','E']] = np.nan
            store.append('df1', df[:10])
            store.append('df1', df[10:])
            tm.assert_frame_equal(store['df1'], df)

            # first column
            df1 = df.copy()
            df1.ix[:,'A1'] = np.nan
            _maybe_remove(store, 'df1')
            store.append('df1', df1[:10])
            store.append('df1', df1[10:])
            tm.assert_frame_equal(store['df1'], df1)

            # 2nd column
            df2 = df.copy()
            df2.ix[:,'A2'] = np.nan
            _maybe_remove(store, 'df2')
            store.append('df2', df2[:10])
            store.append('df2', df2[10:])
            tm.assert_frame_equal(store['df2'], df2)

            # datetimes
            df3 = df.copy()
            df3.ix[:,'E'] = np.nan
            _maybe_remove(store, 'df3')
            store.append('df3', df3[:10])
            store.append('df3', df3[10:])
            tm.assert_frame_equal(store['df3'], df3)

            ##### THIS IS A BUG, should not drop these all-nan rows
            ##### BUT need to store the index which we don't want to do....
            # nan some entire rows
            df = DataFrame({'A1' : np.random.randn(20),
                            'A2' : np.random.randn(20)},
                           index=np.arange(20))

            _maybe_remove(store, 'df4')
            df.ix[0:15,:] = np.nan
            store.append('df4', df[:10])
            store.append('df4', df[10:])
            tm.assert_frame_equal(store['df4'], df[-4:])
            self.assert_(store.get_storer('df4').nrows == 4)

            # nan some entire rows (string are still written!)
            df = DataFrame({'A1' : np.random.randn(20),
                            'A2' : np.random.randn(20),
                            'B' : 'foo', 'C' : 'bar'},
                           index=np.arange(20))

            _maybe_remove(store, 'df5')
            df.ix[0:15,:] = np.nan
            store.append('df5', df[:10])
            store.append('df5', df[10:])
            tm.assert_frame_equal(store['df5'], df)
            self.assert_(store.get_storer('df5').nrows == 20)

            # nan some entire rows (but since we have dates they are still written!)
            df = DataFrame({'A1' : np.random.randn(20),
                            'A2' : np.random.randn(20),
                            'B' : 'foo', 'C' : 'bar', 'D' : Timestamp("20010101"), 'E' : datetime.datetime(2001,1,2,0,0) },
                           index=np.arange(20))

            _maybe_remove(store, 'df6')
            df.ix[0:15,:] = np.nan
            store.append('df6', df[:10])
            store.append('df6', df[10:])
            tm.assert_frame_equal(store['df6'], df)
            self.assert_(store.get_storer('df6').nrows == 20)

    def test_append_frame_column_oriented(self):

        with ensure_clean(self.path) as store:

            # column oriented
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df1')
            store.append('df1', df.ix[:, :2], axes=['columns'])
            store.append('df1', df.ix[:, 2:])
            tm.assert_frame_equal(store['df1'], df)

            result = store.select('df1', 'columns=A')
            expected = df.reindex(columns=['A'])
            tm.assert_frame_equal(expected, result)

            # this isn't supported
            self.assertRaises(TypeError, store.select, 'df1', (
                    'columns=A', Term('index', '>', df.index[4])))

            # selection on the non-indexable
            result = store.select(
                'df1', ('columns=A', Term('index', '=', df.index[0:4])))
            expected = df.reindex(columns=['A'], index=df.index[0:4])
            tm.assert_frame_equal(expected, result)

    def test_append_with_different_block_ordering(self):

        #GH 4096; using same frames, but different block orderings
        with ensure_clean(self.path) as store:

            for i in range(10):

                df = DataFrame(np.random.randn(10,2),columns=list('AB'))
                df['index'] = range(10)
                df['index'] += i*10
                df['int64'] = Series([1]*len(df),dtype='int64')
                df['int16'] = Series([1]*len(df),dtype='int16')

                if i % 2 == 0:
                    del df['int64']
                    df['int64'] = Series([1]*len(df),dtype='int64')
                if i % 3 == 0:
                    a = df.pop('A')
                    df['A'] = a

                df.set_index('index',inplace=True)

                store.append('df',df)


    def test_ndim_indexables(self):
        """ test using ndim tables in new ways"""

        with ensure_clean(self.path) as store:

            p4d = tm.makePanel4D()

            def check_indexers(key, indexers):
                for i, idx in enumerate(indexers):
                    self.assert_(getattr(getattr(
                                store.root, key).table.description, idx)._v_pos == i)

            # append then change (will take existing schema)
            indexers = ['items', 'major_axis', 'minor_axis']

            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=indexers)
            store.append('p4d', p4d.ix[:, :, 10:, :])
            tm.assert_panel4d_equal(store.select('p4d'), p4d)
            check_indexers('p4d', indexers)

            # same as above, but try to append with differnt axes
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=indexers)
            store.append('p4d', p4d.ix[:, :, 10:, :], axes=[
                    'labels', 'items', 'major_axis'])
            tm.assert_panel4d_equal(store.select('p4d'), p4d)
            check_indexers('p4d', indexers)

            # pass incorrect number of axes
            _maybe_remove(store, 'p4d')
            self.assertRaises(ValueError, store.append, 'p4d', p4d.ix[
                    :, :, :10, :], axes=['major_axis', 'minor_axis'])

            # different than default indexables #1
            indexers = ['labels', 'major_axis', 'minor_axis']
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=indexers)
            store.append('p4d', p4d.ix[:, :, 10:, :])
            tm.assert_panel4d_equal(store['p4d'], p4d)
            check_indexers('p4d', indexers)

            # different than default indexables #2
            indexers = ['major_axis', 'labels', 'minor_axis']
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=indexers)
            store.append('p4d', p4d.ix[:, :, 10:, :])
            tm.assert_panel4d_equal(store['p4d'], p4d)
            check_indexers('p4d', indexers)

            # partial selection
            result = store.select('p4d', ['labels=l1'])
            expected = p4d.reindex(labels=['l1'])
            tm.assert_panel4d_equal(result, expected)

            # partial selection2
            result = store.select('p4d', [Term(
                        'labels=l1'), Term('items=ItemA'), Term('minor_axis=B')])
            expected = p4d.reindex(
                labels=['l1'], items=['ItemA'], minor_axis=['B'])
            tm.assert_panel4d_equal(result, expected)

            # non-existant partial selection
            result = store.select('p4d', [Term(
                        'labels=l1'), Term('items=Item1'), Term('minor_axis=B')])
            expected = p4d.reindex(labels=['l1'], items=[], minor_axis=['B'])
            tm.assert_panel4d_equal(result, expected)

    def test_append_with_strings(self):

        with ensure_clean(self.path) as store:
            wp = tm.makePanel()
            wp2 = wp.rename_axis(
                dict([(x, "%s_extra" % x) for x in wp.minor_axis]), axis=2)

            def check_col(key,name,size):
                self.assert_(getattr(store.get_storer(key).table.description,name).itemsize == size)

            store.append('s1', wp, min_itemsize=20)
            store.append('s1', wp2)
            expected = concat([wp, wp2], axis=2)
            expected = expected.reindex(minor_axis=sorted(expected.minor_axis))
            tm.assert_panel_equal(store['s1'], expected)
            check_col('s1', 'minor_axis', 20)

            # test dict format
            store.append('s2', wp, min_itemsize={'minor_axis': 20})
            store.append('s2', wp2)
            expected = concat([wp, wp2], axis=2)
            expected = expected.reindex(minor_axis=sorted(expected.minor_axis))
            tm.assert_panel_equal(store['s2'], expected)
            check_col('s2', 'minor_axis', 20)

            # apply the wrong field (similar to #1)
            store.append('s3', wp, min_itemsize={'major_axis': 20})
            self.assertRaises(ValueError, store.append, 's3', wp2)

            # test truncation of bigger strings
            store.append('s4', wp)
            self.assertRaises(ValueError, store.append, 's4', wp2)

            # avoid truncation on elements
            df = DataFrame([[123, 'asdqwerty'], [345, 'dggnhebbsdfbdfb']])
            store.append('df_big', df)
            tm.assert_frame_equal(store.select('df_big'), df)
            check_col('df_big', 'values_block_1', 15)

            # appending smaller string ok
            df2 = DataFrame([[124, 'asdqy'], [346, 'dggnhefbdfb']])
            store.append('df_big', df2)
            expected = concat([df, df2])
            tm.assert_frame_equal(store.select('df_big'), expected)
            check_col('df_big', 'values_block_1', 15)

            # avoid truncation on elements
            df = DataFrame([[123, 'asdqwerty'], [345, 'dggnhebbsdfbdfb']])
            store.append('df_big2', df, min_itemsize={'values': 50})
            tm.assert_frame_equal(store.select('df_big2'), df)
            check_col('df_big2', 'values_block_1', 50)

            # bigger string on next append
            store.append('df_new', df)
            df_new = DataFrame(
                [[124, 'abcdefqhij'], [346, 'abcdefghijklmnopqrtsuvwxyz']])
            self.assertRaises(ValueError, store.append, 'df_new', df_new)

            # with nans
            _maybe_remove(store, 'df')
            df = tm.makeTimeDataFrame()
            df['string'] = 'foo'
            df.ix[1:4, 'string'] = np.nan
            df['string2'] = 'bar'
            df.ix[4:8, 'string2'] = np.nan
            df['string3'] = 'bah'
            df.ix[1:, 'string3'] = np.nan
            store.append('df', df)
            result = store.select('df')
            tm.assert_frame_equal(result, df)

        with ensure_clean(self.path) as store:

            def check_col(key,name,size):
                self.assert_(getattr(store.get_storer(key).table.description,name).itemsize == size)

            df = DataFrame(dict(A = 'foo', B = 'bar'),index=range(10))

            # a min_itemsize that creates a data_column
            _maybe_remove(store, 'df')
            store.append('df', df, min_itemsize={'A' : 200 })
            check_col('df', 'A', 200)
            self.assert_(store.get_storer('df').data_columns == ['A'])

            # a min_itemsize that creates a data_column2
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns = ['B'], min_itemsize={'A' : 200 })
            check_col('df', 'A', 200)
            self.assert_(store.get_storer('df').data_columns == ['B','A'])

            # a min_itemsize that creates a data_column2
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns = ['B'], min_itemsize={'values' : 200 })
            check_col('df', 'B', 200)
            check_col('df', 'values_block_0', 200)
            self.assert_(store.get_storer('df').data_columns == ['B'])

            # infer the .typ on subsequent appends
            _maybe_remove(store, 'df')
            store.append('df', df[:5], min_itemsize=200)
            store.append('df', df[5:], min_itemsize=200)
            tm.assert_frame_equal(store['df'], df)

            # invalid min_itemsize keys
            df = DataFrame(['foo','foo','foo','barh','barh','barh'],columns=['A'])
            _maybe_remove(store, 'df')
            self.assertRaises(ValueError, store.append, 'df', df, min_itemsize={'foo' : 20, 'foobar' : 20})

    def test_append_with_data_columns(self):

        with ensure_clean(self.path) as store:
            df = tm.makeTimeDataFrame()
            df.loc[:,'B'].iloc[0] = 1.
            _maybe_remove(store, 'df')
            store.append('df', df[:2], data_columns=['B'])
            store.append('df', df[2:])
            tm.assert_frame_equal(store['df'], df)

            # check that we have indicies created
            assert(store._handle.root.df.table.cols.index.is_indexed is True)
            assert(store._handle.root.df.table.cols.B.is_indexed is True)

            # data column searching
            result = store.select('df', [Term('B>0')])
            expected = df[df.B > 0]
            tm.assert_frame_equal(result, expected)

            # data column searching (with an indexable and a data_columns)
            result = store.select(
                'df', [Term('B>0'), Term('index', '>', df.index[3])])
            df_new = df.reindex(index=df.index[4:])
            expected = df_new[df_new.B > 0]
            tm.assert_frame_equal(result, expected)

            # data column selection with a string data_column
            df_new = df.copy()
            df_new['string'] = 'foo'
            df_new['string'][1:4] = np.nan
            df_new['string'][5:6] = 'bar'
            _maybe_remove(store, 'df')
            store.append('df', df_new, data_columns=['string'])
            result = store.select('df', [Term('string', '=', 'foo')])
            expected = df_new[df_new.string == 'foo']
            tm.assert_frame_equal(result, expected)

            # using min_itemsize and a data column
            def check_col(key,name,size):
                self.assert_(getattr(store.get_storer(key).table.description,name).itemsize == size)

        with ensure_clean(self.path) as store:
            _maybe_remove(store, 'df')
            store.append('df', df_new, data_columns=['string'],
                         min_itemsize={'string': 30})
            check_col('df', 'string', 30)
            _maybe_remove(store, 'df')
            store.append(
                'df', df_new, data_columns=['string'], min_itemsize=30)
            check_col('df', 'string', 30)
            _maybe_remove(store, 'df')
            store.append('df', df_new, data_columns=['string'],
                         min_itemsize={'values': 30})
            check_col('df', 'string', 30)

        with ensure_clean(self.path) as store:
            df_new['string2'] = 'foobarbah'
            df_new['string_block1'] = 'foobarbah1'
            df_new['string_block2'] = 'foobarbah2'
            _maybe_remove(store, 'df')
            store.append('df', df_new, data_columns=['string', 'string2'], min_itemsize={'string': 30, 'string2': 40, 'values': 50})
            check_col('df', 'string', 30)
            check_col('df', 'string2', 40)
            check_col('df', 'values_block_1', 50)

        with ensure_clean(self.path) as store:
            # multiple data columns
            df_new = df.copy()
            df_new.loc[:,'A'].iloc[0] = 1.
            df_new.loc[:,'B'].iloc[0] = -1.
            df_new['string'] = 'foo'
            df_new['string'][1:4] = np.nan
            df_new['string'][5:6] = 'bar'
            df_new['string2'] = 'foo'
            df_new['string2'][2:5] = np.nan
            df_new['string2'][7:8] = 'bar'
            _maybe_remove(store, 'df')
            store.append(
                'df', df_new, data_columns=['A', 'B', 'string', 'string2'])
            result = store.select('df', [Term('string', '=', 'foo'), Term(
                        'string2=foo'), Term('A>0'), Term('B<0')])
            expected = df_new[(df_new.string == 'foo') & (
                    df_new.string2 == 'foo') & (df_new.A > 0) & (df_new.B < 0)]
            tm.assert_frame_equal(result, expected)

            # yield an empty frame
            result = store.select('df', [Term('string', '=', 'foo'), Term(
                        'string2=cool')])
            expected = df_new[(df_new.string == 'foo') & (
                    df_new.string2 == 'cool')]
            tm.assert_frame_equal(result, expected)

        with ensure_clean(self.path) as store:
            # doc example
            df_dc = df.copy()
            df_dc['string'] = 'foo'
            df_dc.ix[4:6, 'string'] = np.nan
            df_dc.ix[7:9, 'string'] = 'bar'
            df_dc['string2'] = 'cool'
            df_dc['datetime'] = Timestamp('20010102')
            df_dc = df_dc.convert_objects()
            df_dc.ix[3:5, ['A', 'B', 'datetime']] = np.nan

            _maybe_remove(store, 'df_dc')
            store.append('df_dc', df_dc, data_columns=['B', 'C',
                                                       'string', 'string2', 'datetime'])
            result = store.select('df_dc', [Term('B>0')])

            expected = df_dc[df_dc.B > 0]
            tm.assert_frame_equal(result, expected)

            result = store.select(
                'df_dc', ['B > 0', 'C > 0', 'string == foo'])
            expected = df_dc[(df_dc.B > 0) & (df_dc.C > 0) & (
                    df_dc.string == 'foo')]
            tm.assert_frame_equal(result, expected)

    def test_create_table_index(self):

        with ensure_clean(self.path) as store:

            def col(t,column):
                return getattr(store.get_storer(t).table.cols,column)

            # index=False
            wp = tm.makePanel()
            store.append('p5', wp, index=False)
            store.create_table_index('p5', columns=['major_axis'])
            assert(col('p5', 'major_axis').is_indexed is True)
            assert(col('p5', 'minor_axis').is_indexed is False)

            # index=True
            store.append('p5i', wp, index=True)
            assert(col('p5i', 'major_axis').is_indexed is True)
            assert(col('p5i', 'minor_axis').is_indexed is True)

            # default optlevels
            store.get_storer('p5').create_index()
            assert(col('p5', 'major_axis').index.optlevel == 6)
            assert(col('p5', 'minor_axis').index.kind == 'medium')

            # let's change the indexing scheme
            store.create_table_index('p5')
            assert(col('p5', 'major_axis').index.optlevel == 6)
            assert(col('p5', 'minor_axis').index.kind == 'medium')
            store.create_table_index('p5', optlevel=9)
            assert(col('p5', 'major_axis').index.optlevel == 9)
            assert(col('p5', 'minor_axis').index.kind == 'medium')
            store.create_table_index('p5', kind='full')
            assert(col('p5', 'major_axis').index.optlevel == 9)
            assert(col('p5', 'minor_axis').index.kind == 'full')
            store.create_table_index('p5', optlevel=1, kind='light')
            assert(col('p5', 'major_axis').index.optlevel == 1)
            assert(col('p5', 'minor_axis').index.kind == 'light')

            # data columns
            df = tm.makeTimeDataFrame()
            df['string'] = 'foo'
            df['string2'] = 'bar'
            store.append('f', df, data_columns=['string', 'string2'])
            assert(col('f', 'index').is_indexed is True)
            assert(col('f', 'string').is_indexed is True)
            assert(col('f', 'string2').is_indexed is True)

            # specify index=columns
            store.append(
                'f2', df, index=['string'], data_columns=['string', 'string2'])
            assert(col('f2', 'index').is_indexed is False)
            assert(col('f2', 'string').is_indexed is True)
            assert(col('f2', 'string2').is_indexed is False)

            # try to index a non-table
            _maybe_remove(store, 'f2')
            store.put('f2', df)
            self.assertRaises(TypeError, store.create_table_index, 'f2')

            # try to change the version supports flag
            from pandas.io import pytables
            pytables._table_supports_index = False
            self.assertRaises(Exception, store.create_table_index, 'f')

            # test out some versions
            original = tables.__version__

            for v in ['2.2', '2.2b']:
                pytables._table_mod = None
                pytables._table_supports_index = False
                tables.__version__ = v
                self.assertRaises(Exception, store.create_table_index, 'f')

            for v in ['2.3.1', '2.3.1b', '2.4dev', '2.4', original]:
                pytables._table_mod = None
                pytables._table_supports_index = False
                tables.__version__ = v
                store.create_table_index('f')
                pytables._table_mod = None
                pytables._table_supports_index = False
                tables.__version__ = original

    def test_big_table_frame(self):
        raise nose.SkipTest('no big table frame')

        # create and write a big table
        df = DataFrame(np.random.randn(2000 * 100, 100), index=range(
            2000 * 100), columns=['E%03d' % i for i in xrange(100)])
        for x in range(20):
            df['String%03d' % x] = 'string%03d' % x

        import time
        x = time.time()
        with ensure_clean(self.path,mode='w') as store:
            store.append('df', df)
            rows = store.root.df.table.nrows
            recons = store.select('df')

        print ("\nbig_table frame [%s] -> %5.2f" % (rows, time.time() - x))

    def test_big_table2_frame(self):
        # this is a really big table: 1m rows x 60 float columns, 20 string, 20 datetime
        # columns
        raise nose.SkipTest('no big table2 frame')

        # create and write a big table
        print ("\nbig_table2 start")
        import time
        start_time = time.time()
        df = DataFrame(np.random.randn(1000 * 1000, 60), index=xrange(int(
            1000 * 1000)), columns=['E%03d' % i for i in xrange(60)])
        for x in xrange(20):
            df['String%03d' % x] = 'string%03d' % x
        for x in xrange(20):
            df['datetime%03d' % x] = datetime.datetime(2001, 1, 2, 0, 0)

        print ("\nbig_table2 frame (creation of df) [rows->%s] -> %5.2f"
                    % (len(df.index), time.time() - start_time))

        def f(chunksize):
            with ensure_clean(self.path,mode='w') as store:
                store.append('df', df, chunksize=chunksize)
                r = store.root.df.table.nrows
                return r

        for c in [10000, 50000, 250000]:
            start_time = time.time()
            print ("big_table2 frame [chunk->%s]" % c)
            rows = f(c)
            print ("big_table2 frame [rows->%s,chunk->%s] -> %5.2f"
                    % (rows, c, time.time() - start_time))

    def test_big_put_frame(self):
        raise nose.SkipTest('no big put frame')

        print ("\nbig_put start")
        import time
        start_time = time.time()
        df = DataFrame(np.random.randn(1000 * 1000, 60), index=xrange(int(
            1000 * 1000)), columns=['E%03d' % i for i in xrange(60)])
        for x in xrange(20):
            df['String%03d' % x] = 'string%03d' % x
        for x in xrange(20):
            df['datetime%03d' % x] = datetime.datetime(2001, 1, 2, 0, 0)

        print ("\nbig_put frame (creation of df) [rows->%s] -> %5.2f"
                % (len(df.index), time.time() - start_time))

        with ensure_clean(self.path, mode='w') as store:
            start_time = time.time()
            store = HDFStore(fn, mode='w')
            store.put('df', df)

            print (df.get_dtype_counts())
            print ("big_put frame [shape->%s] -> %5.2f"
                    % (df.shape, time.time() - start_time))

    def test_big_table_panel(self):
        raise nose.SkipTest('no big table panel')

        # create and write a big table
        wp = Panel(
            np.random.randn(20, 1000, 1000), items=['Item%03d' % i for i in xrange(20)],
            major_axis=date_range('1/1/2000', periods=1000), minor_axis=['E%03d' % i for i in xrange(1000)])

        wp.ix[:, 100:200, 300:400] = np.nan

        for x in range(100):
            wp['String%03d'] = 'string%03d' % x

        import time
        x = time.time()


        with ensure_clean(self.path, mode='w') as store:
            store.append('wp', wp)
            rows = store.root.wp.table.nrows
            recons = store.select('wp')

        print ("\nbig_table panel [%s] -> %5.2f" % (rows, time.time() - x))

    def test_append_diff_item_order(self):

        wp = tm.makePanel()
        wp1 = wp.ix[:, :10, :]
        wp2 = wp.ix[['ItemC', 'ItemB', 'ItemA'], 10:, :]

        with ensure_clean(self.path) as store:
            store.put('panel', wp1, table=True)
            self.assertRaises(ValueError, store.put, 'panel', wp2,
                              append=True)

    def test_append_hierarchical(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['foo', 'bar'])
        df = DataFrame(np.random.randn(10, 3), index=index,
                       columns=['A', 'B', 'C'])

        with ensure_clean(self.path) as store:
            store.append('mi', df)
            result = store.select('mi')
            tm.assert_frame_equal(result, df)

            # GH 3748
            result = store.select('mi',columns=['A','B'])
            expected = df.reindex(columns=['A','B'])
            tm.assert_frame_equal(result,expected)

        with tm.ensure_clean('test.hdf') as path:
            df.to_hdf(path,'df',table=True)
            result = read_hdf(path,'df',columns=['A','B'])
            expected = df.reindex(columns=['A','B'])
            tm.assert_frame_equal(result,expected)

    def test_pass_spec_to_storer(self):

        df = tm.makeDataFrame()

        with ensure_clean(self.path) as store:
            store.put('df',df)
            self.assertRaises(TypeError, store.select, 'df', columns=['A'])
            self.assertRaises(TypeError, store.select, 'df',where=[('columns=A')])

    def test_append_misc(self):

        with ensure_clean(self.path) as store:

            # unsuported data types for non-tables
            p4d = tm.makePanel4D()
            self.assertRaises(TypeError, store.put,'p4d',p4d)

            # unsupported data type for table
            s = tm.makeStringSeries()
            self.assertRaises(TypeError, store.append,'s',s)

            # unsuported data types
            self.assertRaises(TypeError, store.put,'abc',None)
            self.assertRaises(TypeError, store.put,'abc','123')
            self.assertRaises(TypeError, store.put,'abc',123)
            self.assertRaises(TypeError, store.put,'abc',np.arange(5))

            df = tm.makeDataFrame()
            store.append('df', df, chunksize=1)
            result = store.select('df')
            tm.assert_frame_equal(result, df)

            store.append('df1', df, expectedrows=10)
            result = store.select('df1')
            tm.assert_frame_equal(result, df)

    def test_append_raise(self):

        with ensure_clean(self.path) as store:

            # test append with invalid input to get good error messages

            # list in column
            df = tm.makeDataFrame()
            df['invalid'] = [['a']] * len(df)
            self.assert_(df.dtypes['invalid'] == np.object_)
            self.assertRaises(TypeError, store.append,'df',df)

            # multiple invalid columns
            df['invalid2'] = [['a']] * len(df)
            df['invalid3'] = [['a']] * len(df)
            self.assertRaises(TypeError, store.append,'df',df)

            # datetime with embedded nans as object
            df = tm.makeDataFrame()
            s = Series(datetime.datetime(2001,1,2),index=df.index,dtype=object)
            s[0:5] = np.nan
            df['invalid'] = s
            self.assert_(df.dtypes['invalid'] == np.object_)
            self.assertRaises(TypeError, store.append,'df', df)

            # directy ndarray
            self.assertRaises(TypeError, store.append,'df',np.arange(10))

            # series directly
            self.assertRaises(TypeError, store.append,'df',Series(np.arange(10)))

            # appending an incompatbile table
            df = tm.makeDataFrame()
            store.append('df',df)

            df['foo'] = 'foo'
            self.assertRaises(ValueError, store.append,'df',df)

    def test_table_index_incompatible_dtypes(self):
        df1 = DataFrame({'a': [1, 2, 3]})
        df2 = DataFrame({'a': [4, 5, 6]},
                        index=date_range('1/1/2000', periods=3))

        with ensure_clean(self.path) as store:
            store.put('frame', df1, table=True)
            self.assertRaises(TypeError, store.put, 'frame', df2,
                              table=True, append=True)

    def test_table_values_dtypes_roundtrip(self):

        with ensure_clean(self.path) as store:
            df1 = DataFrame({'a': [1, 2, 3]}, dtype='f8')
            store.append('df_f8', df1)
            assert df1.dtypes == store['df_f8'].dtypes

            df2 = DataFrame({'a': [1, 2, 3]}, dtype='i8')
            store.append('df_i8', df2)
            assert df2.dtypes == store['df_i8'].dtypes

            # incompatible dtype
            self.assertRaises(ValueError, store.append, 'df_i8', df1)

            # check creation/storage/retrieval of float32 (a bit hacky to actually create them thought)
            df1 = DataFrame(np.array([[1],[2],[3]],dtype='f4'),columns = ['A'])
            store.append('df_f4', df1)
            assert df1.dtypes == store['df_f4'].dtypes
            assert df1.dtypes[0] == 'float32'

            # check with mixed dtypes
            df1 = DataFrame(dict([ (c,Series(np.random.randn(5),dtype=c)) for c in
                                   ['float32','float64','int32','int64','int16','int8'] ]))
            df1['string'] = 'foo'
            df1['float322'] = 1.
            df1['float322'] = df1['float322'].astype('float32')
            df1['bool']     = df1['float32'] > 0
            df1['time1']    = Timestamp('20130101')
            df1['time2']    = Timestamp('20130102')

            store.append('df_mixed_dtypes1', df1)
            result = store.select('df_mixed_dtypes1').get_dtype_counts()
            expected = Series({ 'float32' : 2, 'float64' : 1,'int32' : 1, 'bool' : 1,
                                'int16' : 1, 'int8' : 1, 'int64' : 1, 'object' : 1,
                                'datetime64[ns]' : 2})
            result.sort()
            expected.sort()
            tm.assert_series_equal(result,expected)


    def test_table_mixed_dtypes(self):

        # frame
        df = tm.makeDataFrame()
        df['obj1'] = 'foo'
        df['obj2'] = 'bar'
        df['bool1'] = df['A'] > 0
        df['bool2'] = df['B'] > 0
        df['bool3'] = True
        df['int1'] = 1
        df['int2'] = 2
        df['timestamp1'] = Timestamp('20010102')
        df['timestamp2'] = Timestamp('20010103')
        df['datetime1'] = datetime.datetime(2001, 1, 2, 0, 0)
        df['datetime2'] = datetime.datetime(2001, 1, 3, 0, 0)
        df.ix[3:6, ['obj1']] = np.nan
        df = df.consolidate().convert_objects()

        with ensure_clean(self.path) as store:
            store.append('df1_mixed', df)
            tm.assert_frame_equal(store.select('df1_mixed'), df)

        # panel
        wp = tm.makePanel()
        wp['obj1'] = 'foo'
        wp['obj2'] = 'bar'
        wp['bool1'] = wp['ItemA'] > 0
        wp['bool2'] = wp['ItemB'] > 0
        wp['int1'] = 1
        wp['int2'] = 2
        wp = wp.consolidate()

        with ensure_clean(self.path) as store:
            store.append('p1_mixed', wp)
            tm.assert_panel_equal(store.select('p1_mixed'), wp)

        # ndim
        wp = tm.makePanel4D()
        wp['obj1'] = 'foo'
        wp['obj2'] = 'bar'
        wp['bool1'] = wp['l1'] > 0
        wp['bool2'] = wp['l2'] > 0
        wp['int1'] = 1
        wp['int2'] = 2
        wp = wp.consolidate()

        with ensure_clean(self.path) as store:
            store.append('p4d_mixed', wp)
            tm.assert_panel4d_equal(store.select('p4d_mixed'), wp)

    def test_unimplemented_dtypes_table_columns(self):

        with ensure_clean(self.path) as store:

            l = [('date', datetime.date(2001, 1, 2))]

            # py3 ok for unicode
            if not py3compat.PY3:
                l.append(('unicode', u'\u03c3'))

            ### currently not supported dtypes ####
            for n, f in l:
                df = tm.makeDataFrame()
                df[n] = f
                self.assertRaises(
                    TypeError, store.append, 'df1_%s' % n, df)

        # frame
        df = tm.makeDataFrame()
        df['obj1'] = 'foo'
        df['obj2'] = 'bar'
        df['datetime1'] = datetime.date(2001, 1, 2)
        df = df.consolidate().convert_objects()

        with ensure_clean(self.path) as store:
            # this fails because we have a date in the object block......
            self.assertRaises(TypeError, store.append, 'df_unimplemented', df)

    def test_table_append_with_timezones(self):

        from datetime import timedelta

        def compare(a,b):
            tm.assert_frame_equal(a,b)

            # compare the zones on each element
            for c in a.columns:
                for i in a.index:
                    a_e = a[c][i]
                    b_e = b[c][i]
                    if not (a_e == b_e and a_e.tz == b_e.tz):
                        raise AssertionError("invalid tz comparsion [%s] [%s]" % (a_e,b_e))

        # as columns
        with ensure_clean(self.path) as store:

            _maybe_remove(store, 'df_tz')
            df = DataFrame(dict(A = [ Timestamp('20130102 2:00:00',tz='US/Eastern') + timedelta(hours=1)*i for i in range(5) ]))
            store.append('df_tz',df,data_columns=['A'])
            result = store['df_tz']
            compare(result,df)
            assert_frame_equal(result,df)

            # select with tz aware
            compare(store.select('df_tz',where=Term('A','>=',df.A[3])),df[df.A>=df.A[3]])

            _maybe_remove(store, 'df_tz')
            df = DataFrame(dict(A = Timestamp('20130102',tz='US/Eastern'), B = Timestamp('20130103',tz='US/Eastern')),index=range(5))
            store.append('df_tz',df)
            result = store['df_tz']
            compare(result,df)
            assert_frame_equal(result,df)

            _maybe_remove(store, 'df_tz')
            df = DataFrame(dict(A = Timestamp('20130102',tz='US/Eastern'), B = Timestamp('20130102',tz='EET')),index=range(5))
            self.assertRaises(TypeError, store.append, 'df_tz', df)

            # this is ok
            _maybe_remove(store, 'df_tz')
            store.append('df_tz',df,data_columns=['A','B'])
            result = store['df_tz']
            compare(result,df)
            assert_frame_equal(result,df)

            # can't append with diff timezone
            df = DataFrame(dict(A = Timestamp('20130102',tz='US/Eastern'), B = Timestamp('20130102',tz='CET')),index=range(5))
            self.assertRaises(ValueError, store.append, 'df_tz', df)

        # as index
        with ensure_clean(self.path) as store:

            # GH 4098 example
            df = DataFrame(dict(A = Series(xrange(3), index=date_range('2000-1-1',periods=3,freq='H', tz='US/Eastern'))))

            _maybe_remove(store, 'df')
            store.put('df',df)
            result = store.select('df')
            assert_frame_equal(result,df)

            _maybe_remove(store, 'df')
            store.append('df',df)
            result = store.select('df')
            assert_frame_equal(result,df)

    def test_remove(self):

        with ensure_clean(self.path) as store:

            ts = tm.makeTimeSeries()
            df = tm.makeDataFrame()
            store['a'] = ts
            store['b'] = df
            _maybe_remove(store, 'a')
            self.assertEquals(len(store), 1)
            tm.assert_frame_equal(df, store['b'])

            _maybe_remove(store, 'b')
            self.assertEquals(len(store), 0)

            # nonexistence
            self.assertRaises(KeyError, store.remove, 'a_nonexistent_store')

            # pathing
            store['a'] = ts
            store['b/foo'] = df
            _maybe_remove(store, 'foo')
            _maybe_remove(store, 'b/foo')
            self.assertEquals(len(store), 1)

            store['a'] = ts
            store['b/foo'] = df
            _maybe_remove(store, 'b')
            self.assertEquals(len(store), 1)

            # __delitem__
            store['a'] = ts
            store['b'] = df
            del store['a']
            del store['b']
            self.assertEquals(len(store), 0)

    def test_remove_where(self):

        with ensure_clean(self.path) as store:

            # non-existance
            crit1 = Term('index', '>', 'foo')
            self.assertRaises(KeyError, store.remove, 'a', [crit1])

            # try to remove non-table (with crit)
            # non-table ok (where = None)
            wp = tm.makePanel()
            store.put('wp', wp, table=True)
            store.remove('wp', [('minor_axis', ['A', 'D'])])
            rs = store.select('wp')
            expected = wp.reindex(minor_axis=['B', 'C'])
            tm.assert_panel_equal(rs, expected)

            # empty where
            _maybe_remove(store, 'wp')
            store.put('wp', wp, table=True)

            # deleted number (entire table)
            n = store.remove('wp', [])
            assert(n == 120)

            # non - empty where
            _maybe_remove(store, 'wp')
            store.put('wp', wp, table=True)
            self.assertRaises(ValueError, store.remove,
                              'wp', ['foo'])

            # selectin non-table with a where
            # store.put('wp2', wp, table=False)
            # self.assertRaises(ValueError, store.remove,
            #                  'wp2', [('column', ['A', 'D'])])

    def test_remove_crit(self):

        with ensure_clean(self.path) as store:

            wp = tm.makePanel()

            # group row removal
            date4 = wp.major_axis.take([0, 1, 2, 4, 5, 6, 8, 9, 10])
            crit4 = Term('major_axis', date4)
            store.put('wp3', wp, table=True)
            n = store.remove('wp3', where=[crit4])
            assert(n == 36)
            result = store.select('wp3')
            expected = wp.reindex(major_axis=wp.major_axis - date4)
            tm.assert_panel_equal(result, expected)

            # upper half
            store.put('wp', wp, table=True)
            date = wp.major_axis[len(wp.major_axis) // 2]

            crit1 = Term('major_axis', '>', date)
            crit2 = Term('minor_axis', ['A', 'D'])
            n = store.remove('wp', where=[crit1])

            assert(n == 56)

            n = store.remove('wp', where=[crit2])
            assert(n == 32)

            result = store['wp']
            expected = wp.truncate(after=date).reindex(minor=['B', 'C'])
            tm.assert_panel_equal(result, expected)

            # individual row elements
            store.put('wp2', wp, table=True)

            date1 = wp.major_axis[1:3]
            crit1 = Term('major_axis', date1)
            store.remove('wp2', where=[crit1])
            result = store.select('wp2')
            expected = wp.reindex(major_axis=wp.major_axis - date1)
            tm.assert_panel_equal(result, expected)

            date2 = wp.major_axis[5]
            crit2 = Term('major_axis', date2)
            store.remove('wp2', where=[crit2])
            result = store['wp2']
            expected = wp.reindex(
                major_axis=wp.major_axis - date1 - Index([date2]))
            tm.assert_panel_equal(result, expected)

            date3 = [wp.major_axis[7], wp.major_axis[9]]
            crit3 = Term('major_axis', date3)
            store.remove('wp2', where=[crit3])
            result = store['wp2']
            expected = wp.reindex(
                major_axis=wp.major_axis - date1 - Index([date2]) - Index(date3))
            tm.assert_panel_equal(result, expected)

            # corners
            store.put('wp4', wp, table=True)
            n = store.remove(
                'wp4', where=[Term('major_axis', '>', wp.major_axis[-1])])
            result = store.select('wp4')
            tm.assert_panel_equal(result, wp)

    def test_terms(self):

        with ensure_clean(self.path) as store:

            wp = tm.makePanel()
            p4d = tm.makePanel4D()
            store.put('wp', wp, table=True)
            store.put('p4d', p4d, table=True)

            # some invalid terms
            terms = [
                ['minor', ['A', 'B']],
                ['index', ['20121114']],
                ['index', ['20121114', '20121114']],
                ]
            for t in terms:
                self.assertRaises(Exception, store.select, 'wp', t)

            self.assertRaises(Exception, Term.__init__)
            self.assertRaises(Exception, Term.__init__, 'blah')
            self.assertRaises(Exception, Term.__init__, 'index')
            self.assertRaises(Exception, Term.__init__, 'index', '==')
            self.assertRaises(Exception, Term.__init__, 'index', '>', 5)

            # panel
            result = store.select('wp', [Term(
                        'major_axis<20000108'), Term('minor_axis', '=', ['A', 'B'])])
            expected = wp.truncate(after='20000108').reindex(minor=['A', 'B'])
            tm.assert_panel_equal(result, expected)

            # p4d
            result = store.select('p4d', [Term('major_axis<20000108'),
                                          Term('minor_axis', '=', ['A', 'B']),
                                          Term('items', '=', ['ItemA', 'ItemB'])])
            expected = p4d.truncate(after='20000108').reindex(
                minor=['A', 'B'], items=['ItemA', 'ItemB'])
            tm.assert_panel4d_equal(result, expected)

            # valid terms
            terms = [
                dict(field='major_axis', op='>', value='20121114'),
                ('major_axis', '20121114'),
                ('major_axis', '>', '20121114'),
                (('major_axis', ['20121114', '20121114']),),
                ('major_axis', datetime.datetime(2012, 11, 14)),
                'major_axis> 20121114',
                'major_axis >20121114',
                'major_axis > 20121114',
                (('minor_axis', ['A', 'B']),),
                (('minor_axis', ['A', 'B']),),
                ((('minor_axis', ['A', 'B']),),),
                (('items', ['ItemA', 'ItemB']),),
                ('items=ItemA'),
                ]

            for t in terms:
                store.select('wp', t)
                store.select('p4d', t)

            # valid for p4d only
            terms = [
                (('labels', '=', ['l1', 'l2']),),
                Term('labels', '=', ['l1', 'l2']),
                ]

            for t in terms:
                store.select('p4d', t)

    def test_series(self):

        s = tm.makeStringSeries()
        self._check_roundtrip(s, tm.assert_series_equal)

        ts = tm.makeTimeSeries()
        self._check_roundtrip(ts, tm.assert_series_equal)

        ts2 = Series(ts.index, Index(ts.index, dtype=object))
        self._check_roundtrip(ts2, tm.assert_series_equal)

        ts3 = Series(ts.values, Index(np.asarray(ts.index, dtype=object),
                                      dtype=object))
        self._check_roundtrip(ts3, tm.assert_series_equal)

    def test_sparse_series(self):

        s = tm.makeStringSeries()
        s[3:5] = np.nan
        ss = s.to_sparse()
        self._check_roundtrip(ss, tm.assert_series_equal,
                              check_series_type=True)

        ss2 = s.to_sparse(kind='integer')
        self._check_roundtrip(ss2, tm.assert_series_equal,
                              check_series_type=True)

        ss3 = s.to_sparse(fill_value=0)
        self._check_roundtrip(ss3, tm.assert_series_equal,
                              check_series_type=True)

    def test_sparse_frame(self):

        s = tm.makeDataFrame()
        s.ix[3:5, 1:3] = np.nan
        s.ix[8:10, -2] = np.nan
        ss = s.to_sparse()

        self._check_double_roundtrip(ss, tm.assert_frame_equal,
                                     check_frame_type=True)

        ss2 = s.to_sparse(kind='integer')
        self._check_double_roundtrip(ss2, tm.assert_frame_equal,
                                     check_frame_type=True)

        ss3 = s.to_sparse(fill_value=0)
        self._check_double_roundtrip(ss3, tm.assert_frame_equal,
                                     check_frame_type=True)

    def test_sparse_panel(self):

        items = ['x', 'y', 'z']
        p = Panel(dict((i, tm.makeDataFrame().ix[:2, :2]) for i in items))
        sp = p.to_sparse()

        self._check_double_roundtrip(sp, tm.assert_panel_equal,
                                     check_panel_type=True)

        sp2 = p.to_sparse(kind='integer')
        self._check_double_roundtrip(sp2, tm.assert_panel_equal,
                                     check_panel_type=True)

        sp3 = p.to_sparse(fill_value=0)
        self._check_double_roundtrip(sp3, tm.assert_panel_equal,
                                     check_panel_type=True)

    def test_float_index(self):

        # GH #454
        index = np.random.randn(10)
        s = Series(np.random.randn(10), index=index)
        self._check_roundtrip(s, tm.assert_series_equal)

    def test_tuple_index(self):

        # GH #492
        col = np.arange(10)
        idx = [(0., 1.), (2., 3.), (4., 5.)]
        data = np.random.randn(30).reshape((3, 10))
        DF = DataFrame(data, index=idx, columns=col)
        warnings.filterwarnings('ignore', category=PerformanceWarning)
        self._check_roundtrip(DF, tm.assert_frame_equal)
        warnings.filterwarnings('always', category=PerformanceWarning)

    def test_index_types(self):

        values = np.random.randn(2)

        func = lambda l, r: tm.assert_series_equal(l, r, True, True, True)

        warnings.filterwarnings('ignore', category=PerformanceWarning)
        ser = Series(values, [0, 'y'])
        self._check_roundtrip(ser, func)
        warnings.filterwarnings('always', category=PerformanceWarning)

        ser = Series(values, [datetime.datetime.today(), 0])
        self._check_roundtrip(ser, func)

        ser = Series(values, ['y', 0])
        self._check_roundtrip(ser, func)

        warnings.filterwarnings('ignore', category=PerformanceWarning)
        ser = Series(values, [datetime.date.today(), 'a'])
        self._check_roundtrip(ser, func)
        warnings.filterwarnings('always', category=PerformanceWarning)

        warnings.filterwarnings('ignore', category=PerformanceWarning)
        ser = Series(values, [1.23, 'b'])
        self._check_roundtrip(ser, func)
        warnings.filterwarnings('always', category=PerformanceWarning)

        ser = Series(values, [1, 1.53])
        self._check_roundtrip(ser, func)

        ser = Series(values, [1, 5])
        self._check_roundtrip(ser, func)

        ser = Series(values, [datetime.datetime(
                    2012, 1, 1), datetime.datetime(2012, 1, 2)])
        self._check_roundtrip(ser, func)

    def test_timeseries_preepoch(self):

        if sys.version_info[0] == 2 and sys.version_info[1] < 7:
            raise nose.SkipTest

        dr = bdate_range('1/1/1940', '1/1/1960')
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

        with ensure_clean(self.path) as store:
            # not consolidated
            df['foo'] = np.random.randn(len(df))
            store['df'] = df
            recons = store['df']
            self.assert_(recons._data.is_consolidated())

        # empty
        self._check_roundtrip(df[:0], tm.assert_frame_equal)

    def test_empty_series_frame(self):
        s0 = Series()
        s1 = Series(name='myseries')
        df0 = DataFrame()
        df1 = DataFrame(index=['a', 'b', 'c'])
        df2 = DataFrame(columns=['d', 'e', 'f'])

        self._check_roundtrip(s0, tm.assert_series_equal)
        self._check_roundtrip(s1, tm.assert_series_equal)
        self._check_roundtrip(df0, tm.assert_frame_equal)
        self._check_roundtrip(df1, tm.assert_frame_equal)
        self._check_roundtrip(df2, tm.assert_frame_equal)

    def test_can_serialize_dates(self):

        rng = [x.date() for x in bdate_range('1/1/2000', '1/30/2000')]
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        self._check_roundtrip(frame, tm.assert_frame_equal)

    def test_timezones(self):
        rng = date_range('1/1/2000', '1/30/2000', tz='US/Eastern')
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        with ensure_clean(self.path) as store:
            store['frame'] = frame
            recons = store['frame']
            self.assert_(recons.index.equals(rng))
            self.assertEquals(rng.tz, recons.index.tz)

    def test_fixed_offset_tz(self):
        rng = date_range('1/1/2000 00:00:00-07:00', '1/30/2000 00:00:00-07:00')
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        with ensure_clean(self.path) as store:
            store['frame'] = frame
            recons = store['frame']
            self.assert_(recons.index.equals(rng))
            self.assertEquals(rng.tz, recons.index.tz)

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
        with ensure_clean(self.path) as store:
            store['frame'] = frame
            recons = store['frame']
            assert(recons.index.names == ['foo', 'bar'])

    def test_store_index_name(self):
        df = tm.makeDataFrame()
        df.index.name = 'foo'

        with ensure_clean(self.path) as store:
            store['frame'] = df
            recons = store['frame']
            assert(recons.index.name == 'foo')

    def test_store_series_name(self):
        df = tm.makeDataFrame()
        series = df['A']

        with ensure_clean(self.path) as store:
            store['series'] = series
            recons = store['series']
            assert(recons.name == 'A')

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

        with ensure_clean(self.path) as store:
            store['obj'] = df1
            tm.assert_frame_equal(store['obj'], df1)
            store['obj'] = df2
            tm.assert_frame_equal(store['obj'], df2)

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
        with ensure_clean(self.path) as store:
            store._quiet = True
            store.put('panel', wp, table=True)
            store.put('panel', wp, table=True, append=True)
            recons = store['panel']
            tm.assert_panel_equal(recons, wp)

    def test_long(self):
        def _check(left, right):
            tm.assert_panel_equal(left.to_panel(), right.to_panel())

        wp = tm.makePanel()
        self._check_roundtrip(wp.to_frame(), _check)

        # empty
        # self._check_roundtrip(wp.to_frame()[:0], _check)

    def test_longpanel(self):
        pass

    def test_overwrite_node(self):

        with ensure_clean(self.path) as store:
            store['a'] = tm.makeTimeDataFrame()
            ts = tm.makeTimeSeries()
            store['a'] = ts

            tm.assert_series_equal(store['a'], ts)

    def test_sparse_with_compression(self):

        # GH 2931

        # make sparse dataframe
        df = DataFrame(np.random.binomial(n=1, p=.01, size=(1e3, 10))).to_sparse(fill_value=0)

        # case 1: store uncompressed
        self._check_double_roundtrip(df, tm.assert_frame_equal,
                                     compression = False,
                                     check_frame_type=True)

        # case 2: store compressed (works)
        self._check_double_roundtrip(df, tm.assert_frame_equal,
                                     compression = 'zlib',
                                     check_frame_type=True)

        # set one series to be completely sparse
        df[0] = np.zeros(1e3)

        # case 3: store df with completely sparse series uncompressed
        self._check_double_roundtrip(df, tm.assert_frame_equal,
                                     compression = False,
                                     check_frame_type=True)

        # case 4: try storing df with completely sparse series compressed (fails)
        self._check_double_roundtrip(df, tm.assert_frame_equal,
                                     compression = 'zlib',
                                     check_frame_type=True)

    def test_select(self):
        wp = tm.makePanel()

        with ensure_clean(self.path) as store:

            # put/select ok
            _maybe_remove(store, 'wp')
            store.put('wp', wp, table=True)
            store.select('wp')

            # non-table ok (where = None)
            _maybe_remove(store, 'wp')
            store.put('wp2', wp, table=False)
            store.select('wp2')

            # selection on the non-indexable with a large number of columns
            wp = Panel(
                np.random.randn(100, 100, 100), items=['Item%03d' % i for i in xrange(100)],
                major_axis=date_range('1/1/2000', periods=100), minor_axis=['E%03d' % i for i in xrange(100)])

            _maybe_remove(store, 'wp')
            store.append('wp', wp)
            items = ['Item%03d' % i for i in xrange(80)]
            result = store.select('wp', Term('items', items))
            expected = wp.reindex(items=items)
            tm.assert_panel_equal(expected, result)

            # selectin non-table with a where
            # self.assertRaises(ValueError, store.select,
            #                  'wp2', ('column', ['A', 'D']))

            # select with columns=
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df')
            store.append('df', df)
            result = store.select('df', columns=['A', 'B'])
            expected = df.reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # equivalentsly
            result = store.select('df', [('columns', ['A', 'B'])])
            expected = df.reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # with a data column
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['A'])
            result = store.select('df', ['A > 0'], columns=['A', 'B'])
            expected = df[df.A > 0].reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # all a data columns
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=True)
            result = store.select('df', ['A > 0'], columns=['A', 'B'])
            expected = df[df.A > 0].reindex(columns=['A', 'B'])
            tm.assert_frame_equal(expected, result)

            # with a data column, but different columns
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['A'])
            result = store.select('df', ['A > 0'], columns=['C', 'D'])
            expected = df[df.A > 0].reindex(columns=['C', 'D'])
            tm.assert_frame_equal(expected, result)

    def test_select_dtypes(self):

        with ensure_clean(self.path) as store:

            # with a Timestamp data column (GH #2637)
            df = DataFrame(dict(ts=bdate_range('2012-01-01', periods=300), A=np.random.randn(300)))
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['ts', 'A'])
            result = store.select('df', [Term('ts', '>=', Timestamp('2012-02-01'))])
            expected = df[df.ts >= Timestamp('2012-02-01')]
            tm.assert_frame_equal(expected, result)

            # bool columns (GH #2849)
            df = DataFrame(np.random.randn(5,2), columns =['A','B'])
            df['object'] = 'foo'
            df.ix[4:5,'object'] = 'bar'
            df['bool'] = df['A'] > 0
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns = True)

            expected = df[df.bool == True].reindex(columns=['A','bool'])
            for v in [True,'true',1]:
                result = store.select('df', Term('bool == %s' % str(v)), columns = ['A','bool'])
                tm.assert_frame_equal(expected, result)

            expected = df[df.bool == False ].reindex(columns=['A','bool'])
            for v in [False,'false',0]:
                result = store.select('df', Term('bool == %s' % str(v)), columns = ['A','bool'])
                tm.assert_frame_equal(expected, result)

            # integer index
            df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
            _maybe_remove(store, 'df_int')
            store.append('df_int', df)
            result = store.select(
                'df_int', [Term("index<10"), Term("columns", "=", ["A"])])
            expected = df.reindex(index=list(df.index)[0:10],columns=['A'])
            tm.assert_frame_equal(expected, result)

            # float index
            df = DataFrame(dict(A=np.random.rand(
                        20), B=np.random.rand(20), index=np.arange(20, dtype='f8')))
            _maybe_remove(store, 'df_float')
            store.append('df_float', df)
            result = store.select(
                'df_float', [Term("index<10.0"), Term("columns", "=", ["A"])])
            expected = df.reindex(index=list(df.index)[0:10],columns=['A'])
            tm.assert_frame_equal(expected, result)

    def test_select_with_many_inputs(self):

        with ensure_clean(self.path) as store:

            df = DataFrame(dict(ts=bdate_range('2012-01-01', periods=300),
                                A=np.random.randn(300),
                                B=range(300),
                                users = ['a']*50 + ['b']*50 + ['c']*100 + ['a%03d' % i for i in range(100)]))
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['ts', 'A', 'B', 'users'])

            # regular select
            result = store.select('df', [Term('ts', '>=', Timestamp('2012-02-01'))])
            expected = df[df.ts >= Timestamp('2012-02-01')]
            tm.assert_frame_equal(expected, result)

            # small selector
            result = store.select('df', [Term('ts', '>=', Timestamp('2012-02-01')),Term('users',['a','b','c'])])
            expected = df[ (df.ts >= Timestamp('2012-02-01')) & df.users.isin(['a','b','c']) ]
            tm.assert_frame_equal(expected, result)

            # big selector along the columns
            selector = [ 'a','b','c' ] + [ 'a%03d' % i for i in xrange(60) ]
            result = store.select('df', [Term('ts', '>=', Timestamp('2012-02-01')),Term('users',selector)])
            expected = df[ (df.ts >= Timestamp('2012-02-01')) & df.users.isin(selector) ]
            tm.assert_frame_equal(expected, result)

            selector = range(100,200)
            result = store.select('df', [Term('B', selector)])
            expected = df[ df.B.isin(selector) ]
            tm.assert_frame_equal(expected, result)
            self.assert_(len(result) == 100)

            # big selector along the index
            selector = Index(df.ts[0:100].values)
            result  = store.select('df', [Term('ts', selector)])
            expected = df[ df.ts.isin(selector.values) ]
            tm.assert_frame_equal(expected, result)
            self.assert_(len(result) == 100)

    def test_select_iterator(self):

        # single table
        with ensure_clean(self.path) as store:

            df = tm.makeTimeDataFrame(500)
            _maybe_remove(store, 'df')
            store.append('df', df)

            expected = store.select('df')

            results = []
            for s in store.select('df',iterator=True):
                results.append(s)
            result = concat(results)
            tm.assert_frame_equal(expected, result)
            results = []
            for s in store.select('df',chunksize=100):
                results.append(s)
            self.assert_(len(results) == 5)
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            results = []
            for s in store.select('df',chunksize=150):
                results.append(s)
            result = concat(results)
            tm.assert_frame_equal(result, expected)

        with tm.ensure_clean(self.path) as path:

            df = tm.makeTimeDataFrame(500)
            df.to_hdf(path,'df_non_table')
            self.assertRaises(TypeError, read_hdf, path,'df_non_table',chunksize=100)
            self.assertRaises(TypeError, read_hdf, path,'df_non_table',iterator=True)

        with tm.ensure_clean(self.path) as path:

            df = tm.makeTimeDataFrame(500)
            df.to_hdf(path,'df',table=True)

            results = []
            for x in read_hdf(path,'df',chunksize=100):
                results.append(x)

            self.assert_(len(results) == 5)
            result = concat(results)
            tm.assert_frame_equal(result, df)
            tm.assert_frame_equal(result, read_hdf(path,'df'))

        # multiple

        with ensure_clean(self.path) as store:

            df1 = tm.makeTimeDataFrame(500)
            store.append('df1',df1,data_columns=True)
            df2 = tm.makeTimeDataFrame(500).rename(columns=lambda x: "%s_2" % x)
            df2['foo'] = 'bar'
            store.append('df2',df2)

            df = concat([df1, df2], axis=1)

            # full selection
            expected = store.select_as_multiple(
                ['df1', 'df2'], selector='df1')
            results = []
            for s in store.select_as_multiple(
                ['df1', 'df2'], selector='df1', chunksize=150):
                results.append(s)
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            # where selection
            #expected = store.select_as_multiple(
            #    ['df1', 'df2'], where= Term('A>0'), selector='df1')
            #results = []
            #for s in store.select_as_multiple(
            #    ['df1', 'df2'], where= Term('A>0'), selector='df1', chunksize=25):
            #    results.append(s)
            #result = concat(results)
            #tm.assert_frame_equal(expected, result)

    def test_retain_index_attributes(self):

        # GH 3499, losing frequency info on index recreation
        df = DataFrame(dict(A = Series(xrange(3),
                                       index=date_range('2000-1-1',periods=3,freq='H'))))

        with ensure_clean(self.path) as store:
            _maybe_remove(store,'data')
            store.put('data', df, table=True)

            result = store.get('data')
            tm.assert_frame_equal(df,result)

            for attr in ['freq','tz','name']:
                for idx in ['index','columns']:
                    self.assert_(getattr(getattr(df,idx),attr,None) == getattr(getattr(result,idx),attr,None))


            # try to append a table with a different frequency
            warnings.filterwarnings('ignore', category=AttributeConflictWarning)
            df2 = DataFrame(dict(A = Series(xrange(3),
                                            index=date_range('2002-1-1',periods=3,freq='D'))))
            store.append('data',df2)
            warnings.filterwarnings('always', category=AttributeConflictWarning)

            self.assert_(store.get_storer('data').info['index']['freq'] is None)

            # this is ok
            _maybe_remove(store,'df2')
            df2 = DataFrame(dict(A = Series(xrange(3),
                                            index=[Timestamp('20010101'),Timestamp('20010102'),Timestamp('20020101')])))
            store.append('df2',df2)
            df3 = DataFrame(dict(A = Series(xrange(3),index=date_range('2002-1-1',periods=3,freq='D'))))
            store.append('df2',df3)

    def test_retain_index_attributes2(self):

        with tm.ensure_clean(self.path) as path:

            warnings.filterwarnings('ignore', category=AttributeConflictWarning)

            df  = DataFrame(dict(A = Series(xrange(3), index=date_range('2000-1-1',periods=3,freq='H'))))
            df.to_hdf(path,'data',mode='w',append=True)
            df2 = DataFrame(dict(A = Series(xrange(3), index=date_range('2002-1-1',periods=3,freq='D'))))
            df2.to_hdf(path,'data',append=True)

            idx = date_range('2000-1-1',periods=3,freq='H')
            idx.name = 'foo'
            df  = DataFrame(dict(A = Series(xrange(3), index=idx)))
            df.to_hdf(path,'data',mode='w',append=True)
            self.assert_(read_hdf(path,'data').index.name == 'foo')

            idx2 = date_range('2001-1-1',periods=3,freq='H')
            idx2.name = 'bar'
            df2 = DataFrame(dict(A = Series(xrange(3), index=idx2)))
            df2.to_hdf(path,'data',append=True)
            self.assert_(read_hdf(path,'data').index.name is None)

            warnings.filterwarnings('always', category=AttributeConflictWarning)

    def test_panel_select(self):

        wp = tm.makePanel()

        with ensure_clean(self.path) as store:
            store.put('wp', wp, table=True)
            date = wp.major_axis[len(wp.major_axis) // 2]

            crit1 = ('major_axis', '>=', date)
            crit2 = ('minor_axis', '=', ['A', 'D'])

            result = store.select('wp', [crit1, crit2])
            expected = wp.truncate(before=date).reindex(minor=['A', 'D'])
            tm.assert_panel_equal(result, expected)

            result = store.select(
                'wp', ['major_axis>=20000124', ('minor_axis', '=', ['A', 'B'])])
            expected = wp.truncate(before='20000124').reindex(minor=['A', 'B'])
            tm.assert_panel_equal(result, expected)

    def test_frame_select(self):

        df = tm.makeTimeDataFrame()

        with ensure_clean(self.path) as store:
            store.put('frame', df, table=True)
            date = df.index[len(df) // 2]

            crit1 = ('index', '>=', date)
            crit2 = ('columns', ['A', 'D'])
            crit3 = ('columns', 'A')

            result = store.select('frame', [crit1, crit2])
            expected = df.ix[date:, ['A', 'D']]
            tm.assert_frame_equal(result, expected)

            result = store.select('frame', [crit3])
            expected = df.ix[:, ['A']]
            tm.assert_frame_equal(result, expected)

            # invalid terms
            df = tm.makeTimeDataFrame()
            store.append('df_time', df)
            self.assertRaises(
                ValueError, store.select, 'df_time', [Term("index>0")])

            # can't select if not written as table
            # store['frame'] = df
            # self.assertRaises(ValueError, store.select,
            #                  'frame', [crit1, crit2])

    def test_string_select(self):

        # GH 2973

        df = tm.makeTimeDataFrame()

        with ensure_clean(self.path) as store:


            # test string ==/!=

            df['x'] = 'none'
            df.ix[2:7,'x'] = ''

            store.append('df',df,data_columns=['x'])

            result = store.select('df',Term('x=none'))
            expected = df[df.x == 'none']
            assert_frame_equal(result,expected)

            result = store.select('df',Term('x!=none'))
            expected = df[df.x != 'none']
            assert_frame_equal(result,expected)

            df2 = df.copy()
            df2.x[df2.x==''] = np.nan

            from pandas import isnull
            store.append('df2',df2,data_columns=['x'])
            result = store.select('df2',Term('x!=none'))
            expected = df2[isnull(df2.x)]
            assert_frame_equal(result,expected)


            # int ==/!=
            df['int'] = 1
            df.ix[2:7,'int'] = 2

            store.append('df3',df,data_columns=['int'])

            result = store.select('df3',Term('int=2'))
            expected = df[df.int==2]
            assert_frame_equal(result,expected)

            result = store.select('df3',Term('int!=2'))
            expected = df[df.int!=2]
            assert_frame_equal(result,expected)


    def test_read_column(self):

        df = tm.makeTimeDataFrame()

        with ensure_clean(self.path) as store:
            _maybe_remove(store, 'df')
            store.append('df', df)

            # error
            self.assertRaises(KeyError, store.select_column, 'df', 'foo')

            def f():
                store.select_column('df', 'index', where = ['index>5'])
            self.assertRaises(Exception, f)

            # valid
            result = store.select_column('df', 'index')
            tm.assert_almost_equal(result.values, Series(df.index).values)
            self.assert_(isinstance(result,Series))

            # not a data indexable column
            self.assertRaises(
                ValueError, store.select_column, 'df', 'values_block_0')

            # a data column
            df2 = df.copy()
            df2['string'] = 'foo'
            store.append('df2', df2, data_columns=['string'])
            result = store.select_column('df2', 'string')
            tm.assert_almost_equal(result.values, df2['string'].values)

            # a data column with NaNs, result excludes the NaNs
            df3 = df.copy()
            df3['string'] = 'foo'
            df3.ix[4:6, 'string'] = np.nan
            store.append('df3', df3, data_columns=['string'])
            result = store.select_column('df3', 'string')
            tm.assert_almost_equal(result.values, df3['string'].values)

    def test_coordinates(self):
        df = tm.makeTimeDataFrame()

        with ensure_clean(self.path) as store:

            _maybe_remove(store, 'df')
            store.append('df', df)

            # all
            c = store.select_as_coordinates('df')
            assert((c.values == np.arange(len(df.index))).all() == True)

            # get coordinates back & test vs frame
            _maybe_remove(store, 'df')

            df = DataFrame(dict(A=range(5), B=range(5)))
            store.append('df', df)
            c = store.select_as_coordinates('df', ['index<3'])
            assert((c.values == np.arange(3)).all() == True)
            result = store.select('df', where=c)
            expected = df.ix[0:2, :]
            tm.assert_frame_equal(result, expected)

            c = store.select_as_coordinates('df', ['index>=3', 'index<=4'])
            assert((c.values == np.arange(2) + 3).all() == True)
            result = store.select('df', where=c)
            expected = df.ix[3:4, :]
            tm.assert_frame_equal(result, expected)

            # multiple tables
            _maybe_remove(store, 'df1')
            _maybe_remove(store, 'df2')
            df1 = tm.makeTimeDataFrame()
            df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
            store.append('df1', df1, data_columns=['A', 'B'])
            store.append('df2', df2)

            c = store.select_as_coordinates('df1', ['A>0', 'B>0'])
            df1_result = store.select('df1', c)
            df2_result = store.select('df2', c)
            result = concat([df1_result, df2_result], axis=1)

            expected = concat([df1, df2], axis=1)
            expected = expected[(expected.A > 0) & (expected.B > 0)]
            tm.assert_frame_equal(result, expected)

    def test_append_to_multiple(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
        df2['foo'] = 'bar'
        df = concat([df1, df2], axis=1)

        with ensure_clean(self.path) as store:

            # exceptions
            self.assertRaises(ValueError, store.append_to_multiple,
                              {'df1': ['A', 'B'], 'df2': None}, df, selector='df3')
            self.assertRaises(ValueError, store.append_to_multiple,
                              {'df1': None, 'df2': None}, df, selector='df3')
            self.assertRaises(
                ValueError, store.append_to_multiple, 'df1', df, 'df1')

            # regular operation
            store.append_to_multiple(
                {'df1': ['A', 'B'], 'df2': None}, df, selector='df1')
            result = store.select_as_multiple(
                ['df1', 'df2'], where=['A>0', 'B>0'], selector='df1')
            expected = df[(df.A > 0) & (df.B > 0)]
            tm.assert_frame_equal(result, expected)

    def test_select_as_multiple(self):

        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
        df2['foo'] = 'bar'

        with ensure_clean(self.path) as store:

            # no tables stored
            self.assertRaises(Exception, store.select_as_multiple,
                              None, where=['A>0', 'B>0'], selector='df1')

            store.append('df1', df1, data_columns=['A', 'B'])
            store.append('df2', df2)

            # exceptions
            self.assertRaises(Exception, store.select_as_multiple,
                              None, where=['A>0', 'B>0'], selector='df1')
            self.assertRaises(Exception, store.select_as_multiple,
                              [None], where=['A>0', 'B>0'], selector='df1')
            self.assertRaises(TypeError, store.select_as_multiple,
                              ['df1','df3'], where=['A>0', 'B>0'], selector='df1')
            self.assertRaises(KeyError, store.select_as_multiple,
                              ['df3'], where=['A>0', 'B>0'], selector='df1')
            self.assertRaises(ValueError, store.select_as_multiple,
                              ['df1','df2'], where=['A>0', 'B>0'], selector='df4')

            # default select
            result = store.select('df1', ['A>0', 'B>0'])
            expected = store.select_as_multiple(
                ['df1'], where=['A>0', 'B>0'], selector='df1')
            tm.assert_frame_equal(result, expected)
            expected = store.select_as_multiple(
                'df1', where=['A>0', 'B>0'], selector='df1')
            tm.assert_frame_equal(result, expected)

            # multiple
            result = store.select_as_multiple(
                ['df1', 'df2'], where=['A>0', 'B>0'], selector='df1')
            expected = concat([df1, df2], axis=1)
            expected = expected[(expected.A > 0) & (expected.B > 0)]
            tm.assert_frame_equal(result, expected)

            # multiple (diff selector)
            try:
                result = store.select_as_multiple(['df1', 'df2'], where=[Term(
                            'index', '>', df2.index[4])], selector='df2')
                expected = concat([df1, df2], axis=1)
                expected = expected[5:]
                tm.assert_frame_equal(result, expected)
            except (Exception), detail:
                print ("error in select_as_multiple %s" % str(detail))
                print ("store: %s" % store)
                print ("df1: %s" % df1)
                print ("df2: %s" % df2)


            # test excpection for diff rows
            store.append('df3', tm.makeTimeDataFrame(nper=50))
            self.assertRaises(ValueError, store.select_as_multiple,
                              ['df1','df3'], where=['A>0', 'B>0'], selector='df1')

    def test_start_stop(self):

        with ensure_clean(self.path) as store:

            df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
            store.append('df', df)

            result = store.select(
                'df', [Term("columns", "=", ["A"])], start=0, stop=5)
            expected = df.ix[0:4, ['A']]
            tm.assert_frame_equal(result, expected)

            # out of range
            result = store.select(
                'df', [Term("columns", "=", ["A"])], start=30, stop=40)
            assert(len(result) == 0)
            assert(type(result) == DataFrame)

    def test_select_filter_corner(self):

        df = DataFrame(np.random.randn(50, 100))
        df.index = ['%.3d' % c for c in df.index]
        df.columns = ['%.3d' % c for c in df.columns]

        with ensure_clean(self.path) as store:
            store.put('frame', df, table=True)

            crit = Term('columns', df.columns[:75])
            result = store.select('frame', [crit])
            tm.assert_frame_equal(result, df.ix[:, df.columns[:75]])

    def _check_roundtrip(self, obj, comparator, compression=False, **kwargs):

        options = {}
        if compression:
            options['complib'] = _default_compressor

        with ensure_clean(self.path, 'w', **options) as store:
            store['obj'] = obj
            retrieved = store['obj']
            comparator(retrieved, obj, **kwargs)

    def _check_double_roundtrip(self, obj, comparator, compression=False,
                                **kwargs):
        options = {}
        if compression:
            options['complib'] = compression or _default_compressor

        with ensure_clean(self.path, 'w', **options) as store:
            store['obj'] = obj
            retrieved = store['obj']
            comparator(retrieved, obj, **kwargs)
            store['obj'] = retrieved
            again = store['obj']
            comparator(again, obj, **kwargs)


    def _check_roundtrip_table(self, obj, comparator, compression=False):
        options = {}
        if compression:
            options['complib'] = _default_compressor

        with ensure_clean(self.path, 'w', **options) as store:
            store.put('obj', obj, table=True)
            retrieved = store['obj']
            # sorted_obj = _test_sort(obj)
            comparator(retrieved, obj)

    def test_pytables_native_read(self):

        try:
            store = HDFStore(tm.get_data_path('legacy_hdf/pytables_native.h5'), 'r')
            d2 = store['detector/readout']
        finally:
            safe_close(store)

        try:
            store = HDFStore(tm.get_data_path('legacy_hdf/pytables_native2.h5'), 'r')
            str(store)
            d1 = store['detector']
        finally:
            safe_close(store)

    def test_legacy_read(self):
        try:
            store = HDFStore(tm.get_data_path('legacy_hdf/legacy.h5'), 'r')
            store['a']
            store['b']
            store['c']
            store['d']
        finally:
            safe_close(store)

    def test_legacy_table_read(self):
        # legacy table types
        try:
            store = HDFStore(tm.get_data_path('legacy_hdf/legacy_table.h5'), 'r')
            store.select('df1')
            store.select('df2')
            store.select('wp1')

            # force the frame
            store.select('df2', typ='legacy_frame')

            # old version warning
            warnings.filterwarnings('ignore', category=IncompatibilityWarning)
            self.assertRaises(
                Exception, store.select, 'wp1', Term('minor_axis', '=', 'B'))

            df2 = store.select('df2')
            store.select('df2', Term('index', '>', df2.index[2]))
            warnings.filterwarnings('always', category=IncompatibilityWarning)

        finally:
            safe_close(store)

    def test_legacy_0_10_read(self):
        # legacy from 0.10
        try:
            store = HDFStore(tm.get_data_path('legacy_hdf/legacy_0.10.h5'), 'r')
            str(store)
            for k in store.keys():
                store.select(k)
        finally:
            safe_close(store)

    def test_legacy_0_11_read(self):
        # legacy from 0.11
        try:
            store = HDFStore(tm.get_data_path('legacy_hdf/legacy_table_0.11.h5'), 'r')
            str(store)
            df = store.select('df')
            df1 = store.select('df1')
            mi = store.select('mi')
        finally:
            safe_close(store)

    def test_copy(self):

        def do_copy(f = None, new_f = None, keys = None, propindexes = True, **kwargs):
            try:
                import os

                if f is None:
                    f = tm.get_data_path('legacy_hdf/legacy_0.10.h5')


                store = HDFStore(f, 'r')

                if new_f is None:
                    import tempfile
                    new_f = tempfile.mkstemp()[1]

                tstore = store.copy(new_f, keys = keys, propindexes = propindexes, **kwargs)

                # check keys
                if keys is None:
                    keys = store.keys()
                self.assert_(set(keys) == set(tstore.keys()))

                # check indicies & nrows
                for k in tstore.keys():
                    if tstore.get_storer(k).is_table:
                        new_t = tstore.get_storer(k)
                        orig_t = store.get_storer(k)

                        self.assert_(orig_t.nrows == new_t.nrows)

                        # check propindixes
                        if propindexes:
                            for a in orig_t.axes:
                                if a.is_indexed:
                                    self.assert_(new_t[a.name].is_indexed == True)

            finally:
                safe_close(store)
                safe_close(tstore)
                safe_remove(new_f)

        do_copy()
        do_copy(keys = ['/a','/b','/df1_mixed'])
        do_copy(propindexes = False)

        # new table
        df = tm.makeDataFrame()

        try:
            st = HDFStore(self.path)
            st.append('df', df, data_columns = ['A'])
            st.close()
            do_copy(f = self.path)
            do_copy(f = self.path, propindexes = False)
        finally:
            safe_remove(self.path)

    def test_legacy_table_write(self):
        raise nose.SkipTest

        store = HDFStore(tm.get_data_path('legacy_hdf/legacy_table_%s.h5' % pandas.__version__), 'a')

        df = tm.makeDataFrame()
        wp = tm.makePanel()

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['foo', 'bar'])
        df = DataFrame(np.random.randn(10, 3), index=index,
                       columns=['A', 'B', 'C'])
        store.append('mi', df)

        df = DataFrame(dict(A = 'foo', B = 'bar'),index=range(10))
        store.append('df', df, data_columns = ['B'], min_itemsize={'A' : 200 })

        store.close()

    def test_store_datetime_fractional_secs(self):

        with ensure_clean(self.path) as store:
            dt = datetime.datetime(2012, 1, 2, 3, 4, 5, 123456)
            series = Series([0], [dt])
            store['a'] = series
            self.assertEquals(store['a'].index[0], dt)

    def test_tseries_indices_series(self):

        with ensure_clean(self.path) as store:
            idx = tm.makeDateIndex(10)
            ser = Series(np.random.randn(len(idx)), idx)
            store['a'] = ser
            result = store['a']

            assert_series_equal(result, ser)
            self.assertEquals(type(result.index), type(ser.index))
            self.assertEquals(result.index.freq, ser.index.freq)

            idx = tm.makePeriodIndex(10)
            ser = Series(np.random.randn(len(idx)), idx)
            store['a'] = ser
            result = store['a']

            assert_series_equal(result, ser)
            self.assertEquals(type(result.index), type(ser.index))
            self.assertEquals(result.index.freq, ser.index.freq)

    def test_tseries_indices_frame(self):

        with ensure_clean(self.path) as store:
            idx = tm.makeDateIndex(10)
            df = DataFrame(np.random.randn(len(idx), 3), index=idx)
            store['a'] = df
            result = store['a']

            assert_frame_equal(result, df)
            self.assertEquals(type(result.index), type(df.index))
            self.assertEquals(result.index.freq, df.index.freq)

            idx = tm.makePeriodIndex(10)
            df = DataFrame(np.random.randn(len(idx), 3), idx)
            store['a'] = df
            result = store['a']

            assert_frame_equal(result, df)
            self.assertEquals(type(result.index), type(df.index))
            self.assertEquals(result.index.freq, df.index.freq)

    def test_unicode_index(self):

        unicode_values = [u'\u03c3', u'\u03c3\u03c3']
        warnings.filterwarnings('ignore', category=PerformanceWarning)
        s = Series(np.random.randn(len(unicode_values)), unicode_values)
        self._check_roundtrip(s, tm.assert_series_equal)
        warnings.filterwarnings('always', category=PerformanceWarning)

    def test_store_datetime_mixed(self):

        df = DataFrame(
            {'a': [1, 2, 3], 'b': [1., 2., 3.], 'c': ['a', 'b', 'c']})
        ts = tm.makeTimeSeries()
        df['d'] = ts.index[:3]
        self._check_roundtrip(df, tm.assert_frame_equal)

    # def test_cant_write_multiindex_table(self):
    #    # for now, #1848
    #    df = DataFrame(np.random.randn(10, 4),
    #                   index=[np.arange(5).repeat(2),
    #                          np.tile(np.arange(2), 5)])

    #    self.assertRaises(Exception, store.put, 'foo', df, table=True)


def _test_sort(obj):
    if isinstance(obj, DataFrame):
        return obj.reindex(sorted(obj.index))
    elif isinstance(obj, Panel):
        return obj.reindex(major=sorted(obj.major_axis))
    else:
        raise ValueError('type not supported here')

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
