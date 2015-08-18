import nose
import sys
import os
import warnings
import tempfile
from contextlib import contextmanager

import datetime
import numpy as np

import pandas
import pandas as pd
from pandas import (Series, DataFrame, Panel, MultiIndex, Categorical, bdate_range,
                    date_range, timedelta_range, Index, DatetimeIndex, TimedeltaIndex, isnull)

from pandas.compat import is_platform_windows, PY3
from pandas.io.pytables import _tables, TableIterator
try:
    _tables()
except ImportError as e:
    raise nose.SkipTest(e)


from pandas.io.pytables import (HDFStore, get_store, Term, read_hdf,
                                IncompatibilityWarning, PerformanceWarning,
                                AttributeConflictWarning, DuplicateWarning,
                                PossibleDataLossError, ClosedFileError)
from pandas.io import pytables as pytables
import pandas.util.testing as tm
from pandas.util.testing import (assert_panel4d_equal,
                                 assert_panel_equal,
                                 assert_frame_equal,
                                 assert_series_equal)
from pandas import concat, Timestamp
from pandas import compat
from pandas.compat import range, lrange, u
from pandas.util.testing import assert_produces_warning
from numpy.testing.decorators import slow

try:
    import tables
except ImportError:
    raise nose.SkipTest('no pytables')

from distutils.version import LooseVersion

_default_compressor = LooseVersion(tables.__version__) >= '2.2' \
    and 'blosc' or 'zlib'

_multiprocess_can_split_ = False

# testing on windows/py3 seems to fault
# for using compression
skip_compression = PY3 and is_platform_windows()

# contextmanager to ensure the file cleanup
def safe_remove(path):
    if path is not None:
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


def create_tempfile(path):
    """ create an unopened named temporary file """
    return os.path.join(tempfile.gettempdir(),path)


@contextmanager
def ensure_clean_store(path, mode='a', complevel=None, complib=None,
              fletcher32=False):

    try:

        # put in the temporary path if we don't have one already
        if not len(os.path.dirname(path)):
            path = create_tempfile(path)

        store = HDFStore(path, mode=mode, complevel=complevel,
                         complib=complib, fletcher32=False)
        yield store
    finally:
        safe_close(store)
        if mode == 'w' or mode == 'a':
            safe_remove(path)


@contextmanager
def ensure_clean_path(path):
    """
    return essentially a named temporary file that is not opened
    and deleted on existing; if path is a list, then create and
    return list of filenames
    """
    try:
        if isinstance(path, list):
            filenames = [ create_tempfile(p) for p in path ]
            yield filenames
        else:
            filenames = [ create_tempfile(path) ]
            yield filenames[0]
    finally:
        for f in filenames:
            safe_remove(f)


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


def compat_assert_produces_warning(w, f):
    """ don't produce a warning under PY3 """
    if compat.PY3:
        f()
    else:
        with tm.assert_produces_warning(expected_warning=w):
            f()


class Base(tm.TestCase):

    @classmethod
    def setUpClass(cls):
        super(Base, cls).setUpClass()

        # Pytables 3.0.0 deprecates lots of things
        tm.reset_testing_mode()

    @classmethod
    def tearDownClass(cls):
        super(Base, cls).tearDownClass()

        # Pytables 3.0.0 deprecates lots of things
        tm.set_testing_mode()

    def setUp(self):
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        self.path = 'tmp.__%s__.h5' % tm.rands(10)

    def tearDown(self):
        pass


class TestHDFStore(Base, tm.TestCase):

    def test_factory_fun(self):
        path = create_tempfile(self.path)
        try:
            with get_store(path) as tbl:
                raise ValueError('blah')
        except ValueError:
            pass
        finally:
            safe_remove(path)

        try:
            with get_store(path) as tbl:
                tbl['a'] = tm.makeDataFrame()

            with get_store(path) as tbl:
                self.assertEqual(len(tbl), 1)
                self.assertEqual(type(tbl['a']), DataFrame)
        finally:
            safe_remove(self.path)

    def test_context(self):
        path = create_tempfile(self.path)
        try:
            with HDFStore(path) as tbl:
                raise ValueError('blah')
        except ValueError:
            pass
        finally:
            safe_remove(path)

        try:
            with HDFStore(path) as tbl:
                tbl['a'] = tm.makeDataFrame()

            with HDFStore(path) as tbl:
                self.assertEqual(len(tbl), 1)
                self.assertEqual(type(tbl['a']), DataFrame)
        finally:
            safe_remove(path)

    def test_conv_read_write(self):
        path = create_tempfile(self.path)
        try:
            def roundtrip(key, obj,**kwargs):
                obj.to_hdf(path, key,**kwargs)
                return read_hdf(path, key)

            o = tm.makeTimeSeries()
            assert_series_equal(o, roundtrip('series',o))

            o = tm.makeStringSeries()
            assert_series_equal(o, roundtrip('string_series',o))

            o = tm.makeDataFrame()
            assert_frame_equal(o, roundtrip('frame',o))

            o = tm.makePanel()
            assert_panel_equal(o, roundtrip('panel',o))

            # table
            df = DataFrame(dict(A=lrange(5), B=lrange(5)))
            df.to_hdf(path,'table',append=True)
            result = read_hdf(path, 'table', where = ['index>2'])
            assert_frame_equal(df[df.index>2],result)

        finally:
            safe_remove(path)

    def test_long_strings(self):

        # GH6166
        # unconversion of long strings was being chopped in earlier
        # versions of numpy < 1.7.2
        df = DataFrame({'a': tm.rands_array(100, size=10)},
                       index=tm.rands_array(100, size=10))

        with ensure_clean_store(self.path) as store:
            store.append('df', df, data_columns=['a'])

            result = store.select('df')
            assert_frame_equal(df, result)


    def test_api(self):

        # GH4584
        # API issue when to_hdf doesn't acdept append AND format args
        with ensure_clean_path(self.path) as path:

            df = tm.makeDataFrame()
            df.iloc[:10].to_hdf(path,'df',append=True,format='table')
            df.iloc[10:].to_hdf(path,'df',append=True,format='table')
            assert_frame_equal(read_hdf(path,'df'),df)

            # append to False
            df.iloc[:10].to_hdf(path,'df',append=False,format='table')
            df.iloc[10:].to_hdf(path,'df',append=True,format='table')
            assert_frame_equal(read_hdf(path,'df'),df)

        with ensure_clean_path(self.path) as path:

            df = tm.makeDataFrame()
            df.iloc[:10].to_hdf(path,'df',append=True)
            df.iloc[10:].to_hdf(path,'df',append=True,format='table')
            assert_frame_equal(read_hdf(path,'df'),df)

            # append to False
            df.iloc[:10].to_hdf(path,'df',append=False,format='table')
            df.iloc[10:].to_hdf(path,'df',append=True)
            assert_frame_equal(read_hdf(path,'df'),df)

        with ensure_clean_path(self.path) as path:

            df = tm.makeDataFrame()
            df.to_hdf(path,'df',append=False,format='fixed')
            assert_frame_equal(read_hdf(path,'df'),df)

            df.to_hdf(path,'df',append=False,format='f')
            assert_frame_equal(read_hdf(path,'df'),df)

            df.to_hdf(path,'df',append=False)
            assert_frame_equal(read_hdf(path,'df'),df)

            df.to_hdf(path,'df')
            assert_frame_equal(read_hdf(path,'df'),df)

        with ensure_clean_store(self.path) as store:

            path = store._path
            df = tm.makeDataFrame()

            _maybe_remove(store,'df')
            store.append('df',df.iloc[:10],append=True,format='table')
            store.append('df',df.iloc[10:],append=True,format='table')
            assert_frame_equal(store.select('df'),df)

            # append to False
            _maybe_remove(store,'df')
            store.append('df',df.iloc[:10],append=False,format='table')
            store.append('df',df.iloc[10:],append=True,format='table')
            assert_frame_equal(store.select('df'),df)

            # formats
            _maybe_remove(store,'df')
            store.append('df',df.iloc[:10],append=False,format='table')
            store.append('df',df.iloc[10:],append=True,format='table')
            assert_frame_equal(store.select('df'),df)

            _maybe_remove(store,'df')
            store.append('df',df.iloc[:10],append=False,format='table')
            store.append('df',df.iloc[10:],append=True,format=None)
            assert_frame_equal(store.select('df'),df)

        with ensure_clean_path(self.path) as path:

            # invalid
            df = tm.makeDataFrame()
            self.assertRaises(ValueError, df.to_hdf, path,'df',append=True,format='f')
            self.assertRaises(ValueError, df.to_hdf, path,'df',append=True,format='fixed')

            self.assertRaises(TypeError, df.to_hdf, path,'df',append=True,format='foo')
            self.assertRaises(TypeError, df.to_hdf, path,'df',append=False,format='bar')

        #File path doesn't exist
        path = ""
        self.assertRaises(IOError, read_hdf, path, 'df')

    def test_api_default_format(self):

        # default_format option
        with ensure_clean_store(self.path) as store:
            df = tm.makeDataFrame()

            pandas.set_option('io.hdf.default_format','fixed')
            _maybe_remove(store,'df')
            store.put('df',df)
            self.assertFalse(store.get_storer('df').is_table)
            self.assertRaises(ValueError, store.append, 'df2',df)

            pandas.set_option('io.hdf.default_format','table')
            _maybe_remove(store,'df')
            store.put('df',df)
            self.assertTrue(store.get_storer('df').is_table)
            _maybe_remove(store,'df2')
            store.append('df2',df)
            self.assertTrue(store.get_storer('df').is_table)

            pandas.set_option('io.hdf.default_format',None)

        with ensure_clean_path(self.path) as path:

            df = tm.makeDataFrame()

            pandas.set_option('io.hdf.default_format','fixed')
            df.to_hdf(path,'df')
            with get_store(path) as store:
                self.assertFalse(store.get_storer('df').is_table)
            self.assertRaises(ValueError, df.to_hdf, path,'df2', append=True)

            pandas.set_option('io.hdf.default_format','table')
            df.to_hdf(path,'df3')
            with HDFStore(path) as store:
                self.assertTrue(store.get_storer('df3').is_table)
            df.to_hdf(path,'df4',append=True)
            with HDFStore(path) as store:
                self.assertTrue(store.get_storer('df4').is_table)

            pandas.set_option('io.hdf.default_format',None)

    def test_keys(self):

        with ensure_clean_store(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store['b'] = tm.makeStringSeries()
            store['c'] = tm.makeDataFrame()
            store['d'] = tm.makePanel()
            store['foo/bar'] = tm.makePanel()
            self.assertEqual(len(store), 5)
            self.assertTrue(set(
                store.keys()) == set(['/a', '/b', '/c', '/d', '/foo/bar']))

    def test_repr(self):

        with ensure_clean_store(self.path) as store:
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
            df = df.consolidate()._convert(datetime=True)

            warnings.filterwarnings('ignore', category=PerformanceWarning)
            store['df'] = df
            warnings.filterwarnings('always', category=PerformanceWarning)

            # make a random group in hdf space
            store._handle.create_group(store._handle.root,'bah')

            repr(store)
            str(store)

        # storers
        with ensure_clean_store(self.path) as store:

            df = tm.makeDataFrame()
            store.append('df',df)

            s = store.get_storer('df')
            repr(s)
            str(s)

    def test_contains(self):

        with ensure_clean_store(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store['b'] = tm.makeDataFrame()
            store['foo/bar'] = tm.makeDataFrame()
            self.assertIn('a', store)
            self.assertIn('b', store)
            self.assertNotIn('c', store)
            self.assertIn('foo/bar', store)
            self.assertIn('/foo/bar', store)
            self.assertNotIn('/foo/b', store)
            self.assertNotIn('bar', store)

            # GH 2694
            warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)
            store['node())'] = tm.makeDataFrame()
            self.assertIn('node())', store)

    def test_versioning(self):

        with ensure_clean_store(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store['b'] = tm.makeDataFrame()
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df1')
            store.append('df1', df[:10])
            store.append('df1', df[10:])
            self.assertEqual(store.root.a._v_attrs.pandas_version, '0.15.2')
            self.assertEqual(store.root.b._v_attrs.pandas_version, '0.15.2')
            self.assertEqual(store.root.df1._v_attrs.pandas_version, '0.15.2')

            # write a file and wipe its versioning
            _maybe_remove(store, 'df2')
            store.append('df2', df)

            # this is an error because its table_type is appendable, but no version
            # info
            store.get_node('df2')._v_attrs.pandas_version = None
            self.assertRaises(Exception, store.select, 'df2')

    def test_mode(self):

        df = tm.makeTimeDataFrame()

        def check(mode):

            with ensure_clean_path(self.path) as path:

                # constructor
                if mode in ['r','r+']:
                    self.assertRaises(IOError, HDFStore, path, mode=mode)

                else:
                    store = HDFStore(path,mode=mode)
                    self.assertEqual(store._handle.mode, mode)
                    store.close()

            with ensure_clean_path(self.path) as path:

                # context
                if mode in ['r','r+']:
                    def f():
                        with HDFStore(path,mode=mode) as store:
                            pass
                    self.assertRaises(IOError, f)
                else:
                    with HDFStore(path,mode=mode) as store:
                        self.assertEqual(store._handle.mode, mode)

            with ensure_clean_path(self.path) as path:

                # conv write
                if mode in ['r','r+']:
                    self.assertRaises(IOError, df.to_hdf, path, 'df', mode=mode)
                    df.to_hdf(path,'df',mode='w')
                else:
                    df.to_hdf(path,'df',mode=mode)

                # conv read
                if mode in ['w']:
                    self.assertRaises(KeyError, read_hdf, path, 'df', mode=mode)
                else:
                    result = read_hdf(path,'df',mode=mode)
                    assert_frame_equal(result,df)

        check('r')
        check('r+')
        check('a')
        check('w')

    def test_reopen_handle(self):

        with ensure_clean_path(self.path) as path:

            store = HDFStore(path,mode='a')
            store['a'] = tm.makeTimeSeries()

            # invalid mode change
            self.assertRaises(PossibleDataLossError, store.open, 'w')
            store.close()
            self.assertFalse(store.is_open)

            # truncation ok here
            store.open('w')
            self.assertTrue(store.is_open)
            self.assertEqual(len(store), 0)
            store.close()
            self.assertFalse(store.is_open)

            store = HDFStore(path,mode='a')
            store['a'] = tm.makeTimeSeries()

            # reopen as read
            store.open('r')
            self.assertTrue(store.is_open)
            self.assertEqual(len(store), 1)
            self.assertEqual(store._mode, 'r')
            store.close()
            self.assertFalse(store.is_open)

            # reopen as append
            store.open('a')
            self.assertTrue(store.is_open)
            self.assertEqual(len(store), 1)
            self.assertEqual(store._mode, 'a')
            store.close()
            self.assertFalse(store.is_open)

            # reopen as append (again)
            store.open('a')
            self.assertTrue(store.is_open)
            self.assertEqual(len(store), 1)
            self.assertEqual(store._mode, 'a')
            store.close()
            self.assertFalse(store.is_open)

    def test_open_args(self):

        with ensure_clean_path(self.path) as path:

            df = tm.makeDataFrame()

            # create an in memory store
            store = HDFStore(path,mode='a',driver='H5FD_CORE',driver_core_backing_store=0)
            store['df'] = df
            store.append('df2',df)

            tm.assert_frame_equal(store['df'],df)
            tm.assert_frame_equal(store['df2'],df)

            store.close()

            # the file should not have actually been written
            self.assertFalse(os.path.exists(path))

    def test_flush(self):

        with ensure_clean_store(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            store.flush()
            store.flush(fsync=True)

    def test_get(self):

        with ensure_clean_store(self.path) as store:
            store['a'] = tm.makeTimeSeries()
            left = store.get('a')
            right = store['a']
            tm.assert_series_equal(left, right)

            left = store.get('/a')
            right = store['/a']
            tm.assert_series_equal(left, right)

            self.assertRaises(KeyError, store.get, 'b')

    def test_getattr(self):

        with ensure_clean_store(self.path) as store:

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

        with ensure_clean_store(self.path) as store:

            ts = tm.makeTimeSeries()
            df = tm.makeTimeDataFrame()
            store['a'] = ts
            store['b'] = df[:10]
            store['foo/bar/bah'] = df[:10]
            store['foo'] = df[:10]
            store['/foo'] = df[:10]
            store.put('c', df[:10], format='table')

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
            store.put('c', df[:10], format='table', append=False)
            tm.assert_frame_equal(df[:10], store['c'])

    def test_put_string_index(self):

        with ensure_clean_store(self.path) as store:

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

        with ensure_clean_store(self.path) as store:
            df = tm.makeTimeDataFrame()

            store.put('c', df, format='table', complib='zlib')
            tm.assert_frame_equal(store['c'], df)

            # can't compress if format='fixed'
            self.assertRaises(ValueError, store.put, 'b', df,
                              format='fixed', complib='zlib')

    def test_put_compression_blosc(self):
        tm.skip_if_no_package('tables', '2.2', app='blosc support')
        if skip_compression:
            raise nose.SkipTest("skipping on windows/PY3")

        df = tm.makeTimeDataFrame()

        with ensure_clean_store(self.path) as store:

            # can't compress if format='fixed'
            self.assertRaises(ValueError, store.put, 'b', df,
                              format='fixed', complib='blosc')

            store.put('c', df, format='table', complib='blosc')
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
        df = df.consolidate()._convert(datetime=True)

        with ensure_clean_store(self.path) as store:
            _maybe_remove(store, 'df')

            # cannot use assert_produces_warning here for some reason
            # a PendingDeprecationWarning is also raised?
            warnings.filterwarnings('ignore', category=PerformanceWarning)
            store.put('df',df)
            warnings.filterwarnings('always', category=PerformanceWarning)

            expected = store.get('df')
            tm.assert_frame_equal(expected,df)

    def test_append(self):

        with ensure_clean_store(self.path) as store:
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df1')
            store.append('df1', df[:10])
            store.append('df1', df[10:])
            tm.assert_frame_equal(store['df1'], df)

            _maybe_remove(store, 'df2')
            store.put('df2', df[:10], format='table')
            store.append('df2', df[10:])
            tm.assert_frame_equal(store['df2'], df)

            _maybe_remove(store, 'df3')
            store.append('/df3', df[:10])
            store.append('/df3', df[10:])
            tm.assert_frame_equal(store['df3'], df)

            # this is allowed by almost always don't want to do it
            with tm.assert_produces_warning(expected_warning=tables.NaturalNameWarning):
                _maybe_remove(store, '/df3 foo')
                store.append('/df3 foo', df[:10])
                store.append('/df3 foo', df[10:])
                tm.assert_frame_equal(store['df3 foo'], df)

            # panel
            wp = tm.makePanel()
            _maybe_remove(store, 'wp1')
            store.append('wp1', wp.ix[:, :10, :])
            store.append('wp1', wp.ix[:, 10:, :])
            assert_panel_equal(store['wp1'], wp)

            # ndim
            p4d = tm.makePanel4D()
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :])
            store.append('p4d', p4d.ix[:, :, 10:, :])
            assert_panel4d_equal(store['p4d'], p4d)

            # test using axis labels
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=[
                    'items', 'major_axis', 'minor_axis'])
            store.append('p4d', p4d.ix[:, :, 10:, :], axes=[
                    'items', 'major_axis', 'minor_axis'])
            assert_panel4d_equal(store['p4d'], p4d)

            # test using differnt number of items on each axis
            p4d2 = p4d.copy()
            p4d2['l4'] = p4d['l1']
            p4d2['l5'] = p4d['l1']
            _maybe_remove(store, 'p4d2')
            store.append(
                'p4d2', p4d2, axes=['items', 'major_axis', 'minor_axis'])
            assert_panel4d_equal(store['p4d2'], p4d2)

            # test using differt order of items on the non-index axes
            _maybe_remove(store, 'wp1')
            wp_append1 = wp.ix[:, :10, :]
            store.append('wp1', wp_append1)
            wp_append2 = wp.ix[:, 10:, :].reindex(items=wp.items[::-1])
            store.append('wp1', wp_append2)
            assert_panel_equal(store['wp1'], wp)

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

    def test_append_series(self):

        with ensure_clean_store(self.path) as store:

            # basic
            ss = tm.makeStringSeries()
            ts = tm.makeTimeSeries()
            ns = Series(np.arange(100))

            store.append('ss', ss)
            result = store['ss']
            tm.assert_series_equal(result, ss)
            self.assertIsNone(result.name)

            store.append('ts', ts)
            result = store['ts']
            tm.assert_series_equal(result, ts)
            self.assertIsNone(result.name)

            ns.name = 'foo'
            store.append('ns', ns)
            result = store['ns']
            tm.assert_series_equal(result, ns)
            self.assertEqual(result.name, ns.name)

            # select on the values
            expected = ns[ns>60]
            result = store.select('ns',Term('foo>60'))
            tm.assert_series_equal(result,expected)

            # select on the index and values
            expected = ns[(ns>70) & (ns.index<90)]
            result = store.select('ns',[Term('foo>70'), Term('index<90')])
            tm.assert_series_equal(result,expected)

            # multi-index
            mi = DataFrame(np.random.randn(5,1),columns=['A'])
            mi['B'] = np.arange(len(mi))
            mi['C'] = 'foo'
            mi.loc[3:5,'C'] = 'bar'
            mi.set_index(['C','B'],inplace=True)
            s = mi.stack()
            s.index = s.index.droplevel(2)
            store.append('mi', s)
            tm.assert_series_equal(store['mi'], s)

    def test_store_index_types(self):
        # GH5386
        # test storing various index types

        with ensure_clean_store(self.path) as store:

            def check(format,index):
                df = DataFrame(np.random.randn(10,2),columns=list('AB'))
                df.index = index(len(df))

                _maybe_remove(store, 'df')
                store.put('df',df,format=format)
                assert_frame_equal(df,store['df'])

            for index in [ tm.makeFloatIndex, tm.makeStringIndex, tm.makeIntIndex,
                           tm.makeDateIndex ]:

                check('table',index)
                check('fixed',index)

            # period index currently broken for table
            # seee GH7796 FIXME
            check('fixed',tm.makePeriodIndex)
            #check('table',tm.makePeriodIndex)

            # unicode
            index = tm.makeUnicodeIndex
            if compat.PY3:
                check('table',index)
                check('fixed',index)
            else:

                # only support for fixed types (and they have a perf warning)
                self.assertRaises(TypeError, check, 'table', index)
                with tm.assert_produces_warning(expected_warning=PerformanceWarning):
                    check('fixed',index)

    def test_encoding(self):

        if sys.byteorder != 'little':
            raise nose.SkipTest('system byteorder is not little')

        with ensure_clean_store(self.path) as store:
            df = DataFrame(dict(A='foo',B='bar'),index=range(5))
            df.loc[2,'A'] = np.nan
            df.loc[3,'B'] = np.nan
            _maybe_remove(store, 'df')
            store.append('df', df, encoding='ascii')
            tm.assert_frame_equal(store['df'], df)

            expected = df.reindex(columns=['A'])
            result = store.select('df',Term('columns=A',encoding='ascii'))
            tm.assert_frame_equal(result,expected)

    def test_latin_encoding(self):

        if compat.PY2:
            self.assertRaisesRegexp(TypeError, '\[unicode\] is not implemented as a table column')
            return

        values = [[b'E\xc9, 17', b'', b'a', b'b', b'c'],
                  [b'E\xc9, 17', b'a', b'b', b'c'],
                  [b'EE, 17', b'', b'a', b'b', b'c'],
                  [b'E\xc9, 17', b'\xf8\xfc', b'a', b'b', b'c'],
                  [b'', b'a', b'b', b'c'],
                  [b'\xf8\xfc', b'a', b'b', b'c'],
                  [b'A\xf8\xfc', b'', b'a', b'b', b'c'],
                  [np.nan, b'', b'b', b'c'],
                  [b'A\xf8\xfc', np.nan, b'', b'b', b'c']]

        def _try_decode(x, encoding='latin-1'):
            try:
                return x.decode(encoding)
            except AttributeError:
                return x
        # not sure how to remove latin-1 from code in python 2 and 3
        values = [[_try_decode(x) for x in y] for y in values]

        examples = []
        for dtype in ['category', object]:
            for val in values:
                examples.append(pandas.Series(val, dtype=dtype))

        def roundtrip(s, key='data', encoding='latin-1', nan_rep=''):
            with ensure_clean_path(self.path) as store:
                s.to_hdf(store, key, format='table', encoding=encoding,
                        nan_rep=nan_rep)
                retr = read_hdf(store, key)
                s_nan = s.replace(nan_rep, np.nan)
                assert_series_equal(s_nan, retr)

        for s in examples:
            roundtrip(s)

        # fails:
        # for x in examples:
        #     roundtrip(s, nan_rep=b'\xf8\xfc')


    def test_append_some_nans(self):

        with ensure_clean_store(self.path) as store:
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

    def test_append_all_nans(self):

        with ensure_clean_store(self.path) as store:

            df = DataFrame({'A1' : np.random.randn(20),
                            'A2' : np.random.randn(20)},
                           index=np.arange(20))
            df.ix[0:15,:] = np.nan


            # nan some entire rows (dropna=True)
            _maybe_remove(store, 'df')
            store.append('df', df[:10], dropna=True)
            store.append('df', df[10:], dropna=True)
            tm.assert_frame_equal(store['df'], df[-4:])

            # nan some entire rows (dropna=False)
            _maybe_remove(store, 'df2')
            store.append('df2', df[:10], dropna=False)
            store.append('df2', df[10:], dropna=False)
            tm.assert_frame_equal(store['df2'], df)

            # tests the option io.hdf.dropna_table
            pandas.set_option('io.hdf.dropna_table',False)
            _maybe_remove(store, 'df3')
            store.append('df3', df[:10])
            store.append('df3', df[10:])
            tm.assert_frame_equal(store['df3'], df)

            pandas.set_option('io.hdf.dropna_table',True)
            _maybe_remove(store, 'df4')
            store.append('df4', df[:10])
            store.append('df4', df[10:])
            tm.assert_frame_equal(store['df4'], df[-4:])

            # nan some entire rows (string are still written!)
            df = DataFrame({'A1' : np.random.randn(20),
                            'A2' : np.random.randn(20),
                            'B' : 'foo', 'C' : 'bar'},
                           index=np.arange(20))

            df.ix[0:15,:] = np.nan

            _maybe_remove(store, 'df')
            store.append('df', df[:10], dropna=True)
            store.append('df', df[10:], dropna=True)
            tm.assert_frame_equal(store['df'], df)

            _maybe_remove(store, 'df2')
            store.append('df2', df[:10], dropna=False)
            store.append('df2', df[10:], dropna=False)
            tm.assert_frame_equal(store['df2'], df)

            # nan some entire rows (but since we have dates they are still written!)
            df = DataFrame({'A1' : np.random.randn(20),
                            'A2' : np.random.randn(20),
                            'B' : 'foo', 'C' : 'bar', 'D' : Timestamp("20010101"), 'E' : datetime.datetime(2001,1,2,0,0) },
                           index=np.arange(20))

            df.ix[0:15,:] = np.nan

            _maybe_remove(store, 'df')
            store.append('df', df[:10], dropna=True)
            store.append('df', df[10:], dropna=True)
            tm.assert_frame_equal(store['df'], df)

            _maybe_remove(store, 'df2')
            store.append('df2', df[:10], dropna=False)
            store.append('df2', df[10:], dropna=False)
            tm.assert_frame_equal(store['df2'], df)

        # Test to make sure defaults are to not drop.
        # Corresponding to Issue 9382
        df_with_missing = DataFrame({'col1':[0, np.nan, 2], 'col2':[1, np.nan,  np.nan]})

        with ensure_clean_path(self.path) as path:
            df_with_missing.to_hdf(path, 'df_with_missing', format = 'table')
            reloaded = read_hdf(path, 'df_with_missing')
            tm.assert_frame_equal(df_with_missing, reloaded)

        matrix = [[[np.nan, np.nan, np.nan],[1,np.nan,np.nan]],
            [[np.nan, np.nan, np.nan], [np.nan,5,6]],
            [[np.nan, np.nan, np.nan],[np.nan,3,np.nan]]]

        panel_with_missing = Panel(matrix, items=['Item1', 'Item2','Item3'],
                   major_axis=[1,2],
                   minor_axis=['A', 'B', 'C'])

        with ensure_clean_path(self.path) as path:
           panel_with_missing.to_hdf(path, 'panel_with_missing', format='table')
           reloaded_panel = read_hdf(path, 'panel_with_missing')
           tm.assert_panel_equal(panel_with_missing, reloaded_panel)

    def test_append_frame_column_oriented(self):

        with ensure_clean_store(self.path) as store:

            # column oriented
            df = tm.makeTimeDataFrame()
            _maybe_remove(store, 'df1')
            store.append('df1', df.ix[:, :2], axes=['columns'])
            store.append('df1', df.ix[:, 2:])
            tm.assert_frame_equal(store['df1'], df)

            result = store.select('df1', 'columns=A')
            expected = df.reindex(columns=['A'])
            tm.assert_frame_equal(expected, result)

            # selection on the non-indexable
            result = store.select(
                'df1', ('columns=A', Term('index=df.index[0:4]')))
            expected = df.reindex(columns=['A'], index=df.index[0:4])
            tm.assert_frame_equal(expected, result)

            # this isn't supported
            self.assertRaises(TypeError, store.select, 'df1', (
                    'columns=A', Term('index>df.index[4]')))

    def test_append_with_different_block_ordering(self):

        #GH 4096; using same frames, but different block orderings
        with ensure_clean_store(self.path) as store:

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

        # test a different ordering but with more fields (like invalid combinate)
        with ensure_clean_store(self.path) as store:

            df = DataFrame(np.random.randn(10,2),columns=list('AB'), dtype='float64')
            df['int64'] = Series([1]*len(df),dtype='int64')
            df['int16'] = Series([1]*len(df),dtype='int16')
            store.append('df',df)

            # store additonal fields in different blocks
            df['int16_2'] = Series([1]*len(df),dtype='int16')
            self.assertRaises(ValueError, store.append, 'df', df)

            # store multile additonal fields in different blocks
            df['float_3'] = Series([1.]*len(df),dtype='float64')
            self.assertRaises(ValueError, store.append, 'df', df)

    def test_ndim_indexables(self):
        """ test using ndim tables in new ways"""

        with ensure_clean_store(self.path) as store:

            p4d = tm.makePanel4D()

            def check_indexers(key, indexers):
                for i, idx in enumerate(indexers):
                    self.assertTrue(getattr(getattr(
                        store.root, key).table.description, idx)._v_pos == i)

            # append then change (will take existing schema)
            indexers = ['items', 'major_axis', 'minor_axis']

            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=indexers)
            store.append('p4d', p4d.ix[:, :, 10:, :])
            assert_panel4d_equal(store.select('p4d'), p4d)
            check_indexers('p4d', indexers)

            # same as above, but try to append with differnt axes
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=indexers)
            store.append('p4d', p4d.ix[:, :, 10:, :], axes=[
                    'labels', 'items', 'major_axis'])
            assert_panel4d_equal(store.select('p4d'), p4d)
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
            assert_panel4d_equal(store['p4d'], p4d)
            check_indexers('p4d', indexers)

            # different than default indexables #2
            indexers = ['major_axis', 'labels', 'minor_axis']
            _maybe_remove(store, 'p4d')
            store.append('p4d', p4d.ix[:, :, :10, :], axes=indexers)
            store.append('p4d', p4d.ix[:, :, 10:, :])
            assert_panel4d_equal(store['p4d'], p4d)
            check_indexers('p4d', indexers)

            # partial selection
            result = store.select('p4d', ['labels=l1'])
            expected = p4d.reindex(labels=['l1'])
            assert_panel4d_equal(result, expected)

            # partial selection2
            result = store.select('p4d', [Term(
                        'labels=l1'), Term('items=ItemA'), Term('minor_axis=B')])
            expected = p4d.reindex(
                labels=['l1'], items=['ItemA'], minor_axis=['B'])
            assert_panel4d_equal(result, expected)

            # non-existant partial selection
            result = store.select('p4d', [Term(
                        'labels=l1'), Term('items=Item1'), Term('minor_axis=B')])
            expected = p4d.reindex(labels=['l1'], items=[], minor_axis=['B'])
            assert_panel4d_equal(result, expected)

    def test_append_with_strings(self):

        with ensure_clean_store(self.path) as store:
            wp = tm.makePanel()
            wp2 = wp.rename_axis(
                dict([(x, "%s_extra" % x) for x in wp.minor_axis]), axis=2)

            def check_col(key,name,size):
                self.assertEqual(getattr(store.get_storer(key).table.description,name).itemsize, size)

            store.append('s1', wp, min_itemsize=20)
            store.append('s1', wp2)
            expected = concat([wp, wp2], axis=2)
            expected = expected.reindex(minor_axis=sorted(expected.minor_axis))
            assert_panel_equal(store['s1'], expected)
            check_col('s1', 'minor_axis', 20)

            # test dict format
            store.append('s2', wp, min_itemsize={'minor_axis': 20})
            store.append('s2', wp2)
            expected = concat([wp, wp2], axis=2)
            expected = expected.reindex(minor_axis=sorted(expected.minor_axis))
            assert_panel_equal(store['s2'], expected)
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

        with ensure_clean_store(self.path) as store:

            def check_col(key,name,size):
                self.assertEqual(getattr(store.get_storer(key).table.description,name).itemsize, size)

            df = DataFrame(dict(A = 'foo', B = 'bar'),index=range(10))

            # a min_itemsize that creates a data_column
            _maybe_remove(store, 'df')
            store.append('df', df, min_itemsize={'A' : 200 })
            check_col('df', 'A', 200)
            self.assertEqual(store.get_storer('df').data_columns, ['A'])

            # a min_itemsize that creates a data_column2
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns = ['B'], min_itemsize={'A' : 200 })
            check_col('df', 'A', 200)
            self.assertEqual(store.get_storer('df').data_columns, ['B','A'])

            # a min_itemsize that creates a data_column2
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns = ['B'], min_itemsize={'values' : 200 })
            check_col('df', 'B', 200)
            check_col('df', 'values_block_0', 200)
            self.assertEqual(store.get_storer('df').data_columns, ['B'])

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

        with ensure_clean_store(self.path) as store:
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
                'df', [Term('B>0'), Term('index>df.index[3]')])
            df_new = df.reindex(index=df.index[4:])
            expected = df_new[df_new.B > 0]
            tm.assert_frame_equal(result, expected)

            # data column selection with a string data_column
            df_new = df.copy()
            df_new['string'] = 'foo'
            df_new.loc[1:4,'string'] = np.nan
            df_new.loc[5:6,'string'] = 'bar'
            _maybe_remove(store, 'df')
            store.append('df', df_new, data_columns=['string'])
            result = store.select('df', [Term('string=foo')])
            expected = df_new[df_new.string == 'foo']
            tm.assert_frame_equal(result, expected)

            # using min_itemsize and a data column
            def check_col(key,name,size):
                self.assertEqual(getattr(store.get_storer(key).table.description,name).itemsize, size)

        with ensure_clean_store(self.path) as store:
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

        with ensure_clean_store(self.path) as store:
            df_new['string2'] = 'foobarbah'
            df_new['string_block1'] = 'foobarbah1'
            df_new['string_block2'] = 'foobarbah2'
            _maybe_remove(store, 'df')
            store.append('df', df_new, data_columns=['string', 'string2'], min_itemsize={'string': 30, 'string2': 40, 'values': 50})
            check_col('df', 'string', 30)
            check_col('df', 'string2', 40)
            check_col('df', 'values_block_1', 50)

        with ensure_clean_store(self.path) as store:
            # multiple data columns
            df_new = df.copy()
            df_new.ix[0,'A'] = 1.
            df_new.ix[0,'B'] = -1.
            df_new['string'] = 'foo'
            df_new.loc[1:4,'string'] = np.nan
            df_new.loc[5:6,'string'] = 'bar'
            df_new['string2'] = 'foo'
            df_new.loc[2:5,'string2'] = np.nan
            df_new.loc[7:8,'string2'] = 'bar'
            _maybe_remove(store, 'df')
            store.append(
                'df', df_new, data_columns=['A', 'B', 'string', 'string2'])
            result = store.select('df', [Term('string=foo'), Term(
                        'string2=foo'), Term('A>0'), Term('B<0')])
            expected = df_new[(df_new.string == 'foo') & (
                    df_new.string2 == 'foo') & (df_new.A > 0) & (df_new.B < 0)]
            tm.assert_frame_equal(result, expected, check_index_type=False)

            # yield an empty frame
            result = store.select('df', [Term('string=foo'), Term(
                        'string2=cool')])
            expected = df_new[(df_new.string == 'foo') & (
                    df_new.string2 == 'cool')]
            tm.assert_frame_equal(result, expected, check_index_type=False)

        with ensure_clean_store(self.path) as store:
            # doc example
            df_dc = df.copy()
            df_dc['string'] = 'foo'
            df_dc.ix[4:6, 'string'] = np.nan
            df_dc.ix[7:9, 'string'] = 'bar'
            df_dc['string2'] = 'cool'
            df_dc['datetime'] = Timestamp('20010102')
            df_dc = df_dc._convert(datetime=True)
            df_dc.ix[3:5, ['A', 'B', 'datetime']] = np.nan

            _maybe_remove(store, 'df_dc')
            store.append('df_dc', df_dc, data_columns=['B', 'C',
                                                       'string', 'string2', 'datetime'])
            result = store.select('df_dc', [Term('B>0')])

            expected = df_dc[df_dc.B > 0]
            tm.assert_frame_equal(result, expected, check_index_type=False)

            result = store.select(
                'df_dc', ['B > 0', 'C > 0', 'string == foo'])
            expected = df_dc[(df_dc.B > 0) & (df_dc.C > 0) & (
                    df_dc.string == 'foo')]
            tm.assert_frame_equal(result, expected, check_index_type=False)

        with ensure_clean_store(self.path) as store:
            # doc example part 2
            np.random.seed(1234)
            index = date_range('1/1/2000', periods=8)
            df_dc = DataFrame(np.random.randn(8, 3), index=index,
                              columns=['A', 'B', 'C'])
            df_dc['string'] = 'foo'
            df_dc.ix[4:6,'string'] = np.nan
            df_dc.ix[7:9,'string'] = 'bar'
            df_dc.ix[:,['B','C']] = df_dc.ix[:,['B','C']].abs()
            df_dc['string2'] = 'cool'

            # on-disk operations
            store.append('df_dc', df_dc, data_columns = ['B', 'C', 'string', 'string2'])

            result = store.select('df_dc', [ Term('B>0') ])
            expected = df_dc[df_dc.B>0]
            tm.assert_frame_equal(result,expected)

            result = store.select('df_dc', ['B > 0', 'C > 0', 'string == "foo"'])
            expected = df_dc[(df_dc.B > 0) & (df_dc.C > 0) & (df_dc.string == 'foo')]
            tm.assert_frame_equal(result,expected)

        with ensure_clean_store(self.path) as store:
            # panel
            # GH5717 not handling data_columns
            np.random.seed(1234)
            p = tm.makePanel()

            store.append('p1',p)
            tm.assert_panel_equal(store.select('p1'),p)

            store.append('p2',p,data_columns=True)
            tm.assert_panel_equal(store.select('p2'),p)

            result = store.select('p2',where='ItemA>0')
            expected = p.to_frame()
            expected = expected[expected['ItemA']>0]
            tm.assert_frame_equal(result.to_frame(),expected)

            result = store.select('p2',where='ItemA>0 & minor_axis=["A","B"]')
            expected = p.to_frame()
            expected = expected[expected['ItemA']>0]
            expected = expected[expected.reset_index(level=['major']).index.isin(['A','B'])]
            tm.assert_frame_equal(result.to_frame(),expected)

    def test_create_table_index(self):

        with ensure_clean_store(self.path) as store:

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

    def test_append_diff_item_order(self):

        wp = tm.makePanel()
        wp1 = wp.ix[:, :10, :]
        wp2 = wp.ix[['ItemC', 'ItemB', 'ItemA'], 10:, :]

        with ensure_clean_store(self.path) as store:
            store.put('panel', wp1, format='table')
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

        with ensure_clean_store(self.path) as store:
            store.append('mi', df)
            result = store.select('mi')
            tm.assert_frame_equal(result, df)

            # GH 3748
            result = store.select('mi',columns=['A','B'])
            expected = df.reindex(columns=['A','B'])
            tm.assert_frame_equal(result,expected)

        with ensure_clean_path('test.hdf') as path:
            df.to_hdf(path,'df',format='table')
            result = read_hdf(path,'df',columns=['A','B'])
            expected = df.reindex(columns=['A','B'])
            tm.assert_frame_equal(result,expected)

    def test_column_multiindex(self):
        # GH 4710
        # recreate multi-indexes properly

        index = MultiIndex.from_tuples([('A','a'), ('A','b'), ('B','a'), ('B','b')], names=['first','second'])
        df = DataFrame(np.arange(12).reshape(3,4), columns=index)

        with ensure_clean_store(self.path) as store:

            store.put('df',df)
            tm.assert_frame_equal(store['df'],df,check_index_type=True,check_column_type=True)

            store.put('df1',df,format='table')
            tm.assert_frame_equal(store['df1'],df,check_index_type=True,check_column_type=True)

            self.assertRaises(ValueError, store.put, 'df2',df,format='table',data_columns=['A'])
            self.assertRaises(ValueError, store.put, 'df3',df,format='table',data_columns=True)

        # appending multi-column on existing table (see GH 6167)
        with ensure_clean_store(self.path) as store:
            store.append('df2', df)
            store.append('df2', df)

            tm.assert_frame_equal(store['df2'], concat((df,df)))

        # non_index_axes name
        df = DataFrame(np.arange(12).reshape(3,4), columns=Index(list('ABCD'),name='foo'))

        with ensure_clean_store(self.path) as store:

            store.put('df1',df,format='table')
            tm.assert_frame_equal(store['df1'],df,check_index_type=True,check_column_type=True)

    def test_store_multiindex(self):

        # validate multi-index names
        # GH 5527
        with ensure_clean_store(self.path) as store:

            def make_index(names=None):
                return MultiIndex.from_tuples([( datetime.datetime(2013,12,d), s, t) for d in range(1,3) for s in range(2) for t in range(3)],
                                              names=names)


            # no names
            _maybe_remove(store, 'df')
            df = DataFrame(np.zeros((12,2)), columns=['a','b'], index=make_index())
            store.append('df',df)
            tm.assert_frame_equal(store.select('df'),df)

            # partial names
            _maybe_remove(store, 'df')
            df = DataFrame(np.zeros((12,2)), columns=['a','b'], index=make_index(['date',None,None]))
            store.append('df',df)
            tm.assert_frame_equal(store.select('df'),df)

            # series
            _maybe_remove(store, 's')
            s = Series(np.zeros(12), index=make_index(['date', None, None]))
            store.append('s',s)
            xp = Series(np.zeros(12), index=make_index(['date', 'level_1', 'level_2']))
            tm.assert_series_equal(store.select('s'), xp)

            # dup with column
            _maybe_remove(store, 'df')
            df = DataFrame(np.zeros((12,2)), columns=['a','b'], index=make_index(['date','a','t']))
            self.assertRaises(ValueError, store.append, 'df',df)

            # dup within level
            _maybe_remove(store, 'df')
            df = DataFrame(np.zeros((12,2)), columns=['a','b'], index=make_index(['date','date','date']))
            self.assertRaises(ValueError, store.append, 'df',df)

            # fully names
            _maybe_remove(store, 'df')
            df = DataFrame(np.zeros((12,2)), columns=['a','b'], index=make_index(['date','s','t']))
            store.append('df',df)
            tm.assert_frame_equal(store.select('df'),df)

    def test_select_columns_in_where(self):

        # GH 6169
        # recreate multi-indexes when columns is passed
        # in the `where` argument
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['foo_name', 'bar_name'])

        # With a DataFrame
        df = DataFrame(np.random.randn(10, 3), index=index,
                       columns=['A', 'B', 'C'])

        with ensure_clean_store(self.path) as store:
            store.put('df', df, format='table')
            expected = df[['A']]

            tm.assert_frame_equal(store.select('df', columns=['A']), expected)

            tm.assert_frame_equal(store.select('df', where="columns=['A']"), expected)

        # With a Series
        s = Series(np.random.randn(10), index=index,
                   name='A')
        with ensure_clean_store(self.path) as store:
            store.put('s', s, format='table')
            tm.assert_series_equal(store.select('s', where="columns=['A']"),s)

    def test_pass_spec_to_storer(self):

        df = tm.makeDataFrame()

        with ensure_clean_store(self.path) as store:
            store.put('df',df)
            self.assertRaises(TypeError, store.select, 'df', columns=['A'])
            self.assertRaises(TypeError, store.select, 'df',where=[('columns=A')])

    def test_append_misc(self):

        with ensure_clean_store(self.path) as store:

            # unsuported data types for non-tables
            p4d = tm.makePanel4D()
            self.assertRaises(TypeError, store.put,'p4d',p4d)

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

        # more chunksize in append tests
        def check(obj, comparator):
            for c in [10, 200, 1000]:
                with ensure_clean_store(self.path,mode='w') as store:
                    store.append('obj', obj, chunksize=c)
                    result = store.select('obj')
                    comparator(result,obj)

        df = tm.makeDataFrame()
        df['string'] = 'foo'
        df['float322'] = 1.
        df['float322'] = df['float322'].astype('float32')
        df['bool']     = df['float322'] > 0
        df['time1']    = Timestamp('20130101')
        df['time2']    = Timestamp('20130102')
        check(df, tm.assert_frame_equal)

        p = tm.makePanel()
        check(p, assert_panel_equal)

        p4d = tm.makePanel4D()
        check(p4d, assert_panel4d_equal)

        # empty frame, GH4273
        with ensure_clean_store(self.path) as store:

            # 0 len
            df_empty = DataFrame(columns=list('ABC'))
            store.append('df',df_empty)
            self.assertRaises(KeyError,store.select, 'df')

            # repeated append of 0/non-zero frames
            df = DataFrame(np.random.rand(10,3),columns=list('ABC'))
            store.append('df',df)
            assert_frame_equal(store.select('df'),df)
            store.append('df',df_empty)
            assert_frame_equal(store.select('df'),df)

            # store
            df = DataFrame(columns=list('ABC'))
            store.put('df2',df)
            assert_frame_equal(store.select('df2'),df)

            # 0 len
            p_empty = Panel(items=list('ABC'))
            store.append('p',p_empty)
            self.assertRaises(KeyError,store.select, 'p')

            # repeated append of 0/non-zero frames
            p = Panel(np.random.randn(3,4,5),items=list('ABC'))
            store.append('p',p)
            assert_panel_equal(store.select('p'),p)
            store.append('p',p_empty)
            assert_panel_equal(store.select('p'),p)

            # store
            store.put('p2',p_empty)
            assert_panel_equal(store.select('p2'),p_empty)

    def test_append_raise(self):

        with ensure_clean_store(self.path) as store:

            # test append with invalid input to get good error messages

            # list in column
            df = tm.makeDataFrame()
            df['invalid'] = [['a']] * len(df)
            self.assertEqual(df.dtypes['invalid'], np.object_)
            self.assertRaises(TypeError, store.append,'df',df)

            # multiple invalid columns
            df['invalid2'] = [['a']] * len(df)
            df['invalid3'] = [['a']] * len(df)
            self.assertRaises(TypeError, store.append,'df',df)

            # datetime with embedded nans as object
            df = tm.makeDataFrame()
            s = Series(datetime.datetime(2001,1,2),index=df.index)
            s = s.astype(object)
            s[0:5] = np.nan
            df['invalid'] = s
            self.assertEqual(df.dtypes['invalid'], np.object_)
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

        with ensure_clean_store(self.path) as store:
            store.put('frame', df1, format='table')
            self.assertRaises(TypeError, store.put, 'frame', df2,
                              format='table', append=True)

    def test_table_values_dtypes_roundtrip(self):

        with ensure_clean_store(self.path) as store:
            df1 = DataFrame({'a': [1, 2, 3]}, dtype='f8')
            store.append('df_f8', df1)
            assert_series_equal(df1.dtypes,store['df_f8'].dtypes)

            df2 = DataFrame({'a': [1, 2, 3]}, dtype='i8')
            store.append('df_i8', df2)
            assert_series_equal(df2.dtypes,store['df_i8'].dtypes)

            # incompatible dtype
            self.assertRaises(ValueError, store.append, 'df_i8', df1)

            # check creation/storage/retrieval of float32 (a bit hacky to actually create them thought)
            df1 = DataFrame(np.array([[1],[2],[3]],dtype='f4'),columns = ['A'])
            store.append('df_f4', df1)
            assert_series_equal(df1.dtypes,store['df_f4'].dtypes)
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
        df = df.consolidate()._convert(datetime=True)

        with ensure_clean_store(self.path) as store:
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

        with ensure_clean_store(self.path) as store:
            store.append('p1_mixed', wp)
            assert_panel_equal(store.select('p1_mixed'), wp)

        # ndim
        wp = tm.makePanel4D()
        wp['obj1'] = 'foo'
        wp['obj2'] = 'bar'
        wp['bool1'] = wp['l1'] > 0
        wp['bool2'] = wp['l2'] > 0
        wp['int1'] = 1
        wp['int2'] = 2
        wp = wp.consolidate()

        with ensure_clean_store(self.path) as store:
            store.append('p4d_mixed', wp)
            assert_panel4d_equal(store.select('p4d_mixed'), wp)

    def test_unimplemented_dtypes_table_columns(self):

        with ensure_clean_store(self.path) as store:

            l = [('date', datetime.date(2001, 1, 2))]

            # py3 ok for unicode
            if not compat.PY3:
                l.append(('unicode', u('\\u03c3')))

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
        df = df.consolidate()._convert(datetime=True)

        with ensure_clean_store(self.path) as store:
            # this fails because we have a date in the object block......
            self.assertRaises(TypeError, store.append, 'df_unimplemented', df)

    def test_calendar_roundtrip_issue(self):

        # 8591
        # doc example from tseries holiday section
        weekmask_egypt = 'Sun Mon Tue Wed Thu'
        holidays = ['2012-05-01', datetime.datetime(2013, 5, 1), np.datetime64('2014-05-01')]
        bday_egypt = pandas.offsets.CustomBusinessDay(holidays=holidays, weekmask=weekmask_egypt)
        dt = datetime.datetime(2013, 4, 30)
        dts = date_range(dt, periods=5, freq=bday_egypt)

        s = (Series(dts.weekday, dts).map(Series('Mon Tue Wed Thu Fri Sat Sun'.split())))

        with ensure_clean_store(self.path) as store:

            store.put('fixed',s)
            result = store.select('fixed')
            assert_series_equal(result, s)

            store.append('table',s)
            result = store.select('table')
            assert_series_equal(result, s)

    def test_append_with_timedelta(self):
        # GH 3577
        # append timedelta

        from datetime import timedelta
        df = DataFrame(dict(A = Timestamp('20130101'), B = [ Timestamp('20130101') + timedelta(days=i,seconds=10) for i in range(10) ]))
        df['C'] = df['A']-df['B']
        df.ix[3:5,'C'] = np.nan

        with ensure_clean_store(self.path) as store:

            # table
            _maybe_remove(store, 'df')
            store.append('df',df,data_columns=True)
            result = store.select('df')
            assert_frame_equal(result,df)

            result = store.select('df',Term("C<100000"))
            assert_frame_equal(result,df)

            result = store.select('df',Term("C","<",-3*86400))
            assert_frame_equal(result,df.iloc[3:])

            result = store.select('df',"C<'-3D'")
            assert_frame_equal(result,df.iloc[3:])

            # a bit hacky here as we don't really deal with the NaT properly

            result = store.select('df',"C<'-500000s'")
            result = result.dropna(subset=['C'])
            assert_frame_equal(result,df.iloc[6:])

            result = store.select('df',"C<'-3.5D'")
            result = result.iloc[1:]
            assert_frame_equal(result,df.iloc[4:])

            # fixed
            _maybe_remove(store, 'df2')
            store.put('df2',df)
            result = store.select('df2')
            assert_frame_equal(result,df)

    def test_remove(self):

        with ensure_clean_store(self.path) as store:

            ts = tm.makeTimeSeries()
            df = tm.makeDataFrame()
            store['a'] = ts
            store['b'] = df
            _maybe_remove(store, 'a')
            self.assertEqual(len(store), 1)
            tm.assert_frame_equal(df, store['b'])

            _maybe_remove(store, 'b')
            self.assertEqual(len(store), 0)

            # nonexistence
            self.assertRaises(KeyError, store.remove, 'a_nonexistent_store')

            # pathing
            store['a'] = ts
            store['b/foo'] = df
            _maybe_remove(store, 'foo')
            _maybe_remove(store, 'b/foo')
            self.assertEqual(len(store), 1)

            store['a'] = ts
            store['b/foo'] = df
            _maybe_remove(store, 'b')
            self.assertEqual(len(store), 1)

            # __delitem__
            store['a'] = ts
            store['b'] = df
            del store['a']
            del store['b']
            self.assertEqual(len(store), 0)

    def test_remove_where(self):

        with ensure_clean_store(self.path) as store:

            # non-existance
            crit1 = Term('index>foo')
            self.assertRaises(KeyError, store.remove, 'a', [crit1])

            # try to remove non-table (with crit)
            # non-table ok (where = None)
            wp = tm.makePanel(30)
            store.put('wp', wp, format='table')
            store.remove('wp', ["minor_axis=['A', 'D']"])
            rs = store.select('wp')
            expected = wp.reindex(minor_axis=['B', 'C'])
            assert_panel_equal(rs, expected)

            # empty where
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')

            # deleted number (entire table)
            n = store.remove('wp', [])
            self.assertTrue(n == 120)

            # non - empty where
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')
            self.assertRaises(ValueError, store.remove,
                              'wp', ['foo'])

            # selectin non-table with a where
            # store.put('wp2', wp, format='f')
            # self.assertRaises(ValueError, store.remove,
            #                  'wp2', [('column', ['A', 'D'])])

    def test_remove_startstop(self):
        # GH #4835 and #6177

        with ensure_clean_store(self.path) as store:

            wp = tm.makePanel(30)

            # start
            _maybe_remove(store, 'wp1')
            store.put('wp1', wp, format='t')
            n = store.remove('wp1', start=32)
            self.assertTrue(n == 120-32)
            result = store.select('wp1')
            expected = wp.reindex(major_axis=wp.major_axis[:32//4])
            assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp2')
            store.put('wp2', wp, format='t')
            n = store.remove('wp2', start=-32)
            self.assertTrue(n == 32)
            result = store.select('wp2')
            expected = wp.reindex(major_axis=wp.major_axis[:-32//4])
            assert_panel_equal(result, expected)

            # stop
            _maybe_remove(store, 'wp3')
            store.put('wp3', wp, format='t')
            n = store.remove('wp3', stop=32)
            self.assertTrue(n == 32)
            result = store.select('wp3')
            expected = wp.reindex(major_axis=wp.major_axis[32//4:])
            assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp4')
            store.put('wp4', wp, format='t')
            n = store.remove('wp4', stop=-32)
            self.assertTrue(n == 120-32)
            result = store.select('wp4')
            expected = wp.reindex(major_axis=wp.major_axis[-32//4:])
            assert_panel_equal(result, expected)

            # start n stop
            _maybe_remove(store, 'wp5')
            store.put('wp5', wp, format='t')
            n = store.remove('wp5', start=16, stop=-16)
            self.assertTrue(n == 120-32)
            result = store.select('wp5')
            expected = wp.reindex(major_axis=wp.major_axis[:16//4].union(wp.major_axis[-16//4:]))
            assert_panel_equal(result, expected)

            _maybe_remove(store, 'wp6')
            store.put('wp6', wp, format='t')
            n = store.remove('wp6', start=16, stop=16)
            self.assertTrue(n == 0)
            result = store.select('wp6')
            expected = wp.reindex(major_axis=wp.major_axis)
            assert_panel_equal(result, expected)

            # with where
            _maybe_remove(store, 'wp7')
            date = wp.major_axis.take(np.arange(0,30,3))
            crit = Term('major_axis=date')
            store.put('wp7', wp, format='t')
            n = store.remove('wp7', where=[crit], stop=80)
            self.assertTrue(n == 28)
            result = store.select('wp7')
            expected = wp.reindex(major_axis=wp.major_axis.difference(wp.major_axis[np.arange(0,20,3)]))
            assert_panel_equal(result, expected)

    def test_remove_crit(self):

        with ensure_clean_store(self.path) as store:

            wp = tm.makePanel(30)

            # group row removal
            _maybe_remove(store, 'wp3')
            date4 = wp.major_axis.take([0, 1, 2, 4, 5, 6, 8, 9, 10])
            crit4 = Term('major_axis=date4')
            store.put('wp3', wp, format='t')
            n = store.remove('wp3', where=[crit4])
            self.assertTrue(n == 36)

            result = store.select('wp3')
            expected = wp.reindex(major_axis=wp.major_axis.difference(date4))
            assert_panel_equal(result, expected)

            # upper half
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')
            date = wp.major_axis[len(wp.major_axis) // 2]

            crit1 = Term('major_axis>date')
            crit2 = Term("minor_axis=['A', 'D']")
            n = store.remove('wp', where=[crit1])
            self.assertTrue(n == 56)

            n = store.remove('wp', where=[crit2])
            self.assertTrue(n == 32)

            result = store['wp']
            expected = wp.truncate(after=date).reindex(minor=['B', 'C'])
            assert_panel_equal(result, expected)

            # individual row elements
            _maybe_remove(store, 'wp2')
            store.put('wp2', wp, format='table')

            date1 = wp.major_axis[1:3]
            crit1 = Term('major_axis=date1')
            store.remove('wp2', where=[crit1])
            result = store.select('wp2')
            expected = wp.reindex(major_axis=wp.major_axis.difference(date1))
            assert_panel_equal(result, expected)

            date2 = wp.major_axis[5]
            crit2 = Term('major_axis=date2')
            store.remove('wp2', where=[crit2])
            result = store['wp2']
            expected = wp.reindex(
                major_axis=wp.major_axis.difference(date1).difference(Index([date2])))
            assert_panel_equal(result, expected)

            date3 = [wp.major_axis[7], wp.major_axis[9]]
            crit3 = Term('major_axis=date3')
            store.remove('wp2', where=[crit3])
            result = store['wp2']
            expected = wp.reindex(
                major_axis=wp.major_axis.difference(date1).difference(Index([date2])).difference(Index(date3)))
            assert_panel_equal(result, expected)

            # corners
            _maybe_remove(store, 'wp4')
            store.put('wp4', wp, format='table')
            n = store.remove(
                'wp4', where=[Term('major_axis>wp.major_axis[-1]')])
            result = store.select('wp4')
            assert_panel_equal(result, wp)

    def test_invalid_terms(self):

        with ensure_clean_store(self.path) as store:

            df = tm.makeTimeDataFrame()
            df['string'] = 'foo'
            df.ix[0:4,'string'] = 'bar'
            wp = tm.makePanel()
            p4d = tm.makePanel4D()
            store.put('df', df, format='table')
            store.put('wp', wp, format='table')
            store.put('p4d', p4d, format='table')

            # some invalid terms
            self.assertRaises(ValueError, store.select, 'wp', "minor=['A', 'B']")
            self.assertRaises(ValueError, store.select, 'wp', ["index=['20121114']"])
            self.assertRaises(ValueError, store.select, 'wp', ["index=['20121114', '20121114']"])
            self.assertRaises(TypeError, Term)

            # more invalid
            self.assertRaises(ValueError,  store.select, 'df','df.index[3]')
            self.assertRaises(SyntaxError, store.select, 'df','index>')
            self.assertRaises(ValueError,  store.select, 'wp', "major_axis<'20000108' & minor_axis['A', 'B']")

        # from the docs
        with ensure_clean_path(self.path) as path:
            dfq = DataFrame(np.random.randn(10,4),columns=list('ABCD'),index=date_range('20130101',periods=10))
            dfq.to_hdf(path,'dfq',format='table',data_columns=True)

            # check ok
            read_hdf(path,'dfq',where="index>Timestamp('20130104') & columns=['A', 'B']")
            read_hdf(path,'dfq',where="A>0 or C>0")

        # catch the invalid reference
        with ensure_clean_path(self.path) as path:
            dfq = DataFrame(np.random.randn(10,4),columns=list('ABCD'),index=date_range('20130101',periods=10))
            dfq.to_hdf(path,'dfq',format='table')

            self.assertRaises(ValueError, read_hdf, path,'dfq',where="A>0 or C>0")

    def test_terms(self):

        with ensure_clean_store(self.path) as store:

            wp = tm.makePanel()
            p4d = tm.makePanel4D()
            wpneg = Panel.fromDict({-1: tm.makeDataFrame(), 0: tm.makeDataFrame(),
                                    1: tm.makeDataFrame()})
            store.put('wp', wp, format='table')
            store.put('p4d', p4d, format='table')
            store.put('wpneg', wpneg, format='table')

            # panel
            result = store.select('wp', [Term(
                        'major_axis<"20000108"'), Term("minor_axis=['A', 'B']")])
            expected = wp.truncate(after='20000108').reindex(minor=['A', 'B'])
            assert_panel_equal(result, expected)

            # with deprecation
            result = store.select('wp', [Term(
                'major_axis','<',"20000108"), Term("minor_axis=['A', 'B']")])
            expected = wp.truncate(after='20000108').reindex(minor=['A', 'B'])
            tm.assert_panel_equal(result, expected)

            # p4d
            result = store.select('p4d', [Term('major_axis<"20000108"'),
                                          Term("minor_axis=['A', 'B']"),
                                          Term("items=['ItemA', 'ItemB']")])
            expected = p4d.truncate(after='20000108').reindex(
                minor=['A', 'B'], items=['ItemA', 'ItemB'])
            assert_panel4d_equal(result, expected)

            # back compat invalid terms
            terms = [
                dict(field='major_axis', op='>', value='20121114'),
                [ dict(field='major_axis', op='>', value='20121114') ],
                [ "minor_axis=['A','B']", dict(field='major_axis', op='>', value='20121114') ]
                ]
            for t in terms:
                with tm.assert_produces_warning(expected_warning=DeprecationWarning,
                                                check_stacklevel=False):
                    Term(t)

            # valid terms
            terms = [
                ('major_axis=20121114'),
                ('major_axis>20121114'),
                (("major_axis=['20121114', '20121114']"),),
                ('major_axis=datetime.datetime(2012, 11, 14)'),
                'major_axis> 20121114',
                'major_axis >20121114',
                'major_axis > 20121114',
                (("minor_axis=['A', 'B']"),),
                (("minor_axis=['A', 'B']"),),
                ((("minor_axis==['A', 'B']"),),),
                (("items=['ItemA', 'ItemB']"),),
                ('items=ItemA'),
                ]

            for t in terms:
                store.select('wp', t)
                store.select('p4d', t)

            # valid for p4d only
            terms = [
                (("labels=['l1', 'l2']"),),
                Term("labels=['l1', 'l2']"),
                ]

            for t in terms:
                store.select('p4d', t)

            with tm.assertRaisesRegexp(TypeError, 'Only named functions are supported'):
                store.select('wp', Term('major_axis == (lambda x: x)("20130101")'))

            # check USub node parsing
            res = store.select('wpneg', Term('items == -1'))
            expected = Panel({-1: wpneg[-1]})
            tm.assert_panel_equal(res, expected)

            with tm.assertRaisesRegexp(NotImplementedError,
                                       'Unary addition not supported'):
                store.select('wpneg', Term('items == +1'))

    def test_term_compat(self):
        with ensure_clean_store(self.path) as store:

            wp = Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                       major_axis=date_range('1/1/2000', periods=5),
                       minor_axis=['A', 'B', 'C', 'D'])
            store.append('wp',wp)

            result = store.select('wp', [Term('major_axis>20000102'),
                                         Term('minor_axis', '=', ['A','B']) ])
            expected = wp.loc[:,wp.major_axis>Timestamp('20000102'),['A','B']]
            assert_panel_equal(result, expected)

            store.remove('wp', Term('major_axis>20000103'))
            result = store.select('wp')
            expected = wp.loc[:,wp.major_axis<=Timestamp('20000103'),:]
            assert_panel_equal(result, expected)

        with ensure_clean_store(self.path) as store:

            wp = Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                       major_axis=date_range('1/1/2000', periods=5),
                       minor_axis=['A', 'B', 'C', 'D'])
            store.append('wp',wp)

            # stringified datetimes
            result = store.select('wp', [Term('major_axis','>',datetime.datetime(2000,1,2))])
            expected = wp.loc[:,wp.major_axis>Timestamp('20000102')]
            assert_panel_equal(result, expected)

            result = store.select('wp', [Term('major_axis','>',datetime.datetime(2000,1,2,0,0))])
            expected = wp.loc[:,wp.major_axis>Timestamp('20000102')]
            assert_panel_equal(result, expected)

            result = store.select('wp', [Term('major_axis','=',[datetime.datetime(2000,1,2,0,0),datetime.datetime(2000,1,3,0,0)])])
            expected = wp.loc[:,[Timestamp('20000102'),Timestamp('20000103')]]
            assert_panel_equal(result, expected)

            result = store.select('wp', [Term('minor_axis','=',['A','B'])])
            expected = wp.loc[:,:,['A','B']]
            assert_panel_equal(result, expected)

    def test_backwards_compat_without_term_object(self):
        with ensure_clean_store(self.path) as store:

            wp = Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                       major_axis=date_range('1/1/2000', periods=5),
                       minor_axis=['A', 'B', 'C', 'D'])
            store.append('wp',wp)
            with tm.assert_produces_warning(expected_warning=DeprecationWarning,
                                            check_stacklevel=not compat.PY3):
                result = store.select('wp', [('major_axis>20000102'),
                                             ('minor_axis', '=', ['A','B']) ])
            expected = wp.loc[:,wp.major_axis>Timestamp('20000102'),['A','B']]
            assert_panel_equal(result, expected)

            store.remove('wp', ('major_axis>20000103'))
            result = store.select('wp')
            expected = wp.loc[:,wp.major_axis<=Timestamp('20000103'),:]
            assert_panel_equal(result, expected)

        with ensure_clean_store(self.path) as store:

            wp = Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                       major_axis=date_range('1/1/2000', periods=5),
                       minor_axis=['A', 'B', 'C', 'D'])
            store.append('wp',wp)

            # stringified datetimes
            with tm.assert_produces_warning(expected_warning=DeprecationWarning,
                                            check_stacklevel=not compat.PY3):
                result = store.select('wp', [('major_axis','>',datetime.datetime(2000,1,2))])
            expected = wp.loc[:,wp.major_axis>Timestamp('20000102')]
            assert_panel_equal(result, expected)
            with tm.assert_produces_warning(expected_warning=DeprecationWarning,
                                            check_stacklevel=not compat.PY3):
                result = store.select('wp', [('major_axis','>',datetime.datetime(2000,1,2,0,0))])
            expected = wp.loc[:,wp.major_axis>Timestamp('20000102')]
            assert_panel_equal(result, expected)
            with tm.assert_produces_warning(expected_warning=DeprecationWarning,
                                            check_stacklevel=not compat.PY3):
                result = store.select('wp', [('major_axis','=',[datetime.datetime(2000,1,2,0,0),
                                                                datetime.datetime(2000,1,3,0,0)])])
            expected = wp.loc[:,[Timestamp('20000102'),Timestamp('20000103')]]
            assert_panel_equal(result, expected)
            with tm.assert_produces_warning(expected_warning=DeprecationWarning,
                                            check_stacklevel=not compat.PY3):
                result = store.select('wp', [('minor_axis','=',['A','B'])])
            expected = wp.loc[:,:,['A','B']]
            assert_panel_equal(result, expected)

    def test_same_name_scoping(self):

        with ensure_clean_store(self.path) as store:

            import pandas as pd
            df  = DataFrame(np.random.randn(20, 2),index=pd.date_range('20130101',periods=20))
            store.put('df', df, format='table')
            expected = df[df.index>pd.Timestamp('20130105')]

            import datetime
            result = store.select('df','index>datetime.datetime(2013,1,5)')
            assert_frame_equal(result,expected)

            from datetime import datetime

            # technically an error, but allow it
            result = store.select('df','index>datetime.datetime(2013,1,5)')
            assert_frame_equal(result,expected)

            result = store.select('df','index>datetime(2013,1,5)')
            assert_frame_equal(result,expected)

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

        self._check_double_roundtrip(sp, assert_panel_equal,
                                     check_panel_type=True)

        sp2 = p.to_sparse(kind='integer')
        self._check_double_roundtrip(sp2, assert_panel_equal,
                                     check_panel_type=True)

        sp3 = p.to_sparse(fill_value=0)
        self._check_double_roundtrip(sp3, assert_panel_equal,
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

        expected_warning = Warning if compat.PY35 else PerformanceWarning
        with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):
            self._check_roundtrip(DF, tm.assert_frame_equal)

    def test_index_types(self):

        values = np.random.randn(2)

        func = lambda l, r: tm.assert_series_equal(l, r,
                                                   check_dtype=True,
                                                   check_index_type=True,
                                                   check_series_type=True)

        # nose has a deprecation warning in 3.5
        expected_warning = Warning if compat.PY35 else PerformanceWarning
        with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):
            ser = Series(values, [0, 'y'])
            self._check_roundtrip(ser, func)

        with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):
            ser = Series(values, [datetime.datetime.today(), 0])
            self._check_roundtrip(ser, func)

        with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):
            ser = Series(values, ['y', 0])
            self._check_roundtrip(ser, func)

        with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):
            ser = Series(values, [datetime.date.today(), 'a'])
            self._check_roundtrip(ser, func)

        with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):
            ser = Series(values, [1.23, 'b'])
            self._check_roundtrip(ser, func)

        ser = Series(values, [1, 1.53])
        self._check_roundtrip(ser, func)

        ser = Series(values, [1, 5])
        self._check_roundtrip(ser, func)

        ser = Series(values, [datetime.datetime(
            2012, 1, 1), datetime.datetime(2012, 1, 2)])
        self._check_roundtrip(ser, func)

    def test_timeseries_preepoch(self):

        if sys.version_info[0] == 2 and sys.version_info[1] < 7:
            raise nose.SkipTest("won't work on Python < 2.7")

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

        if not skip_compression:
            self._check_roundtrip_table(df, tm.assert_frame_equal,
                                        compression=True)
            self._check_roundtrip(df, tm.assert_frame_equal,
                                  compression=True)

        tdf = tm.makeTimeDataFrame()
        self._check_roundtrip(tdf, tm.assert_frame_equal)

        if not skip_compression:
            self._check_roundtrip(tdf, tm.assert_frame_equal,
                                  compression=True)

        with ensure_clean_store(self.path) as store:
            # not consolidated
            df['foo'] = np.random.randn(len(df))
            store['df'] = df
            recons = store['df']
            self.assertTrue(recons._data.is_consolidated())

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

    def test_empty_series(self):
        for dtype in [np.int64, np.float64, np.object, 'm8[ns]', 'M8[ns]']:
            s = Series(dtype=dtype)
            self._check_roundtrip(s, tm.assert_series_equal)

    def test_can_serialize_dates(self):

        rng = [x.date() for x in bdate_range('1/1/2000', '1/30/2000')]
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
        with ensure_clean_store(self.path) as store:
            store['frame'] = frame
            recons = store['frame']
            assert(recons.index.names == ('foo', 'bar'))

    def test_store_index_name(self):
        df = tm.makeDataFrame()
        df.index.name = 'foo'

        with ensure_clean_store(self.path) as store:
            store['frame'] = df
            recons = store['frame']
            assert(recons.index.name == 'foo')

    def test_store_series_name(self):
        df = tm.makeDataFrame()
        series = df['A']

        with ensure_clean_store(self.path) as store:
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

        with ensure_clean_store(self.path) as store:
            store['obj'] = df1
            tm.assert_frame_equal(store['obj'], df1)
            store['obj'] = df2
            tm.assert_frame_equal(store['obj'], df2)

        # check that can store Series of all of these types
        self._check_roundtrip(df1['obj1'], tm.assert_series_equal)
        self._check_roundtrip(df1['bool1'], tm.assert_series_equal)
        self._check_roundtrip(df1['int1'], tm.assert_series_equal)

        if not skip_compression:
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
        self._check_roundtrip(wp, assert_panel_equal)

    def test_wide_table(self):

        wp = tm.makePanel()
        self._check_roundtrip_table(wp, assert_panel_equal)

    def test_select_with_dups(self):

        # single dtypes
        df = DataFrame(np.random.randn(10,4),columns=['A','A','B','B'])
        df.index = date_range('20130101 9:30',periods=10,freq='T')

        with ensure_clean_store(self.path) as store:
            store.append('df',df)

            result = store.select('df')
            expected = df
            assert_frame_equal(result,expected,by_blocks=True)

            result = store.select('df',columns=df.columns)
            expected = df
            assert_frame_equal(result,expected,by_blocks=True)

            result = store.select('df',columns=['A'])
            expected = df.loc[:,['A']]
            assert_frame_equal(result,expected)

        # dups accross dtypes
        df = concat([DataFrame(np.random.randn(10,4),columns=['A','A','B','B']),
                     DataFrame(np.random.randint(0,10,size=20).reshape(10,2),columns=['A','C'])],
                    axis=1)
        df.index = date_range('20130101 9:30',periods=10,freq='T')

        with ensure_clean_store(self.path) as store:
            store.append('df',df)

            result = store.select('df')
            expected = df
            assert_frame_equal(result,expected,by_blocks=True)

            result = store.select('df',columns=df.columns)
            expected = df
            assert_frame_equal(result,expected,by_blocks=True)

            expected = df.loc[:,['A']]
            result = store.select('df',columns=['A'])
            assert_frame_equal(result,expected,by_blocks=True)

            expected = df.loc[:,['B','A']]
            result = store.select('df',columns=['B','A'])
            assert_frame_equal(result,expected,by_blocks=True)

        # duplicates on both index and columns
        with ensure_clean_store(self.path) as store:
            store.append('df',df)
            store.append('df',df)

            expected = df.loc[:,['B','A']]
            expected = concat([expected, expected])
            result = store.select('df',columns=['B','A'])
            assert_frame_equal(result,expected,by_blocks=True)

    def test_wide_table_dups(self):
        wp = tm.makePanel()
        with ensure_clean_store(self.path) as store:
            store.put('panel', wp, format='table')
            store.put('panel', wp, format='table', append=True)

            with tm.assert_produces_warning(expected_warning=DuplicateWarning):
                recons = store['panel']

            assert_panel_equal(recons, wp)

    def test_long(self):
        def _check(left, right):
            assert_panel_equal(left.to_panel(), right.to_panel())

        wp = tm.makePanel()
        self._check_roundtrip(wp.to_frame(), _check)

        # empty
        # self._check_roundtrip(wp.to_frame()[:0], _check)

    def test_longpanel(self):
        pass

    def test_overwrite_node(self):

        with ensure_clean_store(self.path) as store:
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

        with ensure_clean_store(self.path) as store:

            # put/select ok
            _maybe_remove(store, 'wp')
            store.put('wp', wp, format='table')
            store.select('wp')

            # non-table ok (where = None)
            _maybe_remove(store, 'wp')
            store.put('wp2', wp)
            store.select('wp2')

            # selection on the non-indexable with a large number of columns
            wp = Panel(
                np.random.randn(100, 100, 100), items=['Item%03d' % i for i in range(100)],
                major_axis=date_range('1/1/2000', periods=100), minor_axis=['E%03d' % i for i in range(100)])

            _maybe_remove(store, 'wp')
            store.append('wp', wp)
            items = ['Item%03d' % i for i in range(80)]
            result = store.select('wp', Term('items=items'))
            expected = wp.reindex(items=items)
            assert_panel_equal(expected, result)

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
            result = store.select('df', [("columns=['A', 'B']")])
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

        with ensure_clean_store(self.path) as store:

            # with a Timestamp data column (GH #2637)
            df = DataFrame(dict(ts=bdate_range('2012-01-01', periods=300), A=np.random.randn(300)))
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['ts', 'A'])

            result = store.select('df', [Term("ts>=Timestamp('2012-02-01')")])
            expected = df[df.ts >= Timestamp('2012-02-01')]
            tm.assert_frame_equal(expected, result)

            # bool columns (GH #2849)
            df = DataFrame(np.random.randn(5,2), columns =['A','B'])
            df['object'] = 'foo'
            df.ix[4:5,'object'] = 'bar'
            df['boolv'] = df['A'] > 0
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns = True)

            expected = df[df.boolv == True].reindex(columns=['A','boolv'])
            for v in [True,'true',1]:
                result = store.select('df', Term('boolv == %s' % str(v)), columns = ['A','boolv'])
                tm.assert_frame_equal(expected, result)

            expected = df[df.boolv == False ].reindex(columns=['A','boolv'])
            for v in [False,'false',0]:
                result = store.select('df', Term('boolv == %s' % str(v)), columns = ['A','boolv'])
                tm.assert_frame_equal(expected, result)

            # integer index
            df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
            _maybe_remove(store, 'df_int')
            store.append('df_int', df)
            result = store.select(
                'df_int', [Term("index<10"), Term("columns=['A']")])
            expected = df.reindex(index=list(df.index)[0:10],columns=['A'])
            tm.assert_frame_equal(expected, result)

            # float index
            df = DataFrame(dict(A=np.random.rand(
                        20), B=np.random.rand(20), index=np.arange(20, dtype='f8')))
            _maybe_remove(store, 'df_float')
            store.append('df_float', df)
            result = store.select(
                'df_float', [Term("index<10.0"), Term("columns=['A']")])
            expected = df.reindex(index=list(df.index)[0:10],columns=['A'])
            tm.assert_frame_equal(expected, result)

        with ensure_clean_store(self.path) as store:

            # floats w/o NaN
            df = DataFrame(dict(cols = range(11), values = range(11)),dtype='float64')
            df['cols'] = (df['cols']+10).apply(str)

            store.append('df1',df,data_columns=True)
            result = store.select(
                'df1', where='values>2.0')
            expected = df[df['values']>2.0]
            tm.assert_frame_equal(expected, result)

            # floats with NaN
            df.iloc[0] = np.nan
            expected = df[df['values']>2.0]

            store.append('df2',df,data_columns=True,index=False)
            result = store.select(
                'df2', where='values>2.0')
            tm.assert_frame_equal(expected, result)

            # https://github.com/PyTables/PyTables/issues/282
            # bug in selection when 0th row has a np.nan and an index
            #store.append('df3',df,data_columns=True)
            #result = store.select(
            #    'df3', where='values>2.0')
            #tm.assert_frame_equal(expected, result)

            # not in first position float with NaN ok too
            df = DataFrame(dict(cols = range(11), values = range(11)),dtype='float64')
            df['cols'] = (df['cols']+10).apply(str)

            df.iloc[1] = np.nan
            expected = df[df['values']>2.0]

            store.append('df4',df,data_columns=True)
            result = store.select(
                'df4', where='values>2.0')
            tm.assert_frame_equal(expected, result)

    def test_select_with_many_inputs(self):

        with ensure_clean_store(self.path) as store:

            df = DataFrame(dict(ts=bdate_range('2012-01-01', periods=300),
                                A=np.random.randn(300),
                                B=range(300),
                                users = ['a']*50 + ['b']*50 + ['c']*100 + ['a%03d' % i for i in range(100)]))
            _maybe_remove(store, 'df')
            store.append('df', df, data_columns=['ts', 'A', 'B', 'users'])

            # regular select
            result = store.select('df', [Term("ts>=Timestamp('2012-02-01')")])
            expected = df[df.ts >= Timestamp('2012-02-01')]
            tm.assert_frame_equal(expected, result)

            # small selector
            result = store.select('df', [Term("ts>=Timestamp('2012-02-01') & users=['a','b','c']")])
            expected = df[ (df.ts >= Timestamp('2012-02-01')) & df.users.isin(['a','b','c']) ]
            tm.assert_frame_equal(expected, result)

            # big selector along the columns
            selector = [ 'a','b','c' ] + [ 'a%03d' % i for i in range(60) ]
            result = store.select('df', [Term("ts>=Timestamp('2012-02-01')"),Term('users=selector')])
            expected = df[ (df.ts >= Timestamp('2012-02-01')) & df.users.isin(selector) ]
            tm.assert_frame_equal(expected, result)

            selector = range(100,200)
            result = store.select('df', [Term('B=selector')])
            expected = df[ df.B.isin(selector) ]
            tm.assert_frame_equal(expected, result)
            self.assertEqual(len(result), 100)

            # big selector along the index
            selector = Index(df.ts[0:100].values)
            result  = store.select('df', [Term('ts=selector')])
            expected = df[ df.ts.isin(selector.values) ]
            tm.assert_frame_equal(expected, result)
            self.assertEqual(len(result), 100)

    def test_select_iterator(self):

        # single table
        with ensure_clean_store(self.path) as store:

            df = tm.makeTimeDataFrame(500)
            _maybe_remove(store, 'df')
            store.append('df', df)

            expected = store.select('df')

            results = [ s for s in store.select('df',iterator=True) ]
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            results = [ s for s in store.select('df',chunksize=100) ]
            self.assertEqual(len(results), 5)
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            results = [ s for s in store.select('df',chunksize=150) ]
            result = concat(results)
            tm.assert_frame_equal(result, expected)

        with ensure_clean_path(self.path) as path:

            df = tm.makeTimeDataFrame(500)
            df.to_hdf(path,'df_non_table')
            self.assertRaises(TypeError, read_hdf, path,'df_non_table',chunksize=100)
            self.assertRaises(TypeError, read_hdf, path,'df_non_table',iterator=True)

        with ensure_clean_path(self.path) as path:

            df = tm.makeTimeDataFrame(500)
            df.to_hdf(path,'df',format='table')

            results = [ s for s in read_hdf(path,'df',chunksize=100) ]
            result = concat(results)

            self.assertEqual(len(results), 5)
            tm.assert_frame_equal(result, df)
            tm.assert_frame_equal(result, read_hdf(path,'df'))

        # multiple

        with ensure_clean_store(self.path) as store:

            df1 = tm.makeTimeDataFrame(500)
            store.append('df1',df1,data_columns=True)
            df2 = tm.makeTimeDataFrame(500).rename(columns=lambda x: "%s_2" % x)
            df2['foo'] = 'bar'
            store.append('df2',df2)

            df = concat([df1, df2], axis=1)

            # full selection
            expected = store.select_as_multiple(
                ['df1', 'df2'], selector='df1')
            results = [ s for s in store.select_as_multiple(
                ['df1', 'df2'], selector='df1', chunksize=150) ]
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

    def test_select_iterator_complete_8014(self):

        # GH 8014
        # using iterator and where clause
        chunksize=1e4

        # no iterator
        with ensure_clean_store(self.path) as store:

            expected = tm.makeTimeDataFrame(100064, 'S')
            _maybe_remove(store, 'df')
            store.append('df',expected)

            beg_dt = expected.index[0]
            end_dt = expected.index[-1]

            # select w/o iteration and no where clause works
            result = store.select('df')
            tm.assert_frame_equal(expected, result)

            # select w/o iterator and where clause, single term, begin
            # of range, works
            where = "index >= '%s'" % beg_dt
            result = store.select('df',where=where)
            tm.assert_frame_equal(expected, result)

            # select w/o iterator and where clause, single term, end
            # of range, works
            where = "index <= '%s'" % end_dt
            result = store.select('df',where=where)
            tm.assert_frame_equal(expected, result)

            # select w/o iterator and where clause, inclusive range,
            # works
            where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
            result = store.select('df',where=where)
            tm.assert_frame_equal(expected, result)

        # with iterator, full range
        with ensure_clean_store(self.path) as store:

            expected = tm.makeTimeDataFrame(100064, 'S')
            _maybe_remove(store, 'df')
            store.append('df',expected)

            beg_dt = expected.index[0]
            end_dt = expected.index[-1]

            # select w/iterator and no where clause works
            results = [ s for s in store.select('df',chunksize=chunksize) ]
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            # select w/iterator and where clause, single term, begin of range
            where = "index >= '%s'" % beg_dt
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            # select w/iterator and where clause, single term, end of range
            where = "index <= '%s'" % end_dt
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            result = concat(results)
            tm.assert_frame_equal(expected, result)

            # select w/iterator and where clause, inclusive range
            where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            result = concat(results)
            tm.assert_frame_equal(expected, result)

    def test_select_iterator_non_complete_8014(self):

        # GH 8014
        # using iterator and where clause
        chunksize=1e4

        # with iterator, non complete range
        with ensure_clean_store(self.path) as store:

            expected = tm.makeTimeDataFrame(100064, 'S')
            _maybe_remove(store, 'df')
            store.append('df',expected)

            beg_dt = expected.index[1]
            end_dt = expected.index[-2]

            # select w/iterator and where clause, single term, begin of range
            where = "index >= '%s'" % beg_dt
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            result = concat(results)
            rexpected = expected[expected.index >= beg_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, single term, end of range
            where = "index <= '%s'" % end_dt
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            result = concat(results)
            rexpected = expected[expected.index <= end_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, inclusive range
            where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            result = concat(results)
            rexpected = expected[(expected.index >= beg_dt) & (expected.index <= end_dt)]
            tm.assert_frame_equal(rexpected, result)

        # with iterator, empty where
        with ensure_clean_store(self.path) as store:

            expected = tm.makeTimeDataFrame(100064, 'S')
            _maybe_remove(store, 'df')
            store.append('df',expected)

            end_dt = expected.index[-1]

            # select w/iterator and where clause, single term, begin of range
            where = "index > '%s'" % end_dt
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            self.assertEqual(0, len(results))

    def test_select_iterator_many_empty_frames(self):

        # GH 8014
        # using iterator and where clause can return many empty
        # frames.
        chunksize=int(1e4)

        # with iterator, range limited to the first chunk
        with ensure_clean_store(self.path) as store:

            expected = tm.makeTimeDataFrame(100000, 'S')
            _maybe_remove(store, 'df')
            store.append('df',expected)

            beg_dt = expected.index[0]
            end_dt = expected.index[chunksize-1]

            # select w/iterator and where clause, single term, begin of range
            where = "index >= '%s'" % beg_dt
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]
            result = concat(results)
            rexpected = expected[expected.index >= beg_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, single term, end of range
            where = "index <= '%s'" % end_dt
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]

            tm.assert_equal(1, len(results))
            result = concat(results)
            rexpected = expected[expected.index <= end_dt]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause, inclusive range
            where = "index >= '%s' & index <= '%s'" % (beg_dt, end_dt)
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]

            # should be 1, is 10
            tm.assert_equal(1, len(results))
            result = concat(results)
            rexpected = expected[(expected.index >= beg_dt) & (expected.index <= end_dt)]
            tm.assert_frame_equal(rexpected, result)

            # select w/iterator and where clause which selects
            # *nothing*.
            #
            # To be consistent with Python idiom I suggest this should
            # return [] e.g. `for e in []: print True` never prints
            # True.

            where = "index <= '%s' & index >= '%s'" % (beg_dt, end_dt)
            results = [ s for s in store.select('df',where=where,chunksize=chunksize) ]

            # should be []
            tm.assert_equal(0, len(results))


    def test_retain_index_attributes(self):

        # GH 3499, losing frequency info on index recreation
        df = DataFrame(dict(A = Series(lrange(3),
                                       index=date_range('2000-1-1',periods=3,freq='H'))))

        with ensure_clean_store(self.path) as store:
            _maybe_remove(store,'data')
            store.put('data', df, format='table')

            result = store.get('data')
            tm.assert_frame_equal(df,result)

            for attr in ['freq','tz','name']:
                for idx in ['index','columns']:
                    self.assertEqual(getattr(getattr(df,idx),attr,None),
                                     getattr(getattr(result,idx),attr,None))


            # try to append a table with a different frequency
            with tm.assert_produces_warning(expected_warning=AttributeConflictWarning):
                df2 = DataFrame(dict(A = Series(lrange(3),
                                                index=date_range('2002-1-1',periods=3,freq='D'))))
                store.append('data',df2)

            self.assertIsNone(store.get_storer('data').info['index']['freq'])

            # this is ok
            _maybe_remove(store,'df2')
            df2 = DataFrame(dict(A = Series(lrange(3),
                                            index=[Timestamp('20010101'),Timestamp('20010102'),Timestamp('20020101')])))
            store.append('df2',df2)
            df3 = DataFrame(dict(A = Series(lrange(3),index=date_range('2002-1-1',periods=3,freq='D'))))
            store.append('df2',df3)

    def test_retain_index_attributes2(self):

        with ensure_clean_path(self.path) as path:

            expected_warning = Warning if compat.PY35 else AttributeConflictWarning
            with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):

                df  = DataFrame(dict(A = Series(lrange(3), index=date_range('2000-1-1',periods=3,freq='H'))))
                df.to_hdf(path,'data',mode='w',append=True)
                df2 = DataFrame(dict(A = Series(lrange(3), index=date_range('2002-1-1',periods=3,freq='D'))))
                df2.to_hdf(path,'data',append=True)

                idx = date_range('2000-1-1',periods=3,freq='H')
                idx.name = 'foo'
                df  = DataFrame(dict(A = Series(lrange(3), index=idx)))
                df.to_hdf(path,'data',mode='w',append=True)

            self.assertEqual(read_hdf(path,'data').index.name, 'foo')

            with tm.assert_produces_warning(expected_warning=expected_warning, check_stacklevel=False):

                idx2 = date_range('2001-1-1',periods=3,freq='H')
                idx2.name = 'bar'
                df2 = DataFrame(dict(A = Series(lrange(3), index=idx2)))
                df2.to_hdf(path,'data',append=True)

            self.assertIsNone(read_hdf(path,'data').index.name)

    def test_panel_select(self):

        wp = tm.makePanel()

        with ensure_clean_store(self.path) as store:
            store.put('wp', wp, format='table')
            date = wp.major_axis[len(wp.major_axis) // 2]

            crit1 = ('major_axis>=date')
            crit2 = ("minor_axis=['A', 'D']")

            result = store.select('wp', [crit1, crit2])
            expected = wp.truncate(before=date).reindex(minor=['A', 'D'])
            assert_panel_equal(result, expected)

            result = store.select(
                'wp', ['major_axis>="20000124"', ("minor_axis=['A', 'B']")])
            expected = wp.truncate(before='20000124').reindex(minor=['A', 'B'])
            assert_panel_equal(result, expected)

    def test_frame_select(self):

        df = tm.makeTimeDataFrame()

        with ensure_clean_store(self.path) as store:
            store.put('frame', df,format='table')
            date = df.index[len(df) // 2]

            crit1 = Term('index>=date')
            self.assertEqual(crit1.env.scope['date'], date)

            crit2 = ("columns=['A', 'D']")
            crit3 = ('columns=A')

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

    def test_frame_select_complex(self):
        # select via complex criteria

        df = tm.makeTimeDataFrame()
        df['string'] = 'foo'
        df.loc[df.index[0:4],'string'] = 'bar'

        with ensure_clean_store(self.path) as store:
            store.put('df', df, format='table', data_columns=['string'])

            # empty
            result = store.select('df', 'index>df.index[3] & string="bar"')
            expected = df.loc[(df.index>df.index[3]) & (df.string=='bar')]
            tm.assert_frame_equal(result, expected)

            result = store.select('df', 'index>df.index[3] & string="foo"')
            expected = df.loc[(df.index>df.index[3]) & (df.string=='foo')]
            tm.assert_frame_equal(result, expected)

            # or
            result = store.select('df', 'index>df.index[3] | string="bar"')
            expected = df.loc[(df.index>df.index[3]) | (df.string=='bar')]
            tm.assert_frame_equal(result, expected)

            result = store.select('df', '(index>df.index[3] & index<=df.index[6]) | string="bar"')
            expected = df.loc[((df.index>df.index[3]) & (df.index<=df.index[6])) | (df.string=='bar')]
            tm.assert_frame_equal(result, expected)

            # invert
            result = store.select('df', 'string!="bar"')
            expected = df.loc[df.string!='bar']
            tm.assert_frame_equal(result, expected)

            # invert not implemented in numexpr :(
            self.assertRaises(NotImplementedError, store.select, 'df', '~(string="bar")')

            # invert ok for filters
            result = store.select('df', "~(columns=['A','B'])")
            expected = df.loc[:,df.columns.difference(['A','B'])]
            tm.assert_frame_equal(result, expected)

            # in
            result = store.select('df', "index>df.index[3] & columns in ['A','B']")
            expected = df.loc[df.index>df.index[3]].reindex(columns=['A','B'])
            tm.assert_frame_equal(result, expected)

    def test_frame_select_complex2(self):

        with ensure_clean_path(['parms.hdf','hist.hdf']) as paths:

            pp, hh = paths

            # use non-trivial selection criteria
            parms = DataFrame({ 'A' : [1,1,2,2,3] })
            parms.to_hdf(pp,'df',mode='w',format='table',data_columns=['A'])

            selection = read_hdf(pp,'df',where='A=[2,3]')
            hist = DataFrame(np.random.randn(25,1),columns=['data'],
                             index=MultiIndex.from_tuples([ (i,j) for i in range(5) for j in range(5) ],
                                                          names=['l1','l2']))

            hist.to_hdf(hh,'df',mode='w',format='table')

            expected = read_hdf(hh,'df',where=Term('l1','=',[2,3,4]))

            # list like
            result = read_hdf(hh,'df',where=Term('l1','=',selection.index.tolist()))
            assert_frame_equal(result, expected)
            l = selection.index.tolist()

            # sccope with list like
            store = HDFStore(hh)
            result = store.select('df',where='l1=l')
            assert_frame_equal(result, expected)
            store.close()

            result = read_hdf(hh,'df',where='l1=l')
            assert_frame_equal(result, expected)

            # index
            index = selection.index
            result = read_hdf(hh,'df',where='l1=index')
            assert_frame_equal(result, expected)

            result = read_hdf(hh,'df',where='l1=selection.index')
            assert_frame_equal(result, expected)

            result = read_hdf(hh,'df',where='l1=selection.index.tolist()')
            assert_frame_equal(result, expected)

            result = read_hdf(hh,'df',where='l1=list(selection.index)')
            assert_frame_equal(result, expected)

            # sccope with index
            store = HDFStore(hh)

            result = store.select('df',where='l1=index')
            assert_frame_equal(result, expected)

            result = store.select('df',where='l1=selection.index')
            assert_frame_equal(result, expected)

            result = store.select('df',where='l1=selection.index.tolist()')
            assert_frame_equal(result, expected)

            result = store.select('df',where='l1=list(selection.index)')
            assert_frame_equal(result, expected)

            store.close()

    def test_invalid_filtering(self):

        # can't use more than one filter (atm)

        df = tm.makeTimeDataFrame()

        with ensure_clean_store(self.path) as store:
            store.put('df', df, format='table')

            # not implemented
            self.assertRaises(NotImplementedError, store.select, 'df', "columns=['A'] | columns=['B']")

            # in theory we could deal with this
            self.assertRaises(NotImplementedError, store.select, 'df', "columns=['A','B'] & columns=['C']")

    def test_string_select(self):
        # GH 2973
        with ensure_clean_store(self.path) as store:

            df = tm.makeTimeDataFrame()

            # test string ==/!=
            df['x'] = 'none'
            df.ix[2:7,'x'] = ''

            store.append('df',df,data_columns=['x'])

            result = store.select('df',Term('x=none'))
            expected = df[df.x == 'none']
            assert_frame_equal(result,expected)

            try:
                result = store.select('df',Term('x!=none'))
                expected = df[df.x != 'none']
                assert_frame_equal(result,expected)
            except Exception as detail:
                com.pprint_thing("[{0}]".format(detail))
                com.pprint_thing(store)
                com.pprint_thing(expected)

            df2 = df.copy()
            df2.loc[df2.x=='','x'] = np.nan

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

        with ensure_clean_store(self.path) as store:
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
            self.assertIsInstance(result, Series)

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

            # start/stop
            result = store.select_column('df3', 'string', start=2)
            tm.assert_almost_equal(result.values, df3['string'].values[2:])

            result = store.select_column('df3', 'string', start=-2)
            tm.assert_almost_equal(result.values, df3['string'].values[-2:])

            result = store.select_column('df3', 'string', stop=2)
            tm.assert_almost_equal(result.values, df3['string'].values[:2])

            result = store.select_column('df3', 'string', stop=-2)
            tm.assert_almost_equal(result.values, df3['string'].values[:-2])

            result = store.select_column('df3', 'string', start=2, stop=-2)
            tm.assert_almost_equal(result.values, df3['string'].values[2:-2])

            result = store.select_column('df3', 'string', start=-2, stop=2)
            tm.assert_almost_equal(result.values, df3['string'].values[-2:2])

            # GH 10392 - make sure column name is preserved
            df4 = DataFrame({'A': np.random.randn(10), 'B': 'foo'})
            store.append('df4', df4, data_columns=True)
            expected = df4['B']
            result = store.select_column('df4', 'B')
            tm.assert_series_equal(result, expected)


    def test_coordinates(self):
        df = tm.makeTimeDataFrame()

        with ensure_clean_store(self.path) as store:

            _maybe_remove(store, 'df')
            store.append('df', df)

            # all
            c = store.select_as_coordinates('df')
            assert((c.values == np.arange(len(df.index))).all() == True)

            # get coordinates back & test vs frame
            _maybe_remove(store, 'df')

            df = DataFrame(dict(A=lrange(5), B=lrange(5)))
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
            self.assertIsInstance(c, Index)

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

        # pass array/mask as the coordinates
        with ensure_clean_store(self.path) as store:

            df = DataFrame(np.random.randn(1000,2),index=date_range('20000101',periods=1000))
            store.append('df',df)
            c = store.select_column('df','index')
            where = c[DatetimeIndex(c).month==5].index
            expected = df.iloc[where]

            # locations
            result = store.select('df',where=where)
            tm.assert_frame_equal(result,expected)

            # boolean
            result = store.select('df',where=where)
            tm.assert_frame_equal(result,expected)

            # invalid
            self.assertRaises(ValueError, store.select, 'df',where=np.arange(len(df),dtype='float64'))
            self.assertRaises(ValueError, store.select, 'df',where=np.arange(len(df)+1))
            self.assertRaises(ValueError, store.select, 'df',where=np.arange(len(df)),start=5)
            self.assertRaises(ValueError, store.select, 'df',where=np.arange(len(df)),start=5,stop=10)

            # selection with filter
            selection = date_range('20000101',periods=500)
            result = store.select('df', where='index in selection')
            expected = df[df.index.isin(selection)]
            tm.assert_frame_equal(result,expected)

            # list
            df = DataFrame(np.random.randn(10,2))
            store.append('df2',df)
            result = store.select('df2',where=[0,3,5])
            expected = df.iloc[[0,3,5]]
            tm.assert_frame_equal(result,expected)

            # boolean
            where = [True] * 10
            where[-2] = False
            result = store.select('df2',where=where)
            expected = df.loc[where]
            tm.assert_frame_equal(result,expected)

            # start/stop
            result = store.select('df2', start=5, stop=10)
            expected = df[5:10]
            tm.assert_frame_equal(result,expected)

    def test_append_to_multiple(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
        df2['foo'] = 'bar'
        df = concat([df1, df2], axis=1)

        with ensure_clean_store(self.path) as store:

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

    def test_append_to_multiple_dropna(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
        df1.ix[1, ['A', 'B']] = np.nan
        df = concat([df1, df2], axis=1)

        with ensure_clean_store(self.path) as store:
            # dropna=True should guarantee rows are synchronized
            store.append_to_multiple(
                {'df1': ['A', 'B'], 'df2': None}, df, selector='df1',
                dropna=True)
            result = store.select_as_multiple(['df1', 'df2'])
            expected = df.dropna()
            tm.assert_frame_equal(result, expected)
            tm.assert_index_equal(store.select('df1').index,
                                  store.select('df2').index)

            # dropna=False shouldn't synchronize row indexes
            store.append_to_multiple(
                {'df1': ['A', 'B'], 'df2': None}, df, selector='df1',
                dropna=False)
            self.assertRaises(
                ValueError, store.select_as_multiple, ['df1', 'df2'])
            assert not store.select('df1').index.equals(
                store.select('df2').index)

    def test_select_as_multiple(self):

        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame().rename(columns=lambda x: "%s_2" % x)
        df2['foo'] = 'bar'

        with ensure_clean_store(self.path) as store:

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
            self.assertRaises(KeyError, store.select_as_multiple,
                              ['df1','df3'], where=['A>0', 'B>0'], selector='df1')
            self.assertRaises(KeyError, store.select_as_multiple,
                              ['df3'], where=['A>0', 'B>0'], selector='df1')
            self.assertRaises(KeyError, store.select_as_multiple,
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
            result = store.select_as_multiple(['df1', 'df2'], where=[Term(
                'index>df2.index[4]')], selector='df2')
            expected = concat([df1, df2], axis=1)
            expected = expected[5:]
            tm.assert_frame_equal(result, expected)

            # test excpection for diff rows
            store.append('df3', tm.makeTimeDataFrame(nper=50))
            self.assertRaises(ValueError, store.select_as_multiple,
                              ['df1','df3'], where=['A>0', 'B>0'], selector='df1')

    def test_nan_selection_bug_4858(self):

        # GH 4858; nan selection bug, only works for pytables >= 3.1
        if LooseVersion(tables.__version__) < '3.1.0':
            raise nose.SkipTest('tables version does not support fix for nan selection bug: GH 4858')

        with ensure_clean_store(self.path) as store:

            df = DataFrame(dict(cols = range(6), values = range(6)), dtype='float64')
            df['cols'] = (df['cols']+10).apply(str)
            df.iloc[0] = np.nan

            expected = DataFrame(dict(cols = ['13.0','14.0','15.0'], values = [3.,4.,5.]), index=[3,4,5])

            # write w/o the index on that particular column
            store.append('df',df, data_columns=True,index=['cols'])
            result = store.select('df',where='values>2.0')
            assert_frame_equal(result,expected)

    def test_start_stop(self):

        with ensure_clean_store(self.path) as store:

            df = DataFrame(dict(A=np.random.rand(20), B=np.random.rand(20)))
            store.append('df', df)

            result = store.select(
                'df', [Term("columns=['A']")], start=0, stop=5)
            expected = df.ix[0:4, ['A']]
            tm.assert_frame_equal(result, expected)

            # out of range
            result = store.select(
                'df', [Term("columns=['A']")], start=30, stop=40)
            assert(len(result) == 0)
            assert(type(result) == DataFrame)

    def test_select_filter_corner(self):

        df = DataFrame(np.random.randn(50, 100))
        df.index = ['%.3d' % c for c in df.index]
        df.columns = ['%.3d' % c for c in df.columns]

        with ensure_clean_store(self.path) as store:
            store.put('frame', df, format='table')

            crit = Term('columns=df.columns[:75]')
            result = store.select('frame', [crit])
            tm.assert_frame_equal(result, df.ix[:, df.columns[:75]])

            crit = Term('columns=df.columns[:75:2]')
            result = store.select('frame', [crit])
            tm.assert_frame_equal(result, df.ix[:, df.columns[:75:2]])

    def _check_roundtrip(self, obj, comparator, compression=False, **kwargs):

        options = {}
        if compression:
            options['complib'] = _default_compressor

        with ensure_clean_store(self.path, 'w', **options) as store:
            store['obj'] = obj
            retrieved = store['obj']
            comparator(retrieved, obj, **kwargs)

    def _check_double_roundtrip(self, obj, comparator, compression=False,
                                **kwargs):
        options = {}
        if compression:
            options['complib'] = compression or _default_compressor

        with ensure_clean_store(self.path, 'w', **options) as store:
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

        with ensure_clean_store(self.path, 'w', **options) as store:
            store.put('obj', obj, format='table')
            retrieved = store['obj']
            # sorted_obj = _test_sort(obj)
            comparator(retrieved, obj)

    def test_multiple_open_close(self):
        # GH 4409, open & close multiple times

        with ensure_clean_path(self.path) as path:

            df = tm.makeDataFrame()
            df.to_hdf(path,'df',mode='w',format='table')

            # single
            store = HDFStore(path)
            self.assertNotIn('CLOSED', str(store))
            self.assertTrue(store.is_open)
            store.close()
            self.assertIn('CLOSED', str(store))
            self.assertFalse(store.is_open)

        with ensure_clean_path(self.path) as path:

            if pytables._table_file_open_policy_is_strict:

                # multiples
                store1 = HDFStore(path)
                def f():
                    HDFStore(path)
                self.assertRaises(ValueError, f)
                store1.close()

            else:

                # multiples
                store1 = HDFStore(path)
                store2 = HDFStore(path)

                self.assertNotIn('CLOSED', str(store1))
                self.assertNotIn('CLOSED', str(store2))
                self.assertTrue(store1.is_open)
                self.assertTrue(store2.is_open)

                store1.close()
                self.assertIn('CLOSED', str(store1))
                self.assertFalse(store1.is_open)
                self.assertNotIn('CLOSED', str(store2))
                self.assertTrue(store2.is_open)

                store2.close()
                self.assertIn('CLOSED', str(store1))
                self.assertIn('CLOSED', str(store2))
                self.assertFalse(store1.is_open)
                self.assertFalse(store2.is_open)

                # nested close
                store = HDFStore(path,mode='w')
                store.append('df',df)

                store2 = HDFStore(path)
                store2.append('df2',df)
                store2.close()
                self.assertIn('CLOSED', str(store2))
                self.assertFalse(store2.is_open)

                store.close()
                self.assertIn('CLOSED', str(store))
                self.assertFalse(store.is_open)

                # double closing
                store = HDFStore(path,mode='w')
                store.append('df', df)

                store2 = HDFStore(path)
                store.close()
                self.assertIn('CLOSED', str(store))
                self.assertFalse(store.is_open)

                store2.close()
                self.assertIn('CLOSED', str(store2))
                self.assertFalse(store2.is_open)

        # ops on a closed store
        with ensure_clean_path(self.path) as path:

            df = tm.makeDataFrame()
            df.to_hdf(path,'df',mode='w',format='table')

            store = HDFStore(path)
            store.close()

            self.assertRaises(ClosedFileError, store.keys)
            self.assertRaises(ClosedFileError, lambda : 'df' in store)
            self.assertRaises(ClosedFileError, lambda : len(store))
            self.assertRaises(ClosedFileError, lambda : store['df'])
            self.assertRaises(ClosedFileError, lambda : store.df)
            self.assertRaises(ClosedFileError, store.select, 'df')
            self.assertRaises(ClosedFileError, store.get, 'df')
            self.assertRaises(ClosedFileError, store.append, 'df2', df)
            self.assertRaises(ClosedFileError, store.put, 'df3', df)
            self.assertRaises(ClosedFileError, store.get_storer, 'df2')
            self.assertRaises(ClosedFileError, store.remove, 'df2')

            def f():
                store.select('df')
            tm.assertRaisesRegexp(ClosedFileError, 'file is not open', f)

    def test_pytables_native_read(self):

        with ensure_clean_store(tm.get_data_path('legacy_hdf/pytables_native.h5'), mode='r') as store:
            d2 = store['detector/readout']
            self.assertIsInstance(d2, DataFrame)

        with ensure_clean_store(tm.get_data_path('legacy_hdf/pytables_native2.h5'), mode='r') as store:
            str(store)
            d1 = store['detector']
            self.assertIsInstance(d1, DataFrame)

    def test_legacy_read(self):
        with ensure_clean_store(tm.get_data_path('legacy_hdf/legacy.h5'), mode='r') as store:
            store['a']
            store['b']
            store['c']
            store['d']

    def test_legacy_table_read(self):
        # legacy table types
        with ensure_clean_store(tm.get_data_path('legacy_hdf/legacy_table.h5'), mode='r') as store:
            store.select('df1')
            store.select('df2')
            store.select('wp1')

            # force the frame
            store.select('df2', typ='legacy_frame')

            # old version warning
            with tm.assert_produces_warning(expected_warning=IncompatibilityWarning):
                self.assertRaises(
                    Exception, store.select, 'wp1', Term('minor_axis=B'))

                df2 = store.select('df2')
                result = store.select('df2', Term('index>df2.index[2]'))
                expected = df2[df2.index > df2.index[2]]
                assert_frame_equal(expected, result)

    def test_legacy_0_10_read(self):
        # legacy from 0.10
        with ensure_clean_store(tm.get_data_path('legacy_hdf/legacy_0.10.h5'), mode='r') as store:
            str(store)
            for k in store.keys():
                store.select(k)

    def test_legacy_0_11_read(self):
        # legacy from 0.11
        path = os.path.join('legacy_hdf', 'legacy_table_0.11.h5')
        with ensure_clean_store(tm.get_data_path(path), mode='r') as store:
            str(store)
            assert 'df' in store
            assert 'df1' in store
            assert 'mi' in store
            df = store.select('df')
            df1 = store.select('df1')
            mi = store.select('mi')
            assert isinstance(df, DataFrame)
            assert isinstance(df1, DataFrame)
            assert isinstance(mi, DataFrame)

    def test_copy(self):

        def do_copy(f = None, new_f = None, keys = None, propindexes = True, **kwargs):
            try:
                if f is None:
                    f = tm.get_data_path(os.path.join('legacy_hdf',
                                                      'legacy_0.10.h5'))


                store = HDFStore(f, 'r')

                if new_f is None:
                    import tempfile
                    fd, new_f = tempfile.mkstemp()

                tstore = store.copy(new_f, keys = keys, propindexes = propindexes, **kwargs)

                # check keys
                if keys is None:
                    keys = store.keys()
                self.assertEqual(set(keys), set(tstore.keys()))

                # check indicies & nrows
                for k in tstore.keys():
                    if tstore.get_storer(k).is_table:
                        new_t = tstore.get_storer(k)
                        orig_t = store.get_storer(k)

                        self.assertEqual(orig_t.nrows, new_t.nrows)

                        # check propindixes
                        if propindexes:
                            for a in orig_t.axes:
                                if a.is_indexed:
                                    self.assertTrue(new_t[a.name].is_indexed)

            finally:
                safe_close(store)
                safe_close(tstore)
                try:
                    os.close(fd)
                except:
                    pass
                safe_remove(new_f)

        do_copy()
        do_copy(keys = ['/a','/b','/df1_mixed'])
        do_copy(propindexes = False)

        # new table
        df = tm.makeDataFrame()

        try:
            path = create_tempfile(self.path)
            st = HDFStore(path)
            st.append('df', df, data_columns = ['A'])
            st.close()
            do_copy(f = path)
            do_copy(f = path, propindexes = False)
        finally:
            safe_remove(path)

    def test_legacy_table_write(self):
        raise nose.SkipTest("cannot write legacy tables")

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

        df = DataFrame(dict(A = 'foo', B = 'bar'),index=lrange(10))
        store.append('df', df, data_columns = ['B'], min_itemsize={'A' : 200 })
        store.append('wp', wp)

        store.close()

    def test_store_datetime_fractional_secs(self):

        with ensure_clean_store(self.path) as store:
            dt = datetime.datetime(2012, 1, 2, 3, 4, 5, 123456)
            series = Series([0], [dt])
            store['a'] = series
            self.assertEqual(store['a'].index[0], dt)

    def test_tseries_indices_series(self):

        with ensure_clean_store(self.path) as store:
            idx = tm.makeDateIndex(10)
            ser = Series(np.random.randn(len(idx)), idx)
            store['a'] = ser
            result = store['a']

            assert_series_equal(result, ser)
            self.assertEqual(type(result.index), type(ser.index))
            self.assertEqual(result.index.freq, ser.index.freq)

            idx = tm.makePeriodIndex(10)
            ser = Series(np.random.randn(len(idx)), idx)
            store['a'] = ser
            result = store['a']

            assert_series_equal(result, ser)
            self.assertEqual(type(result.index), type(ser.index))
            self.assertEqual(result.index.freq, ser.index.freq)

    def test_tseries_indices_frame(self):

        with ensure_clean_store(self.path) as store:
            idx = tm.makeDateIndex(10)
            df = DataFrame(np.random.randn(len(idx), 3), index=idx)
            store['a'] = df
            result = store['a']

            assert_frame_equal(result, df)
            self.assertEqual(type(result.index), type(df.index))
            self.assertEqual(result.index.freq, df.index.freq)

            idx = tm.makePeriodIndex(10)
            df = DataFrame(np.random.randn(len(idx), 3), idx)
            store['a'] = df
            result = store['a']

            assert_frame_equal(result, df)
            self.assertEqual(type(result.index), type(df.index))
            self.assertEqual(result.index.freq, df.index.freq)

    def test_unicode_index(self):

        unicode_values = [u('\u03c3'), u('\u03c3\u03c3')]
        def f():
            s = Series(np.random.randn(len(unicode_values)), unicode_values)
            self._check_roundtrip(s, tm.assert_series_equal)

        compat_assert_produces_warning(PerformanceWarning, f)

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

    #    self.assertRaises(Exception, store.put, 'foo', df, format='table')

    def test_append_with_diff_col_name_types_raises_value_error(self):
        df = DataFrame(np.random.randn(10, 1))
        df2 = DataFrame({'a': np.random.randn(10)})
        df3 = DataFrame({(1, 2): np.random.randn(10)})
        df4 = DataFrame({('1', 2): np.random.randn(10)})
        df5 = DataFrame({('1', 2, object): np.random.randn(10)})

        with ensure_clean_store(self.path) as store:
            name = 'df_%s' % tm.rands(10)
            store.append(name, df)

            for d in (df2, df3, df4, df5):
                with tm.assertRaises(ValueError):
                    store.append(name, d)

    def test_query_with_nested_special_character(self):
        df = DataFrame({'a': ['a', 'a', 'c', 'b', 'test & test', 'c' , 'b', 'e'],
                        'b': [1, 2, 3, 4, 5, 6, 7, 8]})
        expected = df[df.a == 'test & test']
        with ensure_clean_store(self.path) as store:
            store.append('test', df, format='table', data_columns=True)
            result = store.select('test', 'a = "test & test"')
        tm.assert_frame_equal(expected, result)

    def test_categorical(self):

        with ensure_clean_store(self.path) as store:

            # basic
            _maybe_remove(store, 's')
            s = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c'], categories=['a','b','c','d'], ordered=False))
            store.append('s', s, format='table')
            result = store.select('s')
            tm.assert_series_equal(s, result)

            _maybe_remove(store, 's_ordered')
            s = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c'], categories=['a','b','c','d'], ordered=True))
            store.append('s_ordered', s, format='table')
            result = store.select('s_ordered')
            tm.assert_series_equal(s, result)

            _maybe_remove(store, 'df')
            df = DataFrame({"s":s, "vals":[1,2,3,4,5,6]})
            store.append('df', df, format='table')
            result = store.select('df')
            tm.assert_frame_equal(result, df)

            # dtypes
            s = Series([1,1,2,2,3,4,5]).astype('category')
            store.append('si',s)
            result = store.select('si')
            tm.assert_series_equal(result, s)

            s = Series([1,1,np.nan,2,3,4,5]).astype('category')
            store.append('si2',s)
            result = store.select('si2')
            tm.assert_series_equal(result, s)

            # multiple
            df2 = df.copy()
            df2['s2'] = Series(list('abcdefg')).astype('category')
            store.append('df2',df2)
            result = store.select('df2')
            tm.assert_frame_equal(result, df2)

            # make sure the metadata is ok
            self.assertTrue('/df2   ' in str(store))
            self.assertTrue('/df2/meta/values_block_0/meta' in str(store))
            self.assertTrue('/df2/meta/values_block_1/meta' in str(store))

            # unordered
            s = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c'], categories=['a','b','c','d'],ordered=False))
            store.append('s2', s, format='table')
            result = store.select('s2')
            tm.assert_series_equal(result, s)

            # query
            store.append('df3', df, data_columns=['s'])
            expected = df[df.s.isin(['b','c'])]
            result = store.select('df3', where = ['s in ["b","c"]'])
            tm.assert_frame_equal(result, expected)

            expected = df[df.s.isin(['b','c'])]
            result = store.select('df3', where = ['s = ["b","c"]'])
            tm.assert_frame_equal(result, expected)

            expected = df[df.s.isin(['d'])]
            result = store.select('df3', where = ['s in ["d"]'])
            tm.assert_frame_equal(result, expected)

            expected = df[df.s.isin(['f'])]
            result = store.select('df3', where = ['s in ["f"]'])
            tm.assert_frame_equal(result, expected)

            # appending with same categories is ok
            store.append('df3', df)

            df = concat([df,df])
            expected = df[df.s.isin(['b','c'])]
            result = store.select('df3', where = ['s in ["b","c"]'])
            tm.assert_frame_equal(result, expected)

            # appending must have the same categories
            df3 = df.copy()
            df3['s'].cat.remove_unused_categories(inplace=True)

            self.assertRaises(ValueError, lambda : store.append('df3', df3))

            # remove
            # make sure meta data is removed (its a recursive removal so should be)
            result = store.select('df3/meta/s/meta')
            self.assertIsNotNone(result)
            store.remove('df3')
            self.assertRaises(KeyError, lambda : store.select('df3/meta/s/meta'))

    def test_duplicate_column_name(self):
        df = DataFrame(columns=["a", "a"], data=[[0, 0]])

        with ensure_clean_path(self.path) as path:
            self.assertRaises(ValueError, df.to_hdf, path, 'df', format='fixed')

            df.to_hdf(path, 'df', format='table')
            other = read_hdf(path, 'df')

            tm.assert_frame_equal(df, other)
            self.assertTrue(df.equals(other))
            self.assertTrue(other.equals(df))

    def test_round_trip_equals(self):
        # GH 9330
        df = DataFrame({"B": [1,2], "A": ["x","y"]})

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', format='table')
            other = read_hdf(path, 'df')
            tm.assert_frame_equal(df, other)
            self.assertTrue(df.equals(other))
            self.assertTrue(other.equals(df))

    def test_preserve_timedeltaindex_type(self):
        # GH9635
        # Storing TimedeltaIndexed DataFrames in fixed stores did not preserve
        # the type of the index.
        df = DataFrame(np.random.normal(size=(10,5)))
        df.index = timedelta_range(start='0s',periods=10,freq='1s',name='example')

        with ensure_clean_store(self.path) as store:

            store['df'] = df
            assert_frame_equal(store['df'], df)

    def test_colums_multiindex_modified(self):
        # BUG: 7212
        # read_hdf store.select modified the passed columns parameters
        # when multi-indexed.

        df = DataFrame(np.random.rand(4, 5),
                       index=list('abcd'),
                       columns=list('ABCDE'))
        df.index.name = 'letters'
        df = df.set_index(keys='E', append=True)

        data_columns = df.index.names+df.columns.tolist()
        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df',
                      mode='a',
                      append=True,
                      data_columns=data_columns,
                      index=False)
            cols2load = list('BCD')
            cols2load_original = list(cols2load)
            df_loaded = read_hdf(path, 'df', columns=cols2load)
            self.assertTrue(cols2load_original == cols2load)

    def test_to_hdf_with_object_column_names(self):
        # GH9057
        # Writing HDF5 table format should only work for string-like
        # column types

        types_should_fail = [ tm.makeIntIndex, tm.makeFloatIndex,
                                tm.makeDateIndex, tm.makeTimedeltaIndex,
                                tm.makePeriodIndex ]
        types_should_run = [ tm.makeStringIndex, tm.makeCategoricalIndex ]

        if compat.PY3:
            types_should_run.append(tm.makeUnicodeIndex)
        else:
            types_should_fail.append(tm.makeUnicodeIndex)

        for index in types_should_fail:
            df = DataFrame(np.random.randn(10, 2), columns=index(2))
            with ensure_clean_path(self.path) as path:
                with self.assertRaises(ValueError,
                        msg="cannot have non-object label DataIndexableCol"):
                    df.to_hdf(path, 'df', format='table', data_columns=True)

        for index in types_should_run:
            df = DataFrame(np.random.randn(10, 2), columns=index(2))
            with ensure_clean_path(self.path) as path:
                df.to_hdf(path, 'df', format='table', data_columns=True)
                result = pd.read_hdf(path, 'df', where="index = [{0}]".format(df.index[0]))
                assert(len(result))


    def test_read_hdf_open_store(self):
        # GH10330
        # No check for non-string path_or-buf, and no test of open store
        df = DataFrame(np.random.rand(4, 5),
                       index=list('abcd'),
                       columns=list('ABCDE'))
        df.index.name = 'letters'
        df = df.set_index(keys='E', append=True)

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', mode='w')
            direct = read_hdf(path, 'df')
            store = HDFStore(path, mode='r')
            indirect = read_hdf(store, 'df')
            tm.assert_frame_equal(direct, indirect)
            self.assertTrue(store.is_open)
            store.close()

    def test_read_hdf_iterator(self):
        df = DataFrame(np.random.rand(4, 5),
                       index=list('abcd'),
                       columns=list('ABCDE'))
        df.index.name = 'letters'
        df = df.set_index(keys='E', append=True)

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', mode='w', format='t')
            direct = read_hdf(path, 'df')
            iterator = read_hdf(path, 'df', iterator=True)
            self.assertTrue(isinstance(iterator, TableIterator))
            indirect = next(iterator.__iter__())
            tm.assert_frame_equal(direct, indirect)
            iterator.store.close()

    def test_read_hdf_errors(self):
        df = DataFrame(np.random.rand(4, 5),
                       index=list('abcd'),
                       columns=list('ABCDE'))

        with ensure_clean_path(self.path) as path:
            self.assertRaises(IOError, read_hdf, path, 'key')
            df.to_hdf(path, 'df')
            store = HDFStore(path, mode='r')
            store.close()
            self.assertRaises(IOError, read_hdf, store, 'df')
            with open(path, mode='r') as store:
                self.assertRaises(NotImplementedError, read_hdf, store, 'df')

    def test_invalid_complib(self):
        df = DataFrame(np.random.rand(4, 5),
                       index=list('abcd'),
                       columns=list('ABCDE'))
        with ensure_clean_path(self.path) as path:
            self.assertRaises(ValueError, df.to_hdf, path, 'df', complib='blosc:zlib')
    # GH10443
    def test_read_nokey(self):
        df = DataFrame(np.random.rand(4, 5),
                       index=list('abcd'),
                       columns=list('ABCDE'))
        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', mode='a')
            reread = read_hdf(path)
            assert_frame_equal(df, reread)
            df.to_hdf(path, 'df2', mode='a')
            self.assertRaises(ValueError, read_hdf, path)


class TestHDFComplexValues(Base):
    # GH10447
    def test_complex_fixed(self):
        df = DataFrame(np.random.rand(4, 5).astype(np.complex64),
                       index=list('abcd'),
                       columns=list('ABCDE'))

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

        df = DataFrame(np.random.rand(4, 5).astype(np.complex128),
                       index=list('abcd'),
                       columns=list('ABCDE'))
        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    def test_complex_table(self):
        df = DataFrame(np.random.rand(4, 5).astype(np.complex64),
                       index=list('abcd'),
                       columns=list('ABCDE'))

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', format='table')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

        df = DataFrame(np.random.rand(4, 5).astype(np.complex128),
                       index=list('abcd'),
                       columns=list('ABCDE'))

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', format='table', mode='w')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    def test_complex_mixed_fixed(self):
        complex64 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex64)
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j],
                              dtype=np.complex128)
        df = DataFrame({'A': [1, 2, 3, 4],
                        'B': ['a', 'b', 'c', 'd'],
                        'C': complex64,
                        'D': complex128,
                        'E': [1.0, 2.0, 3.0, 4.0]},
                       index=list('abcd'))
        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    def test_complex_mixed_table(self):
        complex64 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j], dtype=np.complex64)
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j],
                              dtype=np.complex128)
        df = DataFrame({'A': [1, 2, 3, 4],
                        'B': ['a', 'b', 'c', 'd'],
                        'C': complex64,
                        'D': complex128,
                        'E': [1.0, 2.0, 3.0, 4.0]},
                       index=list('abcd'))

        with ensure_clean_store(self.path) as store:
            store.append('df', df, data_columns=['A', 'B'])
            result = store.select('df', where=Term('A>2'))
            assert_frame_equal(df.loc[df.A > 2], result)

        with ensure_clean_path(self.path) as path:
            df.to_hdf(path, 'df', format='table')
            reread = read_hdf(path, 'df')
            assert_frame_equal(df, reread)

    def test_complex_across_dimensions_fixed(self):
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
        s = Series(complex128, index=list('abcd'))
        df = DataFrame({'A': s, 'B': s})
        p = Panel({'One': df, 'Two': df})

        objs = [s, df, p]
        comps = [tm.assert_series_equal, tm.assert_frame_equal,
                 tm.assert_panel_equal]
        for obj, comp in zip(objs, comps):
            with ensure_clean_path(self.path) as path:
                obj.to_hdf(path, 'obj', format='fixed')
                reread = read_hdf(path, 'obj')
                comp(obj, reread)

    def test_complex_across_dimensions(self):
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
        s = Series(complex128, index=list('abcd'))
        df = DataFrame({'A': s, 'B': s})
        p = Panel({'One': df, 'Two': df})
        p4d = pd.Panel4D({'i': p, 'ii': p})

        objs = [df, p, p4d]
        comps = [tm.assert_frame_equal, tm.assert_panel_equal,
                 tm.assert_panel4d_equal]
        for obj, comp in zip(objs, comps):
            with ensure_clean_path(self.path) as path:
                obj.to_hdf(path, 'obj', format='table')
                reread = read_hdf(path, 'obj')
                comp(obj, reread)

    def test_complex_indexing_error(self):
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j],
                              dtype=np.complex128)
        df = DataFrame({'A': [1, 2, 3, 4],
                        'B': ['a', 'b', 'c', 'd'],
                        'C': complex128},
                       index=list('abcd'))
        with ensure_clean_store(self.path) as store:
            self.assertRaises(TypeError, store.append, 'df', df, data_columns=['C'])

    def test_complex_series_error(self):
        complex128 = np.array([1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j, 1.0 + 1.0j])
        s = Series(complex128, index=list('abcd'))

        with ensure_clean_path(self.path) as path:
            self.assertRaises(TypeError, s.to_hdf, path, 'obj', format='t')

        with ensure_clean_path(self.path) as path:
            s.to_hdf(path, 'obj', format='t', index=False)
            reread = read_hdf(path, 'obj')
            tm.assert_series_equal(s, reread)

    def test_complex_append(self):
        df = DataFrame({'a': np.random.randn(100).astype(np.complex128),
                        'b': np.random.randn(100)})

        with ensure_clean_store(self.path) as store:
            store.append('df', df, data_columns=['b'])
            store.append('df', df)
            result = store.select('df')
            assert_frame_equal(pd.concat([df, df], 0), result)

class TestTimezones(Base, tm.TestCase):


    def _compare_with_tz(self, a, b):
        tm.assert_frame_equal(a, b)

        # compare the zones on each element
        for c in a.columns:
            for i in a.index:
                a_e = a.loc[i,c]
                b_e = b.loc[i,c]
                if not (a_e == b_e and a_e.tz == b_e.tz):
                    raise AssertionError("invalid tz comparsion [%s] [%s]" % (a_e, b_e))

    def test_append_with_timezones_dateutil(self):

        from datetime import timedelta
        tm._skip_if_no_dateutil()

        # use maybe_get_tz instead of dateutil.tz.gettz to handle the windows filename issues.
        from pandas.tslib import maybe_get_tz
        gettz = lambda x: maybe_get_tz('dateutil/' + x)

        # as columns
        with ensure_clean_store(self.path) as store:

            _maybe_remove(store, 'df_tz')
            df = DataFrame(dict(A=[ Timestamp('20130102 2:00:00', tz=gettz('US/Eastern')) + timedelta(hours=1) * i for i in range(5) ]))

            store.append('df_tz', df, data_columns=['A'])
            result = store['df_tz']
            self._compare_with_tz(result, df)
            assert_frame_equal(result, df)

            # select with tz aware
            expected = df[df.A >= df.A[3]]
            result = store.select('df_tz', where=Term('A>=df.A[3]'))
            self._compare_with_tz(result, expected)

            # ensure we include dates in DST and STD time here.
            _maybe_remove(store, 'df_tz')
            df = DataFrame(dict(A=Timestamp('20130102', tz=gettz('US/Eastern')), B=Timestamp('20130603', tz=gettz('US/Eastern'))), index=range(5))
            store.append('df_tz', df)
            result = store['df_tz']
            self._compare_with_tz(result, df)
            assert_frame_equal(result, df)

            df = DataFrame(dict(A=Timestamp('20130102', tz=gettz('US/Eastern')), B=Timestamp('20130102', tz=gettz('EET'))), index=range(5))
            self.assertRaises(ValueError, store.append, 'df_tz', df)

            # this is ok
            _maybe_remove(store, 'df_tz')
            store.append('df_tz', df, data_columns=['A', 'B'])
            result = store['df_tz']
            self._compare_with_tz(result, df)
            assert_frame_equal(result, df)

            # can't append with diff timezone
            df = DataFrame(dict(A=Timestamp('20130102', tz=gettz('US/Eastern')), B=Timestamp('20130102', tz=gettz('CET'))), index=range(5))
            self.assertRaises(ValueError, store.append, 'df_tz', df)

        # as index
        with ensure_clean_store(self.path) as store:

            # GH 4098 example
            df = DataFrame(dict(A=Series(lrange(3), index=date_range('2000-1-1', periods=3, freq='H', tz=gettz('US/Eastern')))))

            _maybe_remove(store, 'df')
            store.put('df', df)
            result = store.select('df')
            assert_frame_equal(result, df)

            _maybe_remove(store, 'df')
            store.append('df', df)
            result = store.select('df')
            assert_frame_equal(result, df)

    def test_append_with_timezones_pytz(self):

        from datetime import timedelta

        # as columns
        with ensure_clean_store(self.path) as store:

            _maybe_remove(store, 'df_tz')
            df = DataFrame(dict(A = [ Timestamp('20130102 2:00:00',tz='US/Eastern') + timedelta(hours=1)*i for i in range(5) ]))
            store.append('df_tz',df,data_columns=['A'])
            result = store['df_tz']
            self._compare_with_tz(result,df)
            assert_frame_equal(result,df)

            # select with tz aware
            self._compare_with_tz(store.select('df_tz',where=Term('A>=df.A[3]')),df[df.A>=df.A[3]])

            _maybe_remove(store, 'df_tz')
            # ensure we include dates in DST and STD time here.
            df = DataFrame(dict(A = Timestamp('20130102',tz='US/Eastern'), B = Timestamp('20130603',tz='US/Eastern')),index=range(5))
            store.append('df_tz',df)
            result = store['df_tz']
            self._compare_with_tz(result,df)
            assert_frame_equal(result,df)

            df = DataFrame(dict(A = Timestamp('20130102',tz='US/Eastern'), B = Timestamp('20130102',tz='EET')),index=range(5))
            self.assertRaises(ValueError, store.append, 'df_tz', df)

            # this is ok
            _maybe_remove(store, 'df_tz')
            store.append('df_tz',df,data_columns=['A','B'])
            result = store['df_tz']
            self._compare_with_tz(result,df)
            assert_frame_equal(result,df)

            # can't append with diff timezone
            df = DataFrame(dict(A = Timestamp('20130102',tz='US/Eastern'), B = Timestamp('20130102',tz='CET')),index=range(5))
            self.assertRaises(ValueError, store.append, 'df_tz', df)

        # as index
        with ensure_clean_store(self.path) as store:

            # GH 4098 example
            df = DataFrame(dict(A = Series(lrange(3), index=date_range('2000-1-1',periods=3,freq='H', tz='US/Eastern'))))

            _maybe_remove(store, 'df')
            store.put('df',df)
            result = store.select('df')
            assert_frame_equal(result,df)

            _maybe_remove(store, 'df')
            store.append('df',df)
            result = store.select('df')
            assert_frame_equal(result,df)

    def test_tseries_select_index_column(self):
        # GH7777
        # selecting a UTC datetimeindex column did
        # not preserve UTC tzinfo set before storing

        # check that no tz still works
        rng = date_range('1/1/2000', '1/30/2000')
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        with ensure_clean_store(self.path) as store:
            store.append('frame', frame)
            result = store.select_column('frame', 'index')
            self.assertEqual(rng.tz, DatetimeIndex(result.values).tz)

        # check utc
        rng = date_range('1/1/2000', '1/30/2000', tz='UTC')
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        with ensure_clean_store(self.path) as store:
            store.append('frame', frame)
            result = store.select_column('frame', 'index')
            self.assertEqual(rng.tz, result.dt.tz)

        # double check non-utc
        rng = date_range('1/1/2000', '1/30/2000', tz='US/Eastern')
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        with ensure_clean_store(self.path) as store:
            store.append('frame', frame)
            result = store.select_column('frame', 'index')
            self.assertEqual(rng.tz, result.dt.tz)

    def test_timezones(self):
        rng = date_range('1/1/2000', '1/30/2000', tz='US/Eastern')
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        with ensure_clean_store(self.path) as store:
            store['frame'] = frame
            recons = store['frame']
            self.assertTrue(recons.index.equals(rng))
            self.assertEqual(rng.tz, recons.index.tz)

    def test_fixed_offset_tz(self):
        rng = date_range('1/1/2000 00:00:00-07:00', '1/30/2000 00:00:00-07:00')
        frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

        with ensure_clean_store(self.path) as store:
            store['frame'] = frame
            recons = store['frame']
            self.assertTrue(recons.index.equals(rng))
            self.assertEqual(rng.tz, recons.index.tz)

    def test_store_timezone(self):
        # GH2852
        # issue storing datetime.date with a timezone as it resets when read back in a new timezone

        import platform
        if platform.system() == "Windows":
            raise nose.SkipTest("timezone setting not supported on windows")

        import datetime
        import time
        import os

        # original method
        with ensure_clean_store(self.path) as store:

            today = datetime.date(2013,9,10)
            df = DataFrame([1,2,3], index = [today, today, today])
            store['obj1'] = df
            result = store['obj1']
            assert_frame_equal(result, df)

        # with tz setting
        orig_tz = os.environ.get('TZ')

        def setTZ(tz):
            if tz is None:
                try:
                    del os.environ['TZ']
                except:
                    pass
            else:
                os.environ['TZ']=tz
                time.tzset()

        try:

            with ensure_clean_store(self.path) as store:

                setTZ('EST5EDT')
                today = datetime.date(2013,9,10)
                df = DataFrame([1,2,3], index = [today, today, today])
                store['obj1'] = df

                setTZ('CST6CDT')
                result = store['obj1']

                assert_frame_equal(result, df)

        finally:
            setTZ(orig_tz)

    def test_legacy_datetimetz_object(self):
        # legacy from < 0.17.0
        # 8260
        expected = DataFrame(dict(A=Timestamp('20130102', tz='US/Eastern'), B=Timestamp('20130603', tz='CET')), index=range(5))
        with ensure_clean_store(tm.get_data_path('legacy_hdf/datetimetz_object.h5'), mode='r') as store:
            result = store['df']
            assert_frame_equal(result, expected)

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
