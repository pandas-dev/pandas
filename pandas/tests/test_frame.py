# pylint: disable-msg=W0612,E1101
from copy import deepcopy
from datetime import datetime, timedelta, time
from StringIO import StringIO
import cPickle as pickle
import operator
import re
import unittest
import nose

from numpy import random, nan
from numpy.random import randn
import numpy as np
import numpy.ma as ma
from numpy.testing import assert_array_equal

import pandas as pan
import pandas.core.nanops as nanops
import pandas.core.common as com
import pandas.core.format as fmt
import pandas.core.datetools as datetools
from pandas.core.api import (DataFrame, Index, Series, notnull, isnull,
                             MultiIndex, DatetimeIndex, Timestamp, Period)
from pandas import date_range
import pandas as pd
from pandas.io.parsers import read_csv
from pandas.parser import CParserError

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 assertRaisesRegexp,
                                 makeCustomDataframe as mkdf,
                                 ensure_clean)
from pandas.util import py3compat
from pandas.util.compat import OrderedDict

import pandas.util.testing as tm
import pandas.lib as lib

from numpy.testing.decorators import slow

def _skip_if_no_scipy():
    try:
        import scipy.stats
    except ImportError:
        raise nose.SkipTest

#---------------------------------------------------------------------
# DataFrame test cases

JOIN_TYPES = ['inner', 'outer', 'left', 'right']
MIXED_FLOAT_DTYPES = ['float16','float32','float64']
MIXED_INT_DTYPES   = ['uint8','uint16','uint32','uint64','int8','int16',
                      'int32','int64']

def _check_mixed_float(df, dtype = None):

    # float16 are most likely to be upcasted to float32
    dtypes = dict(A = 'float32', B = 'float32', C = 'float16', D = 'float64')
    if isinstance(dtype, basestring):
        dtypes = dict([ (k,dtype) for k, v in dtypes.items() ])
    elif isinstance(dtype, dict):
        dtypes.update(dtype)
    if dtypes.get('A'):
        assert(df.dtypes['A'] == dtypes['A'])
    if dtypes.get('B'):
        assert(df.dtypes['B'] == dtypes['B'])
    if dtypes.get('C'):
        assert(df.dtypes['C'] == dtypes['C'])
    if dtypes.get('D'):
        assert(df.dtypes['D'] == dtypes['D'])

def _check_mixed_int(df, dtype = None):
    dtypes = dict(A = 'int32', B = 'uint64', C = 'uint8', D = 'int64')
    if isinstance(dtype, basestring):
        dtypes = dict([ (k,dtype) for k, v in dtypes.items() ])
    elif isinstance(dtype, dict):
        dtypes.update(dtype)
    if dtypes.get('A'):
        assert(df.dtypes['A'] == dtypes['A'])
    if dtypes.get('B'):
        assert(df.dtypes['B'] == dtypes['B'])
    if dtypes.get('C'):
        assert(df.dtypes['C'] == dtypes['C'])
    if dtypes.get('D'):
        assert(df.dtypes['D'] == dtypes['D'])




class CheckIndexing(object):

    _multiprocess_can_split_ = True

    def test_getitem(self):
        # slicing

        sl = self.frame[:20]
        self.assertEqual(20, len(sl.index))

        # column access

        for _, series in sl.iteritems():
            self.assertEqual(20, len(series.index))
            self.assert_(tm.equalContents(series.index, sl.index))

        for key, _ in self.frame._series.iteritems():
            self.assert_(self.frame[key] is not None)

        self.assert_('random' not in self.frame)
        self.assertRaises(Exception, self.frame.__getitem__, 'random')

    def test_getitem_dupe_cols(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=['a', 'a', 'b'])
        try:
            df[['baf']]
        except KeyError:
            pass
        else:
            self.fail("Dataframe failed to raise KeyError")

    def test_get(self):
        b = self.frame.get('B')
        assert_series_equal(b, self.frame['B'])

        self.assert_(self.frame.get('foo') is None)
        assert_series_equal(self.frame.get('foo', self.frame['B']),
                            self.frame['B'])

    def test_getitem_iterator(self):
        idx = iter(['A', 'B', 'C'])
        result = self.frame.ix[:, idx]
        expected = self.frame.ix[:, ['A', 'B', 'C']]
        assert_frame_equal(result, expected)

    def test_getitem_list(self):
        self.frame.columns.name = 'foo'

        result = self.frame[['B', 'A']]
        result2 = self.frame[Index(['B', 'A'])]

        expected = self.frame.ix[:, ['B', 'A']]
        expected.columns.name = 'foo'

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        self.assertEqual(result.columns.name, 'foo')

        self.assertRaises(Exception, self.frame.__getitem__,
                          ['B', 'A', 'foo'])
        self.assertRaises(Exception, self.frame.__getitem__,
                          Index(['B', 'A', 'foo']))

        # tuples
        df = DataFrame(randn(8, 3),
                       columns=Index([('foo', 'bar'), ('baz', 'qux'),
                                      ('peek', 'aboo')], name='sth'))

        result = df[[('foo', 'bar'), ('baz', 'qux')]]
        expected = df.ix[:, :2]
        assert_frame_equal(result, expected)
        self.assertEqual(result.columns.name, 'sth')

    def test_setitem_list(self):

        self.frame['E'] = 'foo'
        data = self.frame[['A', 'B']]
        self.frame[['B', 'A']] = data

        assert_series_equal(self.frame['B'], data['A'])
        assert_series_equal(self.frame['A'], data['B'])

        df = DataFrame(0, range(3), ['tt1', 'tt2'], dtype=np.int_)
        df.ix[1, ['tt1', 'tt2']] = [1, 2]

        result = df.ix[1, ['tt1', 'tt2']]
        expected = Series([1, 2], df.columns, dtype=np.int_)
        assert_series_equal(result, expected)

        df['tt1'] = df['tt2'] = '0'
        df.ix[1, ['tt1', 'tt2']] = ['1', '2']
        result = df.ix[1, ['tt1', 'tt2']]
        expected = Series(['1', '2'], df.columns)
        assert_series_equal(result, expected)

    def test_setitem_list_not_dataframe(self):
        data = np.random.randn(len(self.frame), 2)
        self.frame[['A', 'B']] = data
        assert_almost_equal(self.frame[['A', 'B']].values, data)

    def test_setitem_list_of_tuples(self):
        tuples = zip(self.frame['A'], self.frame['B'])
        self.frame['tuples'] = tuples

        result = self.frame['tuples']
        expected = Series(tuples, index=self.frame.index)
        assert_series_equal(result, expected)

    def test_getitem_boolean(self):
        # boolean indexing
        d = self.tsframe.index[10]
        indexer = self.tsframe.index > d
        indexer_obj = indexer.astype(object)

        subindex = self.tsframe.index[indexer]
        subframe = self.tsframe[indexer]

        self.assert_(np.array_equal(subindex, subframe.index))
        self.assertRaises(Exception, self.tsframe.__getitem__, indexer[:-1])

        subframe_obj = self.tsframe[indexer_obj]
        assert_frame_equal(subframe_obj, subframe)

        self.assertRaises(ValueError, self.tsframe.__getitem__, self.tsframe)

        # test that Series work
        indexer_obj = Series(indexer_obj, self.tsframe.index)

        subframe_obj = self.tsframe[indexer_obj]
        assert_frame_equal(subframe_obj, subframe)

        # test that Series indexers reindex
        import warnings
        warnings.filterwarnings(action='ignore', category=UserWarning)

        indexer_obj = indexer_obj.reindex(self.tsframe.index[::-1])

        subframe_obj = self.tsframe[indexer_obj]
        assert_frame_equal(subframe_obj, subframe)

        warnings.filterwarnings(action='default', category=UserWarning)

        # test df[df > 0]
        for df in [ self.tsframe, self.mixed_frame, self.mixed_float, self.mixed_int ]:

            data = df._get_numeric_data()
            bif = df[df > 0]
            bifw = DataFrame(dict([ (c,np.where(data[c] > 0, data[c], np.nan)) for c in data.columns ]),
                             index=data.index, columns=data.columns)

            # add back other columns to compare
            for c in df.columns:
                if c not in bifw:
                    bifw[c] = df[c]
            bifw = bifw.reindex(columns = df.columns)

            assert_frame_equal(bif, bifw, check_dtype=False)
            for c in df.columns:
                if bif[c].dtype != bifw[c].dtype:
                    self.assert_(bif[c].dtype == df[c].dtype)

    def test_getitem_boolean_casting(self):

        # don't upcast if we don't need to
        df = self.tsframe.copy()
        df['E'] = 1
        df['E'] = df['E'].astype('int32')
        df['E1'] = df['E'].copy()
        df['F'] = 1
        df['F'] = df['F'].astype('int64')
        df['F1'] = df['F'].copy()

        casted = df[df>0]
        result = casted.get_dtype_counts()
        expected = Series({'float64': 4, 'int32' : 2, 'int64' : 2})
        assert_series_equal(result, expected)

        # int block splitting
        df.ix[1:3,['E1','F1']] = 0
        casted = df[df>0]
        result = casted.get_dtype_counts()
        expected = Series({'float64': 6, 'int32' : 1, 'int64' : 1})
        assert_series_equal(result, expected)

        # where dtype conversions
        # GH 3733
        df = DataFrame(data = np.random.randn(100, 50))
        df = df.where(df > 0) # create nans
        bools = df > 0
        mask = isnull(df)
        expected = bools.astype(float).mask(mask)
        result = bools.mask(mask)
        assert_frame_equal(result,expected)

    def test_getitem_boolean_list(self):
        df = DataFrame(np.arange(12).reshape(3, 4))

        def _checkit(lst):
            result = df[lst]
            expected = df.ix[df.index[lst]]
            assert_frame_equal(result, expected)

        _checkit([True, False, True])
        _checkit([True, True, True])
        _checkit([False, False, False])

    def test_getitem_boolean_iadd(self):
        arr = randn(5, 5)

        df = DataFrame(arr.copy(), columns = ['A','B','C','D','E'])

        df[df < 0] += 1
        arr[arr < 0] += 1

        assert_almost_equal(df.values, arr)

    def test_boolean_index_empty_corner(self):
        # #2096
        blah = DataFrame(np.empty([0, 1]), columns=['A'],
                         index=DatetimeIndex([]))

        # both of these should succeed trivially
        k = np.array([], bool)

        blah[k]
        blah[k] = 0

    def test_getitem_ix_mixed_integer(self):
        df = DataFrame(np.random.randn(4, 3),
                       index=[1, 10, 'C', 'E'], columns=[1, 2, 3])

        result = df.ix[:-1]
        expected = df.ix[df.index[:-1]]
        assert_frame_equal(result, expected)

        result = df.ix[[1, 10]]
        expected = df.ix[Index([1, 10], dtype=object)]
        assert_frame_equal(result, expected)

    def test_getitem_setitem_ix_negative_integers(self):
        result = self.frame.ix[:, -1]
        assert_series_equal(result, self.frame['D'])

        result = self.frame.ix[:, [-1]]
        assert_frame_equal(result, self.frame[['D']])

        result = self.frame.ix[:, [-1, -2]]
        assert_frame_equal(result, self.frame[['D', 'C']])

        self.frame.ix[:, [-1]] = 0
        self.assert_((self.frame['D'] == 0).all())

        df = DataFrame(np.random.randn(8, 4))
        self.assert_(isnull(df.ix[:, [-1]].values).all())

        # #1942
        a = DataFrame(randn(20, 2), index=[chr(x + 65) for x in range(20)])
        a.ix[-1] = a.ix[-2]

        assert_series_equal(a.ix[-1], a.ix[-2])

    def test_getattr(self):
        tm.assert_series_equal(self.frame.A, self.frame['A'])
        self.assertRaises(AttributeError, getattr, self.frame,
                          'NONEXISTENT_NAME')

    def test_setattr_column(self):
        df = DataFrame({'foobar': 1}, index=range(10))

        df.foobar = 5
        self.assert_((df.foobar == 5).all())

    def test_setitem(self):
        # not sure what else to do here
        series = self.frame['A'][::2]
        self.frame['col5'] = series
        self.assert_('col5' in self.frame)
        tm.assert_dict_equal(series, self.frame['col5'],
                             compare_keys=False)

        series = self.frame['A']
        self.frame['col6'] = series
        tm.assert_dict_equal(series, self.frame['col6'],
                             compare_keys=False)

        self.assertRaises(Exception, self.frame.__setitem__,
                          randn(len(self.frame) + 1))

        # set ndarray
        arr = randn(len(self.frame))
        self.frame['col9'] = arr
        self.assert_((self.frame['col9'] == arr).all())

        self.frame['col7'] = 5
        assert((self.frame['col7'] == 5).all())

        self.frame['col0'] = 3.14
        assert((self.frame['col0'] == 3.14).all())

        self.frame['col8'] = 'foo'
        assert((self.frame['col8'] == 'foo').all())

        smaller = self.frame[:2]
        smaller['col10'] = ['1', '2']
        self.assertEqual(smaller['col10'].dtype, np.object_)
        self.assert_((smaller['col10'] == ['1', '2']).all())

        # with a dtype
        for dtype in ['int32','int64','float32','float64']:
            self.frame[dtype] = np.array(arr,dtype=dtype)
            self.assert_(self.frame[dtype].dtype.name == dtype)

    def test_setitem_tuple(self):
        self.frame['A', 'B'] = self.frame['A']
        assert_series_equal(self.frame['A', 'B'], self.frame['A'])

    def test_setitem_always_copy(self):
        s = self.frame['A'].copy()
        self.frame['E'] = s

        self.frame['E'][5:10] = nan
        self.assert_(notnull(s[5:10]).all())

    def test_setitem_boolean(self):
        df = self.frame.copy()
        values = self.frame.values

        df[df['A'] > 0] = 4
        values[values[:, 0] > 0] = 4
        assert_almost_equal(df.values, values)

        # test that column reindexing works
        series = df['A'] == 4
        series = series.reindex(df.index[::-1])
        df[series] = 1
        values[values[:, 0] == 4] = 1
        assert_almost_equal(df.values, values)

        df[df > 0] = 5
        values[values > 0] = 5
        assert_almost_equal(df.values, values)

        df[df == 5] = 0
        values[values == 5] = 0
        assert_almost_equal(df.values, values)

        # a df that needs alignment first
        df[df[:-1] < 0] = 2
        np.putmask(values[:-1], values[:-1] < 0, 2)
        assert_almost_equal(df.values, values)

        # indexed with same shape but rows-reversed df
        df[df[::-1] == 2] = 3
        values[values == 2] = 3
        assert_almost_equal(df.values, values)

        self.assertRaises(Exception, df.__setitem__, df * 0, 2)

        # index with DataFrame
        mask = df > np.abs(df)
        expected = df.copy()
        df[df > np.abs(df)] = nan
        expected.values[mask.values] = nan
        assert_frame_equal(df, expected)

        # set from DataFrame
        expected = df.copy()
        df[df > np.abs(df)] = df * 2
        np.putmask(expected.values, mask.values, df.values * 2)
        assert_frame_equal(df, expected)

    def test_setitem_cast(self):
        self.frame['D'] = self.frame['D'].astype('i8')
        self.assert_(self.frame['D'].dtype == np.int64)

        # #669, should not cast?
        # this is now set to int64, which means a replacement of the column to
        # the value dtype (and nothing to do with the existing dtype)
        self.frame['B'] = 0
        self.assert_(self.frame['B'].dtype == np.int64)

        # cast if pass array of course
        self.frame['B'] = np.arange(len(self.frame))
        self.assert_(issubclass(self.frame['B'].dtype.type, np.integer))

        self.frame['foo'] = 'bar'
        self.frame['foo'] = 0
        self.assert_(self.frame['foo'].dtype == np.int64)

        self.frame['foo'] = 'bar'
        self.frame['foo'] = 2.5
        self.assert_(self.frame['foo'].dtype == np.float64)

        self.frame['something'] = 0
        self.assert_(self.frame['something'].dtype == np.int64)
        self.frame['something'] = 2
        self.assert_(self.frame['something'].dtype == np.int64)
        self.frame['something'] = 2.5
        self.assert_(self.frame['something'].dtype == np.float64)

    def test_setitem_boolean_column(self):
        expected = self.frame.copy()
        mask = self.frame['A'] > 0

        self.frame.ix[mask, 'B'] = 0
        expected.values[mask, 1] = 0

        assert_frame_equal(self.frame, expected)

    def test_setitem_corner(self):
        # corner case
        df = DataFrame({'B': [1., 2., 3.],
                        'C': ['a', 'b', 'c']},
                       index=np.arange(3))
        del df['B']
        df['B'] = [1., 2., 3.]
        self.assert_('B' in df)
        self.assertEqual(len(df.columns), 2)

        df['A'] = 'beginning'
        df['E'] = 'foo'
        df['D'] = 'bar'
        df[datetime.now()] = 'date'
        df[datetime.now()] = 5.

        # what to do when empty frame with index
        dm = DataFrame(index=self.frame.index)
        dm['A'] = 'foo'
        dm['B'] = 'bar'
        self.assertEqual(len(dm.columns), 2)
        self.assertEqual(dm.values.dtype, np.object_)

        # upcast
        dm['C'] = 1
        self.assertEqual(dm['C'].dtype, np.int64)

        dm['E'] = 1.
        self.assertEqual(dm['E'].dtype, np.float64)

        # set existing column
        dm['A'] = 'bar'
        self.assertEqual('bar', dm['A'][0])

        dm = DataFrame(index=np.arange(3))
        dm['A'] = 1
        dm['foo'] = 'bar'
        del dm['foo']
        dm['foo'] = 'bar'
        self.assertEqual(dm['foo'].dtype, np.object_)

        dm['coercable'] = ['1', '2', '3']
        self.assertEqual(dm['coercable'].dtype, np.object_)

    def test_setitem_corner2(self):
        data = {"title": ['foobar', 'bar', 'foobar'] + ['foobar'] * 17,
                "cruft": np.random.random(20)}

        df = DataFrame(data)
        ix = df[df['title'] == 'bar'].index

        df.ix[ix, ['title']] = 'foobar'
        df.ix[ix, ['cruft']] = 0

        assert(df.ix[1, 'title'] == 'foobar')
        assert(df.ix[1, 'cruft'] == 0)

    def test_setitem_ambig(self):
        # difficulties with mixed-type data
        from decimal import Decimal

        # created as float type
        dm = DataFrame(index=range(3), columns=range(3))

        coercable_series = Series([Decimal(1) for _ in range(3)],
                                  index=range(3))
        uncoercable_series = Series(['foo', 'bzr', 'baz'], index=range(3))

        dm[0] = np.ones(3)
        self.assertEqual(len(dm.columns), 3)
        # self.assert_(dm.objects is None)

        dm[1] = coercable_series
        self.assertEqual(len(dm.columns), 3)
        # self.assert_(dm.objects is None)

        dm[2] = uncoercable_series
        self.assertEqual(len(dm.columns), 3)
        # self.assert_(dm.objects is not None)
        self.assert_(dm[2].dtype == np.object_)

    def test_setitem_clear_caches(self):
        # GH #304
        df = DataFrame({'x': [1.1, 2.1, 3.1, 4.1], 'y': [5.1, 6.1, 7.1, 8.1]},
                       index=[0, 1, 2, 3])
        df.insert(2, 'z', np.nan)

        # cache it
        foo = df['z']

        df.ix[2:, 'z'] = 42

        expected = Series([np.nan, np.nan, 42, 42], index=df.index)
        self.assert_(df['z'] is not foo)
        assert_series_equal(df['z'], expected)

    def test_setitem_None(self):
        # GH #766
        self.frame[None] = self.frame['A']
        assert_series_equal(self.frame[None], self.frame['A'])
        repr(self.frame)

    def test_delitem_corner(self):
        f = self.frame.copy()
        del f['D']
        self.assertEqual(len(f.columns), 3)
        self.assertRaises(KeyError, f.__delitem__, 'D')
        del f['B']
        self.assertEqual(len(f.columns), 2)

    def test_getitem_fancy_2d(self):
        f = self.frame
        ix = f.ix

        assert_frame_equal(ix[:, ['B', 'A']], f.reindex(columns=['B', 'A']))

        subidx = self.frame.index[[5, 4, 1]]
        assert_frame_equal(ix[subidx, ['B', 'A']],
                           f.reindex(index=subidx, columns=['B', 'A']))

        # slicing rows, etc.
        assert_frame_equal(ix[5:10], f[5:10])
        assert_frame_equal(ix[5:10, :], f[5:10])
        assert_frame_equal(ix[:5, ['A', 'B']],
                           f.reindex(index=f.index[:5], columns=['A', 'B']))

        # slice rows with labels, inclusive!
        expected = ix[5:11]
        result = ix[f.index[5]:f.index[10]]
        assert_frame_equal(expected, result)

        # slice columns
        assert_frame_equal(ix[:, :2], f.reindex(columns=['A', 'B']))

        # get view
        exp = f.copy()
        ix[5:10].values[:] = 5
        exp.values[5:10] = 5
        assert_frame_equal(f, exp)

        self.assertRaises(ValueError, ix.__getitem__, f > 0.5)

    def test_slice_floats(self):
        index = [52195.504153, 52196.303147, 52198.369883]
        df = DataFrame(np.random.rand(3, 2), index=index)

        s1 = df.ix[52195.1:52196.5]
        self.assertEquals(len(s1), 2)

        s1 = df.ix[52195.1:52196.6]
        self.assertEquals(len(s1), 2)

        s1 = df.ix[52195.1:52198.9]
        self.assertEquals(len(s1), 3)

    def test_getitem_fancy_slice_integers_step(self):
        df = DataFrame(np.random.randn(10, 5))

        # this is OK
        result = df.ix[:8:2]
        df.ix[:8:2] = np.nan
        self.assert_(isnull(df.ix[:8:2]).values.all())

    def test_getitem_setitem_integer_slice_keyerrors(self):
        df = DataFrame(np.random.randn(10, 5), index=range(0, 20, 2))

        # this is OK
        cp = df.copy()
        cp.ix[4:10] = 0
        self.assert_((cp.ix[4:10] == 0).values.all())

        # so is this
        cp = df.copy()
        cp.ix[3:11] = 0
        self.assert_((cp.ix[3:11] == 0).values.all())

        result = df.ix[4:10]
        result2 = df.ix[3:11]
        expected = df.reindex([4, 6, 8, 10])

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        # non-monotonic, raise KeyError
        df2 = df[::-1]
        self.assertRaises(KeyError, df2.ix.__getitem__, slice(3, 11))
        self.assertRaises(KeyError, df2.ix.__setitem__, slice(3, 11), 0)

    def test_setitem_fancy_2d(self):
        f = self.frame
        ix = f.ix

        # case 1
        frame = self.frame.copy()
        expected = frame.copy()
        frame.ix[:, ['B', 'A']] = 1
        expected['B'] = 1.
        expected['A'] = 1.
        assert_frame_equal(frame, expected)

        # case 2
        frame = self.frame.copy()
        frame2 = self.frame.copy()

        expected = frame.copy()

        subidx = self.frame.index[[5, 4, 1]]
        values = randn(3, 2)

        frame.ix[subidx, ['B', 'A']] = values
        frame2.ix[[5, 4, 1], ['B', 'A']] = values

        expected['B'].ix[subidx] = values[:, 0]
        expected['A'].ix[subidx] = values[:, 1]

        assert_frame_equal(frame, expected)
        assert_frame_equal(frame2, expected)

        # case 3: slicing rows, etc.
        frame = self.frame.copy()

        expected1 = self.frame.copy()
        frame.ix[5:10] = 1.
        expected1.values[5:10] = 1.
        assert_frame_equal(frame, expected1)

        expected2 = self.frame.copy()
        arr = randn(5, len(frame.columns))
        frame.ix[5:10] = arr
        expected2.values[5:10] = arr
        assert_frame_equal(frame, expected2)

        # case 4
        frame = self.frame.copy()
        frame.ix[5:10, :] = 1.
        assert_frame_equal(frame, expected1)
        frame.ix[5:10, :] = arr
        assert_frame_equal(frame, expected2)

        # case 5
        frame = self.frame.copy()
        frame2 = self.frame.copy()

        expected = self.frame.copy()
        values = randn(5, 2)

        frame.ix[:5, ['A', 'B']] = values
        expected['A'][:5] = values[:, 0]
        expected['B'][:5] = values[:, 1]
        assert_frame_equal(frame, expected)

        frame2.ix[:5, [0, 1]] = values
        assert_frame_equal(frame2, expected)

        # case 6: slice rows with labels, inclusive!
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[frame.index[5]:frame.index[10]] = 5.
        expected.values[5:11] = 5
        assert_frame_equal(frame, expected)

        # case 7: slice columns
        frame = self.frame.copy()
        frame2 = self.frame.copy()
        expected = self.frame.copy()

        # slice indices
        frame.ix[:, 1:3] = 4.
        expected.values[:, 1:3] = 4.
        assert_frame_equal(frame, expected)

        # slice with labels
        frame.ix[:, 'B':'C'] = 4.
        assert_frame_equal(frame, expected)

        # new corner case of boolean slicing / setting
        frame = DataFrame(zip([2, 3, 9, 6, 7], [np.nan] * 5),
                          columns=['a', 'b'])
        lst = [100]
        lst.extend([np.nan] * 4)
        expected = DataFrame(zip([100, 3, 9, 6, 7], lst), columns=['a', 'b'])
        frame[frame['a'] == 2] = 100
        assert_frame_equal(frame, expected)

    def test_fancy_getitem_slice_mixed(self):
        sliced = self.mixed_frame.ix[:, -3:]
        self.assert_(sliced['D'].dtype == np.float64)

        # get view with single block
        sliced = self.frame.ix[:, -3:]
        sliced['C'] = 4.
        self.assert_((self.frame['C'] == 4).all())

    def test_fancy_setitem_int_labels(self):
        # integer index defers to label-based indexing

        df = DataFrame(np.random.randn(10, 5), index=np.arange(0, 20, 2))

        tmp = df.copy()
        exp = df.copy()
        tmp.ix[[0, 2, 4]] = 5
        exp.values[:3] = 5
        assert_frame_equal(tmp, exp)

        tmp = df.copy()
        exp = df.copy()
        tmp.ix[6] = 5
        exp.values[3] = 5
        assert_frame_equal(tmp, exp)

        tmp = df.copy()
        exp = df.copy()
        tmp.ix[:, 2] = 5
        exp.values[:, 2] = 5
        assert_frame_equal(tmp, exp)

    def test_fancy_getitem_int_labels(self):
        df = DataFrame(np.random.randn(10, 5), index=np.arange(0, 20, 2))

        result = df.ix[[4, 2, 0], [2, 0]]
        expected = df.reindex(index=[4, 2, 0], columns=[2, 0])
        assert_frame_equal(result, expected)

        result = df.ix[[4, 2, 0]]
        expected = df.reindex(index=[4, 2, 0])
        assert_frame_equal(result, expected)

        result = df.ix[4]
        expected = df.xs(4)
        assert_series_equal(result, expected)

        result = df.ix[:, 3]
        expected = df[3]
        assert_series_equal(result, expected)

    def test_fancy_index_int_labels_exceptions(self):
        df = DataFrame(np.random.randn(10, 5), index=np.arange(0, 20, 2))

        # labels that aren't contained
        self.assertRaises(KeyError, df.ix.__setitem__,
                          ([0, 1, 2], [2, 3, 4]), 5)

        # try to set indices not contained in frame
        self.assertRaises(KeyError,
                          self.frame.ix.__setitem__,
                          ['foo', 'bar', 'baz'], 1)
        self.assertRaises(KeyError,
                          self.frame.ix.__setitem__,
                          (slice(None, None), ['E']), 1)
        self.assertRaises(KeyError,
                          self.frame.ix.__setitem__,
                          (slice(None, None), 'E'), 1)

    def test_setitem_fancy_mixed_2d(self):
        self.mixed_frame.ix[:5, ['C', 'B', 'A']] = 5
        result = self.mixed_frame.ix[:5, ['C', 'B', 'A']]
        self.assert_((result.values == 5).all())

        self.mixed_frame.ix[5] = np.nan
        self.assert_(isnull(self.mixed_frame.ix[5]).all())

        self.mixed_frame.ix[5] = self.mixed_frame.ix[6]
        assert_series_equal(self.mixed_frame.ix[5], self.mixed_frame.ix[6])

        # #1432
        df = DataFrame({1: [1., 2., 3.],
                        2: [3, 4, 5]})
        self.assert_(df._is_mixed_type)

        df.ix[1] = [5, 10]

        expected = DataFrame({1: [1., 5., 3.],
                              2: [3, 10, 5]})

        assert_frame_equal(df, expected)

    def test_ix_align(self):
        b = Series(randn(10))
        b.sort()
        df_orig = DataFrame(randn(10, 4))
        df = df_orig.copy()

        df.ix[:, 0] = b
        assert_series_equal(df.ix[:, 0].reindex(b.index), b)

        dft = df_orig.T
        dft.ix[0, :] = b
        assert_series_equal(dft.ix[0, :].reindex(b.index), b)

        df = df_orig.copy()
        df.ix[:5, 0] = b
        s = df.ix[:5, 0]
        assert_series_equal(s, b.reindex(s.index))

        dft = df_orig.T
        dft.ix[0, :5] = b
        s = dft.ix[0, :5]
        assert_series_equal(s, b.reindex(s.index))

        df = df_orig.copy()
        idx = [0, 1, 3, 5]
        df.ix[idx, 0] = b
        s = df.ix[idx, 0]
        assert_series_equal(s, b.reindex(s.index))

        dft = df_orig.T
        dft.ix[0, idx] = b
        s = dft.ix[0, idx]
        assert_series_equal(s, b.reindex(s.index))

    def test_ix_frame_align(self):
        b = DataFrame(np.random.randn(3, 4))
        df_orig = DataFrame(randn(10, 4))
        df = df_orig.copy()

        df.ix[:3] = b
        out = b.ix[:3]
        assert_frame_equal(out, b)

        b.sort_index(inplace=True)

        df = df_orig.copy()
        df.ix[[0, 1, 2]] = b
        out = df.ix[[0, 1, 2]].reindex(b.index)
        assert_frame_equal(out, b)

        df = df_orig.copy()
        df.ix[:3] = b
        out = df.ix[:3]
        assert_frame_equal(out, b.reindex(out.index))

    def test_getitem_setitem_non_ix_labels(self):
        df = tm.makeTimeDataFrame()

        start, end = df.index[[5, 10]]

        result = df.ix[start:end]
        result2 = df[start:end]
        expected = df[5:11]
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        result = df.copy()
        result.ix[start:end] = 0
        result2 = df.copy()
        result2[start:end] = 0
        expected = df.copy()
        expected[5:11] = 0
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

    def test_ix_multi_take(self):
        df = DataFrame(np.random.randn(3, 2))
        rs = df.ix[df.index == 0, :]
        xp = df.reindex([0])
        assert_frame_equal(rs, xp)

        """ #1321
        df = DataFrame(np.random.randn(3, 2))
        rs = df.ix[df.index==0, df.columns==1]
        xp = df.reindex([0], [1])
        assert_frame_equal(rs, xp)
        """

    def test_ix_multi_take_nonint_index(self):
        df = DataFrame(np.random.randn(3, 2), index=['x', 'y', 'z'],
                       columns=['a', 'b'])
        rs = df.ix[[0], [0]]
        xp = df.reindex(['x'], columns=['a'])
        assert_frame_equal(rs, xp)

    def test_ix_multi_take_multiindex(self):
        df = DataFrame(np.random.randn(3, 2), index=['x', 'y', 'z'],
                       columns=[['a', 'b'], ['1', '2']])
        rs = df.ix[[0], [0]]
        xp = df.reindex(['x'], columns=[('a', '1')])
        assert_frame_equal(rs, xp)

    def test_ix_dup(self):
        idx = Index(['a', 'a', 'b', 'c', 'd', 'd'])
        df = DataFrame(np.random.randn(len(idx), 3), idx)

        sub = df.ix[:'d']
        assert_frame_equal(sub, df)

        sub = df.ix['a':'c']
        assert_frame_equal(sub, df.ix[0:4])

        sub = df.ix['b':'d']
        assert_frame_equal(sub, df.ix[2:])

    def test_getitem_fancy_1d(self):
        f = self.frame
        ix = f.ix

        # return self if no slicing...for now
        self.assert_(ix[:, :] is f)

        # low dimensional slice
        xs1 = ix[2, ['C', 'B', 'A']]
        xs2 = f.xs(f.index[2]).reindex(['C', 'B', 'A'])
        assert_series_equal(xs1, xs2)

        ts1 = ix[5:10, 2]
        ts2 = f[f.columns[2]][5:10]
        assert_series_equal(ts1, ts2)

        # positional xs
        xs1 = ix[0]
        xs2 = f.xs(f.index[0])
        assert_series_equal(xs1, xs2)

        xs1 = ix[f.index[5]]
        xs2 = f.xs(f.index[5])
        assert_series_equal(xs1, xs2)

        # single column
        assert_series_equal(ix[:, 'A'], f['A'])

        # return view
        exp = f.copy()
        exp.values[5] = 4
        ix[5][:] = 4
        assert_frame_equal(exp, f)

        exp.values[:, 1] = 6
        ix[:, 1][:] = 6
        assert_frame_equal(exp, f)

        # slice of mixed-frame
        xs = self.mixed_frame.ix[5]
        exp = self.mixed_frame.xs(self.mixed_frame.index[5])
        assert_series_equal(xs, exp)

    def test_setitem_fancy_1d(self):
        # case 1: set cross-section for indices
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[2, ['C', 'B', 'A']] = [1., 2., 3.]
        expected['C'][2] = 1.
        expected['B'][2] = 2.
        expected['A'][2] = 3.
        assert_frame_equal(frame, expected)

        frame2 = self.frame.copy()
        frame2.ix[2, [3, 2, 1]] = [1., 2., 3.]
        assert_frame_equal(frame, expected)

        # case 2, set a section of a column
        frame = self.frame.copy()
        expected = self.frame.copy()

        vals = randn(5)
        expected.values[5:10, 2] = vals
        frame.ix[5:10, 2] = vals
        assert_frame_equal(frame, expected)

        frame2 = self.frame.copy()
        frame2.ix[5:10, 'B'] = vals
        assert_frame_equal(frame, expected)

        # case 3: full xs
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[4] = 5.
        expected.values[4] = 5.
        assert_frame_equal(frame, expected)

        frame.ix[frame.index[4]] = 6.
        expected.values[4] = 6.
        assert_frame_equal(frame, expected)

        # single column
        frame = self.frame.copy()
        expected = self.frame.copy()

        frame.ix[:, 'A'] = 7.
        expected['A'] = 7.
        assert_frame_equal(frame, expected)

    def test_getitem_fancy_scalar(self):
        f = self.frame
        ix = f.ix
        # individual value
        for col in f.columns:
            ts = f[col]
            for idx in f.index[::5]:
                assert_almost_equal(ix[idx, col], ts[idx])

    def test_setitem_fancy_scalar(self):
        f = self.frame
        expected = self.frame.copy()
        ix = f.ix
        # individual value
        for j, col in enumerate(f.columns):
            ts = f[col]
            for idx in f.index[::5]:
                i = f.index.get_loc(idx)
                val = randn()
                expected.values[i, j] = val
                ix[idx, col] = val
                assert_frame_equal(f, expected)

    def test_getitem_fancy_boolean(self):
        f = self.frame
        ix = f.ix

        expected = f.reindex(columns=['B', 'D'])
        result = ix[:, [False, True, False, True]]
        assert_frame_equal(result, expected)

        expected = f.reindex(index=f.index[5:10], columns=['B', 'D'])
        result = ix[5:10, [False, True, False, True]]
        assert_frame_equal(result, expected)

        boolvec = f.index > f.index[7]
        expected = f.reindex(index=f.index[boolvec])
        result = ix[boolvec]
        assert_frame_equal(result, expected)
        result = ix[boolvec, :]
        assert_frame_equal(result, expected)

        result = ix[boolvec, 2:]
        expected = f.reindex(index=f.index[boolvec],
                             columns=['C', 'D'])
        assert_frame_equal(result, expected)

    def test_setitem_fancy_boolean(self):
        # from 2d, set with booleans
        frame = self.frame.copy()
        expected = self.frame.copy()

        mask = frame['A'] > 0
        frame.ix[mask] = 0.
        expected.values[mask] = 0.
        assert_frame_equal(frame, expected)

        frame = self.frame.copy()
        expected = self.frame.copy()
        frame.ix[mask, ['A', 'B']] = 0.
        expected.values[mask, :2] = 0.
        assert_frame_equal(frame, expected)

    def test_getitem_fancy_ints(self):
        result = self.frame.ix[[1, 4, 7]]
        expected = self.frame.ix[self.frame.index[[1, 4, 7]]]
        assert_frame_equal(result, expected)

        result = self.frame.ix[:, [2, 0, 1]]
        expected = self.frame.ix[:, self.frame.columns[[2, 0, 1]]]
        assert_frame_equal(result, expected)

    def test_getitem_setitem_fancy_exceptions(self):
        ix = self.frame.ix
        self.assertRaises(Exception, ix.__getitem__,
                          (slice(None, None, None),
                           slice(None, None, None),
                           slice(None, None, None)))
        self.assertRaises(Exception, ix.__setitem__,
                          (slice(None, None, None),
                           slice(None, None, None),
                           slice(None, None, None)), 1)

    def test_getitem_setitem_boolean_misaligned(self):
        # boolean index misaligned labels
        mask = self.frame['A'][::-1] > 1

        result = self.frame.ix[mask]
        expected = self.frame.ix[mask[::-1]]
        assert_frame_equal(result, expected)

        cp = self.frame.copy()
        expected = self.frame.copy()
        cp.ix[mask] = 0
        expected.ix[mask] = 0
        assert_frame_equal(cp, expected)

    def test_getitem_setitem_boolean_multi(self):
        df = DataFrame(np.random.randn(3, 2))

        # get
        k1 = np.array([True, False, True])
        k2 = np.array([False, True])
        result = df.ix[k1, k2]
        expected = df.ix[[0, 2], [1]]
        assert_frame_equal(result, expected)

        expected = df.copy()
        df.ix[np.array([True, False, True]),
              np.array([False, True])] = 5
        expected.ix[[0, 2], [1]] = 5
        assert_frame_equal(df, expected)

    def test_getitem_setitem_float_labels(self):
        index = Index([1.5, 2, 3, 4, 5])
        df = DataFrame(np.random.randn(5, 5), index=index)

        result = df.ix[1.5:4]
        expected = df.reindex([1.5, 2, 3, 4])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 4)

        result = df.ix[4:5]
        expected = df.reindex([4, 5])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 2)

        # this should raise an exception
        self.assertRaises(Exception, df.ix.__getitem__, slice(1, 2))
        self.assertRaises(Exception, df.ix.__setitem__, slice(1, 2), 0)

        # #2727
        index = Index([1.0, 2.5, 3.5, 4.5, 5.0])
        df = DataFrame(np.random.randn(5, 5), index=index)

        # positional slicing!
        result = df.ix[1.0:5]
        expected = df.reindex([2.5, 3.5, 4.5, 5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 4)

        # positional again
        result = df.ix[4:5]
        expected = df.reindex([5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 1)

        # label-based
        result = df.ix[1.0:5.0]
        expected = df.reindex([1.0, 2.5, 3.5, 4.5, 5.0])
        assert_frame_equal(result, expected)
        self.assertEqual(len(result), 5)

        cp = df.copy()
        # positional slicing!
        cp.ix[1.0:5] = 0
        self.assert_((cp.ix[1.0:5] == 0).values.all())
        self.assert_((cp.ix[0:1] == df.ix[0:1]).values.all())

        cp = df.copy()
        # positional again
        cp.ix[4:5] = 0
        self.assert_((cp.ix[4:5] == 0).values.all())
        self.assert_((cp.ix[0:4] == df.ix[0:4]).values.all())

        cp = df.copy()
        # label-based
        cp.ix[1.0:5.0] = 0
        self.assert_((cp.ix[1.0:5.0] == 0).values.all())

    def test_setitem_single_column_mixed(self):
        df = DataFrame(randn(5, 3), index=['a', 'b', 'c', 'd', 'e'],
                       columns=['foo', 'bar', 'baz'])
        df['str'] = 'qux'
        df.ix[::2, 'str'] = nan
        expected = [nan, 'qux', nan, 'qux', nan]
        assert_almost_equal(df['str'].values, expected)

    def test_setitem_single_column_mixed_datetime(self):
        df = DataFrame(randn(5, 3), index=['a', 'b', 'c', 'd', 'e'],
                       columns=['foo', 'bar', 'baz'])

        df['timestamp'] = Timestamp('20010102')

        # check our dtypes
        result = df.get_dtype_counts()
        expected = Series({'float64': 3, 'datetime64[ns]': 1})
        assert_series_equal(result, expected)

        # set an allowable datetime64 type
        from pandas import tslib
        df.ix['b', 'timestamp'] = tslib.iNaT
        self.assert_(com.isnull(df.ix['b', 'timestamp']))

        # allow this syntax
        df.ix['c', 'timestamp'] = nan
        self.assert_(com.isnull(df.ix['c', 'timestamp']))

        # allow this syntax
        df.ix['d', :] = nan
        self.assert_(com.isnull(df.ix['c', :]).all() == False)

        # as of GH 3216 this will now work!
        # try to set with a list like item
        #self.assertRaises(
        #    Exception, df.ix.__setitem__, ('d', 'timestamp'), [nan])

    def test_setitem_frame(self):
        piece = self.frame.ix[:2, ['A', 'B']]
        self.frame.ix[-2:, ['A', 'B']] = piece.values
        assert_almost_equal(self.frame.ix[-2:, ['A', 'B']].values,
                            piece.values)

        # GH 3216

        # already aligned
        f = self.mixed_frame.copy()
        piece = DataFrame([[ 1, 2], [3, 4]], index=f.index[0:2],columns=['A', 'B'])
        key = (slice(None,2), ['A', 'B'])
        f.ix[key] = piece
        assert_almost_equal(f.ix[0:2, ['A', 'B']].values,
                            piece.values)

        # rows unaligned
        f = self.mixed_frame.copy()
        piece = DataFrame([[ 1, 2 ], [3, 4], [5, 6], [7, 8]], index=list(f.index[0:2]) + ['foo','bar'],columns=['A', 'B'])
        key = (slice(None,2), ['A', 'B'])
        f.ix[key] = piece
        assert_almost_equal(f.ix[0:2:, ['A', 'B']].values,
                            piece.values[0:2])

        # key is unaligned with values
        f = self.mixed_frame.copy()
        piece = f.ix[:2, ['A']]
        key = (slice(-2, None), ['A', 'B'])
        f.ix[key] = piece
        piece['B'] = np.nan
        assert_almost_equal(f.ix[-2:, ['A', 'B']].values,
                            piece.values)

        # ndarray
        f = self.mixed_frame.copy()
        piece = self.mixed_frame.ix[:2, ['A', 'B']]
        key = (slice(-2, None), ['A', 'B'])
        f.ix[key] = piece.values
        assert_almost_equal(f.ix[-2:, ['A', 'B']].values,
                            piece.values)


        # needs upcasting
        df = DataFrame([[1,2,'foo'],[3,4,'bar']],columns=['A','B','C'])
        df2 = df.copy()
        df2.ix[:,['A','B']] = df.ix[:,['A','B']]+0.5
        expected = df.reindex(columns=['A','B'])
        expected += 0.5
        expected['C'] = df['C']
        assert_frame_equal(df2, expected)

    def test_setitem_frame_align(self):
        piece = self.frame.ix[:2, ['A', 'B']]
        piece.index = self.frame.index[-2:]
        piece.columns = ['A', 'B']
        self.frame.ix[-2:, ['A', 'B']] = piece
        assert_almost_equal(self.frame.ix[-2:, ['A', 'B']].values,
                            piece.values)

    def test_setitem_fancy_exceptions(self):
        pass

    def test_getitem_boolean_missing(self):
        pass

    def test_setitem_boolean_missing(self):
        pass

    def test_getitem_setitem_ix_duplicates(self):
        # #1201
        df = DataFrame(np.random.randn(5, 3),
                       index=['foo', 'foo', 'bar', 'baz', 'bar'])

        result = df.ix['foo']
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.ix['bar']
        expected = df.ix[[2, 4]]
        assert_frame_equal(result, expected)

        result = df.ix['baz']
        expected = df.ix[3]
        assert_series_equal(result, expected)

    def test_getitem_ix_boolean_duplicates_multiple(self):
        # #1201
        df = DataFrame(np.random.randn(5, 3),
                       index=['foo', 'foo', 'bar', 'baz', 'bar'])

        result = df.ix[['bar']]
        exp = df.ix[[2, 4]]
        assert_frame_equal(result, exp)

        result = df.ix[df[1] > 0]
        exp = df[df[1] > 0]
        assert_frame_equal(result, exp)

        result = df.ix[df[0] > 0]
        exp = df[df[0] > 0]
        assert_frame_equal(result, exp)

    def test_getitem_setitem_ix_bool_keyerror(self):
        # #2199
        df = DataFrame({'a': [1, 2, 3]})

        self.assertRaises(KeyError, df.ix.__getitem__, False)
        self.assertRaises(KeyError, df.ix.__getitem__, True)

        self.assertRaises(KeyError, df.ix.__setitem__, False, 0)
        self.assertRaises(KeyError, df.ix.__setitem__, True, 0)

    def test_getitem_list_duplicates(self):
        # #1943
        df = DataFrame(np.random.randn(4, 4), columns=list('AABC'))
        df.columns.name = 'foo'

        result = df[['B', 'C']]
        self.assert_(result.columns.name == 'foo')

        expected = df.ix[:, 2:]
        assert_frame_equal(result, expected)

    def test_get_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                result = self.frame.get_value(idx, col)
                expected = self.frame[col][idx]
                assert_almost_equal(result, expected)

    def test_iteritems(self):
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=['a', 'a', 'b'])
        for k, v in df.iteritems():
            self.assertEqual(type(v), Series)

    def test_lookup(self):
        def alt(df, rows, cols):
            result = []
            for r, c in zip(rows, cols):
                result.append(df.get_value(r, c))
            return result

        def testit(df):
            rows = list(df.index) * len(df.columns)
            cols = list(df.columns) * len(df.index)
            result = df.lookup(rows, cols)
            expected = alt(df, rows, cols)
            assert_almost_equal(result, expected)

        testit(self.mixed_frame)
        testit(self.frame)

        df = DataFrame({'label': ['a', 'b', 'a', 'c'],
                        'mask_a': [True, True, False, True],
                        'mask_b': [True, False, False, False],
                        'mask_c': [False, True, False, True]})
        df['mask'] = df.lookup(df.index, 'mask_' + df['label'])
        exp_mask = alt(df, df.index, 'mask_' + df['label'])
        assert_almost_equal(df['mask'], exp_mask)
        self.assert_(df['mask'].dtype == np.bool_)

        self.assertRaises(ValueError, self.frame.lookup,
                          ['xyz'], ['A'])

        self.assertRaises(ValueError, self.frame.lookup,
                          [self.frame.index[0]], ['xyz'])

    def test_set_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                self.frame.set_value(idx, col, 1)
                assert_almost_equal(self.frame[col][idx], 1)

    def test_set_value_resize(self):
        res = self.frame.set_value('foobar', 'B', 0)
        self.assert_(res is not self.frame)
        self.assert_(res.index[-1] == 'foobar')
        self.assertEqual(res.get_value('foobar', 'B'), 0)

        res2 = res.set_value('foobar', 'qux', 0)
        self.assert_(res2 is not res)
        self.assert_(np.array_equal(res2.columns,
                                    list(self.frame.columns) + ['qux']))
        self.assertEqual(res2.get_value('foobar', 'qux'), 0)

        res3 = res.set_value('foobar', 'baz', 'sam')
        self.assert_(res3['baz'].dtype == np.object_)

        res3 = res.set_value('foobar', 'baz', True)
        self.assert_(res3['baz'].dtype == np.object_)

        res3 = res.set_value('foobar', 'baz', 5)
        self.assert_(com.is_float_dtype(res3['baz']))
        self.assert_(isnull(res3['baz'].drop(['foobar'])).values.all())
        self.assertRaises(ValueError, res3.set_value, 'foobar', 'baz', 'sam')

    def test_set_value_with_index_dtype_change(self):
        df = DataFrame(randn(3, 3), index=range(3), columns=list('ABC'))
        res = df.set_value('C', 2, 1.0)
        self.assert_(list(res.index) == list(df.index) + ['C'])
        self.assert_(list(res.columns) == list(df.columns) + [2])

    def test_get_set_value_no_partial_indexing(self):
        # partial w/ MultiIndex raise exception
        index = MultiIndex.from_tuples([(0, 1), (0, 2), (1, 1), (1, 2)])
        df = DataFrame(index=index, columns=range(4))
        self.assertRaises(KeyError, df.get_value, 0, 1)
        # self.assertRaises(KeyError, df.set_value, 0, 1, 0)

    def test_single_element_ix_dont_upcast(self):
        self.frame['E'] = 1
        self.assert_(issubclass(self.frame['E'].dtype.type,
                                (int, np.integer)))

        result = self.frame.ix[self.frame.index[5], 'E']
        self.assert_(com.is_integer(result))

    def test_irow(self):
        df = DataFrame(np.random.randn(10, 4), index=range(0, 20, 2))

        result = df.irow(1)
        exp = df.ix[2]
        assert_series_equal(result, exp)

        result = df.irow(2)
        exp = df.ix[4]
        assert_series_equal(result, exp)

        # slice
        result = df.irow(slice(4, 8))
        expected = df.ix[8:14]
        assert_frame_equal(result, expected)

        # verify slice is view
        result[2] = 0.
        exp_col = df[2].copy()
        exp_col[4:8] = 0.
        assert_series_equal(df[2], exp_col)

        # list of integers
        result = df.irow([1, 2, 4, 6])
        expected = df.reindex(df.index[[1, 2, 4, 6]])
        assert_frame_equal(result, expected)

    def test_icol(self):
        df = DataFrame(np.random.randn(4, 10), columns=range(0, 20, 2))

        result = df.icol(1)
        exp = df.ix[:, 2]
        assert_series_equal(result, exp)

        result = df.icol(2)
        exp = df.ix[:, 4]
        assert_series_equal(result, exp)

        # slice
        result = df.icol(slice(4, 8))
        expected = df.ix[:, 8:14]
        assert_frame_equal(result, expected)

        # verify slice is view
        result[8] = 0.
        self.assert_((df[8] == 0).all())

        # list of integers
        result = df.icol([1, 2, 4, 6])
        expected = df.reindex(columns=df.columns[[1, 2, 4, 6]])
        assert_frame_equal(result, expected)

    def test_irow_icol_duplicates(self):
        df = DataFrame(np.random.rand(3, 3), columns=list('ABC'),
                       index=list('aab'))

        result = df.irow(0)
        result2 = df.ix[0]
        self.assert_(isinstance(result, Series))
        assert_almost_equal(result.values, df.values[0])
        assert_series_equal(result, result2)

        result = df.T.icol(0)
        result2 = df.T.ix[:, 0]
        self.assert_(isinstance(result, Series))
        assert_almost_equal(result.values, df.values[0])
        assert_series_equal(result, result2)

        # multiindex
        df = DataFrame(np.random.randn(3, 3), columns=[['i', 'i', 'j'],
                                                       ['A', 'A', 'B']],
                       index=[['i', 'i', 'j'], ['X', 'X', 'Y']])
        rs = df.irow(0)
        xp = df.ix[0]
        assert_series_equal(rs, xp)

        rs = df.icol(0)
        xp = df.T.ix[0]
        assert_series_equal(rs, xp)

        rs = df.icol([0])
        xp = df.ix[:, [0]]
        assert_frame_equal(rs, xp)

        # #2259
        df = DataFrame([[1, 2, 3], [4, 5, 6]], columns=[1, 1, 2])
        result = df.icol([0])
        expected = df.take([0], axis=1)
        assert_frame_equal(result, expected)

    def test_icol_sparse_propegate_fill_value(self):
        from pandas.sparse.api import SparseDataFrame
        df = SparseDataFrame({'A': [999, 1]}, default_fill_value=999)
        self.assertTrue(len(df['A'].sp_values) == len(df.icol(0).sp_values))

    def test_iget_value(self):
        for i, row in enumerate(self.frame.index):
            for j, col in enumerate(self.frame.columns):
                result = self.frame.iget_value(i, j)
                expected = self.frame.get_value(row, col)
                assert_almost_equal(result, expected)

    def test_nested_exception(self):
        # Ignore the strange way of triggering the problem
        # (which may get fixed), it's just a way to trigger
        # the issue or reraising an outer exception without
        # a named argument
        df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8,
                       9]}).set_index(["a", "b"])
        l = list(df.index)
        l[0] = ["a", "b"]
        df.index = l

        try:
            repr(df)
        except Exception, e:
            self.assertNotEqual(type(e), UnboundLocalError)

_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()

_frame = DataFrame(_seriesd)
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])
_intframe = DataFrame(dict((k, v.astype(int))
                           for k, v in _seriesd.iteritems()))

_tsframe = DataFrame(_tsd)

_mixed_frame = _frame.copy()
_mixed_frame['foo'] = 'bar'


class SafeForSparse(object):

    _multiprocess_can_split_ = True

    def test_getitem_pop_assign_name(self):
        s = self.frame['A']
        self.assertEqual(s.name, 'A')

        s = self.frame.pop('A')
        self.assertEqual(s.name, 'A')

        s = self.frame.ix[:, 'B']
        self.assertEqual(s.name, 'B')

        s2 = s.ix[:]
        self.assertEqual(s2.name, 'B')

    def test_get_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                result = self.frame.get_value(idx, col)
                expected = self.frame[col][idx]
                assert_almost_equal(result, expected)

    def test_join_index(self):
        # left / right

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2)
        self.assert_(f.index.equals(joined.index))
        self.assertEqual(len(joined.columns), 4)

        joined = f.join(f2, how='left')
        self.assert_(joined.index.equals(f.index))
        self.assertEqual(len(joined.columns), 4)

        joined = f.join(f2, how='right')
        self.assert_(joined.index.equals(f2.index))
        self.assertEqual(len(joined.columns), 4)

        # corner case
        self.assertRaises(Exception, self.frame.join, self.frame,
                          how='left')

        # inner

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='inner')
        self.assert_(joined.index.equals(f.index.intersection(f2.index)))
        self.assertEqual(len(joined.columns), 4)

        # corner case
        self.assertRaises(Exception, self.frame.join, self.frame,
                          how='inner')

        # outer

        f = self.frame.reindex(columns=['A', 'B'])[:10]
        f2 = self.frame.reindex(columns=['C', 'D'])

        joined = f.join(f2, how='outer')
        self.assert_(tm.equalContents(self.frame.index, joined.index))
        self.assertEqual(len(joined.columns), 4)

        # corner case
        self.assertRaises(Exception, self.frame.join, self.frame,
                          how='outer')

        self.assertRaises(Exception, f.join, f2, how='foo')

    def test_join_index_more(self):
        af = self.frame.ix[:, ['A', 'B']]
        bf = self.frame.ix[::2, ['C', 'D']]

        expected = af.copy()
        expected['C'] = self.frame['C'][::2]
        expected['D'] = self.frame['D'][::2]

        result = af.join(bf)
        assert_frame_equal(result, expected)

        result = af.join(bf, how='right')
        assert_frame_equal(result, expected[::2])

        result = bf.join(af, how='right')
        assert_frame_equal(result, expected.ix[:, result.columns])

    def test_join_index_series(self):
        df = self.frame.copy()
        s = df.pop(self.frame.columns[-1])
        joined = df.join(s)

        assert_frame_equal(joined, self.frame, check_names=False) # TODO should this check_names ?

        s.name = None
        self.assertRaises(Exception, df.join, s)

    def test_join_overlap(self):
        df1 = self.frame.ix[:, ['A', 'B', 'C']]
        df2 = self.frame.ix[:, ['B', 'C', 'D']]

        joined = df1.join(df2, lsuffix='_df1', rsuffix='_df2')
        df1_suf = df1.ix[:, ['B', 'C']].add_suffix('_df1')
        df2_suf = df2.ix[:, ['B', 'C']].add_suffix('_df2')
        no_overlap = self.frame.ix[:, ['A', 'D']]
        expected = df1_suf.join(df2_suf).join(no_overlap)

        # column order not necessarily sorted
        assert_frame_equal(joined, expected.ix[:, joined.columns])

    def test_add_prefix_suffix(self):
        with_prefix = self.frame.add_prefix('foo#')
        expected = ['foo#%s' % c for c in self.frame.columns]
        self.assert_(np.array_equal(with_prefix.columns, expected))

        with_suffix = self.frame.add_suffix('#foo')
        expected = ['%s#foo' % c for c in self.frame.columns]
        self.assert_(np.array_equal(with_suffix.columns, expected))


class TestDataFrame(unittest.TestCase, CheckIndexing,
                    SafeForSparse):
    klass = DataFrame

    _multiprocess_can_split_ = True

    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()

        # force these all to int64 to avoid platform testing issues
        self.intframe = DataFrame(dict([ (c,s) for c,s in _intframe.iteritems() ]), dtype = np.int64)
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()
        self.mixed_float  = DataFrame({ 'A': _frame['A'].copy().astype('float32'),
                                        'B': _frame['B'].copy().astype('float32'),
                                        'C': _frame['C'].copy().astype('float16'),
                                        'D': _frame['D'].copy().astype('float64') })
        self.mixed_float2 = DataFrame({ 'A': _frame2['A'].copy().astype('float32'),
                                        'B': _frame2['B'].copy().astype('float32'),
                                        'C': _frame2['C'].copy().astype('float16'),
                                        'D': _frame2['D'].copy().astype('float64') })
        self.mixed_int    = DataFrame({ 'A': _intframe['A'].copy().astype('int32'),
                                        'B': np.ones(len(_intframe['B']),dtype='uint64'),
                                        'C': _intframe['C'].copy().astype('uint8'),
                                        'D': _intframe['D'].copy().astype('int64') })
        self.all_mixed    = DataFrame({'a': 1., 'b': 2, 'c': 'foo', 'float32' : np.array([1.]*10,dtype='float32'),
                                       'int32' : np.array([1]*10,dtype='int32'),
                                       }, index=np.arange(10))

        self.ts1 = tm.makeTimeSeries()
        self.ts2 = tm.makeTimeSeries()[5:]
        self.ts3 = tm.makeTimeSeries()[-5:]
        self.ts4 = tm.makeTimeSeries()[1:-1]

        self.ts_dict = {
            'col1': self.ts1,
            'col2': self.ts2,
            'col3': self.ts3,
            'col4': self.ts4,
        }
        self.empty = DataFrame({})

        arr = np.array([[1., 2., 3.],
                        [4., 5., 6.],
                        [7., 8., 9.]])

        self.simple = DataFrame(arr, columns=['one', 'two', 'three'],
                                index=['a', 'b', 'c'])

    def test_get_axis(self):
        f = self.frame
        self.assert_(f._get_axis_name(0) == 'index')
        self.assert_(f._get_axis_name(1) == 'columns')
        self.assert_(f._get_axis_name('index') == 'index')
        self.assert_(f._get_axis_name('columns') == 'columns')
        self.assertRaises(Exception, f._get_axis_name, 'foo')
        self.assertRaises(Exception, f._get_axis_name, None)

        self.assert_(f._get_axis_number(0) == 0)
        self.assert_(f._get_axis_number(1) == 1)
        self.assert_(f._get_axis_number('index') == 0)
        self.assert_(f._get_axis_number('columns') == 1)
        self.assertRaises(Exception, f._get_axis_number, 2)
        self.assertRaises(Exception, f._get_axis_number, None)

        self.assert_(self.frame._get_axis(0) is self.frame.index)
        self.assert_(self.frame._get_axis(1) is self.frame.columns)

    def test_set_index(self):
        idx = Index(np.arange(len(self.mixed_frame)))

        # cache it
        _ = self.mixed_frame['foo']
        self.mixed_frame.index = idx
        self.assert_(self.mixed_frame['foo'].index is idx)
        self.assertRaises(Exception, setattr, self.mixed_frame, 'index',
                          idx[::2])

    def test_set_index_cast(self):

        # issue casting an index then set_index
        df = DataFrame({'A' : [1.1,2.2,3.3], 'B' : [5.0,6.1,7.2]},
                       index = [2010,2011,2012])
        expected = df.ix[2010]
        new_index = df.index.astype(np.int32)
        df.index = new_index
        result = df.ix[2010]
        assert_series_equal(result,expected)

    def test_set_index2(self):
        df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                        'B': ['one', 'two', 'three', 'one', 'two'],
                        'C': ['a', 'b', 'c', 'd', 'e'],
                        'D': np.random.randn(5),
                        'E': np.random.randn(5)})

        # new object, single-column
        result = df.set_index('C')
        result_nodrop = df.set_index('C', drop=False)

        index = Index(df['C'], name='C')

        expected = df.ix[:, ['A', 'B', 'D', 'E']]
        expected.index = index

        expected_nodrop = df.copy()
        expected_nodrop.index = index

        assert_frame_equal(result, expected)
        assert_frame_equal(result_nodrop, expected_nodrop)
        self.assertEqual(result.index.name, index.name)

        # inplace, single
        df2 = df.copy()

        df2.set_index('C', inplace=True)

        assert_frame_equal(df2, expected)

        df3 = df.copy()
        df3.set_index('C', drop=False, inplace=True)

        assert_frame_equal(df3, expected_nodrop)

        # create new object, multi-column
        result = df.set_index(['A', 'B'])
        result_nodrop = df.set_index(['A', 'B'], drop=False)

        index = MultiIndex.from_arrays([df['A'], df['B']], names=['A', 'B'])

        expected = df.ix[:, ['C', 'D', 'E']]
        expected.index = index

        expected_nodrop = df.copy()
        expected_nodrop.index = index

        assert_frame_equal(result, expected)
        assert_frame_equal(result_nodrop, expected_nodrop)
        self.assertEqual(result.index.names, index.names)

        # inplace
        df2 = df.copy()
        df2.set_index(['A', 'B'], inplace=True)
        assert_frame_equal(df2, expected)

        df3 = df.copy()
        df3.set_index(['A', 'B'], drop=False, inplace=True)
        assert_frame_equal(df3, expected_nodrop)

        # corner case
        self.assertRaises(Exception, df.set_index, 'A', verify_integrity=True)

        # append
        result = df.set_index(['A', 'B'], append=True)
        xp = df.reset_index().set_index(['index', 'A', 'B'])
        xp.index.names = [None, 'A', 'B']
        assert_frame_equal(result, xp)

        # append to existing multiindex
        rdf = df.set_index(['A'], append=True)
        rdf = rdf.set_index(['B', 'C'], append=True)
        expected = df.set_index(['A', 'B', 'C'], append=True)
        assert_frame_equal(rdf, expected)

        # Series
        result = df.set_index(df.C)
        self.assertEqual(result.index.name, 'C')

    def test_set_index_nonuniq(self):
        df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                        'B': ['one', 'two', 'three', 'one', 'two'],
                        'C': ['a', 'b', 'c', 'd', 'e'],
                        'D': np.random.randn(5),
                        'E': np.random.randn(5)})
        self.assertRaises(Exception, df.set_index, 'A', verify_integrity=True,
                          inplace=True)
        self.assert_('A' in df)

    def test_set_index_bug(self):
        # GH1590
        df = DataFrame({'val': [0, 1, 2], 'key': ['a', 'b', 'c']})
        df2 = df.select(lambda indx: indx >= 1)
        rs = df2.set_index('key')
        xp = DataFrame({'val': [1, 2]},
                       Index(['b', 'c'], name='key'))
        assert_frame_equal(rs, xp)

    def test_set_index_pass_arrays(self):
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'foo', 'foo'],
                        'B': ['one', 'one', 'two', 'three',
                              'two', 'two', 'one', 'three'],
                        'C': np.random.randn(8),
                        'D': np.random.randn(8)})

        # multiple columns
        result = df.set_index(['A', df['B'].values], drop=False)
        expected = df.set_index(['A', 'B'], drop=False)
        assert_frame_equal(result, expected, check_names=False) # TODO should set_index check_names ?

    def test_set_index_cast_datetimeindex(self):
        df = DataFrame({'A': [datetime(2000, 1, 1) + timedelta(i)
                              for i in range(1000)],
                        'B': np.random.randn(1000)})

        idf = df.set_index('A')
        self.assert_(isinstance(idf.index, DatetimeIndex))

    def test_set_index_multiindexcolumns(self):
        columns = MultiIndex.from_tuples([('foo', 1), ('foo', 2), ('bar', 1)])
        df = DataFrame(np.random.randn(3, 3), columns=columns)
        rs = df.set_index(df.columns[0])
        xp = df.ix[:, 1:]
        xp.index = df.ix[:, 0].values
        xp.index.names = [df.columns[0]]
        assert_frame_equal(rs, xp)

    def test_set_index_empty_column(self):
        # #1971
        df = DataFrame([
            dict(a=1, p=0),
            dict(a=2, m=10),
            dict(a=3, m=11, p=20),
            dict(a=4, m=12, p=21)
        ], columns=('a', 'm', 'p', 'x'))

        # it works!
        result = df.set_index(['a', 'x'])
        repr(result)

    def test_set_columns(self):
        cols = Index(np.arange(len(self.mixed_frame.columns)))
        self.mixed_frame.columns = cols
        self.assertRaises(Exception, setattr, self.mixed_frame, 'columns',
                          cols[::2])

    def test_keys(self):
        getkeys = self.frame.keys
        self.assert_(getkeys() is self.frame.columns)

    def test_column_contains_typeerror(self):
        try:
            self.frame.columns in self.frame
        except TypeError:
            pass

    def test_constructor(self):
        df = DataFrame()
        self.assert_(len(df.index) == 0)

        df = DataFrame(data={})
        self.assert_(len(df.index) == 0)

    def test_constructor_mixed(self):
        index, data = tm.getMixedTypeDict()

        indexed_frame = DataFrame(data, index=index)
        unindexed_frame = DataFrame(data)

        self.assertEqual(self.mixed_frame['foo'].dtype, np.object_)

    def test_constructor_cast_failure(self):
        foo = DataFrame({'a': ['a', 'b', 'c']}, dtype=np.float64)
        self.assert_(foo['a'].dtype == object)

        # GH 3010, constructing with odd arrays
        df = DataFrame(np.ones((4,2)))

        # this is ok
        df['foo'] = np.ones((4,2)).tolist()

        # this is not ok
        self.assertRaises(ValueError, df.__setitem__, tuple(['test']), np.ones((4,2)))

        # this is ok
        df['foo2'] = np.ones((4,2)).tolist()


    def test_constructor_dtype_nocast_view(self):
        df = DataFrame([[1, 2]])
        should_be_view = DataFrame(df, dtype=df[0].dtype)
        should_be_view[0][0] = 99
        self.assertEqual(df.values[0, 0], 99)

        should_be_view = DataFrame(df.values, dtype=df[0].dtype)
        should_be_view[0][0] = 97
        self.assertEqual(df.values[0, 0], 97)

    def test_constructor_dtype_list_data(self):
        df = DataFrame([[1, '2'],
                        [None, 'a']], dtype=object)
        self.assert_(df.ix[1, 0] is None)
        self.assert_(df.ix[0, 1] == '2')

    def test_constructor_list_frames(self):

        # GH 3243
        result = DataFrame([DataFrame([])])
        self.assert_(result.shape == (1,0))

        result = DataFrame([DataFrame(dict(A = range(5)))])
        self.assert_(type(result.iloc[0,0]) == DataFrame)

    def test_constructor_mixed_dtypes(self):

        def _make_mixed_dtypes_df(typ, ad = None):

            if typ == 'int':
                dtypes = MIXED_INT_DTYPES
                arrays = [ np.array(np.random.rand(10), dtype = d) for d in dtypes ]
            elif typ == 'float':
                dtypes = MIXED_FLOAT_DTYPES
                arrays = [ np.array(np.random.randint(10, size=10), dtype = d) for d in dtypes ]

            zipper = zip(dtypes,arrays)
            for d,a in zipper:
                assert(a.dtype == d)
            if ad is None:
                ad = dict()
            ad.update(dict([ (d,a) for d,a in zipper ]))
            return DataFrame(ad)

        def _check_mixed_dtypes(df, dtypes = None):
            if dtypes is None:
                dtypes = MIXED_FLOAT_DTYPES + MIXED_INT_DTYPES
            for d in dtypes:
                if d in df:
                    assert(df.dtypes[d] == d)

        # mixed floating and integer coexinst in the same frame
        df     = _make_mixed_dtypes_df('float')
        _check_mixed_dtypes(df)

        # add lots of types
        df     = _make_mixed_dtypes_df('float', dict(A = 1, B = 'foo', C = 'bar'))
        _check_mixed_dtypes(df)

        # GH 622
        df     = _make_mixed_dtypes_df('int')
        _check_mixed_dtypes(df)

    def test_constructor_rec(self):
        rec = self.frame.to_records(index=False)

        # Assigning causes segfault in NumPy < 1.5.1
        # rec.dtype.names = list(rec.dtype.names)[::-1]

        index = self.frame.index

        df = DataFrame(rec)
        self.assert_(np.array_equal(df.columns, rec.dtype.names))

        df2 = DataFrame(rec, index=index)
        self.assert_(np.array_equal(df2.columns, rec.dtype.names))
        self.assert_(df2.index.equals(index))

        rng = np.arange(len(rec))[::-1]
        df3 = DataFrame(rec, index=rng, columns=['C', 'B'])
        expected = DataFrame(rec, index=rng).reindex(columns=['C', 'B'])
        assert_frame_equal(df3, expected)

    def test_constructor_bool(self):
        df = DataFrame({0: np.ones(10, dtype=bool),
                        1: np.zeros(10, dtype=bool)})
        self.assertEqual(df.values.dtype, np.bool_)

    def test_constructor_overflow_int64(self):
        values = np.array([2 ** 64 - i for i in range(1, 10)],
                          dtype=np.uint64)

        result = DataFrame({'a': values})
        self.assert_(result['a'].dtype == object)

        # #2355
        data_scores = [(6311132704823138710, 273), (2685045978526272070, 23),
                       (8921811264899370420, 45), (17019687244989530680L, 270),
                       (9930107427299601010L, 273)]
        dtype = [('uid', 'u8'), ('score', 'u8')]
        data = np.zeros((len(data_scores),), dtype=dtype)
        data[:] = data_scores
        df_crawls = DataFrame(data)
        self.assert_(df_crawls['uid'].dtype == object)

    def test_is_mixed_type(self):
        self.assert_(not self.frame._is_mixed_type)
        self.assert_(self.mixed_frame._is_mixed_type)

    def test_constructor_ordereddict(self):
        import random
        nitems = 100
        nums = range(nitems)
        random.shuffle(nums)
        expected = ['A%d' % i for i in nums]
        df = DataFrame(OrderedDict(zip(expected, [[0]] * nitems)))
        self.assertEqual(expected, list(df.columns))

    def test_constructor_dict(self):
        frame = DataFrame({'col1': self.ts1,
                           'col2': self.ts2})

        tm.assert_dict_equal(self.ts1, frame['col1'], compare_keys=False)
        tm.assert_dict_equal(self.ts2, frame['col2'], compare_keys=False)

        frame = DataFrame({'col1': self.ts1,
                           'col2': self.ts2},
                          columns=['col2', 'col3', 'col4'])

        self.assertEqual(len(frame), len(self.ts2))
        self.assert_('col1' not in frame)
        self.assert_(isnull(frame['col3']).all())

        # Corner cases
        self.assertEqual(len(DataFrame({})), 0)
        self.assertRaises(Exception, lambda x: DataFrame([self.ts1, self.ts2]))

        # mix dict and array, wrong size
        self.assertRaises(Exception, DataFrame,
                          {'A': {'a': 'a', 'b': 'b'},
                           'B': ['a', 'b', 'c']})

        # Length-one dict micro-optimization
        frame = DataFrame({'A': {'1': 1, '2': 2}})
        self.assert_(np.array_equal(frame.index, ['1', '2']))

        # empty dict plus index
        idx = Index([0, 1, 2])
        frame = DataFrame({}, index=idx)
        self.assert_(frame.index is idx)

        # empty with index and columns
        idx = Index([0, 1, 2])
        frame = DataFrame({}, index=idx, columns=idx)
        self.assert_(frame.index is idx)
        self.assert_(frame.columns is idx)
        self.assertEqual(len(frame._series), 3)

        # with dict of empty list and Series
        frame = DataFrame({'A': [], 'B': []}, columns=['A', 'B'])
        self.assert_(frame.index.equals(Index([])))

    def test_constructor_error_msgs(self):

        # mix dict and array, wrong size
        def testit():
            DataFrame({'A': {'a': 'a', 'b': 'b'},
                       'B': ['a', 'b', 'c']})
        assertRaisesRegexp(ValueError, "Mixing dicts with non-Series may lead to ambiguous ordering.", testit)

        # wrong size ndarray, GH 3105
        def testit():
            DataFrame(np.arange(12).reshape((4, 3)), columns=['foo', 'bar', 'baz'],
                      index=date_range('2000-01-01', periods=3))
        assertRaisesRegexp(ValueError, "Shape of passed values is \(3, 4\), indices imply \(3, 3\)", testit)

        # higher dim raise exception
        def testit():
            DataFrame(np.zeros((3, 3, 3)), columns=['A', 'B', 'C'], index=[1])
        assertRaisesRegexp(ValueError, "Must pass 2-d input", testit)

        # wrong size axis labels
        def testit():
            DataFrame(np.random.rand(2,3), columns=['A', 'B', 'C'], index=[1])
        assertRaisesRegexp(ValueError, "Shape of passed values is \(3, 2\), indices imply \(3, 1\)", testit)

        def testit():
            DataFrame(np.random.rand(2,3), columns=['A', 'B'], index=[1, 2])
        assertRaisesRegexp(ValueError, "Shape of passed values is \(3, 2\), indices imply \(2, 2\)", testit)

        def testit():
            DataFrame({'a': False, 'b': True})
        assertRaisesRegexp(ValueError, 'If using all scalar values, you must must pass an index', testit)

    def test_insert_error_msmgs(self):

        # GH 4107, more descriptive error message
        df = DataFrame(np.random.randint(0,2,(4,4)),
                       columns=['a', 'b', 'c', 'd'])

        def testit():
            df['gr'] = df.groupby(['b', 'c']).count()

        assertRaisesRegexp(TypeError, 'incompatible index of inserted column with frame index', testit)

    def test_constructor_subclass_dict(self):
        # Test for passing dict subclass to constructor
        data = {'col1': tm.TestSubDict((x, 10.0 * x) for x in xrange(10)),
                'col2': tm.TestSubDict((x, 20.0 * x) for x in xrange(10))}
        df = DataFrame(data)
        refdf = DataFrame(dict((col, dict(val.iteritems()))
                               for col, val in data.iteritems()))
        assert_frame_equal(refdf, df)

        data = tm.TestSubDict(data.iteritems())
        df = DataFrame(data)
        assert_frame_equal(refdf, df)

        # try with defaultdict
        from collections import defaultdict
        data = {}
        self.frame['B'][:10] = np.nan
        for k, v in self.frame.iterkv():
            dct = defaultdict(dict)
            dct.update(v.to_dict())
            data[k] = dct
        frame = DataFrame(data)
        assert_frame_equal(self.frame.sort_index(), frame)

    def test_constructor_dict_block(self):
        expected = [[4., 3., 2., 1.]]
        df = DataFrame({'d': [4.], 'c': [3.], 'b': [2.], 'a': [1.]},
                       columns=['d', 'c', 'b', 'a'])
        assert_almost_equal(df.values, expected)

    def test_constructor_dict_cast(self):
        # cast float tests
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }
        frame = DataFrame(test_data, dtype=float)
        self.assertEqual(len(frame), 3)
        self.assert_(frame['B'].dtype == np.float64)
        self.assert_(frame['A'].dtype == np.float64)

        frame = DataFrame(test_data)
        self.assertEqual(len(frame), 3)
        self.assert_(frame['B'].dtype == np.object_)
        self.assert_(frame['A'].dtype == np.float64)

        # can't cast to float
        test_data = {
            'A': dict(zip(range(20), tm.makeStringIndex(20))),
            'B': dict(zip(range(15), randn(15)))
        }
        frame = DataFrame(test_data, dtype=float)
        self.assertEqual(len(frame), 20)
        self.assert_(frame['A'].dtype == np.object_)
        self.assert_(frame['B'].dtype == np.float64)

    def test_constructor_dict_dont_upcast(self):
        d = {'Col1': {'Row1': 'A String', 'Row2': np.nan}}
        df = DataFrame(d)
        self.assert_(isinstance(df['Col1']['Row2'], float))

        dm = DataFrame([[1, 2], ['a', 'b']], index=[1, 2], columns=[1, 2])
        self.assert_(isinstance(dm[1][1], int))

    def test_constructor_dict_of_tuples(self):
        # GH #1491
        data = {'a': (1, 2, 3), 'b': (4, 5, 6)}

        result = DataFrame(data)
        expected = DataFrame(dict((k, list(v)) for k, v in data.iteritems()))
        assert_frame_equal(result, expected, check_dtype=False)

    def test_constructor_ndarray(self):
        mat = np.zeros((2, 3), dtype=float)

        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=np.int64)
        self.assert_(frame.values.dtype == np.int64)

        # 1-D input
        frame = DataFrame(np.zeros(3), columns=['A'], index=[1, 2, 3])
        self.assertEqual(len(frame.index), 3)
        self.assertEqual(len(frame.columns), 1)

        frame = DataFrame(['foo', 'bar'], index=[0, 1], columns=['A'])
        self.assertEqual(len(frame), 2)

        # higher dim raise exception
        self.assertRaises(Exception, DataFrame, np.zeros((3, 3, 3)),
                          columns=['A', 'B', 'C'], index=[1])

        # wrong size axis labels
        self.assertRaises(Exception, DataFrame, mat,
                          columns=['A', 'B', 'C'], index=[1])

        self.assertRaises(Exception, DataFrame, mat,
                          columns=['A', 'B'], index=[1, 2])

        # automatic labeling
        frame = DataFrame(mat)
        self.assert_(np.array_equal(frame.index, range(2)))
        self.assert_(np.array_equal(frame.columns, range(3)))

        frame = DataFrame(mat, index=[1, 2])
        self.assert_(np.array_equal(frame.columns, range(3)))

        frame = DataFrame(mat, columns=['A', 'B', 'C'])
        self.assert_(np.array_equal(frame.index, range(2)))

        # 0-length axis
        frame = DataFrame(np.empty((0, 3)))
        self.assert_(len(frame.index) == 0)

        frame = DataFrame(np.empty((3, 0)))
        self.assert_(len(frame.columns) == 0)

    def test_constructor_maskedarray(self):
        mat = ma.masked_all((2, 3), dtype=float)

        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)
        self.assertTrue(np.all(~np.asarray(frame == frame)))

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=np.int64)
        self.assert_(frame.values.dtype == np.int64)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = 1.0
        mat2[1, 2] = 2.0
        frame = DataFrame(mat2, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(1.0, frame['A'][1])
        self.assertEqual(2.0, frame['C'][2])

        # 1-D input
        frame = DataFrame(ma.masked_all((3,)), columns=['A'], index=[1, 2, 3])
        self.assertEqual(len(frame.index), 3)
        self.assertEqual(len(frame.columns), 1)
        self.assertTrue(np.all(~np.asarray(frame == frame)))

        # higher dim raise exception
        self.assertRaises(Exception, DataFrame, ma.masked_all((3, 3, 3)),
                          columns=['A', 'B', 'C'], index=[1])

        # wrong size axis labels
        self.assertRaises(Exception, DataFrame, mat,
                          columns=['A', 'B', 'C'], index=[1])

        self.assertRaises(Exception, DataFrame, mat,
                          columns=['A', 'B'], index=[1, 2])

        # automatic labeling
        frame = DataFrame(mat)
        self.assert_(np.array_equal(frame.index, range(2)))
        self.assert_(np.array_equal(frame.columns, range(3)))

        frame = DataFrame(mat, index=[1, 2])
        self.assert_(np.array_equal(frame.columns, range(3)))

        frame = DataFrame(mat, columns=['A', 'B', 'C'])
        self.assert_(np.array_equal(frame.index, range(2)))

        # 0-length axis
        frame = DataFrame(ma.masked_all((0, 3)))
        self.assert_(len(frame.index) == 0)

        frame = DataFrame(ma.masked_all((3, 0)))
        self.assert_(len(frame.columns) == 0)

    def test_constructor_maskedarray_nonfloat(self):
        # masked int promoted to float
        mat = ma.masked_all((2, 3), dtype=int)
        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)
        self.assertTrue(np.all(~np.asarray(frame == frame)))

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=np.float64)
        self.assert_(frame.values.dtype == np.float64)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = 1
        mat2[1, 2] = 2
        frame = DataFrame(mat2, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(1, frame['A'][1])
        self.assertEqual(2, frame['C'][2])

        # masked np.datetime64 stays (use lib.NaT as null)
        mat = ma.masked_all((2, 3), dtype='M8[ns]')
        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)
        self.assertTrue(isnull(frame).values.all())

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=np.int64)
        self.assert_(frame.values.dtype == np.int64)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = 1
        mat2[1, 2] = 2
        frame = DataFrame(mat2, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(1, frame['A'].view('i8')[1])
        self.assertEqual(2, frame['C'].view('i8')[2])

        # masked bool promoted to object
        mat = ma.masked_all((2, 3), dtype=bool)
        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)
        self.assertTrue(np.all(~np.asarray(frame == frame)))

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                          index=[1, 2], dtype=object)
        self.assert_(frame.values.dtype == object)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0, 0] = True
        mat2[1, 2] = False
        frame = DataFrame(mat2, columns=['A', 'B', 'C'], index=[1, 2])
        self.assertEqual(True, frame['A'][1])
        self.assertEqual(False, frame['C'][2])

    def test_constructor_corner(self):
        df = DataFrame(index=[])
        self.assertEqual(df.values.shape, (0, 0))

        # empty but with specified dtype
        df = DataFrame(index=range(10), columns=['a', 'b'], dtype=object)
        self.assert_(df.values.dtype == np.object_)

        # does not error but ends up float
        df = DataFrame(index=range(10), columns=['a', 'b'], dtype=int)
        self.assert_(df.values.dtype == np.object_)

        # #1783 empty dtype object
        df = DataFrame({}, columns=['foo', 'bar'])
        self.assert_(df.values.dtype == np.object_)

    def test_constructor_scalar_inference(self):
        data = {'int': 1, 'bool': True,
                'float': 3., 'complex': 4j, 'object': 'foo'}
        df = DataFrame(data, index=np.arange(10))

        self.assert_(df['int'].dtype == np.int64)
        self.assert_(df['bool'].dtype == np.bool_)
        self.assert_(df['float'].dtype == np.float64)
        self.assert_(df['complex'].dtype == np.complex128)
        self.assert_(df['object'].dtype == np.object_)

    def test_constructor_arrays_and_scalars(self):
        df = DataFrame({'a': randn(10), 'b': True})
        exp = DataFrame({'a': df['a'].values, 'b': [True] * 10})

        assert_frame_equal(df, exp)

        self.assertRaises(ValueError, DataFrame, {'a': False, 'b': True})

    def test_constructor_DataFrame(self):
        df = DataFrame(self.frame)
        assert_frame_equal(df, self.frame)

        df_casted = DataFrame(self.frame, dtype=np.int64)
        self.assert_(df_casted.values.dtype == np.int64)

    def test_constructor_more(self):
        # used to be in test_matrix.py
        arr = randn(10)
        dm = DataFrame(arr, columns=['A'], index=np.arange(10))
        self.assertEqual(dm.values.ndim, 2)

        arr = randn(0)
        dm = DataFrame(arr)
        self.assertEqual(dm.values.ndim, 2)
        self.assertEqual(dm.values.ndim, 2)

        # no data specified
        dm = DataFrame(columns=['A', 'B'], index=np.arange(10))
        self.assertEqual(dm.values.shape, (10, 2))

        dm = DataFrame(columns=['A', 'B'])
        self.assertEqual(dm.values.shape, (0, 2))

        dm = DataFrame(index=np.arange(10))
        self.assertEqual(dm.values.shape, (10, 0))

        # corner, silly
        self.assertRaises(Exception, DataFrame, (1, 2, 3))

        # can't cast
        mat = np.array(['foo', 'bar'], dtype=object).reshape(2, 1)
        self.assertRaises(ValueError, DataFrame, mat, index=[0, 1],
                          columns=[0], dtype=float)

        dm = DataFrame(DataFrame(self.frame._series))
        tm.assert_frame_equal(dm, self.frame)

        # int cast
        dm = DataFrame({'A': np.ones(10, dtype=int),
                        'B': np.ones(10, dtype=np.float64)},
                       index=np.arange(10))

        self.assertEqual(len(dm.columns), 2)
        self.assert_(dm.values.dtype == np.float64)

    def test_constructor_empty_list(self):
        df = DataFrame([], index=[])
        expected = DataFrame(index=[])
        assert_frame_equal(df, expected)

    def test_constructor_list_of_lists(self):
        # GH #484
        l = [[1, 'a'], [2, 'b']]
        df = DataFrame(data=l, columns=["num", "str"])
        self.assert_(com.is_integer_dtype(df['num']))
        self.assert_(df['str'].dtype == np.object_)

    def test_constructor_list_of_dicts(self):
        data = [OrderedDict([['a', 1.5], ['b', 3], ['c', 4], ['d', 6]]),
                OrderedDict([['a', 1.5], ['b', 3], ['d', 6]]),
                OrderedDict([['a', 1.5], ['d', 6]]),
                OrderedDict(),
                OrderedDict([['a', 1.5], ['b', 3], ['c', 4]]),
                OrderedDict([['b', 3], ['c', 4], ['d', 6]])]

        result = DataFrame(data)
        expected = DataFrame.from_dict(dict(zip(range(len(data)), data)),
                                       orient='index')
        assert_frame_equal(result, expected.reindex(result.index))

        result = DataFrame([{}])
        expected = DataFrame(index=[0])
        assert_frame_equal(result, expected)

    def test_constructor_list_of_series(self):
        data = [OrderedDict([['a', 1.5], ['b', 3.0], ['c', 4.0]]),
                OrderedDict([['a', 1.5], ['b', 3.0], ['c', 6.0]])]
        sdict = OrderedDict(zip(['x', 'y'], data))
        idx = Index(['a', 'b', 'c'])

        # all named
        data2 = [Series([1.5, 3, 4], idx, dtype='O', name='x'),
                 Series([1.5, 3, 6], idx, name='y')]
        result = DataFrame(data2)
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result, expected)

        # some unnamed
        data2 = [Series([1.5, 3, 4], idx, dtype='O', name='x'),
                 Series([1.5, 3, 6], idx)]
        result = DataFrame(data2)

        sdict = OrderedDict(zip(['x', 'Unnamed 0'], data))
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result.sort_index(), expected)

        # none named
        data = [OrderedDict([['a', 1.5], ['b', 3], ['c', 4], ['d', 6]]),
                OrderedDict([['a', 1.5], ['b', 3], ['d', 6]]),
                OrderedDict([['a', 1.5], ['d', 6]]),
                OrderedDict(),
                OrderedDict([['a', 1.5], ['b', 3], ['c', 4]]),
                OrderedDict([['b', 3], ['c', 4], ['d', 6]])]
        data = [Series(d) for d in data]

        result = DataFrame(data)
        sdict = OrderedDict(zip(range(len(data)), data))
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result, expected.reindex(result.index))

        result2 = DataFrame(data, index=np.arange(6))
        assert_frame_equal(result, result2)

        result = DataFrame([Series({})])
        expected = DataFrame(index=[0])
        assert_frame_equal(result, expected)

        data = [OrderedDict([['a', 1.5], ['b', 3.0], ['c', 4.0]]),
                OrderedDict([['a', 1.5], ['b', 3.0], ['c', 6.0]])]
        sdict = OrderedDict(zip(range(len(data)), data))

        idx = Index(['a', 'b', 'c'])
        data2 = [Series([1.5, 3, 4], idx, dtype='O'),
                 Series([1.5, 3, 6], idx)]
        result = DataFrame(data2)
        expected = DataFrame.from_dict(sdict, orient='index')
        assert_frame_equal(result, expected)

    def test_constructor_list_of_derived_dicts(self):
        class CustomDict(dict):
            pass
        d = {'a': 1.5, 'b': 3}

        data_custom = [CustomDict(d)]
        data = [d]

        result_custom = DataFrame(data_custom)
        result = DataFrame(data)
        assert_frame_equal(result, result_custom)

    def test_constructor_ragged(self):
        data = {'A': randn(10),
                'B': randn(8)}
        self.assertRaises(Exception, DataFrame, data)

    def test_constructor_scalar(self):
        idx = Index(range(3))
        df = DataFrame({"a": 0}, index=idx)
        expected = DataFrame({"a": [0, 0, 0]}, index=idx)
        assert_frame_equal(df, expected, check_dtype=False)

    def test_constructor_Series_copy_bug(self):
        df = DataFrame(self.frame['A'], index=self.frame.index, columns=['A'])
        df.copy()

    def test_constructor_mixed_dict_and_Series(self):
        data = {}
        data['A'] = {'foo': 1, 'bar': 2, 'baz': 3}
        data['B'] = Series([4, 3, 2, 1], index=['bar', 'qux', 'baz', 'foo'])

        result = DataFrame(data)
        self.assert_(result.index.is_monotonic)

        # ordering ambiguous, raise exception
        self.assertRaises(Exception, DataFrame,
                          {'A': ['a', 'b'], 'B': {'a': 'a', 'b': 'b'}})

        # this is OK though
        result = DataFrame({'A': ['a', 'b'],
                            'B': Series(['a', 'b'], index=['a', 'b'])})
        expected = DataFrame({'A': ['a', 'b'], 'B': ['a', 'b']},
                             index=['a', 'b'])
        assert_frame_equal(result, expected)

    def test_constructor_tuples(self):
        result = DataFrame({'A': [(1, 2), (3, 4)]})
        expected = DataFrame({'A': Series([(1, 2), (3, 4)])})
        assert_frame_equal(result, expected)

    def test_constructor_orient(self):
        data_dict = self.mixed_frame.T._series
        recons = DataFrame.from_dict(data_dict, orient='index')
        expected = self.mixed_frame.sort_index()
        assert_frame_equal(recons, expected)

        # dict of sequence
        a = {'hi': [32, 3, 3],
             'there': [3, 5, 3]}
        rs = DataFrame.from_dict(a, orient='index')
        xp = DataFrame.from_dict(a).T.reindex(a.keys())
        assert_frame_equal(rs, xp)

    def test_constructor_Series_named(self):
        a = Series([1, 2, 3], index=['a', 'b', 'c'], name='x')
        df = DataFrame(a)
        self.assert_(df.columns[0] == 'x')
        self.assert_(df.index.equals(a.index))

        # #2234
        a = Series([], name='x')
        df = DataFrame(a)
        self.assert_(df.columns[0] == 'x')

    def test_constructor_Series_differently_indexed(self):
        # name
        s1 = Series([1, 2, 3], index=['a', 'b', 'c'], name='x')

        # no name
        s2 = Series([1, 2, 3], index=['a', 'b', 'c'])

        other_index = Index(['a', 'b'])

        df1 = DataFrame(s1, index=other_index)
        exp1 = DataFrame(s1.reindex(other_index))
        self.assert_(df1.columns[0] == 'x')
        assert_frame_equal(df1, exp1)

        df2 = DataFrame(s2, index=other_index)
        exp2 = DataFrame(s2.reindex(other_index))
        self.assert_(df2.columns[0] == 0)
        self.assert_(df2.index.equals(other_index))
        assert_frame_equal(df2, exp2)

    def test_constructor_manager_resize(self):
        index = list(self.frame.index[:5])
        columns = list(self.frame.columns[:3])

        result = DataFrame(self.frame._data, index=index,
                           columns=columns)
        self.assert_(np.array_equal(result.index, index))
        self.assert_(np.array_equal(result.columns, columns))

    def test_constructor_from_items(self):
        items = [(c, self.frame[c]) for c in self.frame.columns]
        recons = DataFrame.from_items(items)
        assert_frame_equal(recons, self.frame)

        # pass some columns
        recons = DataFrame.from_items(items, columns=['C', 'B', 'A'])
        assert_frame_equal(recons, self.frame.ix[:, ['C', 'B', 'A']])

        # orient='index'

        row_items = [(idx, self.mixed_frame.xs(idx))
                     for idx in self.mixed_frame.index]

        recons = DataFrame.from_items(row_items,
                                      columns=self.mixed_frame.columns,
                                      orient='index')
        assert_frame_equal(recons, self.mixed_frame)
        self.assert_(recons['A'].dtype == np.float64)

        self.assertRaises(ValueError, DataFrame.from_items, row_items,
                          orient='index')

        # orient='index', but thar be tuples
        arr = lib.list_to_object_array(
            [('bar', 'baz')] * len(self.mixed_frame))
        self.mixed_frame['foo'] = arr
        row_items = [(idx, list(self.mixed_frame.xs(idx)))
                     for idx in self.mixed_frame.index]
        recons = DataFrame.from_items(row_items,
                                      columns=self.mixed_frame.columns,
                                      orient='index')
        assert_frame_equal(recons, self.mixed_frame)
        self.assert_(isinstance(recons['foo'][0], tuple))

        rs = DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6])],
                                  orient='index', columns=['one', 'two', 'three'])
        xp = DataFrame([[1, 2, 3], [4, 5, 6]], index=['A', 'B'],
                       columns=['one', 'two', 'three'])
        assert_frame_equal(rs, xp)

    def test_constructor_mix_series_nonseries(self):
        df = DataFrame({'A': self.frame['A'],
                        'B': list(self.frame['B'])}, columns=['A', 'B'])
        assert_frame_equal(df, self.frame.ix[:, ['A', 'B']])

        self.assertRaises(ValueError, DataFrame,
                          {'A': self.frame['A'],
                           'B': list(self.frame['B'])[:-2]})

    def test_constructor_miscast_na_int_dtype(self):
        df = DataFrame([[np.nan, 1], [1, 0]], dtype=np.int64)
        expected = DataFrame([[np.nan, 1], [1, 0]])
        assert_frame_equal(df, expected)

    def test_constructor_column_duplicates(self):
        # it works! #2079
        df = DataFrame([[8, 5]], columns=['a', 'a'])
        edf = DataFrame([[8, 5]])
        edf.columns = ['a', 'a']

        assert_frame_equal(df, edf)

        idf = DataFrame.from_items(
            [('a', [8]), ('a', [5])], columns=['a', 'a'])
        assert_frame_equal(idf, edf)

        self.assertRaises(ValueError, DataFrame.from_items,
                          [('a', [8]), ('a', [5]), ('b', [6])],
                          columns=['b', 'a', 'a'])


    def test_column_duplicates_operations(self):

        def check(result, expected=None):
            if expected is not None:
                assert_frame_equal(result,expected)
            result.dtypes
            str(result)

        # assignment
        # GH 3687
        arr = np.random.randn(3, 2)
        idx = range(2)
        df = DataFrame(arr, columns=['A', 'A'])
        df.columns = idx
        expected = DataFrame(arr,columns=idx)
        check(df,expected)

        idx = date_range('20130101',periods=4,freq='Q-NOV')
        df = DataFrame([[1,1,1,5],[1,1,2,5],[2,1,3,5]],columns=['a','a','a','a'])
        df.columns = idx
        expected = DataFrame([[1,1,1,5],[1,1,2,5],[2,1,3,5]],columns=idx)
        check(df,expected)

        # insert
        df = DataFrame([[1,1,1,5],[1,1,2,5],[2,1,3,5]],columns=['foo','bar','foo','hello'])
        df['string'] = 'bah'
        expected = DataFrame([[1,1,1,5,'bah'],[1,1,2,5,'bah'],[2,1,3,5,'bah']],columns=['foo','bar','foo','hello','string'])
        check(df,expected)

        # insert same dtype
        df['foo2'] = 3
        expected = DataFrame([[1,1,1,5,'bah',3],[1,1,2,5,'bah',3],[2,1,3,5,'bah',3]],columns=['foo','bar','foo','hello','string','foo2'])
        check(df,expected)

        # set (non-dup)
        df['foo2'] = 4
        expected = DataFrame([[1,1,1,5,'bah',4],[1,1,2,5,'bah',4],[2,1,3,5,'bah',4]],columns=['foo','bar','foo','hello','string','foo2'])
        check(df,expected)
        df['foo2'] = 3

        # delete (non dup)
        del df['bar']
        expected = DataFrame([[1,1,5,'bah',3],[1,2,5,'bah',3],[2,3,5,'bah',3]],columns=['foo','foo','hello','string','foo2'])
        check(df,expected)

        # try to delete again (its not consolidated)
        del df['hello']
        expected = DataFrame([[1,1,'bah',3],[1,2,'bah',3],[2,3,'bah',3]],columns=['foo','foo','string','foo2'])
        check(df,expected)

        # consolidate
        df = df.consolidate()
        expected = DataFrame([[1,1,'bah',3],[1,2,'bah',3],[2,3,'bah',3]],columns=['foo','foo','string','foo2'])
        check(df,expected)

        # insert
        df.insert(2,'new_col',5.)
        expected = DataFrame([[1,1,5.,'bah',3],[1,2,5.,'bah',3],[2,3,5.,'bah',3]],columns=['foo','foo','new_col','string','foo2'])
        check(df,expected)

        # insert a dup
        self.assertRaises(Exception, df.insert, 2, 'new_col', 4.)
        df.insert(2,'new_col',4.,allow_duplicates=True)
        expected = DataFrame([[1,1,4.,5.,'bah',3],[1,2,4.,5.,'bah',3],[2,3,4.,5.,'bah',3]],columns=['foo','foo','new_col','new_col','string','foo2'])
        check(df,expected)

        # delete (dup)
        del df['foo']
        expected = DataFrame([[4.,5.,'bah',3],[4.,5.,'bah',3],[4.,5.,'bah',3]],columns=['new_col','new_col','string','foo2'])
        assert_frame_equal(df,expected)

        # dup across dtypes
        df = DataFrame([[1,1,1.,5],[1,1,2.,5],[2,1,3.,5]],columns=['foo','bar','foo','hello'])
        check(df)

        df['foo2'] = 7.
        expected = DataFrame([[1,1,1.,5,7.],[1,1,2.,5,7.],[2,1,3.,5,7.]],columns=['foo','bar','foo','hello','foo2'])
        check(df,expected)

        result = df['foo']
        expected = DataFrame([[1,1.],[1,2.],[2,3.]],columns=['foo','foo'])
        check(result,expected)

        # multiple replacements
        df['foo'] = 'string'
        expected = DataFrame([['string',1,'string',5,7.],['string',1,'string',5,7.],['string',1,'string',5,7.]],columns=['foo','bar','foo','hello','foo2'])
        check(df,expected)

        del df['foo']
        expected = DataFrame([[1,5,7.],[1,5,7.],[1,5,7.]],columns=['bar','hello','foo2'])
        check(df,expected)

        # reindex
        df = DataFrame([[1,5,7.],[1,5,7.],[1,5,7.]],columns=['bar','a','a'])
        expected = DataFrame([[1],[1],[1]],columns=['bar'])
        result = df.reindex(columns=['bar'])
        check(result,expected)

        result1 = DataFrame([[1],[1],[1]],columns=['bar']).reindex(columns=['bar','foo'])
        result2 = df.reindex(columns=['bar','foo'])
        check(result2,result1)

        # drop
        df = DataFrame([[1,5,7.],[1,5,7.],[1,5,7.]],columns=['bar','a','a'])
        df = df.drop(['a'],axis=1)
        expected = DataFrame([[1],[1],[1]],columns=['bar'])
        check(df,expected)

    def test_insert_benchmark(self):
        # from the vb_suite/frame_methods/frame_insert_columns
        N = 10
        K = 5
        df = DataFrame(index=range(N))
        new_col = np.random.randn(N)
        for i in range(K):
            df[i] = new_col
        expected = DataFrame(np.repeat(new_col,K).reshape(N,K),index=range(N))
        assert_frame_equal(df,expected)

    def test_constructor_single_value(self):

        # expecting single value upcasting here
        df = DataFrame(0., index=[1, 2, 3], columns=['a', 'b', 'c'])
        assert_frame_equal(df, DataFrame(np.zeros(df.shape).astype('float64'), df.index,
                                         df.columns))

        df = DataFrame(0, index=[1, 2, 3], columns=['a', 'b', 'c'])
        assert_frame_equal(df, DataFrame(np.zeros(df.shape).astype('int64'), df.index,
                                         df.columns))


        df = DataFrame('a', index=[1, 2], columns=['a', 'c'])
        assert_frame_equal(df, DataFrame(np.array([['a', 'a'],
                                                   ['a', 'a']],
                                                  dtype=object),
                                         index=[1, 2],
                                         columns=['a', 'c']))

        self.assertRaises(com.PandasError, DataFrame, 'a', [1, 2])
        self.assertRaises(com.PandasError, DataFrame, 'a', columns=['a', 'c'])
        self.assertRaises(
            com.PandasError, DataFrame, 'a', [1, 2], ['a', 'c'], float)


    def test_constructor_with_datetimes(self):
        intname = np.dtype(np.int_).name
        floatname = np.dtype(np.float_).name
        datetime64name = np.dtype('M8[ns]').name
        objectname = np.dtype(np.object_).name

        # single item
        df = DataFrame({'A' : 1, 'B' : 'foo', 'C' : 'bar', 'D' : Timestamp("20010101"), 'E' : datetime(2001,1,2,0,0) },
                       index=np.arange(10))
        result = df.get_dtype_counts()
        expected = Series({'int64': 1, datetime64name: 2, objectname : 2})
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

        # check with ndarray construction ndim==0 (e.g. we are passing a ndim 0 ndarray with a dtype specified)
        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo', floatname : np.array(1.,dtype=floatname),
                        intname : np.array(1,dtype=intname)}, index=np.arange(10))
        result = df.get_dtype_counts()
        expected = { objectname : 1 }
        if intname == 'int64':
            expected['int64'] = 2
        else:
            expected['int64'] = 1
            expected[intname] = 1
        if floatname == 'float64':
            expected['float64'] = 2
        else:
            expected['float64'] = 1
            expected[floatname] = 1

        result.sort()
        expected = Series(expected)
        expected.sort()
        assert_series_equal(result, expected)

        # check with ndarray construction ndim>0
        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo', floatname : np.array([1.]*10,dtype=floatname),
                        intname : np.array([1]*10,dtype=intname)}, index=np.arange(10))
        result = df.get_dtype_counts()
        result.sort()
        assert_series_equal(result, expected)

        # GH 2809
        ind = date_range(start="2000-01-01", freq="D", periods=10)
        datetimes = [ts.to_pydatetime() for ts in ind]
        datetime_s = Series(datetimes)
        self.assert_(datetime_s.dtype == 'M8[ns]')
        df = DataFrame({'datetime_s':datetime_s})
        result = df.get_dtype_counts()
        expected = Series({ datetime64name : 1 })
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

        # GH 2810
        ind = date_range(start="2000-01-01", freq="D", periods=10)
        datetimes = [ts.to_pydatetime() for ts in ind]
        dates = [ts.date() for ts in ind]
        df = DataFrame({'datetimes': datetimes, 'dates':dates})
        result = df.get_dtype_counts()
        expected = Series({ datetime64name : 1, objectname : 1 })
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

    def test_constructor_for_list_with_dtypes(self):
        intname = np.dtype(np.int_).name
        floatname = np.dtype(np.float_).name
        datetime64name = np.dtype('M8[ns]').name
        objectname = np.dtype(np.object_).name

        # test list of lists/ndarrays
        df = DataFrame([np.arange(5) for x in range(5)])
        result = df.get_dtype_counts()
        expected = Series({'int64' : 5})

        df = DataFrame([np.array(np.arange(5),dtype='int32') for x in range(5)])
        result = df.get_dtype_counts()
        expected = Series({'int32' : 5})

        # overflow issue? (we always expecte int64 upcasting here)
        df = DataFrame({'a' : [2**31,2**31+1]})
        result = df.get_dtype_counts()
        expected = Series({'int64' : 1 })
        assert_series_equal(result, expected)

        # GH #2751 (construction with no index specified), make sure we cast to platform values
        df = DataFrame([1, 2])
        result = df.get_dtype_counts()
        expected = Series({'int64': 1 })
        assert_series_equal(result, expected)

        df = DataFrame([1.,2.])
        result = df.get_dtype_counts()
        expected = Series({'float64' : 1 })
        assert_series_equal(result, expected)

        df = DataFrame({'a' : [1, 2]})
        result = df.get_dtype_counts()
        expected = Series({'int64' : 1})
        assert_series_equal(result, expected)

        df = DataFrame({'a' : [1., 2.]})
        result = df.get_dtype_counts()
        expected = Series({'float64' : 1})
        assert_series_equal(result, expected)

        df = DataFrame({'a' : 1 }, index=range(3))
        result = df.get_dtype_counts()
        expected = Series({'int64': 1})
        assert_series_equal(result, expected)

        df = DataFrame({'a' : 1. }, index=range(3))
        result = df.get_dtype_counts()
        expected = Series({'float64': 1 })
        assert_series_equal(result, expected)

        # with object list
        df = DataFrame({'a':[1,2,4,7], 'b':[1.2, 2.3, 5.1, 6.3],
                        'c':list('abcd'), 'd':[datetime(2000,1,1) for i in range(4)],
                        'e' : [1.,2,4.,7]})
        result = df.get_dtype_counts()
        expected = Series({'int64': 1, 'float64' : 2, datetime64name: 1, objectname : 1})
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

    def test_not_hashable(self):
        df = pd.DataFrame([1])
        self.assertRaises(TypeError, hash, df)
        self.assertRaises(TypeError, hash, self.empty)

    def test_timedeltas(self):

        df = DataFrame(dict(A = Series(date_range('2012-1-1', periods=3, freq='D')),
                            B = Series([ timedelta(days=i) for i in range(3) ])))
        result = df.get_dtype_counts()
        expected = Series({'datetime64[ns]': 1, 'timedelta64[ns]' : 1 })
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

        df['C'] = df['A'] + df['B']
        expected = Series({'datetime64[ns]': 2, 'timedelta64[ns]' : 1 })
        result = df.get_dtype_counts()
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

        # mixed int types
        df['D'] = 1
        expected = Series({'datetime64[ns]': 2, 'timedelta64[ns]' : 1, 'int64' : 1 })
        result = df.get_dtype_counts()
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

    def test_operators_timedelta64(self):

        from datetime import datetime, timedelta
        df = DataFrame(dict(A = date_range('2012-1-1', periods=3, freq='D'),
                            B = date_range('2012-1-2', periods=3, freq='D'),
                            C = Timestamp('20120101')-timedelta(minutes=5,seconds=5)))

        diffs = DataFrame(dict(A = df['A']-df['C'],
                               B = df['A']-df['B']))


        # min
        result = diffs.min()
        self.assert_(result[0] == diffs.ix[0,'A'])
        self.assert_(result[1] == diffs.ix[0,'B'])

        result = diffs.min(axis=1)
        self.assert_((result == diffs.ix[0,'B']).all() == True)

        # max
        result = diffs.max()
        self.assert_(result[0] == diffs.ix[2,'A'])
        self.assert_(result[1] == diffs.ix[2,'B'])

        result = diffs.max(axis=1)
        self.assert_((result == diffs['A']).all() == True)

        # abs ###### THIS IS BROKEN NOW ###### (results are dtype=timedelta64[us]
        # even though fixed in series
        #result = np.abs(df['A']-df['B'])
        #result = diffs.abs()
        #expected = DataFrame(dict(A = df['A']-df['C'],
        #                          B = df['B']-df['A']))
        #assert_frame_equal(result,expected)

        # mixed frame
        mixed = diffs.copy()
        mixed['C'] = 'foo'
        mixed['D'] = 1
        mixed['E'] = 1.

        # this is ok
        result = mixed.min()

        # this is not
        result = mixed.min(axis=1)

        # GH 3106
        df = DataFrame({'time' : date_range('20130102',periods=5),
                        'time2' : date_range('20130105',periods=5) })
        df['off1'] = df['time2']-df['time']
        self.assert_(df['off1'].dtype == 'timedelta64[ns]')

        df['off2'] = df['time']-df['time2']
        df._consolidate_inplace()
        self.assertTrue(df['off1'].dtype == 'timedelta64[ns]')
        self.assertTrue(df['off2'].dtype == 'timedelta64[ns]')

    def test__slice_consolidate_invalidate_item_cache(self):
        # #3970
        df = DataFrame({ "aa":range(5), "bb":[2.2]*5})

        # Creates a second float block
        df["cc"] = 0.0

        # caches a reference to the 'bb' series
        df["bb"]

        # repr machinery triggers consolidation
        repr(df)

        # Assignment to wrong series
        df['bb'].iloc[0] = 0.17
        df._clear_item_cache()
        self.assertAlmostEqual(df['bb'][0], 0.17)

    def test_new_empty_index(self):
        df1 = DataFrame(randn(0, 3))
        df2 = DataFrame(randn(0, 3))
        df1.index.name = 'foo'
        self.assert_(df2.index.name is None)

    def test_astype(self):
        casted = self.frame.astype(int)
        expected = DataFrame(self.frame.values.astype(int),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

        casted = self.frame.astype(np.int32)
        expected = DataFrame(self.frame.values.astype(np.int32),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

        self.frame['foo'] = '5'
        casted = self.frame.astype(int)
        expected = DataFrame(self.frame.values.astype(int),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

        # mixed casting
        def _check_cast(df, v):
            self.assert_(list(set([ s.dtype.name for _, s in df.iteritems() ]))[0] == v)

        mn = self.all_mixed._get_numeric_data().copy()
        mn['little_float'] = np.array(12345.,dtype='float16')
        mn['big_float']    = np.array(123456789101112.,dtype='float64')

        casted = mn.astype('float64')
        _check_cast(casted, 'float64')

        casted = mn.astype('int64')
        _check_cast(casted, 'int64')

        casted = self.mixed_float.reindex(columns = ['A','B']).astype('float32')
        _check_cast(casted, 'float32')

        casted = mn.reindex(columns = ['little_float']).astype('float16')
        _check_cast(casted, 'float16')

        casted = self.mixed_float.reindex(columns = ['A','B']).astype('float16')
        _check_cast(casted, 'float16')

        casted = mn.astype('float32')
        _check_cast(casted, 'float32')

        casted = mn.astype('int32')
        _check_cast(casted, 'int32')

        # to object
        casted = mn.astype('O')
        _check_cast(casted, 'object')

    def test_astype_with_exclude_string(self):
        df = self.frame.copy()
        expected = self.frame.astype(int)
        df['string'] = 'foo'
        casted = df.astype(int, raise_on_error = False)

        expected['string'] = 'foo'
        assert_frame_equal(casted, expected)

        df = self.frame.copy()
        expected = self.frame.astype(np.int32)
        df['string'] = 'foo'
        casted = df.astype(np.int32, raise_on_error = False)

        expected['string'] = 'foo'
        assert_frame_equal(casted, expected)

    def test_astype_with_view(self):

        tf = self.mixed_float.reindex(columns = ['A','B','C'])
        self.assertRaises(TypeError, self.frame.astype, np.int32, copy = False)

        self.assertRaises(TypeError, tf, np.int32, copy = False)

        self.assertRaises(TypeError, tf, np.int64, copy = False)
        casted = tf.astype(np.int64)

        self.assertRaises(TypeError, tf, np.float32, copy = False)
        casted = tf.astype(np.float32)

        # this is the only real reason to do it this way
        tf = np.round(self.frame).astype(np.int32)
        casted = tf.astype(np.float32, copy = False)
        #self.assert_(casted.values.data == tf.values.data)

        tf = self.frame.astype(np.float64)
        casted = tf.astype(np.int64, copy = False)
        #self.assert_(casted.values.data == tf.values.data)

        # can't view to an object array
        self.assertRaises(Exception, self.frame.astype, 'O', copy = False)

    def test_astype_cast_nan_int(self):
        df = DataFrame(data={"Values": [1.0, 2.0, 3.0, np.nan]})
        self.assertRaises(ValueError, df.astype, np.int64)

    def test_array_interface(self):
        result = np.sqrt(self.frame)
        self.assert_(type(result) is type(self.frame))
        self.assert_(result.index is self.frame.index)
        self.assert_(result.columns is self.frame.columns)

        assert_frame_equal(result, self.frame.apply(np.sqrt))

    def test_pickle(self):
        unpickled = pickle.loads(pickle.dumps(self.mixed_frame))
        assert_frame_equal(self.mixed_frame, unpickled)

        # buglet
        self.mixed_frame._data.ndim

        # empty
        unpickled = pickle.loads(pickle.dumps(self.empty))
        repr(unpickled)

    def test_to_dict(self):
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }
        recons_data = DataFrame(test_data).to_dict()

        for k, v in test_data.iteritems():
            for k2, v2 in v.iteritems():
                self.assertEqual(v2, recons_data[k][k2])

        recons_data = DataFrame(test_data).to_dict("l")

        for k, v in test_data.iteritems():
            for k2, v2 in v.iteritems():
                self.assertEqual(v2, recons_data[k][int(k2) - 1])

        recons_data = DataFrame(test_data).to_dict("s")

        for k, v in test_data.iteritems():
            for k2, v2 in v.iteritems():
                self.assertEqual(v2, recons_data[k][k2])

    def test_to_records_dt64(self):
        df = DataFrame([["one", "two", "three"],
                        ["four", "five", "six"]],
                       index=pan.date_range("2012-01-01", "2012-01-02"))
        self.assert_(df.to_records()['index'][0] == df.index[0])

        rs = df.to_records(convert_datetime64=False)
        self.assert_(rs['index'][0] == df.index.values[0])

    def test_to_records_with_multindex(self):
        # GH3189
        index = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                 ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        data = np.zeros((8, 4))
        df = DataFrame(data, index=index)
        r = df.to_records(index=True)['level_0']
        self.assertTrue('bar' in r)
        self.assertTrue('one' not in r)

    def test_to_records_with_Mapping_type(self):
        import email
        from email.parser import Parser
        import collections

        collections.Mapping.register(email.message.Message)

        headers = Parser().parsestr('From: <user@example.com>\n'
                'To: <someone_else@example.com>\n'
                'Subject: Test message\n'
                '\n'
                'Body would go here\n')

        frame = DataFrame.from_records([headers])
        all( x in frame for x in ['Type','Subject','From'])

    def test_from_records_to_records(self):
        # from numpy documentation
        arr = np.zeros((2,), dtype=('i4,f4,a10'))
        arr[:] = [(1, 2., 'Hello'), (2, 3., "World")]

        frame = DataFrame.from_records(arr)

        index = np.arange(len(arr))[::-1]
        indexed_frame = DataFrame.from_records(arr, index=index)
        self.assert_(np.array_equal(indexed_frame.index, index))

        # wrong length
        self.assertRaises(Exception, DataFrame.from_records, arr,
                          index=index[:-1])

        indexed_frame = DataFrame.from_records(arr, index='f1')
        self.assertRaises(Exception, DataFrame.from_records, np.zeros((2, 3)))

        # what to do?
        records = indexed_frame.to_records()
        self.assertEqual(len(records.dtype.names), 3)

        records = indexed_frame.to_records(index=False)
        self.assertEqual(len(records.dtype.names), 2)
        self.assert_('index' not in records.dtype.names)

    def test_from_records_nones(self):
        tuples = [(1, 2, None, 3),
                  (1, 2, None, 3),
                  (None, 2, 5, 3)]

        df = DataFrame.from_records(tuples, columns=['a', 'b', 'c', 'd'])
        self.assert_(np.isnan(df['c'][0]))

    def test_from_records_iterator(self):
        arr = np.array([(1.0, 1.0, 2, 2), (3.0, 3.0, 4, 4), (5., 5., 6, 6), (7., 7., 8, 8)],
                       dtype=[('x', np.float64), ('u', np.float32), ('y', np.int64), ('z', np.int32) ])
        df = DataFrame.from_records(iter(arr), nrows=2)
        xp = DataFrame({'x': np.array([1.0, 3.0], dtype=np.float64),
                        'u': np.array([1.0, 3.0], dtype=np.float32),
                        'y': np.array([2, 4], dtype=np.int64),
                        'z': np.array([2, 4], dtype=np.int32)})
        assert_frame_equal(df.reindex_like(xp), xp)

        # no dtypes specified here, so just compare with the default
        arr = [(1.0, 2), (3.0, 4), (5., 6), (7., 8)]
        df = DataFrame.from_records(iter(arr), columns=['x', 'y'],
                                    nrows=2)
        assert_frame_equal(df, xp.reindex(columns=['x','y']), check_dtype=False)

    def test_from_records_columns_not_modified(self):
        tuples = [(1, 2, 3),
                  (1, 2, 3),
                  (2, 5, 3)]

        columns = ['a', 'b', 'c']
        original_columns = list(columns)
        df = DataFrame.from_records(tuples, columns=columns, index='a')
        self.assertEqual(columns, original_columns)

    def test_from_records_decimal(self):
        from decimal import Decimal

        tuples = [(Decimal('1.5'),), (Decimal('2.5'),), (None,)]

        df = DataFrame.from_records(tuples, columns=['a'])
        self.assert_(df['a'].dtype == object)

        df = DataFrame.from_records(tuples, columns=['a'], coerce_float=True)
        self.assert_(df['a'].dtype == np.float64)
        self.assert_(np.isnan(df['a'].values[-1]))

    def test_from_records_duplicates(self):
        result = DataFrame.from_records([(1, 2, 3), (4, 5, 6)],
                                        columns=['a', 'b', 'a'])

        expected = DataFrame([(1, 2, 3), (4, 5, 6)],
                             columns=['a', 'b', 'a'])

        assert_frame_equal(result, expected)

    def test_from_records_set_index_name(self):
        def create_dict(order_id):
            return {'order_id': order_id, 'quantity': np.random.randint(1, 10),
                    'price': np.random.randint(1, 10)}
        documents = [create_dict(i) for i in range(10)]
        # demo missing data
        documents.append({'order_id': 10, 'quantity': 5})

        result = DataFrame.from_records(documents, index='order_id')
        self.assert_(result.index.name == 'order_id')

        # MultiIndex
        result = DataFrame.from_records(documents,
                                        index=['order_id', 'quantity'])
        self.assert_(result.index.names == ['order_id', 'quantity'])

    def test_from_records_misc_brokenness(self):
        # #2179

        data = {1: ['foo'], 2: ['bar']}

        result = DataFrame.from_records(data, columns=['a', 'b'])
        exp = DataFrame(data, columns=['a', 'b'])
        assert_frame_equal(result, exp)

        # overlap in index/index_names

        data = {'a': [1, 2, 3], 'b': [4, 5, 6]}

        result = DataFrame.from_records(data, index=['a', 'b', 'c'])
        exp = DataFrame(data, index=['a', 'b', 'c'])
        assert_frame_equal(result, exp)


        # GH 2623
        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), 'hi']) # test col upconverts to obj
        df2_obj = DataFrame.from_records(rows, columns=['date', 'test'])
        results = df2_obj.get_dtype_counts()
        expected = Series({ 'datetime64[ns]' : 1, 'object' : 1 })

        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), 1])
        df2_obj = DataFrame.from_records(rows, columns=['date', 'test'])
        results = df2_obj.get_dtype_counts()
        expected = Series({ 'datetime64[ns]' : 1, 'int64' : 1 })

    def test_from_records_empty(self):
        # 3562
        result = DataFrame.from_records([], columns=['a','b','c'])
        expected = DataFrame(columns=['a','b','c'])
        assert_frame_equal(result, expected)

        result = DataFrame.from_records([], columns=['a','b','b'])
        expected = DataFrame(columns=['a','b','b'])
        assert_frame_equal(result, expected)

    def test_from_records_empty_with_nonempty_fields_gh3682(self):
        a = np.array([(1, 2)], dtype=[('id', np.int64), ('value', np.int64)])
        df = DataFrame.from_records(a, index='id')
        assert_array_equal(df.index, Index([1], name='id'))
        self.assertEqual(df.index.name, 'id')
        assert_array_equal(df.columns, Index(['value']))

        b = np.array([], dtype=[('id', np.int64), ('value', np.int64)])
        df = DataFrame.from_records(b, index='id')
        assert_array_equal(df.index, Index([], name='id'))
        self.assertEqual(df.index.name, 'id')

    def test_to_records_floats(self):
        df = DataFrame(np.random.rand(10, 10))
        df.to_records()

    def test_to_recods_index_name(self):
        df = DataFrame(np.random.randn(3, 3))
        df.index.name = 'X'
        rs = df.to_records()
        self.assert_('X' in rs.dtype.fields)

        df = DataFrame(np.random.randn(3, 3))
        rs = df.to_records()
        self.assert_('index' in rs.dtype.fields)

        df.index = MultiIndex.from_tuples([('a', 'x'), ('a', 'y'), ('b', 'z')])
        df.index.names = ['A', None]
        rs = df.to_records()
        self.assert_('level_0' in rs.dtype.fields)

    def test_join_str_datetime(self):
        str_dates = ['20120209', '20120222']
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        A = DataFrame(str_dates, index=range(2), columns=['aa'])
        C = DataFrame([[1, 2], [3, 4]], index=str_dates, columns=dt_dates)

        tst = A.join(C, on='aa')

        self.assert_(len(tst.columns) == 3)

    def test_from_records_sequencelike(self):
        df = DataFrame({'A' : np.array(np.random.randn(6), dtype = np.float64),
                        'A1': np.array(np.random.randn(6), dtype = np.float64),
                        'B' : np.array(np.arange(6), dtype = np.int64),
                        'C' : ['foo'] * 6,
                        'D' : np.array([True, False] * 3, dtype=bool),
                        'E' : np.array(np.random.randn(6), dtype = np.float32),
                        'E1': np.array(np.random.randn(6), dtype = np.float32),
                        'F' : np.array(np.arange(6), dtype = np.int32) })

        # this is actually tricky to create the recordlike arrays and have the dtypes be intact
        blocks = df.blocks
        tuples = []
        columns = []
        dtypes  = []
        for dtype, b in blocks.iteritems():
            columns.extend(b.columns)
            dtypes.extend([ (c,np.dtype(dtype).descr[0][1]) for c in b.columns ])
        for i in xrange(len(df.index)):
            tup = []
            for _, b in blocks.iteritems():
                tup.extend(b.irow(i).values)
            tuples.append(tuple(tup))

        recarray  = np.array(tuples, dtype=dtypes).view(np.recarray)
        recarray2 = df.to_records()
        lists     = [list(x) for x in tuples]

        # tuples (lose the dtype info)
        result  = DataFrame.from_records(tuples,    columns=columns).reindex(columns=df.columns)

        # created recarray and with to_records recarray (have dtype info)
        result2 = DataFrame.from_records(recarray,  columns=columns).reindex(columns=df.columns)
        result3 = DataFrame.from_records(recarray2, columns=columns).reindex(columns=df.columns)

        # list of tupels (no dtype info)
        result4 = DataFrame.from_records(lists,     columns=columns).reindex(columns=df.columns)

        assert_frame_equal(result, df, check_dtype=False)
        assert_frame_equal(result2, df)
        assert_frame_equal(result3, df)
        assert_frame_equal(result4, df, check_dtype=False)

        # tuples is in the order of the columns
        result = DataFrame.from_records(tuples)
        self.assert_(np.array_equal(result.columns, range(8)))

        # test exclude parameter & we are casting the results here (as we don't have dtype info to recover)
        columns_to_test = [ columns.index('C'), columns.index('E1') ]

        exclude = list(set(xrange(8))-set(columns_to_test))
        result = DataFrame.from_records(tuples, exclude=exclude)
        result.columns = [ columns[i] for i in sorted(columns_to_test) ]
        assert_series_equal(result['C'], df['C'])
        assert_series_equal(result['E1'], df['E1'].astype('float64'))

        # empty case
        result = DataFrame.from_records([], columns=['foo', 'bar', 'baz'])
        self.assertEqual(len(result), 0)
        self.assert_(np.array_equal(result.columns, ['foo', 'bar', 'baz']))

        result = DataFrame.from_records([])
        self.assertEqual(len(result), 0)
        self.assertEqual(len(result.columns), 0)

    def test_from_records_dictlike(self):

        # test the dict methods
        df = DataFrame({'A' : np.array(np.random.randn(6), dtype = np.float64),
                        'A1': np.array(np.random.randn(6), dtype = np.float64),
                        'B' : np.array(np.arange(6), dtype = np.int64),
                        'C' : ['foo'] * 6,
                        'D' : np.array([True, False] * 3, dtype=bool),
                        'E' : np.array(np.random.randn(6), dtype = np.float32),
                        'E1': np.array(np.random.randn(6), dtype = np.float32),
                        'F' : np.array(np.arange(6), dtype = np.int32) })

        # columns is in a different order here than the actual items iterated from the dict
        columns = []
        for dtype, b in df.blocks.iteritems():
            columns.extend(b.columns)

        asdict    = dict((x, y) for x, y in df.iteritems())
        asdict2   = dict((x, y.values) for x, y in df.iteritems())

        # dict of series & dict of ndarrays (have dtype info)
        results = []
        results.append(DataFrame.from_records(asdict).reindex(columns=df.columns))
        results.append(DataFrame.from_records(asdict,    columns=columns).reindex(columns=df.columns))
        results.append(DataFrame.from_records(asdict2,   columns=columns).reindex(columns=df.columns))

        for r in results:
            assert_frame_equal(r, df)

    def test_from_records_with_index_data(self):
        df = DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])

        data = np.random.randn(10)
        df1 = DataFrame.from_records(df, index=data)
        assert(df1.index.equals(Index(data)))

    def test_from_records_bad_index_column(self):
        df = DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])

        # should pass
        df1 = DataFrame.from_records(df, index=['C'])
        assert(df1.index.equals(Index(df.C)))

        df1 = DataFrame.from_records(df, index='C')
        assert(df1.index.equals(Index(df.C)))

        # should fail
        self.assertRaises(ValueError, DataFrame.from_records, df, index=[2])
        self.assertRaises(KeyError, DataFrame.from_records, df, index=2)

    def test_from_records_non_tuple(self):
        class Record(object):

            def __init__(self, *args):
                self.args = args

            def __getitem__(self, i):
                return self.args[i]

            def __iter__(self):
                return iter(self.args)

        recs = [Record(1, 2, 3), Record(4, 5, 6), Record(7, 8, 9)]
        tups = map(tuple, recs)

        result = DataFrame.from_records(recs)
        expected = DataFrame.from_records(tups)
        assert_frame_equal(result, expected)

    def test_from_records_len0_with_columns(self):
        # #2633
        result = DataFrame.from_records([], index='foo',
                                        columns=['foo', 'bar'])

        self.assertTrue(np.array_equal(result.columns, ['bar']))
        self.assertEqual(len(result), 0)
        self.assertEqual(result.index.name, 'foo')

    def test_get_agg_axis(self):
        cols = self.frame._get_agg_axis(0)
        self.assert_(cols is self.frame.columns)

        idx = self.frame._get_agg_axis(1)
        self.assert_(idx is self.frame.index)

        self.assertRaises(Exception, self.frame._get_agg_axis, 2)

    def test_nonzero(self):
        self.assertTrue(self.empty.empty)

        self.assertFalse(self.frame.empty)
        self.assertFalse(self.mixed_frame.empty)

        # corner case
        df = DataFrame({'A': [1., 2., 3.],
                        'B': ['a', 'b', 'c']},
                       index=np.arange(3))
        del df['A']
        self.assertFalse(df.empty)

    def test_repr_empty(self):
        buf = StringIO()

        # empty
        foo = repr(self.empty)

        # empty with index
        frame = DataFrame(index=np.arange(1000))
        foo = repr(frame)

    def test_repr_mixed(self):
        buf = StringIO()

        # mixed
        foo = repr(self.mixed_frame)
        self.mixed_frame.info(verbose=False, buf=buf)

    @slow
    def test_repr_mixed_big(self):
        # big mixed
        biggie = DataFrame({'A': randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=range(200))
        biggie['A'][:20] = nan
        biggie['B'][:20] = nan

        foo = repr(biggie)

    def test_repr(self):
        buf = StringIO()

        # small one
        foo = repr(self.frame)
        self.frame.info(verbose=False, buf=buf)

        # even smaller
        self.frame.reindex(columns=['A']).info(verbose=False, buf=buf)
        self.frame.reindex(columns=['A', 'B']).info(verbose=False, buf=buf)

        # exhausting cases in DataFrame.info

        # columns but no index
        no_index = DataFrame(columns=[0, 1, 3])
        foo = repr(no_index)

        # no columns or index
        self.empty.info(buf=buf)

        df = DataFrame(["a\n\r\tb"], columns=["a\n\r\td"], index=["a\n\r\tf"])
        self.assertFalse("\t" in repr(df))
        self.assertFalse("\r" in repr(df))
        self.assertFalse("a\n" in repr(df))

    @slow
    def test_repr_big(self):
        buf = StringIO()

        # big one
        biggie = DataFrame(np.zeros((200, 4)), columns=range(4),
                           index=range(200))
        foo = repr(biggie)

    def test_repr_unsortable(self):
        # columns are not sortable
        import warnings
        warn_filters = warnings.filters
        warnings.filterwarnings('ignore',
                                category=FutureWarning,
                                module=".*format")

        unsortable = DataFrame({'foo': [1] * 50,
                                datetime.today(): [1] * 50,
                                'bar': ['bar'] * 50,
                                datetime.today(
                                ) + timedelta(1): ['bar'] * 50},
                               index=np.arange(50))
        foo = repr(unsortable)

        fmt.set_option('display.precision', 3, 'display.column_space', 10)
        repr(self.frame)

        fmt.set_option('display.max_rows', 10, 'display.max_columns', 2)
        repr(self.frame)

        fmt.set_option('display.max_rows', 1000, 'display.max_columns', 1000)
        repr(self.frame)

        fmt.reset_option('^display\.')

        warnings.filters = warn_filters

    def test_repr_unicode(self):
        uval = u'\u03c3\u03c3\u03c3\u03c3'
        bval = uval.encode('utf-8')
        df = DataFrame({'A': [uval, uval]})

        result = repr(df)
        ex_top = '      A'
        self.assertEqual(result.split('\n')[0].rstrip(), ex_top)

        df = DataFrame({'A': [uval, uval]})
        result = repr(df)
        self.assertEqual(result.split('\n')[0].rstrip(), ex_top)

    def test_unicode_string_with_unicode(self):
        df = DataFrame({'A': [u"\u05d0"]})

        if py3compat.PY3:
            str(df)
        else:
            unicode(df)

    def test_bytestring_with_unicode(self):
        df = DataFrame({'A': [u"\u05d0"]})
        if py3compat.PY3:
            bytes(df)
        else:
            str(df)

    def test_very_wide_info_repr(self):
        df = DataFrame(np.random.randn(10, 20),
                       columns=[tm.rands(10) for _ in xrange(20)])
        repr(df)

    def test_repr_column_name_unicode_truncation_bug(self):
        # #1906
        df = DataFrame({'Id': [7117434],
                        'StringCol': ('Is it possible to modify drop plot code'
                                      ' so that the output graph is displayed '
                                      'in iphone simulator, Is it possible to '
                                      'modify drop plot code so that the '
                                      'output graph is \xe2\x80\xa8displayed '
                                      'in iphone simulator.Now we are adding '
                                      'the CSV file externally. I want to Call'
                                      ' the File through the code..')})

        result = repr(df)
        self.assert_('StringCol' in result)

    def test_head_tail(self):
        assert_frame_equal(self.frame.head(), self.frame[:5])
        assert_frame_equal(self.frame.tail(), self.frame[-5:])

    def test_insert(self):
        df = DataFrame(np.random.randn(5, 3), index=np.arange(5),
                       columns=['c', 'b', 'a'])

        df.insert(0, 'foo', df['a'])
        self.assert_(np.array_equal(df.columns, ['foo', 'c', 'b', 'a']))
        assert_almost_equal(df['a'], df['foo'])

        df.insert(2, 'bar', df['c'])
        self.assert_(np.array_equal(df.columns, ['foo', 'c', 'bar', 'b', 'a']))
        assert_almost_equal(df['c'], df['bar'])

        # diff dtype

        # new item
        df['x'] = df['a'].astype('float32')
        result = Series(dict(float64 = 5, float32 = 1))
        self.assert_((df.get_dtype_counts() == result).all() == True)

        # replacing current (in different block)
        df['a'] = df['a'].astype('float32')
        result = Series(dict(float64 = 4, float32 = 2))
        self.assert_((df.get_dtype_counts() == result).all() == True)

        df['y'] = df['a'].astype('int32')
        result = Series(dict(float64 = 4, float32 = 2, int32 = 1))
        self.assert_((df.get_dtype_counts() == result).all() == True)

        self.assertRaises(Exception, df.insert, 1, 'a', df['b'])
        self.assertRaises(Exception, df.insert, 1, 'c', df['b'])

        df.columns.name = 'some_name'
        # preserve columns name field
        df.insert(0, 'baz', df['c'])
        self.assertEqual(df.columns.name, 'some_name')

    def test_delitem(self):
        del self.frame['A']
        self.assert_('A' not in self.frame)

    def test_pop(self):
        self.frame.columns.name = 'baz'

        A = self.frame.pop('A')
        self.assert_('A' not in self.frame)

        self.frame['foo'] = 'bar'
        foo = self.frame.pop('foo')
        self.assert_('foo' not in self.frame)
        # TODO self.assert_(self.frame.columns.name == 'baz')

    def test_pop_non_unique_cols(self):
        df = DataFrame({0: [0, 1], 1: [0, 1], 2: [4, 5]})
        df.columns = ["a", "b", "a"]

        res = df.pop("a")
        self.assertEqual(type(res), DataFrame)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(df.columns), 1)
        self.assertTrue("b" in df.columns)
        self.assertFalse("a" in df.columns)
        self.assertEqual(len(df.index), 2)

    def test_iter(self):
        self.assert_(tm.equalContents(list(self.frame), self.frame.columns))

    def test_iterrows(self):
        for i, (k, v) in enumerate(self.frame.iterrows()):
            exp = self.frame.xs(self.frame.index[i])
            assert_series_equal(v, exp)

        for i, (k, v) in enumerate(self.mixed_frame.iterrows()):
            exp = self.mixed_frame.xs(self.mixed_frame.index[i])
            assert_series_equal(v, exp)

    def test_itertuples(self):
        for i, tup in enumerate(self.frame.itertuples()):
            s = Series(tup[1:])
            s.name = tup[0]
            expected = self.frame.ix[i, :].reset_index(drop=True)
            assert_series_equal(s, expected)

        df = DataFrame({'floats': np.random.randn(5),
                        'ints': range(5)}, columns=['floats', 'ints'])

        for tup in df.itertuples(index=False):
            self.assert_(isinstance(tup[1], np.integer))

        df = DataFrame(data={"a": [1, 2, 3], "b": [4, 5, 6]})
        dfaa = df[['a', 'a']]
        self.assertEqual(list(dfaa.itertuples()), [(0, 1, 1), (1, 2, 2), (2, 3, 3)])

    def test_len(self):
        self.assertEqual(len(self.frame), len(self.frame.index))

    def test_operators(self):
        garbage = random.random(4)
        colSeries = Series(garbage, index=np.array(self.frame.columns))

        idSum = self.frame + self.frame
        seriesSum = self.frame + colSeries

        for col, series in idSum.iteritems():
            for idx, val in series.iteritems():
                origVal = self.frame[col][idx] * 2
                if not np.isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assert_(np.isnan(origVal))

        for col, series in seriesSum.iteritems():
            for idx, val in series.iteritems():
                origVal = self.frame[col][idx] + colSeries[col]
                if not np.isnan(val):
                    self.assertEqual(val, origVal)
                else:
                    self.assert_(np.isnan(origVal))

        added = self.frame2 + self.frame2
        expected = self.frame2 * 2
        assert_frame_equal(added, expected)

        df = DataFrame({'a': ['a', None, 'b']})
        assert_frame_equal(df + df, DataFrame({'a': ['aa', np.nan, 'bb']}))

    def test_operators_none_as_na(self):
        df = DataFrame({"col1": [2, 5.0, 123, None],
                        "col2": [1, 2, 3, 4]}, dtype=object)

        ops = [operator.add, operator.sub, operator.mul, operator.truediv]

        for op in ops:
            filled = df.fillna(np.nan)
            result = op(df, 3)
            expected = op(filled, 3)
            expected[com.isnull(expected)] = None
            assert_frame_equal(result, expected)

            result = op(df, df)
            expected = op(filled, filled)
            expected[com.isnull(expected)] = None
            assert_frame_equal(result, expected)

            result = op(df, df.fillna(7))
            assert_frame_equal(result, expected)

            result = op(df.fillna(7), df)
            assert_frame_equal(result, expected)

    def test_modulo(self):

        # GH3590, modulo as ints
        p = DataFrame({ 'first' : [3,4,5,8], 'second' : [0,0,0,3] })

        ### this is technically wrong as the integer portion is coerced to float ###
        expected = DataFrame({ 'first' : Series([0,0,0,0],dtype='float64'), 'second' : Series([np.nan,np.nan,np.nan,0]) })
        result = p % p
        assert_frame_equal(result,expected)

        # numpy has a slightly different (wrong) treatement
        result2 = DataFrame(p.values % p.values,index=p.index,columns=p.columns,dtype='float64')
        result2.iloc[0:3,1] = np.nan
        assert_frame_equal(result2,expected)

        result = p % 0
        expected = DataFrame(np.nan,index=p.index,columns=p.columns)
        assert_frame_equal(result,expected)

        # numpy has a slightly different (wrong) treatement
        result2 = DataFrame(p.values.astype('float64') % 0,index=p.index,columns=p.columns)
        assert_frame_equal(result2,expected)

        # not commutative with series
        p = DataFrame(np.random.randn(10, 5))
        s = p[0]
        res = s % p
        res2 = p % s
        self.assertFalse(np.array_equal(res.fillna(0), res2.fillna(0)))

    def test_div(self):

        # integer div, but deal with the 0's
        p = DataFrame({ 'first' : [3,4,5,8], 'second' : [0,0,0,3] })
        result = p / p

        ### this is technically wrong as the integer portion is coerced to float ###
        expected = DataFrame({ 'first' : Series([1,1,1,1],dtype='float64'), 'second' : Series([np.inf,np.inf,np.inf,1]) })
        assert_frame_equal(result,expected)

        result2 = DataFrame(p.values.astype('float64')/p.values,index=p.index,columns=p.columns).fillna(np.inf)
        assert_frame_equal(result2,expected)

        result = p / 0
        expected = DataFrame(np.inf,index=p.index,columns=p.columns)
        assert_frame_equal(result,expected)

        # numpy has a slightly different (wrong) treatement
        result2 = DataFrame(p.values.astype('float64')/0,index=p.index,columns=p.columns).fillna(np.inf)
        assert_frame_equal(result2,expected)

        p = DataFrame(np.random.randn(10, 5))
        s = p[0]
        res = s / p
        res2 = p / s
        self.assertFalse(np.array_equal(res.fillna(0), res2.fillna(0)))

    def test_logical_operators(self):
        import operator

        def _check_bin_op(op):
            result = op(df1, df2)
            expected = DataFrame(op(df1.values, df2.values), index=df1.index,
                                 columns=df1.columns)
            self.assert_(result.values.dtype == np.bool_)
            assert_frame_equal(result, expected)

        def _check_unary_op(op):
            result = op(df1)
            expected = DataFrame(op(df1.values), index=df1.index,
                                 columns=df1.columns)
            self.assert_(result.values.dtype == np.bool_)
            assert_frame_equal(result, expected)

        df1 = {'a': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True},
               'b': {'a': False, 'b': True, 'c': False,
                     'd': False, 'e': False},
               'c': {'a': False, 'b': False, 'c': True,
                     'd': False, 'e': False},
               'd': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True},
               'e': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True}}

        df2 = {'a': {'a': True, 'b': False, 'c': True, 'd': False, 'e': False},
               'b': {'a': False, 'b': True, 'c': False,
                     'd': False, 'e': False},
               'c': {'a': True, 'b': False, 'c': True, 'd': False, 'e': False},
               'd': {'a': False, 'b': False, 'c': False,
                     'd': True, 'e': False},
               'e': {'a': False, 'b': False, 'c': False,
                     'd': False, 'e': True}}

        df1 = DataFrame(df1)
        df2 = DataFrame(df2)

        _check_bin_op(operator.and_)
        _check_bin_op(operator.or_)
        _check_bin_op(operator.xor)

        _check_unary_op(operator.neg)

    def test_logical_typeerror(self):
        if py3compat.PY3:
            pass
        else:
            self.assertRaises(TypeError, self.frame.__eq__, 'foo')
            self.assertRaises(TypeError, self.frame.__lt__, 'foo')
            self.assertRaises(TypeError, self.frame.__gt__, 'foo')
            self.assertRaises(TypeError, self.frame.__ne__, 'foo')

    def test_constructor_lists_to_object_dtype(self):
        # from #1074
        d = DataFrame({'a': [np.nan, False]})
        self.assert_(d['a'].dtype == np.object_)
        self.assert_(d['a'][1] is False)

    def test_logical_with_nas(self):
        d = DataFrame({'a': [np.nan, False], 'b': [True, True]})

        result = d['a'] | d['b']
        expected = Series([np.nan, True])
        assert_series_equal(result, expected)

        result = d['a'].fillna(False) | d['b']
        expected = Series([True, True], dtype=object)
        assert_series_equal(result, expected)

    def test_neg(self):
        # what to do?
        assert_frame_equal(-self.frame, -1 * self.frame)

    def test_invert(self):
        assert_frame_equal(-(self.frame < 0), ~(self.frame < 0))

    def test_first_last_valid(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan
        mat[-5:] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        index = frame.first_valid_index()

        self.assert_(index == frame.index[5])

        index = frame.last_valid_index()
        self.assert_(index == frame.index[-6])

    def test_arith_flex_frame(self):
        ops = ['add', 'sub', 'mul', 'div', 'pow']
        aliases = {'div': 'truediv'}

        for op in ops:
            alias = aliases.get(op, op)
            f = getattr(operator, alias)
            result = getattr(self.frame, op)(2 * self.frame)
            exp = f(self.frame, 2 * self.frame)
            assert_frame_equal(result, exp)

            # vs mix float
            result = getattr(self.mixed_float, op)(2 * self.mixed_float)
            exp = f(self.mixed_float, 2 * self.mixed_float)
            assert_frame_equal(result, exp)
            _check_mixed_float(result, dtype = dict(C = None))

            # vs mix int
            if op in ['add','sub','mul']:
                result = getattr(self.mixed_int, op)(2 + self.mixed_int)
                exp = f(self.mixed_int, 2 + self.mixed_int)

                # overflow in the uint
                dtype = None
                if op in ['sub']:
                    dtype = dict(B = 'object', C = None)
                elif op in ['add','mul']:
                    dtype = dict(C = None)
                assert_frame_equal(result, exp)
                _check_mixed_int(result, dtype = dtype)

        # res_add = self.frame.add(self.frame)
        # res_sub = self.frame.sub(self.frame)
        # res_mul = self.frame.mul(self.frame)
        # res_div = self.frame.div(2 * self.frame)

        # assert_frame_equal(res_add, self.frame + self.frame)
        # assert_frame_equal(res_sub, self.frame - self.frame)
        # assert_frame_equal(res_mul, self.frame * self.frame)
        # assert_frame_equal(res_div, self.frame / (2 * self.frame))

        const_add = self.frame.add(1)
        assert_frame_equal(const_add, self.frame + 1)

        # corner cases
        result = self.frame.add(self.frame[:0])
        assert_frame_equal(result, self.frame * np.nan)

        result = self.frame[:0].add(self.frame)
        assert_frame_equal(result, self.frame * np.nan)

    def test_arith_mixed(self):

        left = DataFrame({'A': ['a', 'b', 'c'],
                          'B': [1, 2, 3]})

        result = left + left
        expected = DataFrame({'A': ['aa', 'bb', 'cc'],
                              'B': [2, 4, 6]})
        assert_frame_equal(result, expected)

    def test_arith_getitem_commute(self):
        df = DataFrame({'A': [1.1, 3.3], 'B': [2.5, -3.9]})

        self._test_op(df, operator.add)
        self._test_op(df, operator.sub)
        self._test_op(df, operator.mul)
        self._test_op(df, operator.truediv)
        self._test_op(df, operator.floordiv)
        self._test_op(df, operator.pow)

        self._test_op(df, lambda x, y: y + x)
        self._test_op(df, lambda x, y: y - x)
        self._test_op(df, lambda x, y: y * x)
        self._test_op(df, lambda x, y: y / x)
        self._test_op(df, lambda x, y: y ** x)

        self._test_op(df, lambda x, y: x + y)
        self._test_op(df, lambda x, y: x - y)
        self._test_op(df, lambda x, y: x * y)
        self._test_op(df, lambda x, y: x / y)
        self._test_op(df, lambda x, y: x ** y)

    @staticmethod
    def _test_op(df, op):
        result = op(df, 1)

        if not df.columns.is_unique:
            raise ValueError("Only unique columns supported by this test")

        for col in result.columns:
            assert_series_equal(result[col], op(df[col], 1))

    def test_bool_flex_frame(self):
        data = np.random.randn(5, 3)
        other_data = np.random.randn(5, 3)
        df = DataFrame(data)
        other = DataFrame(other_data)

        # No NAs

        # DataFrame
        self.assert_(df.eq(df).values.all())
        self.assert_(not df.ne(df).values.any())

        assert_frame_equal((df == other), df.eq(other))
        assert_frame_equal((df != other), df.ne(other))
        assert_frame_equal((df > other), df.gt(other))
        assert_frame_equal((df < other), df.lt(other))
        assert_frame_equal((df >= other), df.ge(other))
        assert_frame_equal((df <= other), df.le(other))

        # Unaligned
        def _check_unaligned_frame(meth, op, df, other, default=False):
            part_o = other.ix[3:, 1:].copy()
            rs = meth(df, part_o)
            xp = op(df, part_o.reindex(index=df.index, columns=df.columns))
            assert_frame_equal(rs, xp)

        _check_unaligned_frame(DataFrame.eq, operator.eq, df, other)
        _check_unaligned_frame(DataFrame.ne, operator.ne, df, other,
                               default=True)
        _check_unaligned_frame(DataFrame.gt, operator.gt, df, other)
        _check_unaligned_frame(DataFrame.lt, operator.lt, df, other)
        _check_unaligned_frame(DataFrame.ge, operator.ge, df, other)
        _check_unaligned_frame(DataFrame.le, operator.le, df, other)

        # Series
        def _test_seq(df, idx_ser, col_ser):
            idx_eq = df.eq(idx_ser, axis=0)
            col_eq = df.eq(col_ser)
            idx_ne = df.ne(idx_ser, axis=0)
            col_ne = df.ne(col_ser)
            assert_frame_equal(col_eq, df == Series(col_ser))
            assert_frame_equal(col_eq, -col_ne)
            assert_frame_equal(idx_eq, -idx_ne)
            assert_frame_equal(idx_eq, df.T.eq(idx_ser).T)
            assert_frame_equal(col_eq, df.eq(list(col_ser)))
            assert_frame_equal(idx_eq, df.eq(Series(idx_ser), axis=0))
            assert_frame_equal(idx_eq, df.eq(list(idx_ser), axis=0))

            idx_gt = df.gt(idx_ser, axis=0)
            col_gt = df.gt(col_ser)
            idx_le = df.le(idx_ser, axis=0)
            col_le = df.le(col_ser)

            assert_frame_equal(col_gt, df > Series(col_ser))
            assert_frame_equal(col_gt, -col_le)
            assert_frame_equal(idx_gt, -idx_le)
            assert_frame_equal(idx_gt, df.T.gt(idx_ser).T)

            idx_ge = df.ge(idx_ser, axis=0)
            col_ge = df.ge(col_ser)
            idx_lt = df.lt(idx_ser, axis=0)
            col_lt = df.lt(col_ser)
            assert_frame_equal(col_ge, df >= Series(col_ser))
            assert_frame_equal(col_ge, -col_lt)
            assert_frame_equal(idx_ge, -idx_lt)
            assert_frame_equal(idx_ge, df.T.ge(idx_ser).T)

        idx_ser = Series(np.random.randn(5))
        col_ser = Series(np.random.randn(3))
        _test_seq(df, idx_ser, col_ser)

        # ndarray

        assert_frame_equal((df == other.values), df.eq(other.values))
        assert_frame_equal((df != other.values), df.ne(other.values))
        assert_frame_equal((df > other.values), df.gt(other.values))
        assert_frame_equal((df < other.values), df.lt(other.values))
        assert_frame_equal((df >= other.values), df.ge(other.values))
        assert_frame_equal((df <= other.values), df.le(other.values))

        # list/tuple
        _test_seq(df, idx_ser.values, col_ser.values)

        # NA
        df.ix[0, 0] = np.nan
        rs = df.eq(df)
        self.assert_(not rs.ix[0, 0])
        rs = df.ne(df)
        self.assert_(rs.ix[0, 0])
        rs = df.gt(df)
        self.assert_(not rs.ix[0, 0])
        rs = df.lt(df)
        self.assert_(not rs.ix[0, 0])
        rs = df.ge(df)
        self.assert_(not rs.ix[0, 0])
        rs = df.le(df)
        self.assert_(not rs.ix[0, 0])

        # scalar
        assert_frame_equal(df.eq(0), df == 0)
        assert_frame_equal(df.ne(0), df != 0)
        assert_frame_equal(df.gt(0), df > 0)
        assert_frame_equal(df.lt(0), df < 0)
        assert_frame_equal(df.ge(0), df >= 0)
        assert_frame_equal(df.le(0), df <= 0)

        assert_frame_equal(df.eq(np.nan), df == np.nan)
        assert_frame_equal(df.ne(np.nan), df != np.nan)
        assert_frame_equal(df.gt(np.nan), df > np.nan)
        assert_frame_equal(df.lt(np.nan), df < np.nan)
        assert_frame_equal(df.ge(np.nan), df >= np.nan)
        assert_frame_equal(df.le(np.nan), df <= np.nan)

        # complex
        arr = np.array([np.nan, 1, 6, np.nan])
        arr2 = np.array([2j, np.nan, 7, None])
        df = DataFrame({'a': arr})
        df2 = DataFrame({'a': arr2})
        rs = df.gt(df2)
        self.assert_(not rs.values.any())
        rs = df.ne(df2)
        self.assert_(rs.values.all())

        arr3 = np.array([2j, np.nan, None])
        df3 = DataFrame({'a': arr3})
        rs = df3.gt(2j)
        self.assert_(not rs.values.any())

        # corner, dtype=object
        df1 = DataFrame({'col': ['foo', np.nan, 'bar']})
        df2 = DataFrame({'col': ['foo', datetime.now(), 'bar']})
        result = df1.ne(df2)
        exp = DataFrame({'col': [False, True, False]})
        assert_frame_equal(result, exp)

    def test_arith_flex_series(self):
        df = self.simple

        row = df.xs('a')
        col = df['two']

        assert_frame_equal(df.add(row), df + row)
        assert_frame_equal(df.add(row, axis=None), df + row)
        assert_frame_equal(df.sub(row), df - row)
        assert_frame_equal(df.div(row), df / row)
        assert_frame_equal(df.mul(row), df * row)
        assert_frame_equal(df.mod(row), df % row)

        assert_frame_equal(df.add(col, axis=0), (df.T + col).T)
        assert_frame_equal(df.sub(col, axis=0), (df.T - col).T)
        assert_frame_equal(df.div(col, axis=0), (df.T / col).T)
        assert_frame_equal(df.mul(col, axis=0), (df.T * col).T)
        assert_frame_equal(df.mod(col, axis=0), (df.T % col).T)

    def test_arith_non_pandas_object(self):
        df = self.simple

        val1 = df.xs('a').values
        added = DataFrame(df.values + val1, index=df.index, columns=df.columns)
        assert_frame_equal(df + val1, added)

        added = DataFrame((df.values.T + val1).T,
                          index=df.index, columns=df.columns)
        assert_frame_equal(df.add(val1, axis=0), added)

        val2 = list(df['two'])

        added = DataFrame(df.values + val2, index=df.index, columns=df.columns)
        assert_frame_equal(df + val2, added)

        added = DataFrame((df.values.T + val2).T, index=df.index,
                          columns=df.columns)
        assert_frame_equal(df.add(val2, axis='index'), added)

        val3 = np.random.rand(*df.shape)
        added = DataFrame(df.values + val3, index=df.index, columns=df.columns)
        assert_frame_equal(df.add(val3), added)

    def test_combineFrame(self):
        frame_copy = self.frame.reindex(self.frame.index[::2])

        del frame_copy['D']
        frame_copy['C'][:5] = nan

        added = self.frame + frame_copy
        tm.assert_dict_equal(added['A'].valid(),
                             self.frame['A'] * 2,
                             compare_keys=False)

        self.assert_(np.isnan(added['C'].reindex(frame_copy.index)[:5]).all())

        # assert(False)

        self.assert_(np.isnan(added['D']).all())

        self_added = self.frame + self.frame
        self.assert_(self_added.index.equals(self.frame.index))

        added_rev = frame_copy + self.frame
        self.assert_(np.isnan(added['D']).all())

        # corner cases

        # empty
        plus_empty = self.frame + self.empty
        self.assert_(np.isnan(plus_empty.values).all())

        empty_plus = self.empty + self.frame
        self.assert_(np.isnan(empty_plus.values).all())

        empty_empty = self.empty + self.empty
        self.assertTrue(empty_empty.empty)

        # out of order
        reverse = self.frame.reindex(columns=self.frame.columns[::-1])

        assert_frame_equal(reverse + self.frame, self.frame * 2)

        # mix vs float64, upcast
        added = self.frame + self.mixed_float
        _check_mixed_float(added, dtype = 'float64')
        added = self.mixed_float + self.frame
        _check_mixed_float(added, dtype = 'float64')

        # mix vs mix
        added = self.mixed_float + self.mixed_float2
        _check_mixed_float(added, dtype = dict(C = None))
        added = self.mixed_float2 + self.mixed_float
        _check_mixed_float(added, dtype = dict(C = None))

        # with int
        added = self.frame + self.mixed_int
        _check_mixed_float(added, dtype = 'float64')

    def test_combineSeries(self):

        # Series
        series = self.frame.xs(self.frame.index[0])

        added = self.frame + series

        for key, s in added.iteritems():
            assert_series_equal(s, self.frame[key] + series[key])

        larger_series = series.to_dict()
        larger_series['E'] = 1
        larger_series = Series(larger_series)
        larger_added = self.frame + larger_series

        for key, s in self.frame.iteritems():
            assert_series_equal(larger_added[key], s + series[key])
        self.assert_('E' in larger_added)
        self.assert_(np.isnan(larger_added['E']).all())

        # vs mix (upcast) as needed
        added = self.mixed_float + series
        _check_mixed_float(added, dtype = 'float64')
        added = self.mixed_float + series.astype('float32')
        _check_mixed_float(added, dtype = dict(C = None))
        added = self.mixed_float + series.astype('float16')
        _check_mixed_float(added, dtype = dict(C = None))

        #### these raise with numexpr.....as we are adding an int64 to an uint64....weird
        # vs int
        #added = self.mixed_int + (100*series).astype('int64')
        #_check_mixed_int(added, dtype = dict(A = 'int64', B = 'float64', C = 'int64', D = 'int64'))
        #added = self.mixed_int + (100*series).astype('int32')
        #_check_mixed_int(added, dtype = dict(A = 'int32', B = 'float64', C = 'int32', D = 'int64'))

        # TimeSeries
        import sys

        buf = StringIO()
        tmp = sys.stderr
        sys.stderr = buf

        try:
            ts = self.tsframe['A']
            added = self.tsframe + ts

            for key, col in self.tsframe.iteritems():
                assert_series_equal(added[key], col + ts)

            smaller_frame = self.tsframe[:-5]
            smaller_added = smaller_frame + ts

            self.assert_(smaller_added.index.equals(self.tsframe.index))

            smaller_ts = ts[:-5]
            smaller_added2 = self.tsframe + smaller_ts
            assert_frame_equal(smaller_added, smaller_added2)

            # length 0
            result = self.tsframe + ts[:0]

            # Frame is length 0
            result = self.tsframe[:0] + ts
            self.assertEqual(len(result), 0)

            # empty but with non-empty index
            frame = self.tsframe[:1].reindex(columns=[])
            result = frame * ts
            self.assertEqual(len(result), len(ts))
        finally:
            sys.stderr = tmp

    def test_combineFunc(self):
        result = self.frame * 2
        self.assert_(np.array_equal(result.values, self.frame.values * 2))

        # vs mix
        result = self.mixed_float * 2
        for c, s in result.iteritems():
            self.assert_(np.array_equal(s.values, self.mixed_float[c].values * 2))
        _check_mixed_float(result, dtype = dict(C = None))

        result = self.empty * 2
        self.assert_(result.index is self.empty.index)
        self.assertEqual(len(result.columns), 0)

    def test_comparisons(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame()

        row = self.simple.xs('a')

        def test_comp(func):
            result = func(df1, df2)
            self.assert_(np.array_equal(result.values,
                                        func(df1.values, df2.values)))

            result2 = func(self.simple, row)
            self.assert_(np.array_equal(result2.values,
                                        func(self.simple.values, row.values)))

            result3 = func(self.frame, 0)
            self.assert_(np.array_equal(result3.values,
                                        func(self.frame.values, 0)))

            self.assertRaises(Exception, func, self.simple, self.simple[:2])

        test_comp(operator.eq)
        test_comp(operator.ne)
        test_comp(operator.lt)
        test_comp(operator.gt)
        test_comp(operator.ge)
        test_comp(operator.le)

    def test_string_comparison(self):
        df = DataFrame([{"a": 1, "b": "foo"}, {"a": 2, "b": "bar"}])
        mask_a = df.a > 1
        assert_frame_equal(df[mask_a], df.ix[1:1, :])
        assert_frame_equal(df[-mask_a], df.ix[0:0, :])

        mask_b = df.b == "foo"
        assert_frame_equal(df[mask_b], df.ix[0:0, :])
        assert_frame_equal(df[-mask_b], df.ix[1:1, :])

    def test_float_none_comparison(self):
        df = DataFrame(np.random.randn(8, 3), index=range(8),
                       columns=['A', 'B', 'C'])

        self.assertRaises(TypeError, df.__eq__, None)

    def test_to_csv_deprecated_options(self):

        pname = '__tmp_to_csv_deprecated_options__'
        with ensure_clean(pname) as path:

            self.tsframe[1:3] = np.nan
            self.tsframe.to_csv(path, nanRep='foo')
            recons = read_csv(path,index_col=0,parse_dates=[0],na_values=['foo'])
            assert_frame_equal(self.tsframe, recons)

    def test_to_csv_from_csv(self):

        pname = '__tmp_to_csv_from_csv__'
        with ensure_clean(pname) as path:

             self.frame['A'][:5] = nan

             self.frame.to_csv(path)
             self.frame.to_csv(path, cols=['A', 'B'])
             self.frame.to_csv(path, header=False)
             self.frame.to_csv(path, index=False)

             # test roundtrip
             self.tsframe.to_csv(path)
             recons = DataFrame.from_csv(path)

             assert_frame_equal(self.tsframe, recons)

             self.tsframe.to_csv(path, index_label='index')
             recons = DataFrame.from_csv(path, index_col=None)
             assert(len(recons.columns) == len(self.tsframe.columns) + 1)

             # no index
             self.tsframe.to_csv(path, index=False)
             recons = DataFrame.from_csv(path, index_col=None)
             assert_almost_equal(self.tsframe.values, recons.values)

             # corner case
             dm = DataFrame({'s1': Series(range(3), range(3)),
                             's2': Series(range(2), range(2))})
             dm.to_csv(path)
             recons = DataFrame.from_csv(path)
             assert_frame_equal(dm, recons)

        with ensure_clean(pname) as path:

             # duplicate index
             df = DataFrame(np.random.randn(3, 3), index=['a', 'a', 'b'],
                            columns=['x', 'y', 'z'])
             df.to_csv(path)
             result = DataFrame.from_csv(path)
             assert_frame_equal(result, df)

             midx = MultiIndex.from_tuples([('A', 1, 2), ('A', 1, 2), ('B', 1, 2)])
             df = DataFrame(np.random.randn(3, 3), index=midx,
                            columns=['x', 'y', 'z'])
             df.to_csv(path)
             result = DataFrame.from_csv(path, index_col=[0, 1, 2],
                                         parse_dates=False)
             assert_frame_equal(result, df, check_names=False)  # TODO from_csv names index ['Unnamed: 1', 'Unnamed: 2'] should it ?

             # column aliases
             col_aliases = Index(['AA', 'X', 'Y', 'Z'])
             self.frame2.to_csv(path, header=col_aliases)
             rs = DataFrame.from_csv(path)
             xp = self.frame2.copy()
             xp.columns = col_aliases

             assert_frame_equal(xp, rs)

             self.assertRaises(ValueError, self.frame2.to_csv, path,
                               header=['AA', 'X'])

        with ensure_clean(pname) as path:
            import pandas as pd
            df1 = DataFrame(np.random.randn(3, 1))
            df2 = DataFrame(np.random.randn(3, 1))

            df1.to_csv(path)
            df2.to_csv(path,mode='a',header=False)
            xp = pd.concat([df1,df2])
            rs = pd.read_csv(path,index_col=0)
            rs.columns = map(int,rs.columns)
            xp.columns = map(int,xp.columns)
            assert_frame_equal(xp,rs)

    def test_to_csv_cols_reordering(self):
        # GH3454
        import pandas as pd

        def _check_df(df,cols=None):
            with ensure_clean() as path:
                df.to_csv(path,cols = cols,engine='python')
                rs_p = pd.read_csv(path,index_col=0)
                df.to_csv(path,cols = cols,chunksize=chunksize)
                rs_c = pd.read_csv(path,index_col=0)

            if cols:
                df = df[cols]
            assert (rs_c.columns==rs_p.columns).all()
            assert_frame_equal(df,rs_c,check_names=False)

        chunksize=5
        N = int(chunksize*2.5)

        df= mkdf(N, 3)
        cs = df.columns
        cols = [cs[2],cs[0]]
        _check_df(df,cols)

    def test_to_csv_legacy_raises_on_dupe_cols(self):
        df= mkdf(10, 3)
        df.columns = ['a','a','b']
        with ensure_clean() as path:
            self.assertRaises(NotImplementedError,df.to_csv,path,engine='python')

    def test_to_csv_new_dupe_cols(self):
        import pandas as pd
        def _check_df(df,cols=None):
            with ensure_clean() as path:
                df.to_csv(path,cols = cols,chunksize=chunksize)
                rs_c = pd.read_csv(path,index_col=0)

                # we wrote them in a different order
                # so compare them in that order
                if cols is not None:

                    if df.columns.is_unique:
                        rs_c.columns = cols
                    else:
                        indexer, missing = df.columns.get_indexer_non_unique(cols)
                        rs_c.columns = df.columns.take(indexer)

                    for c in cols:
                       obj_df = df[c]
                       obj_rs = rs_c[c]
                       if isinstance(obj_df,Series):
                           assert_series_equal(obj_df,obj_rs)
                       else:
                           assert_frame_equal(obj_df,obj_rs,check_names=False)

                # wrote in the same order
                else:
                    rs_c.columns = df.columns
                    assert_frame_equal(df,rs_c,check_names=False)

        chunksize=5
        N = int(chunksize*2.5)

        # dupe cols
        df= mkdf(N, 3)
        df.columns = ['a','a','b']
        _check_df(df,None)

        # dupe cols with selection
        cols = ['b','a']
        _check_df(df,cols)

    @slow
    def test_to_csv_moar(self):
        path = '__tmp_to_csv_moar__'

        def _do_test(df,path,r_dtype=None,c_dtype=None,rnlvl=None,cnlvl=None,
                     dupe_col=False):

               if cnlvl:
                   header = range(cnlvl)
                   with ensure_clean(path) as path:
                        df.to_csv(path,encoding='utf8',chunksize=chunksize,tupleize_cols=False)
                        recons = DataFrame.from_csv(path,header=range(cnlvl),tupleize_cols=False,parse_dates=False)
               else:
                   with ensure_clean(path) as path:
                       df.to_csv(path,encoding='utf8',chunksize=chunksize)
                       recons = DataFrame.from_csv(path,header=0,parse_dates=False)

               def _to_uni(x):
                   if not isinstance(x,unicode):
                       return x.decode('utf8')
                   return x
               if dupe_col:
                   # read_Csv disambiguates the columns by
                   # labeling them dupe.1,dupe.2, etc'. monkey patch columns
                   recons.columns = df.columns
               if rnlvl:
                   delta_lvl = [recons.icol(i).values for i in range(rnlvl-1)]
                   ix=MultiIndex.from_arrays([list(recons.index)]+delta_lvl)
                   recons.index = ix
                   recons = recons.iloc[:,rnlvl-1:]

               type_map = dict(i='i',f='f',s='O',u='O',dt='O',p='O')
               if r_dtype:
                    if r_dtype == 'u': # unicode
                        r_dtype='O'
                        recons.index = np.array(map(_to_uni,recons.index),
                                                dtype=r_dtype )
                        df.index = np.array(map(_to_uni,df.index),dtype=r_dtype )
                    if r_dtype == 'dt': # unicode
                        r_dtype='O'
                        recons.index = np.array(map(Timestamp,recons.index),
                                                dtype=r_dtype )
                        df.index = np.array(map(Timestamp,df.index),dtype=r_dtype )
                    elif r_dtype == 'p':
                        r_dtype='O'
                        recons.index = np.array(map(Timestamp,recons.index.to_datetime()),
                                                dtype=r_dtype )
                        df.index = np.array(map(Timestamp,df.index.to_datetime()),dtype=r_dtype )
                    else:
                        r_dtype= type_map.get(r_dtype)
                        recons.index = np.array(recons.index,dtype=r_dtype )
                        df.index = np.array(df.index,dtype=r_dtype )
               if c_dtype:
                    if c_dtype == 'u':
                        c_dtype='O'
                        recons.columns = np.array(map(_to_uni,recons.columns),
                                                dtype=c_dtype )
                        df.columns = np.array(map(_to_uni,df.columns),dtype=c_dtype )
                    elif c_dtype == 'dt':
                        c_dtype='O'
                        recons.columns = np.array(map(Timestamp,recons.columns),
                                                dtype=c_dtype )
                        df.columns = np.array(map(Timestamp,df.columns),dtype=c_dtype )
                    elif c_dtype == 'p':
                        c_dtype='O'
                        recons.columns = np.array(map(Timestamp,recons.columns.to_datetime()),
                                                dtype=c_dtype )
                        df.columns = np.array(map(Timestamp,df.columns.to_datetime()),dtype=c_dtype )
                    else:
                        c_dtype= type_map.get(c_dtype)
                        recons.columns = np.array(recons.columns,dtype=c_dtype )
                        df.columns = np.array(df.columns,dtype=c_dtype )

               assert_frame_equal(df, recons,check_names=False,check_less_precise=True)

        N = 100
        chunksize=1000

        # GH3437
        from pandas import NaT
        def make_dtnat_arr(n,nnat=None):
             if nnat is None:
                 nnat= int(n*0.1) # 10%
             s=list(date_range('2000',freq='5min',periods=n))
             if nnat:
                 for i in np.random.randint(0,len(s),nnat):
                     s[i] = NaT
                 i = np.random.randint(100)
                 s[-i] = NaT
                 s[i] = NaT
             return s

        # N=35000
        s1=make_dtnat_arr(chunksize+5)
        s2=make_dtnat_arr(chunksize+5,0)

        # s3=make_dtnat_arr(chunksize+5,0)
        with ensure_clean('1.csv') as path:
            df=DataFrame(dict(a=s1,b=s2))
            df.to_csv(path,chunksize=chunksize)
            recons = DataFrame.from_csv(path).convert_objects('coerce')
            assert_frame_equal(df, recons,check_names=False,check_less_precise=True)

        for ncols in [4]:
            base = int((chunksize// ncols or 1) or 1)
            for nrows in [2,10,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                  base-1,base,base+1]:
                _do_test(mkdf(nrows, ncols,r_idx_type='dt',
                              c_idx_type='s'),path, 'dt','s')
                pass


        for ncols in [4]:
            base = int((chunksize// ncols or 1) or 1)
            for nrows in [2,10,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                  base-1,base,base+1]:
                _do_test(mkdf(nrows, ncols,r_idx_type='dt',
                              c_idx_type='s'),path, 'dt','s')
                pass

        for r_idx_type,c_idx_type  in [('i','i'),('s','s'),('u','dt'),('p','p')]:
            for ncols in [1,2,3,4]:
                base = int((chunksize// ncols or 1) or 1)
                for nrows in [2,10,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                      base-1,base,base+1]:
                    _do_test(mkdf(nrows, ncols,r_idx_type=r_idx_type,
                                  c_idx_type=c_idx_type),path,r_idx_type,c_idx_type)

        for ncols in [1,2,3,4]:
            base = int((chunksize// ncols or 1) or 1)
            for nrows in [10,N-2,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                      base-1,base,base+1]:
                print( nrows,ncols)
                _do_test(mkdf(nrows, ncols),path)

        for nrows in [10,N-2,N-1,N,N+1,N+2]:
            df = mkdf(nrows, 3)
            cols = list(df.columns)
            cols[:2] = ["dupe","dupe"]
            cols[-2:] = ["dupe","dupe"]
            ix = list(df.index)
            ix[:2] = ["rdupe","rdupe"]
            ix[-2:] = ["rdupe","rdupe"]
            df.index=ix
            df.columns=cols
            _do_test(df,path,dupe_col=True)


        _do_test(DataFrame(index=range(10)),path)
        _do_test(mkdf(chunksize//2+1, 2,r_idx_nlevels=2),path,rnlvl=2)
        for ncols in [2,3,4]:
            base = int(chunksize//ncols)
            for nrows in [10,N-2,N-1,N,N+1,N+2,2*N-2,2*N-1,2*N,2*N+1,2*N+2,
                      base-1,base,base+1]:
                print(nrows, ncols)
                _do_test(mkdf(nrows, ncols,r_idx_nlevels=2),path,rnlvl=2)
                _do_test(mkdf(nrows, ncols,c_idx_nlevels=2),path,cnlvl=2)
                _do_test(mkdf(nrows, ncols,r_idx_nlevels=2,c_idx_nlevels=2),
                         path,rnlvl=2,cnlvl=2)



    def test_to_csv_from_csv_w_some_infs(self):

        # test roundtrip with inf, -inf, nan, as full columns and mix
        self.frame['G'] = np.nan
        f = lambda x: [np.inf, np.nan][np.random.rand() < .5]
        self.frame['H'] = self.frame.index.map(f)

        with ensure_clean() as path:
             self.frame.to_csv(path)
             recons = DataFrame.from_csv(path)

             assert_frame_equal(self.frame, recons, check_names=False)  # TODO to_csv drops column name
             assert_frame_equal(np.isinf(self.frame), np.isinf(recons), check_names=False)

    def test_to_csv_from_csv_w_all_infs(self):

        # test roundtrip with inf, -inf, nan, as full columns and mix
        self.frame['E'] = np.inf
        self.frame['F'] = -np.inf

        with ensure_clean() as path:
            self.frame.to_csv(path)
            recons = DataFrame.from_csv(path)

            assert_frame_equal(self.frame, recons, check_names=False)  # TODO to_csv drops column name
            assert_frame_equal(np.isinf(self.frame), np.isinf(recons), check_names=False)

    def test_to_csv_no_index(self):
        # GH 3624, after appending columns, to_csv fails
        pname = '__tmp_to_csv_no_index__'
        with ensure_clean(pname) as path:
            df = DataFrame({'c1':[1,2,3], 'c2':[4,5,6]})
            df.to_csv(path, index=False)
            result = read_csv(path)
            assert_frame_equal(df,result)
            df['c3'] = Series([7,8,9],dtype='int64')
            df.to_csv(path, index=False)
            result = read_csv(path)
            assert_frame_equal(df,result)

    def test_to_csv_multiindex(self):

        pname = '__tmp_to_csv_multiindex__'
        frame = self.frame
        old_index = frame.index
        arrays = np.arange(len(old_index) * 2).reshape(2, -1)
        new_index = MultiIndex.from_arrays(arrays, names=['first', 'second'])
        frame.index = new_index

        with ensure_clean(pname) as path:

             frame.to_csv(path, header=False)
             frame.to_csv(path, cols=['A', 'B'])

             # round trip
             frame.to_csv(path)
             df = DataFrame.from_csv(path, index_col=[0, 1], parse_dates=False)

             assert_frame_equal(frame, df, check_names=False)  # TODO to_csv drops column name
             self.assertEqual(frame.index.names, df.index.names)
             self.frame.index = old_index  # needed if setUP becomes a classmethod

             # try multiindex with dates
             tsframe = self.tsframe
             old_index = tsframe.index
             new_index = [old_index, np.arange(len(old_index))]
             tsframe.index = MultiIndex.from_arrays(new_index)

             tsframe.to_csv(path, index_label=['time', 'foo'])
             recons = DataFrame.from_csv(path, index_col=[0, 1])
             assert_frame_equal(tsframe, recons, check_names=False)  # TODO to_csv drops column name

             # do not load index
             tsframe.to_csv(path)
             recons = DataFrame.from_csv(path, index_col=None)
             np.testing.assert_equal(len(recons.columns), len(tsframe.columns) + 2)

             # no index
             tsframe.to_csv(path, index=False)
             recons = DataFrame.from_csv(path, index_col=None)
             assert_almost_equal(recons.values, self.tsframe.values)
             self.tsframe.index = old_index  # needed if setUP becomes classmethod

        with ensure_clean(pname) as path:
            # GH3571, GH1651, GH3141

            def _make_frame(names=None):
                if names is True:
                    names = ['first','second']
                return DataFrame(np.random.randint(0,10,size=(3,3)),
                                 columns=MultiIndex.from_tuples([('bah', 'foo'),
                                                                 ('bah', 'bar'),
                                                                 ('ban', 'baz')],
                                                                names=names),
                                 dtype='int64')

            # column & index are multi-index
            df = mkdf(5,3,r_idx_nlevels=2,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1,2,3],index_col=[0,1],tupleize_cols=False)
            assert_frame_equal(df,result)

            # column is mi
            df = mkdf(5,3,r_idx_nlevels=1,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1,2,3],index_col=0,tupleize_cols=False)
            assert_frame_equal(df,result)

            # dup column names?
            df = mkdf(5,3,r_idx_nlevels=3,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1,2,3],index_col=[0,1],tupleize_cols=False)
            result.columns = ['R2','A','B','C']
            new_result = result.reset_index().set_index(['R0','R1','R2'])
            new_result.columns = df.columns
            assert_frame_equal(df,new_result)

            # writing with no index
            df = _make_frame()
            df.to_csv(path,tupleize_cols=False,index=False)
            result = read_csv(path,header=[0,1],tupleize_cols=False)
            assert_frame_equal(df,result)

            # we lose the names here
            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=False,index=False)
            result = read_csv(path,header=[0,1],tupleize_cols=False)
            self.assert_(all([ x is None for x in result.columns.names ]))
            result.columns.names = df.columns.names
            assert_frame_equal(df,result)

            # tupleize_cols=True and index=False
            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=True,index=False)
            result = read_csv(path,header=0,tupleize_cols=True,index_col=None)
            result.columns = df.columns
            assert_frame_equal(df,result)

            # whatsnew example
            df = _make_frame()
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1],index_col=[0],tupleize_cols=False)
            assert_frame_equal(df,result)

            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=False)
            result = read_csv(path,header=[0,1],index_col=[0],tupleize_cols=False)
            assert_frame_equal(df,result)

            # column & index are multi-index (compatibility)
            df = mkdf(5,3,r_idx_nlevels=2,c_idx_nlevels=4)
            df.to_csv(path,tupleize_cols=True)
            result = read_csv(path,header=0,index_col=[0,1],tupleize_cols=True)
            result.columns = df.columns
            assert_frame_equal(df,result)

            # invalid options
            df = _make_frame(True)
            df.to_csv(path,tupleize_cols=False)

            # catch invalid headers
            def testit():
                read_csv(path,tupleize_cols=False,header=range(3),index_col=0)
            assertRaisesRegexp(CParserError, 'Passed header=\[0,1,2\] are too many rows for this multi_index of columns', testit)

            def testit():
                read_csv(path,tupleize_cols=False,header=range(7),index_col=0)
            assertRaisesRegexp(CParserError, 'Passed header=\[0,1,2,3,4,5,6\], len of 7, but only 6 lines in file', testit)

            for i in [3,4,5,6,7]:
                 self.assertRaises(Exception, read_csv, path, tupleize_cols=False, header=range(i), index_col=0)
            self.assertRaises(Exception, read_csv, path, tupleize_cols=False, header=[0,2], index_col=0)

            # write with cols
            self.assertRaises(Exception, df.to_csv, path,tupleize_cols=False,cols=['foo','bar'])

        with ensure_clean(pname) as path:
            # empty
            tsframe[:0].to_csv(path)
            recons = DataFrame.from_csv(path)
            exp = tsframe[:0]
            exp.index = []

            self.assert_(recons.columns.equals(exp.columns))
            self.assert_(len(recons) == 0)

    def test_to_csv_float32_nanrep(self):
        df = DataFrame(np.random.randn(1, 4).astype(np.float32))
        df[1] = np.nan

        with ensure_clean('__tmp_to_csv_float32_nanrep__.csv') as path:
            df.to_csv(path, na_rep=999)

            with open(path) as f:
                 lines = f.readlines()
                 self.assert_(lines[1].split(',')[2] == '999')

    def test_to_csv_withcommas(self):

        # Commas inside fields should be correctly escaped when saving as CSV.
        df = DataFrame({'A': [1, 2, 3], 'B': ['5,6', '7,8', '9,0']})

        with ensure_clean('__tmp_to_csv_withcommas__.csv') as path:
            df.to_csv(path)
            df2 = DataFrame.from_csv(path)
            assert_frame_equal(df2, df)

    def test_to_csv_mixed(self):

        def create_cols(name):
            return [ "%s%03d" % (name,i) for i in xrange(5) ]

        df_float  = DataFrame(np.random.randn(100, 5),dtype='float64',columns=create_cols('float'))
        df_int    = DataFrame(np.random.randn(100, 5),dtype='int64',columns=create_cols('int'))
        df_bool   = DataFrame(True,index=df_float.index,columns=create_cols('bool'))
        df_object = DataFrame('foo',index=df_float.index,columns=create_cols('object'))
        df_dt     = DataFrame(Timestamp('20010101'),index=df_float.index,columns=create_cols('date'))

        # add in some nans
        df_float.ix[30:50,1:3] = np.nan

        #### this is a bug in read_csv right now ####
        #df_dt.ix[30:50,1:3] = np.nan

        df        = pan.concat([ df_float, df_int, df_bool, df_object, df_dt ], axis=1)

        # dtype
        dtypes = dict()
        for n,dtype in [('float',np.float64),('int',np.int64),('bool',np.bool),('object',np.object)]:
            for c in create_cols(n):
                dtypes[c] = dtype

        with ensure_clean() as filename:
            df.to_csv(filename)
            rs = pan.read_csv(filename, index_col=0, dtype=dtypes, parse_dates=create_cols('date'))
            assert_frame_equal(rs, df)

    def test_to_csv_dups_cols(self):

        df        = DataFrame(np.random.randn(1000, 30),columns=range(15)+range(15),dtype='float64')

        with ensure_clean() as filename:
            df.to_csv(filename) # single dtype, fine
            result = read_csv(filename,index_col=0)
            result.columns = df.columns
            assert_frame_equal(result,df)

        df_float  = DataFrame(np.random.randn(1000, 3),dtype='float64')
        df_int    = DataFrame(np.random.randn(1000, 3),dtype='int64')
        df_bool   = DataFrame(True,index=df_float.index,columns=range(3))
        df_object = DataFrame('foo',index=df_float.index,columns=range(3))
        df_dt     = DataFrame(Timestamp('20010101'),index=df_float.index,columns=range(3))
        df        = pan.concat([ df_float, df_int, df_bool, df_object, df_dt ], axis=1, ignore_index=True)

        cols = []
        for i in range(5):
            cols.extend([0,1,2])
        df.columns = cols

        from pandas import to_datetime
        with ensure_clean() as filename:
            df.to_csv(filename)
            result = read_csv(filename,index_col=0)

            # date cols
            for i in ['0.4','1.4','2.4']:
                 result[i] = to_datetime(result[i])

            result.columns = df.columns
            assert_frame_equal(result,df)

        # GH3457
        from pandas.util.testing import makeCustomDataframe as mkdf

        N=10
        df= mkdf(N, 3)
        df.columns = ['a','a','b']

        with ensure_clean() as filename:
            df.to_csv(filename)

            # read_csv will rename the dups columns
            result = read_csv(filename,index_col=0)
            result = result.rename(columns={ 'a.1' : 'a' })
            assert_frame_equal(result,df)

    def test_to_csv_chunking(self):

        aa=DataFrame({'A':range(100000)})
        aa['B'] = aa.A + 1.0
        aa['C'] = aa.A + 2.0
        aa['D'] = aa.A + 3.0

        for chunksize in [10000,50000,100000]:
            with ensure_clean() as filename:
                aa.to_csv(filename,chunksize=chunksize)
                rs = pan.read_csv(filename,index_col=0)
                assert_frame_equal(rs, aa)

    def test_to_csv_bug(self):
        f1 = StringIO('a,1.0\nb,2.0')
        df = DataFrame.from_csv(f1, header=None)
        newdf = DataFrame({'t': df[df.columns[0]]})

        with ensure_clean() as path:
            newdf.to_csv(path)

            recons = pan.read_csv(path, index_col=0)
            assert_frame_equal(recons, newdf, check_names=False)  # don't check_names as t != 1

    def test_to_csv_unicode(self):

        df = DataFrame({u'c/\u03c3': [1, 2, 3]})
        with ensure_clean() as path:

            df.to_csv(path, encoding='UTF-8')
            df2 = pan.read_csv(path, index_col=0, encoding='UTF-8')
            assert_frame_equal(df, df2)

            df.to_csv(path, encoding='UTF-8', index=False)
            df2 = pan.read_csv(path, index_col=None, encoding='UTF-8')
            assert_frame_equal(df, df2)

    def test_to_csv_unicode_index_col(self):
        buf = StringIO('')
        df = DataFrame(
            [[u"\u05d0", "d2", "d3", "d4"], ["a1", "a2", "a3", "a4"]],
            columns=[u"\u05d0",
                     u"\u05d1", u"\u05d2", u"\u05d3"],
            index=[u"\u05d0", u"\u05d1"])

        df.to_csv(buf, encoding='UTF-8')
        buf.seek(0)

        df2 = pan.read_csv(buf, index_col=0, encoding='UTF-8')
        assert_frame_equal(df, df2)

    def test_to_csv_stringio(self):
        buf = StringIO()
        self.frame.to_csv(buf)
        buf.seek(0)
        recons = pan.read_csv(buf, index_col=0)
        assert_frame_equal(recons, self.frame, check_names=False)  # TODO to_csv drops column name

    def test_to_csv_float_format(self):

        df = DataFrame([[0.123456, 0.234567, 0.567567],
                        [12.32112, 123123.2, 321321.2]],
                       index=['A', 'B'], columns=['X', 'Y', 'Z'])

        with ensure_clean() as filename:

            df.to_csv(filename, float_format='%.2f')

            rs = pan.read_csv(filename, index_col=0)
            xp = DataFrame([[0.12, 0.23, 0.57],
                            [12.32, 123123.20, 321321.20]],
                           index=['A', 'B'], columns=['X', 'Y', 'Z'])
            assert_frame_equal(rs, xp)

    def test_to_csv_quoting(self):
        import csv

        df = DataFrame({'A': [1, 2, 3], 'B': ['foo', 'bar', 'baz']})

        buf = StringIO()
        df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC)

        result = buf.getvalue()
        expected = ('"A","B"\n'
                    '1,"foo"\n'
                    '2,"bar"\n'
                    '3,"baz"\n')

        self.assertEqual(result, expected)

        # quoting windows line terminators, presents with encoding?
        # #3503
        text = 'a,b,c\n1,"test \r\n",3\n'
        df = pd.read_csv(StringIO(text))
        buf = StringIO()
        df.to_csv(buf, encoding='utf-8', index=False)
        self.assertEqual(buf.getvalue(), text)

    def test_to_csv_unicodewriter_quoting(self):
        import csv

        df = DataFrame({'A': [1, 2, 3], 'B': ['foo', 'bar', 'baz']})

        buf = StringIO()
        df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC,
                  encoding='utf-8')

        result = buf.getvalue()
        expected = ('"A","B"\n'
                    '1,"foo"\n'
                    '2,"bar"\n'
                    '3,"baz"\n')

        self.assertEqual(result, expected)

    def test_to_csv_index_no_leading_comma(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]},
                       index=['one', 'two', 'three'])

        buf = StringIO()
        df.to_csv(buf, index_label=False)
        expected = ('A,B\n'
                    'one,1,4\n'
                    'two,2,5\n'
                    'three,3,6\n')
        self.assertEqual(buf.getvalue(), expected)

    def test_to_csv_line_terminators(self):
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]},
                       index=['one', 'two', 'three'])

        buf = StringIO()
        df.to_csv(buf, line_terminator='\r\n')
        expected = (',A,B\r\n'
                    'one,1,4\r\n'
                    'two,2,5\r\n'
                    'three,3,6\r\n')
        self.assertEqual(buf.getvalue(), expected)

        buf = StringIO()
        df.to_csv(buf)  # The default line terminator remains \n
        expected = (',A,B\n'
                    'one,1,4\n'
                    'two,2,5\n'
                    'three,3,6\n')
        self.assertEqual(buf.getvalue(), expected)

    def test_info(self):
        io = StringIO()
        self.frame.info(buf=io)
        self.tsframe.info(buf=io)

        frame = DataFrame(np.random.randn(5, 3))

        import sys
        sys.stdout = StringIO()
        frame.info()
        frame.info(verbose=False)
        sys.stdout = sys.__stdout__

    def test_info_wide(self):
        from pandas import set_option, reset_option
        io = StringIO()
        df = DataFrame(np.random.randn(5, 101))
        df.info(buf=io)
        rs = io.getvalue()
        self.assert_(len(rs.splitlines()) == 4)

        io = StringIO()
        df.info(buf=io, max_cols=101)
        rs = io.getvalue()
        self.assert_(len(rs.splitlines()) > 100)
        xp = rs

        set_option('display.max_info_columns', 101)
        io = StringIO()
        df.info(buf=io)
        self.assert_(rs == xp)
        reset_option('display.max_info_columns')

    def test_info_duplicate_columns(self):
        io = StringIO()

        # it works!
        frame = DataFrame(np.random.randn(1500, 4),
                          columns=['a', 'a', 'b', 'b'])
        frame.info(buf=io)

    def test_dtypes(self):
        self.mixed_frame['bool'] = self.mixed_frame['A'] > 0
        result = self.mixed_frame.dtypes
        expected = Series(dict((k, v.dtype)
                               for k, v in self.mixed_frame.iteritems()),
                          index=result.index)
        assert_series_equal(result, expected)

    def test_convert_objects(self):

        oops = self.mixed_frame.T.T
        converted = oops.convert_objects()
        assert_frame_equal(converted, self.mixed_frame)
        self.assert_(converted['A'].dtype == np.float64)

        # force numeric conversion
        self.mixed_frame['H'] = '1.'
        self.mixed_frame['I'] = '1'

        # add in some items that will be nan
        l = len(self.mixed_frame)
        self.mixed_frame['J'] = '1.'
        self.mixed_frame['K'] = '1'
        self.mixed_frame.ix[0:5,['J','K']] = 'garbled'
        converted = self.mixed_frame.convert_objects(convert_numeric=True)
        self.assert_(converted['H'].dtype == 'float64')
        self.assert_(converted['I'].dtype == 'int64')
        self.assert_(converted['J'].dtype == 'float64')
        self.assert_(converted['K'].dtype == 'float64')
        self.assert_(len(converted['J'].dropna()) == l-5)
        self.assert_(len(converted['K'].dropna()) == l-5)

        # via astype
        converted = self.mixed_frame.copy()
        converted['H'] = converted['H'].astype('float64')
        converted['I'] = converted['I'].astype('int64')
        self.assert_(converted['H'].dtype == 'float64')
        self.assert_(converted['I'].dtype == 'int64')

        # via astype, but errors
        converted = self.mixed_frame.copy()
        self.assertRaises(Exception, converted['H'].astype, 'int32')

        # mixed in a single column
        df = DataFrame(dict(s = Series([1, 'na', 3 ,4])))
        result = df.convert_objects(convert_numeric=True)
        expected = DataFrame(dict(s = Series([1, np.nan, 3 ,4])))
        assert_frame_equal(result, expected)

    def test_convert_objects_no_conversion(self):
        mixed1 = DataFrame(
            {'a': [1, 2, 3], 'b': [4.0, 5, 6], 'c': ['x', 'y', 'z']})
        mixed2 = mixed1.convert_objects()
        assert_frame_equal(mixed1, mixed2)

    def test_append_series_dict(self):
        df = DataFrame(np.random.randn(5, 4),
                       columns=['foo', 'bar', 'baz', 'qux'])

        series = df.ix[4]
        self.assertRaises(ValueError, df.append, series, verify_integrity=True)
        series.name = None
        self.assertRaises(Exception, df.append, series, verify_integrity=True)

        result = df.append(series[::-1], ignore_index=True)
        expected = df.append(DataFrame({0: series[::-1]}, index=df.columns).T,
                             ignore_index=True)
        assert_frame_equal(result, expected)

        # dict
        result = df.append(series.to_dict(), ignore_index=True)
        assert_frame_equal(result, expected)

        result = df.append(series[::-1][:3], ignore_index=True)
        expected = df.append(DataFrame({0: series[::-1][:3]}).T,
                             ignore_index=True)
        assert_frame_equal(result, expected.ix[:, result.columns])

        # can append when name set
        row = df.ix[4]
        row.name = 5
        result = df.append(row)
        expected = df.append(df[-1:], ignore_index=True)
        assert_frame_equal(result, expected)

    def test_append_list_of_series_dicts(self):
        df = DataFrame(np.random.randn(5, 4),
                       columns=['foo', 'bar', 'baz', 'qux'])

        dicts = [x.to_dict() for idx, x in df.iterrows()]

        result = df.append(dicts, ignore_index=True)
        expected = df.append(df, ignore_index=True)
        assert_frame_equal(result, expected)

        # different columns
        dicts = [{'foo': 1, 'bar': 2, 'baz': 3, 'peekaboo': 4},
                 {'foo': 5, 'bar': 6, 'baz': 7, 'peekaboo': 8}]
        result = df.append(dicts, ignore_index=True)
        expected = df.append(DataFrame(dicts), ignore_index=True)
        assert_frame_equal(result, expected)

    def test_append_empty_dataframe(self):

        # Empty df append empty df
        df1 = DataFrame([])
        df2 = DataFrame([])
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Non-empty df append empty df
        df1 = DataFrame(np.random.randn(5, 2))
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Empty df with columns append empty df
        df1 = DataFrame(columns=['bar', 'foo'])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # Non-Empty df with columns append empty df
        df1 = DataFrame(np.random.randn(5, 2), columns=['bar', 'foo'])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        assert_frame_equal(result, expected)

    def test_asfreq(self):
        offset_monthly = self.tsframe.asfreq(datetools.bmonthEnd)
        rule_monthly = self.tsframe.asfreq('BM')

        assert_almost_equal(offset_monthly['A'], rule_monthly['A'])

        filled = rule_monthly.asfreq('B', method='pad')
        # TODO: actually check that this worked.

        # don't forget!
        filled_dep = rule_monthly.asfreq('B', method='pad')

        # test does not blow up on length-0 DataFrame
        zero_length = self.tsframe.reindex([])
        result = zero_length.asfreq('BM')
        self.assert_(result is not zero_length)

    def test_asfreq_datetimeindex(self):
        df = DataFrame({'A': [1, 2, 3]},
                       index=[datetime(2011, 11, 01), datetime(2011, 11, 2),
                              datetime(2011, 11, 3)])
        df = df.asfreq('B')
        self.assert_(isinstance(df.index, DatetimeIndex))

        ts = df['A'].asfreq('B')
        self.assert_(isinstance(ts.index, DatetimeIndex))

    def test_at_time_between_time_datetimeindex(self):
        index = pan.date_range("2012-01-01", "2012-01-05", freq='30min')
        df = DataFrame(randn(len(index), 5), index=index)
        akey = time(12, 0, 0)
        bkey = slice(time(13, 0, 0), time(14, 0, 0))
        ainds = [24, 72, 120, 168]
        binds = [26, 27, 28, 74, 75, 76, 122, 123, 124, 170, 171, 172]

        result = df.at_time(akey)
        expected = df.ix[akey]
        expected2 = df.ix[ainds]
        assert_frame_equal(result, expected)
        assert_frame_equal(result, expected2)
        self.assert_(len(result) == 4)

        result = df.between_time(bkey.start, bkey.stop)
        expected = df.ix[bkey]
        expected2 = df.ix[binds]
        assert_frame_equal(result, expected)
        assert_frame_equal(result, expected2)
        self.assert_(len(result) == 12)

        result = df.copy()
        result.ix[akey] = 0
        result = result.ix[akey]
        expected = df.ix[akey].copy()
        expected.ix[:] = 0
        assert_frame_equal(result, expected)

        result = df.copy()
        result.ix[akey] = 0
        result.ix[akey] = df.ix[ainds]
        assert_frame_equal(result, df)

        result = df.copy()
        result.ix[bkey] = 0
        result = result.ix[bkey]
        expected = df.ix[bkey].copy()
        expected.ix[:] = 0
        assert_frame_equal(result, expected)

        result = df.copy()
        result.ix[bkey] = 0
        result.ix[bkey] = df.ix[binds]
        assert_frame_equal(result, df)

    def test_as_matrix(self):
        frame = self.frame
        mat = frame.as_matrix()

        frameCols = frame.columns
        for i, row in enumerate(mat):
            for j, value in enumerate(row):
                col = frameCols[j]
                if np.isnan(value):
                    self.assert_(np.isnan(frame[col][i]))
                else:
                    self.assertEqual(value, frame[col][i])

        # mixed type
        mat = self.mixed_frame.as_matrix(['foo', 'A'])
        self.assertEqual(mat[0, 0], 'bar')

        df = DataFrame({'real': [1, 2, 3], 'complex': [1j, 2j, 3j]})
        mat = df.as_matrix()
        self.assertEqual(mat[0, 0], 1j)

        # single block corner case
        mat = self.frame.as_matrix(['A', 'B'])
        expected = self.frame.reindex(columns=['A', 'B']).values
        assert_almost_equal(mat, expected)

    def test_as_matrix_duplicates(self):
        df = DataFrame([[1, 2, 'a', 'b'],
                        [1, 2, 'a', 'b']],
                       columns=['one', 'one', 'two', 'two'])

        result = df.values
        expected = np.array([[1, 2, 'a', 'b'], [1, 2, 'a', 'b']],
                            dtype=object)

        self.assertTrue(np.array_equal(result, expected))

    def test_as_blocks(self):
        frame = self.mixed_float
        mat = frame.blocks
        self.assert_(set([ x.name for x in frame.dtypes.values ]) == set(mat.keys()))

    def test_values(self):
        self.frame.values[:, 0] = 5.
        self.assert_((self.frame.values[:, 0] == 5).all())

    def test_deepcopy(self):
        cp = deepcopy(self.frame)
        series = cp['A']
        series[:] = 10
        for idx, value in series.iteritems():
            self.assertNotEqual(self.frame['A'][idx], value)

    def test_copy(self):
        cop = self.frame.copy()
        cop['E'] = cop['A']
        self.assert_('E' not in self.frame)

        # copy objects
        copy = self.mixed_frame.copy()
        self.assert_(copy._data is not self.mixed_frame._data)

    # def test_copy_index_name_checking(self):
    #     # don't want to be able to modify the index stored elsewhere after
    #     # making a copy

    #     self.frame.columns.name = None
    #     cp = self.frame.copy()
    #     cp.columns.name = 'foo'

    #     self.assert_(self.frame.columns.name is None)

    def _check_method(self, method='pearson', check_minp=False):
        if not check_minp:
            correls = self.frame.corr(method=method)
            exp = self.frame['A'].corr(self.frame['C'], method=method)
            assert_almost_equal(correls['A']['C'], exp)
        else:
            result = self.frame.corr(min_periods=len(self.frame) - 8)
            expected = self.frame.corr()
            expected.ix['A', 'B'] = expected.ix['B', 'A'] = nan

    def test_corr_pearson(self):
        _skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('pearson')

    def test_corr_kendall(self):
        _skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('kendall')

    def test_corr_spearman(self):
        _skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        self._check_method('spearman')

    def test_corr_non_numeric(self):
        _skip_if_no_scipy()
        self.frame['A'][:5] = nan
        self.frame['B'][5:10] = nan

        # exclude non-numeric types
        result = self.mixed_frame.corr()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].corr()
        assert_frame_equal(result, expected)

    def test_corr_nooverlap(self):
        _skip_if_no_scipy()

        # nothing in common
        for meth in ['pearson', 'kendall', 'spearman']:
            df = DataFrame({'A': [1, 1.5, 1, np.nan, np.nan, np.nan],
                            'B': [np.nan, np.nan, np.nan, 1, 1.5, 1]})
            rs = df.corr(meth)
            self.assert_(isnull(rs.ix['A', 'B']))
            self.assert_(isnull(rs.ix['B', 'A']))
            self.assert_(rs.ix['A', 'A'] == 1)
            self.assert_(rs.ix['B', 'B'] == 1)

    def test_corr_constant(self):
        _skip_if_no_scipy()

        # constant --> all NA

        for meth in ['pearson', 'spearman']:
            df = DataFrame({'A': [1, 1, 1, np.nan, np.nan, np.nan],
                            'B': [np.nan, np.nan, np.nan, 1, 1, 1]})
            rs = df.corr(meth)
            self.assert_(isnull(rs.values).all())

    def test_corr_int(self):
        # dtypes other than float64 #1761
        df3 = DataFrame({"a": [1, 2, 3, 4], "b": [1, 2, 3, 4]})

        # it works!
        df3.cov()
        df3.corr()

    def test_cov(self):
        # min_periods no NAs (corner case)
        expected = self.frame.cov()
        result = self.frame.cov(min_periods=len(self.frame))

        assert_frame_equal(expected, result)

        result = self.frame.cov(min_periods=len(self.frame) + 1)
        self.assert_(isnull(result.values).all())

        # with NAs
        frame = self.frame.copy()
        frame['A'][:5] = nan
        frame['B'][5:10] = nan
        result = self.frame.cov(min_periods=len(self.frame) - 8)
        expected = self.frame.cov()
        expected.ix['A', 'B'] = np.nan
        expected.ix['B', 'A'] = np.nan

        # regular
        self.frame['A'][:5] = nan
        self.frame['B'][:10] = nan
        cov = self.frame.cov()

        assert_almost_equal(cov['A']['C'],
                            self.frame['A'].cov(self.frame['C']))

        # exclude non-numeric types
        result = self.mixed_frame.cov()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].cov()
        assert_frame_equal(result, expected)

        # Single column frame
        df = DataFrame(np.linspace(0.0,1.0,10))
        result = df.cov()
        expected = DataFrame(np.cov(df.values.T).reshape((1,1)),
                             index=df.columns,columns=df.columns)
        assert_frame_equal(result, expected)
        df.ix[0] = np.nan
        result = df.cov()
        expected = DataFrame(np.cov(df.values[1:].T).reshape((1,1)),
                             index=df.columns,columns=df.columns)
        assert_frame_equal(result, expected)

    def test_corrwith(self):
        a = self.tsframe
        noise = Series(randn(len(a)), index=a.index)

        b = self.tsframe + noise

        # make sure order does not matter
        b = b.reindex(columns=b.columns[::-1], index=b.index[::-1][10:])
        del b['B']

        colcorr = a.corrwith(b, axis=0)
        assert_almost_equal(colcorr['A'], a['A'].corr(b['A']))

        rowcorr = a.corrwith(b, axis=1)
        assert_series_equal(rowcorr, a.T.corrwith(b.T, axis=0))

        dropped = a.corrwith(b, axis=0, drop=True)
        assert_almost_equal(dropped['A'], a['A'].corr(b['A']))
        self.assert_('B' not in dropped)

        dropped = a.corrwith(b, axis=1, drop=True)
        self.assert_(a.index[-1] not in dropped.index)

        # non time-series data
        index = ['a', 'b', 'c', 'd', 'e']
        columns = ['one', 'two', 'three', 'four']
        df1 = DataFrame(randn(5, 4), index=index, columns=columns)
        df2 = DataFrame(randn(4, 4), index=index[:4], columns=columns)
        correls = df1.corrwith(df2, axis=1)
        for row in index[:4]:
            assert_almost_equal(correls[row], df1.ix[row].corr(df2.ix[row]))

    def test_corrwith_with_objects(self):
        df1 = tm.makeTimeDataFrame()
        df2 = tm.makeTimeDataFrame()
        cols = ['A', 'B', 'C', 'D']

        df1['obj'] = 'foo'
        df2['obj'] = 'bar'

        result = df1.corrwith(df2)
        expected = df1.ix[:, cols].corrwith(df2.ix[:, cols])
        assert_series_equal(result, expected)

        result = df1.corrwith(df2, axis=1)
        expected = df1.ix[:, cols].corrwith(df2.ix[:, cols], axis=1)
        assert_series_equal(result, expected)

    def test_corrwith_series(self):
        result = self.tsframe.corrwith(self.tsframe['A'])
        expected = self.tsframe.apply(self.tsframe['A'].corr)

        assert_series_equal(result, expected)

    def test_drop_names(self):
        df = DataFrame([[1, 2, 3],[3, 4, 5],[5, 6, 7]], index=['a', 'b', 'c'], columns=['d', 'e', 'f'])
        df.index.name, df.columns.name = 'first', 'second'
        df_dropped_b = df.drop('b')
        df_dropped_e = df.drop('e', axis=1)
        self.assert_(df_dropped_b.index.name == 'first')
        self.assert_(df_dropped_e.index.name == 'first')
        self.assert_(df_dropped_b.columns.name == 'second')
        self.assert_(df_dropped_e.columns.name == 'second')

    def test_dropEmptyRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)

        smaller_frame = frame.dropna(how='all')
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

        smaller_frame = frame.dropna(how='all', subset=['foo'])
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

    def test_dropIncompleteRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo': mat}, index=self.frame.index)
        frame['bar'] = 5

        smaller_frame = frame.dropna()
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

        samesize_frame = frame.dropna(subset=['bar'])
        self.assert_(samesize_frame.index.equals(self.frame.index))

    def test_dropna(self):
        df = DataFrame(np.random.randn(6, 4))
        df[2][:2] = nan

        dropped = df.dropna(axis=1)
        expected = df.ix[:, [0, 1, 3]]
        assert_frame_equal(dropped, expected)

        dropped = df.dropna(axis=0)
        expected = df.ix[range(2, 6)]
        assert_frame_equal(dropped, expected)

        # threshold
        dropped = df.dropna(axis=1, thresh=5)
        expected = df.ix[:, [0, 1, 3]]
        assert_frame_equal(dropped, expected)

        dropped = df.dropna(axis=0, thresh=4)
        expected = df.ix[range(2, 6)]
        assert_frame_equal(dropped, expected)

        dropped = df.dropna(axis=1, thresh=4)
        assert_frame_equal(dropped, df)

        dropped = df.dropna(axis=1, thresh=3)
        assert_frame_equal(dropped, df)

        # subset
        dropped = df.dropna(axis=0, subset=[0, 1, 3])
        assert_frame_equal(dropped, df)

        # all
        dropped = df.dropna(axis=1, how='all')
        assert_frame_equal(dropped, df)

        df[2] = nan
        dropped = df.dropna(axis=1, how='all')
        expected = df.ix[:, [0, 1, 3]]
        assert_frame_equal(dropped, expected)

    def test_dropna_corner(self):
        # bad input
        self.assertRaises(ValueError, self.frame.dropna, how='foo')
        self.assertRaises(ValueError, self.frame.dropna, how=None)

    def test_dropna_multiple_axes(self):
        df = DataFrame([[1, np.nan, 2, 3],
                        [4, np.nan, 5, 6],
                        [np.nan, np.nan, np.nan, np.nan],
                        [7, np.nan, 8, 9]])

        result = df.dropna(how='all', axis=[0, 1])
        result2 = df.dropna(how='all', axis=(0, 1))
        expected = df.dropna(how='all').dropna(how='all', axis=1)

        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

    def test_drop_duplicates(self):
        df = DataFrame({'AAA': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': range(8)})

        # single column
        result = df.drop_duplicates('AAA')
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('AAA', take_last=True)
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        expected = df.ix[[0, 1, 2, 3]]
        result = df.drop_duplicates(np.array(['AAA', 'B']))
        assert_frame_equal(result, expected)
        result = df.drop_duplicates(['AAA', 'B'])
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AAA', 'B'), take_last=True)
        expected = df.ix[[0, 5, 6, 7]]
        assert_frame_equal(result, expected)

        # consider everything
        df2 = df.ix[:, ['AAA', 'B', 'C']]

        result = df2.drop_duplicates()
        # in this case only
        expected = df2.drop_duplicates(['AAA', 'B'])
        assert_frame_equal(result, expected)

        result = df2.drop_duplicates(take_last=True)
        expected = df2.drop_duplicates(['AAA', 'B'], take_last=True)
        assert_frame_equal(result, expected)

    def test_drop_duplicates_tuple(self):
        df = DataFrame({('AA', 'AB'): ['foo', 'bar', 'foo', 'bar',
                                       'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1, 1, 2, 2, 2, 2, 1, 2],
                        'D': range(8)})

        # single column
        result = df.drop_duplicates(('AA', 'AB'))
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(('AA', 'AB'), take_last=True)
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        expected = df.ix[[0, 1, 2, 3]]
        result = df.drop_duplicates((('AA', 'AB'), 'B'))
        assert_frame_equal(result, expected)

    def test_drop_duplicates_NA(self):
        # none
        df = DataFrame({'A': [None, None, 'foo', 'bar',
                              'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                        'D': range(8)})

        # single column
        result = df.drop_duplicates('A')
        expected = df.ix[[0, 2, 3]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', take_last=True)
        expected = df.ix[[1, 6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        result = df.drop_duplicates(['A', 'B'])
        expected = df.ix[[0, 2, 3, 6]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['A', 'B'], take_last=True)
        expected = df.ix[[1, 5, 6, 7]]
        assert_frame_equal(result, expected)

        # nan
        df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                              'foo', 'bar', 'bar', 'foo'],
                        'B': ['one', 'one', 'two', 'two',
                              'two', 'two', 'one', 'two'],
                        'C': [1.0, np.nan, np.nan, np.nan, 1., 1., 1, 1.],
                        'D': range(8)})

        # single column
        result = df.drop_duplicates('C')
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('C', take_last=True)
        expected = df.ix[[3, 7]]
        assert_frame_equal(result, expected)

        # multi column
        result = df.drop_duplicates(['C', 'B'])
        expected = df.ix[[0, 1, 2, 4]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['C', 'B'], take_last=True)
        expected = df.ix[[1, 3, 6, 7]]
        assert_frame_equal(result, expected)

    def test_drop_duplicates_inplace(self):
        orig = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                                'foo', 'bar', 'bar', 'foo'],
                          'B': ['one', 'one', 'two', 'two',
                                'two', 'two', 'one', 'two'],
                          'C': [1, 1, 2, 2, 2, 2, 1, 2],
                          'D': range(8)})

        # single column
        df = orig.copy()
        df.drop_duplicates('A', inplace=True)
        expected = orig[:2]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates('A', take_last=True, inplace=True)
        expected = orig.ix[[6, 7]]
        result = df
        assert_frame_equal(result, expected)

        # multi column
        df = orig.copy()
        df.drop_duplicates(['A', 'B'], inplace=True)
        expected = orig.ix[[0, 1, 2, 3]]
        result = df
        assert_frame_equal(result, expected)

        df = orig.copy()
        df.drop_duplicates(['A', 'B'], take_last=True, inplace=True)
        expected = orig.ix[[0, 5, 6, 7]]
        result = df
        assert_frame_equal(result, expected)

        # consider everything
        orig2 = orig.ix[:, ['A', 'B', 'C']].copy()

        df2 = orig2.copy()
        df2.drop_duplicates(inplace=True)
        # in this case only
        expected = orig2.drop_duplicates(['A', 'B'])
        result = df2
        assert_frame_equal(result, expected)

        df2 = orig2.copy()
        df2.drop_duplicates(take_last=True, inplace=True)
        expected = orig2.drop_duplicates(['A', 'B'], take_last=True)
        result = df2
        assert_frame_equal(result, expected)

    def test_drop_col_still_multiindex(self):
        arrays = [['a', 'b', 'c', 'top'],
                  ['', '', '', 'OD'],
                  ['', '', '', 'wx']]

        tuples = zip(*arrays)
        tuples.sort()
        index = MultiIndex.from_tuples(tuples)

        df = DataFrame(randn(3, 4), columns=index)
        del df[('a', '', '')]
        assert(isinstance(df.columns, MultiIndex))

    def test_fillna(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        zero_filled = self.tsframe.fillna(0)
        self.assert_((zero_filled['A'][:5] == 0).all())

        padded = self.tsframe.fillna(method='pad')
        self.assert_(np.isnan(padded['A'][:5]).all())
        self.assert_((padded['A'][-5:] == padded['A'][-5]).all())

        # mixed type
        self.mixed_frame['foo'][5:20] = nan
        self.mixed_frame['A'][-10:] = nan
        result = self.mixed_frame.fillna(value=0)
        result = self.mixed_frame.fillna(method='pad')

        self.assertRaises(ValueError, self.tsframe.fillna)
        self.assertRaises(ValueError, self.tsframe.fillna, 5, method='ffill')

        # mixed numeric (but no float16)
        mf = self.mixed_float.reindex(columns=['A','B','D'])
        mf['A'][-10:] = nan
        result = mf.fillna(value=0)
        _check_mixed_float(result, dtype = dict(C = None))

        result = mf.fillna(method='pad')
        _check_mixed_float(result, dtype = dict(C = None))

        # empty frame (GH #2778)
        df = DataFrame(columns=['x'])
        for m in ['pad','backfill']:
            df.x.fillna(method=m,inplace=1)
            df.x.fillna(method=m)

    def test_ffill(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        assert_frame_equal(self.tsframe.ffill(),
                           self.tsframe.fillna(method='ffill'))

    def test_bfill(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        assert_frame_equal(self.tsframe.bfill(),
                           self.tsframe.fillna(method='bfill'))

    def test_fillna_skip_certain_blocks(self):
        # don't try to fill boolean, int blocks

        df = DataFrame(np.random.randn(10, 4).astype(int))

        # it works!
        df.fillna(np.nan)

    def test_fillna_inplace(self):
        df = DataFrame(np.random.randn(10, 4))
        df[1][:4] = np.nan
        df[3][-4:] = np.nan

        expected = df.fillna(value=0)
        self.assert_(expected is not df)

        df.fillna(value=0, inplace=True)
        assert_frame_equal(df, expected)

        df[1][:4] = np.nan
        df[3][-4:] = np.nan
        expected = df.fillna(method='ffill')
        self.assert_(expected is not df)

        df.fillna(method='ffill', inplace=True)
        assert_frame_equal(df, expected)

    def test_fillna_dict_series(self):
        df = DataFrame({'a': [nan, 1, 2, nan, nan],
                        'b': [1, 2, 3, nan, nan],
                        'c': [nan, 1, 2, 3, 4]})

        result = df.fillna({'a': 0, 'b': 5})

        expected = df.copy()
        expected['a'] = expected['a'].fillna(0)
        expected['b'] = expected['b'].fillna(5)
        assert_frame_equal(result, expected)

        # it works
        result = df.fillna({'a': 0, 'b': 5, 'd': 7})

        # Series treated same as dict
        result = df.fillna(df.max())
        expected = df.fillna(df.max().to_dict())
        assert_frame_equal(result, expected)

        # disable this for now
        self.assertRaises(Exception, df.fillna, df.max(1), axis=1)

    def test_fillna_columns(self):
        df = DataFrame(np.random.randn(10, 10))
        df.values[:, ::2] = np.nan

        result = df.fillna(method='ffill', axis=1)
        expected = df.T.fillna(method='pad').T
        assert_frame_equal(result, expected)

        df.insert(6, 'foo', 5)
        result = df.fillna(method='ffill', axis=1)
        expected = df.astype(float).fillna(method='ffill', axis=1)
        assert_frame_equal(result, expected)

    def test_fillna_invalid_method(self):
        try:
            self.frame.fillna(method='ffil')
        except ValueError, inst:
            self.assert_('ffil' in str(inst))

    def test_fillna_invalid_value(self):
        # list
        self.assertRaises(TypeError, self.frame.fillna, [1, 2])
        # tuple
        self.assertRaises(TypeError, self.frame.fillna, (1, 2))

    def test_replace_inplace(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        tsframe = self.tsframe.copy()
        res = tsframe.replace(nan, 0, inplace=True)
        assert_frame_equal(tsframe, self.tsframe.fillna(0))

        self.assertRaises(TypeError, self.tsframe.replace, nan, inplace=True)

        # mixed type
        self.mixed_frame['foo'][5:20] = nan
        self.mixed_frame['A'][-10:] = nan

        result = self.mixed_frame.replace(np.nan, 0)
        expected = self.mixed_frame.fillna(value=0)
        assert_frame_equal(result, expected)

        tsframe = self.tsframe.copy()
        tsframe.replace([nan], [0], inplace=True)
        assert_frame_equal(tsframe, self.tsframe.fillna(0))

    def test_regex_replace_scalar(self):
        obj = {'a': list('ab..'), 'b': list('efgh')}
        dfobj = DataFrame(obj)
        mix = {'a': range(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)

        ### simplest cases
        ## regex -> value
        # obj frame
        res = dfobj.replace(r'\s*\.\s*', nan, regex=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.replace(r'\s*\.\s*', nan, regex=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.replace(re.compile(r'\s*\.\s*'), nan, regex=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.replace(re.compile(r'\s*\.\s*'), nan, regex=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1')
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1')
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        res = dfmix.replace(regex=re.compile(r'\s*(\.)\s*'), value=r'\1\1\1')
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        res = dfmix.replace(regex=r'\s*(\.)\s*', value=r'\1\1\1')
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

    def test_regex_replace_scalar_inplace(self):
        obj = {'a': list('ab..'), 'b': list('efgh')}
        dfobj = DataFrame(obj)
        mix = {'a': range(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)

        ### simplest cases
        ## regex -> value
        # obj frame
        res = dfobj.copy()
        res.replace(r'\s*\.\s*', nan, regex=True, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(r'\s*\.\s*', nan, regex=True, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True, inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(r'\s*(\.)\s*', r'\1\1\1', regex=True, inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.copy()
        res.replace(re.compile(r'\s*\.\s*'), nan, regex=True, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(re.compile(r'\s*\.\s*'), nan, regex=True, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1', regex=True,
                    inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(re.compile(r'\s*(\.)\s*'), r'\1\1\1', regex=True,
                    inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        res = dfobj.copy()
        res.replace(regex=r'\s*\.\s*', value=nan, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(regex=r'\s*\.\s*', value=nan, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(regex=r'\s*(\.)\s*', value=r'\1\1\1', inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(regex=r'\s*(\.)\s*', value=r'\1\1\1', inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

        # everything with compiled regexs as well
        res = dfobj.copy()
        res.replace(regex=re.compile(r'\s*\.\s*'), value=nan, inplace=True)
        assert_frame_equal(dfobj, res.fillna('.'))

        # mixed
        res = dfmix.copy()
        res.replace(regex=re.compile(r'\s*\.\s*'), value=nan, inplace=True)
        assert_frame_equal(dfmix, res.fillna('.'))

        ## regex -> regex
        # obj frame
        res = dfobj.copy()
        res.replace(regex=re.compile(r'\s*(\.)\s*'), value=r'\1\1\1',
                    inplace=True)
        objc = obj.copy()
        objc['a'] = ['a', 'b', '...', '...']
        expec = DataFrame(objc)
        assert_frame_equal(res, expec)

        # with mixed
        res = dfmix.copy()
        res.replace(regex=re.compile(r'\s*(\.)\s*'), value=r'\1\1\1',
                    inplace=True)
        mixc = mix.copy()
        mixc['b'] = ['a', 'b', '...', '...']
        expec = DataFrame(mixc)
        assert_frame_equal(res, expec)

    def test_regex_replace_list_obj(self):
        obj = {'a': list('ab..'), 'b': list('efgh'), 'c': list('helo')}
        dfobj = DataFrame(obj)

        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'e|f|g']
        values = [nan, 'crap']
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': ['a', 'b', nan, nan], 'b': ['crap'] * 3 +
                           ['h'], 'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(e|f|g)']
        values = [r'\1\1', r'\1_crap']
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['e_crap',
                                                              'f_crap',
                                                              'g_crap', 'h'],
                           'c': ['h', 'e_crap', 'l', 'o']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.replace(value=values, regex=to_replace_res)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

    def test_regex_replace_list_obj_inplace(self):
        ### same as above with inplace=True
        ## lists of regexes and values
        obj = {'a': list('ab..'), 'b': list('efgh'), 'c': list('helo')}
        dfobj = DataFrame(obj)

        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'e|f|g']
        values = [nan, 'crap']
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': ['a', 'b', nan, nan], 'b': ['crap'] * 3 +
                           ['h'], 'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(e|f|g)']
        values = [r'\1\1', r'\1_crap']
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['e_crap',
                                                              'f_crap',
                                                              'g_crap', 'h'],
                           'c': ['h', 'e_crap', 'l', 'o']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'e']
        values = [r'\1\1', r'crap']
        res = dfobj.copy()
        res.replace(value=values, regex=to_replace_res, inplace=True)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['crap', 'f', 'g',
                                                              'h'],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

    def test_regex_replace_list_mixed(self):
        ## mixed frame to make sure this doesn't break things
        mix = {'a': range(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)

        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'a']
        values = [nan, 'crap']
        mix2 = {'a': range(4), 'b': list('ab..'), 'c': list('halo')}
        dfmix2 = DataFrame(mix2)
        res = dfmix2.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': mix2['a'], 'b': ['crap', 'b', nan, nan],
                           'c': ['h', 'crap', 'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(a|b)']
        values = [r'\1\1', r'\1_crap']
        res = dfmix.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a_crap', 'b_crap', '..',
                                                '..']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.replace(to_replace_res, values, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.replace(regex=to_replace_res, value=values)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

    def test_regex_replace_list_mixed_inplace(self):
        mix = {'a': range(4), 'b': list('ab..')}
        dfmix = DataFrame(mix)
        # the same inplace
        ## lists of regexes and values
        # list of [re1, re2, ..., reN] -> [v1, v2, ..., vN]
        to_replace_res = [r'\s*\.\s*', r'a']
        values = [nan, 'crap']
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b', nan, nan]})
        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [re1, re2, .., reN]
        to_replace_res = [r'\s*(\.)\s*', r'(a|b)']
        values = [r'\1\1', r'\1_crap']
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a_crap', 'b_crap', '..',
                                                '..']})

        assert_frame_equal(res, expec)

        # list of [re1, re2, ..., reN] -> [(re1 or v1), (re2 or v2), ..., (reN
        # or vN)]
        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.copy()
        res.replace(to_replace_res, values, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

        to_replace_res = [r'\s*(\.)\s*', r'a', r'(b)']
        values = [r'\1\1', r'crap', r'\1_crap']
        res = dfmix.copy()
        res.replace(regex=to_replace_res, value=values, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['crap', 'b_crap', '..', '..']})
        assert_frame_equal(res, expec)

    def test_regex_replace_dict_mixed(self):
        mix = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        dfmix = DataFrame(mix)

        ## dicts
        # single dict {re1: v1}, search the whole frame
        # need test for this...

        # list of dicts {re1: v1, re2: v2, ..., re3: v3}, search the whole
        # frame
        res = dfmix.replace({'b': r'\s*\.\s*'}, {'b': nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace({'b': r'\s*\.\s*'}, {'b': nan}, inplace=True, regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', nan, nan], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        # list of dicts {re1: re11, re2: re12, ..., reN: re1N}, search the
        # whole frame
        res = dfmix.replace({'b': r'\s*(\.)\s*'}, {'b': r'\1ty'}, regex=True)
        res2 = dfmix.copy()
        res2.replace({'b': r'\s*(\.)\s*'}, {'b': r'\1ty'}, inplace=True,
                     regex=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', '.ty', '.ty'], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        res = dfmix.replace(regex={'b': r'\s*(\.)\s*'}, value={'b': r'\1ty'})
        res2 = dfmix.copy()
        res2.replace(regex={'b': r'\s*(\.)\s*'}, value={'b': r'\1ty'},
                     inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', '.ty', '.ty'], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        # scalar -> dict
        # to_replace regex, {value: value}
        res = dfmix.replace('a', {'b': nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace('a', {'b': nan}, regex=True, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': [nan, 'b', '.', '.'], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

        res = dfmix.replace('a', {'b': nan}, regex=True)
        res2 = dfmix.copy()
        res2.replace(regex='a', value={'b': nan}, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': [nan, 'b', '.', '.'], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)

    def test_regex_replace_dict_nested(self):
        # nested dicts will not work until this is implemented for Series
        mix = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        dfmix = DataFrame(mix)
        res = dfmix.replace({'b': {r'\s*\.\s*': nan}}, regex=True)
        res2 = dfmix.copy()
        res4 = dfmix.copy()
        res2.replace({'b': {r'\s*\.\s*': nan}}, inplace=True, regex=True)
        res3 = dfmix.replace(regex={'b': {r'\s*\.\s*': nan}})
        res4.replace(regex={'b': {r'\s*\.\s*': nan}}, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', nan, nan], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)
        assert_frame_equal(res4, expec)

    def test_regex_replace_dict_nested_gh4115(self):
        df = pd.DataFrame({'Type':['Q','T','Q','Q','T'], 'tmp':2})
        expected = DataFrame({'Type': [0,1,0,0,1], 'tmp': 2})
        assert_frame_equal(df.replace({'Type': {'Q':0,'T':1}}), expected)

    def test_regex_replace_list_to_scalar(self):
        mix = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace([r'\s*\.\s*', 'a|b'], nan, regex=True)
        res2 = df.copy()
        res3 = df.copy()
        res2.replace([r'\s*\.\s*', 'a|b'], nan, regex=True, inplace=True)
        res3.replace(regex=[r'\s*\.\s*', 'a|b'], value=nan, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': np.array([nan] * 4),
                           'c': [nan, nan, nan, 'd']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_str_to_numeric(self):
        # what happens when you try to replace a numeric value with a regex?
        mix = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace(r'\s*\.\s*', 0, regex=True)
        res2 = df.copy()
        res2.replace(r'\s*\.\s*', 0, inplace=True, regex=True)
        res3 = df.copy()
        res3.replace(regex=r'\s*\.\s*', value=0, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', 0, 0], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_regex_list_to_numeric(self):
        mix = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace([r'\s*\.\s*', 'b'], 0, regex=True)
        res2 = df.copy()
        res2.replace([r'\s*\.\s*', 'b'], 0, regex=True, inplace=True)
        res3 = df.copy()
        res3.replace(regex=[r'\s*\.\s*', 'b'], value=0, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 0, 0, 0], 'c': ['a', 0,
                                                                     nan,
                                                                     'd']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_series_of_regexes(self):
        mix = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        s1 = Series({'b': r'\s*\.\s*'})
        s2 = Series({'b': nan})
        res = df.replace(s1, s2, regex=True)
        res2 = df.copy()
        res2.replace(s1, s2, inplace=True, regex=True)
        res3 = df.copy()
        res3.replace(regex=s1, value=s2, inplace=True)
        expec = DataFrame({'a': mix['a'], 'b': ['a', 'b', nan, nan], 'c':
                           mix['c']})
        assert_frame_equal(res, expec)
        assert_frame_equal(res2, expec)
        assert_frame_equal(res3, expec)

    def test_regex_replace_numeric_to_object_conversion(self):
        mix = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
        df = DataFrame(mix)
        res = df.replace(0, 'a')
        expec = DataFrame({'a': ['a', 1, 2, 3], 'b': mix['b'], 'c': mix['c']})
        assert_frame_equal(res, expec)
        self.assertEqual(res.a.dtype, np.object_)

    def test_replace(self):
        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan

        zero_filled = self.tsframe.replace(nan, -1e8)
        assert_frame_equal(zero_filled, self.tsframe.fillna(-1e8))
        assert_frame_equal(zero_filled.replace(-1e8, nan), self.tsframe)

        self.tsframe['A'][:5] = nan
        self.tsframe['A'][-5:] = nan
        self.tsframe['B'][:5] = -1e8

        # empty
        df = DataFrame(index=['a', 'b'])
        assert_frame_equal(df, df.replace(5, 7))

    def test_replace_list(self):
        obj = {'a': list('ab..'), 'b': list('efgh'), 'c': list('helo')}
        dfobj = DataFrame(obj)

        ## lists of regexes and values
        # list of [v1, v2, ..., vN] -> [v1, v2, ..., vN]
        to_replace_res = [r'.', r'e']
        values = [nan, 'crap']
        res = dfobj.replace(to_replace_res, values)
        expec = DataFrame({'a': ['a', 'b', nan, nan],
                           'b': ['crap', 'f', 'g', 'h'], 'c': ['h', 'crap',
                                                               'l', 'o']})
        assert_frame_equal(res, expec)

        # list of [v1, v2, ..., vN] -> [v1, v2, .., vN]
        to_replace_res = [r'.', r'f']
        values = [r'..', r'crap']
        res = dfobj.replace(to_replace_res, values)
        expec = DataFrame({'a': ['a', 'b', '..', '..'], 'b': ['e', 'crap', 'g',
                                                              'h'],
                           'c': ['h', 'e', 'l', 'o']})

        assert_frame_equal(res, expec)

    def test_replace_series_dict(self):
        # from GH 3064
        df = DataFrame({'zero': {'a': 0.0, 'b': 1}, 'one': {'a': 2.0, 'b': 0}})
        result = df.replace(0, {'zero': 0.5, 'one': 1.0})
        expected = DataFrame({'zero': {'a': 0.5, 'b': 1}, 'one': {'a': 2.0, 'b': 1.0}})
        assert_frame_equal(result, expected)

        result = df.replace(0, df.mean())
        assert_frame_equal(result, expected)

        # series to series/dict
        df = DataFrame({'zero': {'a': 0.0, 'b': 1}, 'one': {'a': 2.0, 'b': 0}})
        s = Series({'zero': 0.0, 'one': 2.0})
        result = df.replace(s, {'zero': 0.5, 'one': 1.0})
        expected = DataFrame({'zero': {'a': 0.5, 'b': 1}, 'one': {'a': 1.0, 'b': 0.0}})
        assert_frame_equal(result, expected)

        result = df.replace(s, df.mean())
        assert_frame_equal(result, expected)

    def test_replace_convert(self):
        # gh 3907
        df = DataFrame([['foo', 'bar', 'bah'], ['bar', 'foo', 'bah']])
        m = {'foo': 1, 'bar': 2, 'bah': 3}
        rep = df.replace(m)
        expec = Series([ np.int64] * 3)
        res = rep.dtypes
        assert_series_equal(expec, res)

    def test_replace_mixed(self):
        self.mixed_frame['foo'][5:20] = nan
        self.mixed_frame['A'][-10:] = nan

        result = self.mixed_frame.replace(np.nan, -18)
        expected = self.mixed_frame.fillna(value=-18)
        assert_frame_equal(result, expected)
        assert_frame_equal(result.replace(-18, nan), self.mixed_frame)

        result = self.mixed_frame.replace(np.nan, -1e8)
        expected = self.mixed_frame.fillna(value=-1e8)
        assert_frame_equal(result, expected)
        assert_frame_equal(result.replace(-1e8, nan), self.mixed_frame)

        # int block upcasting
        df = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0,1],dtype='int64') })
        expected = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0.5,1],dtype='float64') })
        result = df.replace(0, 0.5)
        assert_frame_equal(result,expected)

        df.replace(0, 0.5, inplace=True)
        assert_frame_equal(df,expected)

        # int block splitting
        df = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0,1],dtype='int64'), 'C' : Series([1,2],dtype='int64') })
        expected = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0.5,1],dtype='float64'), 'C' : Series([1,2],dtype='int64') })
        result = df.replace(0, 0.5)
        assert_frame_equal(result,expected)

        # to object block upcasting
        df = DataFrame({ 'A' : Series([1.0,2.0],dtype='float64'), 'B' : Series([0,1],dtype='int64') })
        expected = DataFrame({ 'A' : Series([1,'foo'],dtype='object'), 'B' : Series([0,1],dtype='int64') })
        result = df.replace(2, 'foo')
        assert_frame_equal(result,expected)

        expected = DataFrame({ 'A' : Series(['foo','bar'],dtype='object'), 'B' : Series([0,'foo'],dtype='object') })
        result = df.replace([1,2], ['foo','bar'])
        assert_frame_equal(result,expected)

        # test case from
        from pandas.util.testing import makeCustomDataframe as mkdf
        df = DataFrame({'A' : Series([3,0],dtype='int64'), 'B' : Series([0,3],dtype='int64') })
        result = df.replace(3, df.mean().to_dict())
        expected = df.copy().astype('float64')
        m = df.mean()
        expected.iloc[0,0] = m[0]
        expected.iloc[1,1] = m[1]
        assert_frame_equal(result,expected)

    def test_interpolate(self):
        pass

    def test_replace_value_is_none(self):
        self.assertRaises(TypeError, self.tsframe.replace, nan)
        orig_value = self.tsframe.iloc[0, 0]
        orig2 = self.tsframe.iloc[1, 0]

        self.tsframe.iloc[0, 0] = nan
        self.tsframe.iloc[1, 0] = 1

        result = self.tsframe.replace(to_replace={nan: 0})
        expected = self.tsframe.T.replace(to_replace={nan: 0}).T
        assert_frame_equal(result, expected)

        result = self.tsframe.replace(to_replace={nan: 0, 1: -1e8})
        tsframe = self.tsframe.copy()
        tsframe.iloc[0, 0] = 0
        tsframe.iloc[1, 0] = -1e8
        expected = tsframe
        assert_frame_equal(expected, result)
        self.tsframe.iloc[0, 0] = orig_value
        self.tsframe.iloc[1, 0] = orig2

    def test_replace_for_new_dtypes(self):

        # dtypes
        tsframe = self.tsframe.copy().astype(np.float32)
        tsframe['A'][:5] = nan
        tsframe['A'][-5:] = nan

        zero_filled = tsframe.replace(nan, -1e8)
        assert_frame_equal(zero_filled, tsframe.fillna(-1e8))
        assert_frame_equal(zero_filled.replace(-1e8, nan), tsframe)

        tsframe['A'][:5] = nan
        tsframe['A'][-5:] = nan
        tsframe['B'][:5] = -1e8

        b = tsframe['B']
        b[b == -1e8] = nan
        tsframe['B'] = b
        result = tsframe.fillna(method='bfill')
        assert_frame_equal(result, tsframe.fillna(method='bfill'))

    def test_replace_dtypes(self):
        # int
        df = DataFrame({'ints': [1, 2, 3]})
        result = df.replace(1, 0)
        expected = DataFrame({'ints': [0, 2, 3]})
        assert_frame_equal(result, expected)

        df = DataFrame({'ints': [1, 2, 3]}, dtype=np.int32)
        result = df.replace(1, 0)
        expected = DataFrame({'ints': [0, 2, 3]}, dtype=np.int32)
        assert_frame_equal(result, expected)

        df = DataFrame({'ints': [1, 2, 3]}, dtype=np.int16)
        result = df.replace(1, 0)
        expected = DataFrame({'ints': [0, 2, 3]}, dtype=np.int16)
        assert_frame_equal(result, expected)

        # bools
        df = DataFrame({'bools': [True, False, True]})
        result = df.replace(False, True)
        self.assert_(result.values.all())

        # complex blocks
        df = DataFrame({'complex': [1j, 2j, 3j]})
        result = df.replace(1j, 0j)
        expected = DataFrame({'complex': [0j, 2j, 3j]})
        assert_frame_equal(result, expected)

        # datetime blocks
        prev = datetime.today()
        now = datetime.today()
        df = DataFrame({'datetime64': Index([prev, now, prev])})
        result = df.replace(prev, now)
        expected = DataFrame({'datetime64': Index([now] * 3)})
        assert_frame_equal(result, expected)

    def test_replace_input_formats(self):
        # both dicts
        to_rep = {'A': np.nan, 'B': 0, 'C': ''}
        values = {'A': 0, 'B': -1, 'C': 'missing'}
        df = DataFrame({'A': [np.nan, 0, np.inf], 'B': [0, 2, 5],
                        'C': ['', 'asdf', 'fd']})
        filled = df.replace(to_rep, values)
        expected = {}
        for k, v in df.iteritems():
            expected[k] = v.replace(to_rep[k], values[k])
        assert_frame_equal(filled, DataFrame(expected))

        result = df.replace([0, 2, 5], [5, 2, 0])
        expected = DataFrame({'A': [np.nan, 5, np.inf], 'B': [5, 2, 0],
                              'C': ['', 'asdf', 'fd']})
        assert_frame_equal(result, expected)

        # dict to scalar
        filled = df.replace(to_rep, 0)
        expected = {}
        for k, v in df.iteritems():
            expected[k] = v.replace(to_rep[k], 0)
        assert_frame_equal(filled, DataFrame(expected))

        self.assertRaises(TypeError, df.replace, to_rep, [np.nan, 0, ''])

        # scalar to dict
        values = {'A': 0, 'B': -1, 'C': 'missing'}
        df = DataFrame({'A': [np.nan, 0, np.nan], 'B': [0, 2, 5],
                        'C': ['', 'asdf', 'fd']})
        filled = df.replace(np.nan, values)
        expected = {}
        for k, v in df.iteritems():
            expected[k] = v.replace(np.nan, values[k])
        assert_frame_equal(filled, DataFrame(expected))

        # list to list
        to_rep = [np.nan, 0, '']
        values = [-2, -1, 'missing']
        result = df.replace(to_rep, values)
        expected = df.copy()
        for i in range(len(to_rep)):
            expected.replace(to_rep[i], values[i], inplace=True)
        assert_frame_equal(result, expected)

        self.assertRaises(ValueError, df.replace, to_rep, values[1:])

        # list to scalar
        to_rep = [np.nan, 0, '']
        result = df.replace(to_rep, -1)
        expected = df.copy()
        for i in range(len(to_rep)):
            expected.replace(to_rep[i], -1, inplace=True)
        assert_frame_equal(result, expected)

    def test_replace_limit(self):
        pass

    def test_combine_multiple_frames_dtypes(self):
        from pandas import concat

        # GH 2759
        A = DataFrame(data=np.ones((10, 2)), columns=['foo', 'bar'], dtype=np.float64)
        B = DataFrame(data=np.ones((10, 2)), dtype=np.float32)
        results = concat((A, B), axis=1).get_dtype_counts()
        expected = Series(dict( float64 = 2, float32 = 2 ))
        assert_series_equal(results,expected)

    def test_truncate(self):
        offset = datetools.bday

        ts = self.tsframe[::3]

        start, end = self.tsframe.index[3], self.tsframe.index[6]

        start_missing = self.tsframe.index[2]
        end_missing = self.tsframe.index[7]

        # neither specified
        truncated = ts.truncate()
        assert_frame_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        assert_frame_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        assert_frame_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        assert_frame_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        assert_frame_equal(truncated, expected)

    def test_truncate_copy(self):
        index = self.tsframe.index
        truncated = self.tsframe.truncate(index[5], index[10])
        truncated.values[:] = 5.
        self.assert_(not (self.tsframe.values[5:11] == 5).any())

    def test_xs(self):
        idx = self.frame.index[5]
        xs = self.frame.xs(idx)
        for item, value in xs.iteritems():
            if np.isnan(value):
                self.assert_(np.isnan(self.frame[item][idx]))
            else:
                self.assertEqual(value, self.frame[item][idx])

        # mixed-type xs
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }
        frame = DataFrame(test_data)
        xs = frame.xs('1')
        self.assert_(xs.dtype == np.object_)
        self.assertEqual(xs['A'], 1)
        self.assertEqual(xs['B'], '1')

        self.assertRaises(Exception, self.tsframe.xs,
                          self.tsframe.index[0] - datetools.bday)

        # xs get column
        series = self.frame.xs('A', axis=1)
        expected = self.frame['A']
        assert_series_equal(series, expected)

        # no view by default
        series[:] = 5
        self.assert_((expected != 5).all())

        # view
        series = self.frame.xs('A', axis=1, copy=False)
        series[:] = 5
        self.assert_((expected == 5).all())

    def test_xs_corner(self):
        # pathological mixed-type reordering case
        df = DataFrame(index=[0])
        df['A'] = 1.
        df['B'] = 'foo'
        df['C'] = 2.
        df['D'] = 'bar'
        df['E'] = 3.

        xs = df.xs(0)
        assert_almost_equal(xs, [1., 'foo', 2., 'bar', 3.])

        # no columns but index
        df = DataFrame(index=['a', 'b', 'c'])
        result = df.xs('a')
        expected = Series([])
        assert_series_equal(result, expected)

    def test_xs_duplicates(self):
        df = DataFrame(randn(5, 2), index=['b', 'b', 'c', 'b', 'a'])

        cross = df.xs('c')
        exp = df.irow(2)
        assert_series_equal(cross, exp)

    def test_pivot(self):
        data = {
            'index': ['A', 'B', 'C', 'C', 'B', 'A'],
            'columns': ['One', 'One', 'One', 'Two', 'Two', 'Two'],
            'values': [1., 2., 3., 3., 2., 1.]
        }

        frame = DataFrame(data)
        pivoted = frame.pivot(
            index='index', columns='columns', values='values')

        expected = DataFrame({
            'One': {'A': 1., 'B': 2., 'C': 3.},
            'Two': {'A': 1., 'B': 2., 'C': 3.}
        })
        expected.index.name, expected.columns.name = 'index', 'columns'

        assert_frame_equal(pivoted, expected)

        # name tracking
        self.assertEqual(pivoted.index.name, 'index')
        self.assertEqual(pivoted.columns.name, 'columns')

        # don't specify values
        pivoted = frame.pivot(index='index', columns='columns')
        self.assertEqual(pivoted.index.name, 'index')
        self.assertEqual(pivoted.columns.names, [None, 'columns'])

        # pivot multiple columns
        wp = tm.makePanel()
        lp = wp.to_frame()
        df = lp.reset_index()
        assert_frame_equal(df.pivot('major', 'minor'), lp.unstack())

    def test_pivot_duplicates(self):
        data = DataFrame({'a': ['bar', 'bar', 'foo', 'foo', 'foo'],
                          'b': ['one', 'two', 'one', 'one', 'two'],
                          'c': [1., 2., 3., 3., 4.]})
        self.assertRaises(Exception, data.pivot, 'a', 'b', 'c')

    def test_pivot_empty(self):
        df = DataFrame({}, columns=['a', 'b', 'c'])
        result = df.pivot('a', 'b', 'c')
        expected = DataFrame({})
        assert_frame_equal(result, expected, check_names=False)

    def test_pivot_integer_bug(self):
        df = DataFrame(data=[("A", "1", "A1"), ("B", "2", "B2")])

        result = df.pivot(index=1, columns=0, values=2)
        repr(result)
        self.assert_(np.array_equal(result.columns, ['A', 'B']))

    def test_reindex(self):
        newFrame = self.frame.reindex(self.ts1.index)

        for col in newFrame.columns:
            for idx, val in newFrame[col].iteritems():
                if idx in self.frame.index:
                    if np.isnan(val):
                        self.assert_(np.isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assert_(np.isnan(val))

        for col, series in newFrame.iteritems():
            self.assert_(tm.equalContents(series.index, newFrame.index))
        emptyFrame = self.frame.reindex(Index([]))
        self.assert_(len(emptyFrame.index) == 0)

        # Cython code should be unit-tested directly
        nonContigFrame = self.frame.reindex(self.ts1.index[::2])

        for col in nonContigFrame.columns:
            for idx, val in nonContigFrame[col].iteritems():
                if idx in self.frame.index:
                    if np.isnan(val):
                        self.assert_(np.isnan(self.frame[col][idx]))
                    else:
                        self.assertEqual(val, self.frame[col][idx])
                else:
                    self.assert_(np.isnan(val))

        for col, series in nonContigFrame.iteritems():
            self.assert_(tm.equalContents(series.index,
                                          nonContigFrame.index))

        # corner cases

        # Same index, copies values
        newFrame = self.frame.reindex(self.frame.index)
        self.assert_(newFrame.index is self.frame.index)

        # length zero
        newFrame = self.frame.reindex([])
        self.assert_(newFrame.empty)
        self.assertEqual(len(newFrame.columns), len(self.frame.columns))

        # length zero with columns reindexed with non-empty index
        newFrame = self.frame.reindex([])
        newFrame = newFrame.reindex(self.frame.index)
        self.assertEqual(len(newFrame.index), len(self.frame.index))
        self.assertEqual(len(newFrame.columns), len(self.frame.columns))

        # pass non-Index
        newFrame = self.frame.reindex(list(self.ts1.index))
        self.assert_(newFrame.index.equals(self.ts1.index))

    def test_reindex_name_remains(self):
        s = Series(random.rand(10))
        df = DataFrame(s, index=np.arange(len(s)))
        i = Series(np.arange(10), name='iname')

        df = df.reindex(i)
        self.assert_(df.index.name == 'iname')

        df = df.reindex(Index(np.arange(10), name='tmpname'))
        self.assert_(df.index.name == 'tmpname')

        s = Series(random.rand(10))
        df = DataFrame(s.T, index=np.arange(len(s)))
        i = Series(np.arange(10), name='iname')
        df = df.reindex(columns=i)
        self.assert_(df.columns.name == 'iname')

    def test_reindex_int(self):
        smaller = self.intframe.reindex(self.intframe.index[::2])

        self.assert_(smaller['A'].dtype == np.int64)

        bigger = smaller.reindex(self.intframe.index)
        self.assert_(bigger['A'].dtype == np.float64)

        smaller = self.intframe.reindex(columns=['A', 'B'])
        self.assert_(smaller['A'].dtype == np.int64)

    def test_reindex_like(self):
        other = self.frame.reindex(index=self.frame.index[:10],
                                   columns=['C', 'B'])

        assert_frame_equal(other, self.frame.reindex_like(other))

    def test_reindex_columns(self):
        newFrame = self.frame.reindex(columns=['A', 'B', 'E'])

        assert_series_equal(newFrame['B'], self.frame['B'])
        self.assert_(np.isnan(newFrame['E']).all())
        self.assert_('C' not in newFrame)

        # length zero
        newFrame = self.frame.reindex(columns=[])
        self.assert_(newFrame.empty)

    def test_reindex_fill_value(self):
        df = DataFrame(np.random.randn(10, 4))

        # axis=0
        result = df.reindex(range(15))
        self.assert_(np.isnan(result.values[-5:]).all())

        result = df.reindex(range(15), fill_value=0)
        expected = df.reindex(range(15)).fillna(0)
        assert_frame_equal(result, expected)

        # axis=1
        result = df.reindex(columns=range(5), fill_value=0.)
        expected = df.copy()
        expected[4] = 0.
        assert_frame_equal(result, expected)

        result = df.reindex(columns=range(5), fill_value=0)
        expected = df.copy()
        expected[4] = 0
        assert_frame_equal(result, expected)

        result = df.reindex(columns=range(5), fill_value='foo')
        expected = df.copy()
        expected[4] = 'foo'
        assert_frame_equal(result, expected)

        # reindex_axis
        result = df.reindex_axis(range(15), fill_value=0., axis=0)
        expected = df.reindex(range(15)).fillna(0)
        assert_frame_equal(result, expected)

        result = df.reindex_axis(range(5), fill_value=0., axis=1)
        expected = df.reindex(columns=range(5)).fillna(0)
        assert_frame_equal(result, expected)

        # other dtypes
        df['foo'] = 'foo'
        result = df.reindex(range(15), fill_value=0)
        expected = df.reindex(range(15)).fillna(0)
        assert_frame_equal(result, expected)

    def test_align(self):
        af, bf = self.frame.align(self.frame)
        self.assert_(af._data is not self.frame._data)

        af, bf = self.frame.align(self.frame, copy=False)
        self.assert_(af._data is self.frame._data)

        # axis = 0
        other = self.frame.ix[:-5, :3]
        af, bf = self.frame.align(other, axis=0, fill_value=-1)
        self.assert_(bf.columns.equals(other.columns))
        # test fill value
        join_idx = self.frame.index.join(other.index)
        diff_a = self.frame.index.diff(join_idx)
        diff_b = other.index.diff(join_idx)
        diff_a_vals = af.reindex(diff_a).values
        diff_b_vals = bf.reindex(diff_b).values
        self.assert_((diff_a_vals == -1).all())

        af, bf = self.frame.align(other, join='right', axis=0)
        self.assert_(bf.columns.equals(other.columns))
        self.assert_(bf.index.equals(other.index))
        self.assert_(af.index.equals(other.index))

        # axis = 1
        other = self.frame.ix[:-5, :3].copy()
        af, bf = self.frame.align(other, axis=1)
        self.assert_(bf.columns.equals(self.frame.columns))
        self.assert_(bf.index.equals(other.index))

        # test fill value
        join_idx = self.frame.index.join(other.index)
        diff_a = self.frame.index.diff(join_idx)
        diff_b = other.index.diff(join_idx)
        diff_a_vals = af.reindex(diff_a).values
        diff_b_vals = bf.reindex(diff_b).values
        self.assert_((diff_a_vals == -1).all())

        af, bf = self.frame.align(other, join='inner', axis=1)
        self.assert_(bf.columns.equals(other.columns))

        af, bf = self.frame.align(other, join='inner', axis=1, method='pad')
        self.assert_(bf.columns.equals(other.columns))

        # test other non-float types
        af, bf = self.intframe.align(other, join='inner', axis=1, method='pad')
        self.assert_(bf.columns.equals(other.columns))

        af, bf = self.mixed_frame.align(self.mixed_frame,
                                        join='inner', axis=1, method='pad')
        self.assert_(bf.columns.equals(self.mixed_frame.columns))

        af, bf = self.frame.align(other.ix[:, 0], join='inner', axis=1,
                                  method=None, fill_value=None)
        self.assert_(bf.index.equals(Index([])))

        af, bf = self.frame.align(other.ix[:, 0], join='inner', axis=1,
                                  method=None, fill_value=0)
        self.assert_(bf.index.equals(Index([])))

        # mixed floats/ints
        af, bf = self.mixed_float.align(other.ix[:, 0], join='inner', axis=1,
                                        method=None, fill_value=0)
        self.assert_(bf.index.equals(Index([])))

        af, bf = self.mixed_int.align(other.ix[:, 0], join='inner', axis=1,
                                        method=None, fill_value=0)
        self.assert_(bf.index.equals(Index([])))

        # try to align dataframe to series along bad axis
        self.assertRaises(ValueError, self.frame.align, af.ix[0, :3],
                          join='inner', axis=2)

    def _check_align(self, a, b, axis, fill_axis, how, method, limit=None):
        aa, ab = a.align(b, axis=axis, join=how, method=method, limit=limit,
                         fill_axis=fill_axis)

        join_index, join_columns = None, None

        ea, eb = a, b
        if axis is None or axis == 0:
            join_index = a.index.join(b.index, how=how)
            ea = ea.reindex(index=join_index)
            eb = eb.reindex(index=join_index)

        if axis is None or axis == 1:
            join_columns = a.columns.join(b.columns, how=how)
            ea = ea.reindex(columns=join_columns)
            eb = eb.reindex(columns=join_columns)

        ea = ea.fillna(axis=fill_axis, method=method, limit=limit)
        eb = eb.fillna(axis=fill_axis, method=method, limit=limit)

        assert_frame_equal(aa, ea)
        assert_frame_equal(ab, eb)

    def test_align_fill_method_inner(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('inner', meth, ax, fax)

    def test_align_fill_method_outer(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('outer', meth, ax, fax)

    def test_align_fill_method_left(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('left', meth, ax, fax)

    def test_align_fill_method_right(self):
        for meth in ['pad', 'bfill']:
            for ax in [0, 1, None]:
                for fax in [0, 1]:
                    self._check_align_fill('right', meth, ax, fax)

    def _check_align_fill(self, kind, meth, ax, fax):
        left = self.frame.ix[0:4, :10]
        right = self.frame.ix[2:, 6:]
        empty = self.frame.ix[:0, :0]

        self._check_align(left, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(left, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

        # empty left
        self._check_align(empty, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(empty, right, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

        # empty right
        self._check_align(left, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(left, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

        # both empty
        self._check_align(empty, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth)
        self._check_align(empty, empty, axis=ax, fill_axis=fax,
                          how=kind, method=meth, limit=1)

    def test_align_int_fill_bug(self):
        # GH #910
        X = np.arange(10*10, dtype='float64').reshape(10, 10)
        Y = np.ones((10, 1), dtype=int)

        df1 = DataFrame(X)
        df1['0.X'] = Y.squeeze()

        df2 = df1.astype(float)

        result = df1 - df1.mean()
        expected = df2 - df2.mean()
        assert_frame_equal(result, expected)

    def test_where(self):
        default_frame = DataFrame(np.random.randn(5, 3),columns=['A','B','C'])

        def _safe_add(df):
            # only add to the numeric items
            def is_ok(s):
                return issubclass(s.dtype.type, (np.integer,np.floating)) and s.dtype != 'uint8'
            return DataFrame(dict([ (c,s+1) if is_ok(s) else (c,s) for c, s in df.iteritems() ]))

        def _check_get(df, cond, check_dtypes = True):
            other1 = _safe_add(df)
            rs = df.where(cond, other1)
            rs2 = df.where(cond.values, other1)
            for k, v in rs.iteritems():
                assert_series_equal(v, np.where(cond[k], df[k], other1[k]))
            assert_frame_equal(rs, rs2)

            # dtypes
            if check_dtypes:
                self.assert_((rs.dtypes == df.dtypes).all() == True)


        # check getting
        for df in [ default_frame, self.mixed_frame, self.mixed_float, self.mixed_int ]:
            cond = df > 0
            _check_get(df, cond)


        # upcasting case (GH # 2794)
        df = DataFrame(dict([ (c,Series([1]*3,dtype=c)) for c in ['int64','int32','float32','float64'] ]))
        df.ix[1,:] = 0

        result = df.where(df>=0).get_dtype_counts()

        #### when we don't preserver boolean casts ####
        #expected = Series({ 'float32' : 1, 'float64' : 3 })

        expected = Series({ 'float32' : 1, 'float64' : 1, 'int32' : 1, 'int64' : 1 })
        assert_series_equal(result, expected)

        # aligning
        def _check_align(df, cond, other, check_dtypes = True):
            rs = df.where(cond, other)
            for i, k in enumerate(rs.columns):
                result = rs[k]
                d = df[k].values
                c = cond[k].reindex(df[k].index).fillna(False).values

                if np.isscalar(other):
                    o = other
                else:
                    if isinstance(other,np.ndarray):
                        o = Series(other[:,i],index=result.index).values
                    else:
                        o = other[k].values

                new_values = d if c.all() else np.where(c, d, o)
                expected = Series(new_values,index=result.index)

                # since we can't always have the correct numpy dtype
                # as numpy doesn't know how to downcast, don't check
                assert_series_equal(result, expected, check_dtype=False)

            # dtypes
            # can't check dtype when other is an ndarray

            if check_dtypes and not isinstance(other,np.ndarray):
                self.assert_((rs.dtypes == df.dtypes).all() == True)

        for df in [ self.mixed_frame, self.mixed_float, self.mixed_int ]:

            # other is a frame
            cond = (df > 0)[1:]
            _check_align(df, cond, _safe_add(df))

            # check other is ndarray
            cond = df > 0
            _check_align(df, cond, (_safe_add(df).values))

            # integers are upcast, so don't check the dtypes
            cond = df > 0
            check_dtypes = all([ not issubclass(s.type,np.integer) for s in df.dtypes ])
            _check_align(df, cond, np.nan, check_dtypes = check_dtypes)

        # invalid conditions
        df = default_frame
        err1 = (df + 1).values[0:2, :]
        self.assertRaises(ValueError, df.where, cond, err1)

        err2 = cond.ix[:2, :].values
        other1 = _safe_add(df)
        self.assertRaises(ValueError, df.where, err2, other1)

        self.assertRaises(ValueError, df.mask, True)
        self.assertRaises(ValueError, df.mask, 0)

        # where inplace
        def _check_set(df, cond, check_dtypes = True):
            dfi = df.copy()
            econd = cond.reindex_like(df).fillna(True)
            expected = dfi.mask(~econd)

            dfi.where(cond, np.nan, inplace=True)
            assert_frame_equal(dfi, expected)

            # dtypes (and confirm upcasts)x
            if check_dtypes:
                for k, v in df.dtypes.iteritems():
                    if issubclass(v.type,np.integer) and not cond[k].all():
                        v = np.dtype('float64')
                    self.assert_(dfi[k].dtype == v)

        for df in [ default_frame, self.mixed_frame, self.mixed_float, self.mixed_int ]:

            cond = df > 0
            _check_set(df, cond)

            cond = df >= 0
            _check_set(df, cond)

            # aligining
            cond = (df >= 0)[1:]
            _check_set(df, cond)

    def test_where_bug(self):

        # GH 2793

        df = DataFrame({'a': [1.0, 2.0, 3.0, 4.0], 'b': [4.0, 3.0, 2.0, 1.0]}, dtype = 'float64')
        expected = DataFrame({'a': [np.nan, np.nan, 3.0, 4.0], 'b': [4.0, 3.0, np.nan, np.nan]}, dtype = 'float64')
        result   = df.where(df > 2, np.nan)
        assert_frame_equal(result, expected)

        result = df.copy()
        result.where(result > 2, np.nan, inplace=True)
        assert_frame_equal(result, expected)

        # mixed
        for dtype in ['int16','int8','int32','int64']:
            df = DataFrame({'a': np.array([1, 2, 3, 4],dtype=dtype), 'b': np.array([4.0, 3.0, 2.0, 1.0], dtype = 'float64') })
            expected = DataFrame({'a': [np.nan, np.nan, 3.0, 4.0], 'b': [4.0, 3.0, np.nan, np.nan]}, dtype = 'float64')
            result   = df.where(df > 2, np.nan)
            assert_frame_equal(result, expected)

            result = df.copy()
            result.where(result > 2, np.nan, inplace=True)
            assert_frame_equal(result, expected)

    def test_where_datetime(self):

        # GH 3311
        df = DataFrame(dict(A = date_range('20130102',periods=5),
                            B = date_range('20130104',periods=5),
                            C = np.random.randn(5)))

        stamp = datetime(2013,1,3)
        result = df[df>stamp]
        expected = df.copy()
        expected.loc[[0,1],'A'] = np.nan
        assert_frame_equal(result,expected)

    def test_mask(self):
        df = DataFrame(np.random.randn(5, 3))
        cond = df > 0

        rs = df.where(cond, np.nan)
        assert_frame_equal(rs, df.mask(df <= 0))
        assert_frame_equal(rs, df.mask(~cond))

    def test_mask_edge_case_1xN_frame(self):
        # GH4071
        df = DataFrame([[1, 2]])
        res = df.mask(DataFrame([[True, False]]))
        expec = DataFrame([[nan, 2]])
        assert_frame_equal(res, expec)

    #----------------------------------------------------------------------
    # Transposing
    def test_transpose(self):
        frame = self.frame
        dft = frame.T
        for idx, series in dft.iteritems():
            for col, value in series.iteritems():
                if np.isnan(value):
                    self.assert_(np.isnan(frame[col][idx]))
                else:
                    self.assertEqual(value, frame[col][idx])

        # mixed type
        index, data = tm.getMixedTypeDict()
        mixed = DataFrame(data, index=index)

        mixed_T = mixed.T
        for col, s in mixed_T.iteritems():
            self.assert_(s.dtype == np.object_)

    def test_transpose_get_view(self):
        dft = self.frame.T
        dft.values[:, 5:10] = 5

        self.assert_((self.frame.values[5:10] == 5).all())

    #----------------------------------------------------------------------
    # Renaming

    def test_rename(self):
        mapping = {
            'A': 'a',
            'B': 'b',
            'C': 'c',
            'D': 'd'
        }

        renamed = self.frame.rename(columns=mapping)
        renamed2 = self.frame.rename(columns=str.lower)

        assert_frame_equal(renamed, renamed2)
        assert_frame_equal(renamed2.rename(columns=str.upper),
                           self.frame, check_names=False)

        # index
        data = {
            'A': {'foo': 0, 'bar': 1}
        }

        # gets sorted alphabetical
        df = DataFrame(data)
        renamed = df.rename(index={'foo': 'bar', 'bar': 'foo'})
        self.assert_(np.array_equal(renamed.index, ['foo', 'bar']))

        renamed = df.rename(index=str.upper)
        self.assert_(np.array_equal(renamed.index, ['BAR', 'FOO']))

        # have to pass something
        self.assertRaises(Exception, self.frame.rename)

        # partial columns
        renamed = self.frame.rename(columns={'C': 'foo', 'D': 'bar'})
        self.assert_(np.array_equal(renamed.columns, ['A', 'B', 'foo', 'bar']))

        # other axis
        renamed = self.frame.T.rename(index={'C': 'foo', 'D': 'bar'})
        self.assert_(np.array_equal(renamed.index, ['A', 'B', 'foo', 'bar']))

        # index with name
        index = Index(['foo', 'bar'], name='name')
        renamer = DataFrame(data, index=index)
        renamed = renamer.rename(index={'foo': 'bar', 'bar': 'foo'})
        self.assert_(np.array_equal(renamed.index, ['bar', 'foo']))
        self.assertEquals(renamed.index.name, renamer.index.name)

        # MultiIndex
        tuples_index = [('foo1', 'bar1'), ('foo2', 'bar2')]
        tuples_columns = [('fizz1', 'buzz1'), ('fizz2', 'buzz2')]
        index = MultiIndex.from_tuples(tuples_index, names=['foo', 'bar'])
        columns = MultiIndex.from_tuples(tuples_columns, names=['fizz', 'buzz'])
        renamer = DataFrame([(0,0),(1,1)], index=index, columns=columns)
        renamed = renamer.rename(index={'foo1': 'foo3', 'bar2': 'bar3'},
                                 columns={'fizz1': 'fizz3', 'buzz2': 'buzz3'})
        new_index = MultiIndex.from_tuples([('foo3', 'bar1'), ('foo2', 'bar3')])
        new_columns = MultiIndex.from_tuples([('fizz3', 'buzz1'), ('fizz2', 'buzz3')])
        self.assert_(np.array_equal(renamed.index, new_index))
        self.assert_(np.array_equal(renamed.columns, new_columns))
        self.assertEquals(renamed.index.names, renamer.index.names)
        self.assertEquals(renamed.columns.names, renamer.columns.names)

    def test_rename_nocopy(self):
        renamed = self.frame.rename(columns={'C': 'foo'}, copy=False)
        renamed['foo'] = 1.
        self.assert_((self.frame['C'] == 1.).all())

    def test_rename_inplace(self):
        self.frame.rename(columns={'C': 'foo'})
        self.assert_('C' in self.frame)
        self.assert_('foo' not in self.frame)

        c_id = id(self.frame['C'])
        frame = self.frame.copy()
        frame.rename(columns={'C': 'foo'}, inplace=True)

        self.assert_('C' not in frame)
        self.assert_('foo' in frame)
        self.assert_(id(frame['foo']) != c_id)

    #----------------------------------------------------------------------
    # Time series related
    def test_diff(self):
        the_diff = self.tsframe.diff(1)

        assert_series_equal(the_diff['A'],
                            self.tsframe['A'] - self.tsframe['A'].shift(1))

        # int dtype
        a = 10000000000000000
        b = a + 1
        s = Series([a, b])

        rs = DataFrame({'s': s}).diff()
        self.assertEqual(rs.s[1], 1)

        # mixed numeric
        tf = self.tsframe.astype('float32')
        the_diff = tf.diff(1)
        assert_series_equal(the_diff['A'],
                            tf['A'] - tf['A'].shift(1))


    def test_diff_mixed_dtype(self):
        df = DataFrame(np.random.randn(5, 3))
        df['A'] = np.array([1, 2, 3, 4, 5], dtype=object)

        result = df.diff()
        self.assert_(result[0].dtype == np.float64)

    def test_diff_neg_n(self):
        rs = self.tsframe.diff(-1)
        xp = self.tsframe - self.tsframe.shift(-1)
        assert_frame_equal(rs, xp)

    def test_diff_float_n(self):
        rs = self.tsframe.diff(1.)
        xp = self.tsframe.diff(1)
        assert_frame_equal(rs, xp)

    def test_pct_change(self):
        rs = self.tsframe.pct_change(fill_method=None)
        assert_frame_equal(rs, self.tsframe / self.tsframe.shift(1) - 1)

        rs = self.tsframe.pct_change(2)
        filled = self.tsframe.fillna(method='pad')
        assert_frame_equal(rs, filled / filled.shift(2) - 1)

        rs = self.tsframe.pct_change(fill_method='bfill', limit=1)
        filled = self.tsframe.fillna(method='bfill', limit=1)
        assert_frame_equal(rs, filled / filled.shift(1) - 1)

        rs = self.tsframe.pct_change(freq='5D')
        filled = self.tsframe.fillna(method='pad')
        assert_frame_equal(rs, filled / filled.shift(freq='5D') - 1)

    def test_pct_change_shift_over_nas(self):
        s = Series([1., 1.5, np.nan, 2.5, 3.])

        df = DataFrame({'a': s, 'b': s})

        chg = df.pct_change()
        expected = Series([np.nan, 0.5, np.nan, 2.5 / 1.5 - 1, .2])
        edf = DataFrame({'a': expected, 'b': expected})
        assert_frame_equal(chg, edf)

    def test_shift(self):
        # naive shift
        shiftedFrame = self.tsframe.shift(5)
        self.assert_(shiftedFrame.index.equals(self.tsframe.index))

        shiftedSeries = self.tsframe['A'].shift(5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        shiftedFrame = self.tsframe.shift(-5)
        self.assert_(shiftedFrame.index.equals(self.tsframe.index))

        shiftedSeries = self.tsframe['A'].shift(-5)
        assert_series_equal(shiftedFrame['A'], shiftedSeries)

        # shift by 0
        unshifted = self.tsframe.shift(0)
        assert_frame_equal(unshifted, self.tsframe)

        # shift by DateOffset
        shiftedFrame = self.tsframe.shift(5, freq=datetools.BDay())
        self.assert_(len(shiftedFrame) == len(self.tsframe))

        shiftedFrame2 = self.tsframe.shift(5, freq='B')
        assert_frame_equal(shiftedFrame, shiftedFrame2)

        d = self.tsframe.index[0]
        shifted_d = d + datetools.BDay(5)
        assert_series_equal(self.tsframe.xs(d),
                            shiftedFrame.xs(shifted_d))

        # shift int frame
        int_shifted = self.intframe.shift(1)

        # Shifting with PeriodIndex
        ps = tm.makePeriodFrame()
        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        self.assert_(shifted.index.equals(ps.index))

        tm.assert_dict_equal(unshifted.ix[:, 0].valid(), ps.ix[:, 0],
                             compare_keys=False)

        shifted2 = ps.shift(1, 'B')
        shifted3 = ps.shift(1, datetools.bday)
        assert_frame_equal(shifted2, shifted3)
        assert_frame_equal(ps, shifted2.shift(-1, 'B'))

        self.assertRaises(ValueError, ps.shift, freq='D')

    def test_shift_bool(self):
        df = DataFrame({'high': [True, False],
                        'low': [False, False]})
        rs = df.shift(1)
        xp = DataFrame(np.array([[np.nan, np.nan],
                                 [True, False]], dtype=object),
                       columns=['high', 'low'])
        assert_frame_equal(rs, xp)

    def test_tshift(self):
        # PeriodIndex
        ps = tm.makePeriodFrame()
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_frame_equal(unshifted, ps)

        shifted2 = ps.tshift(freq='B')
        assert_frame_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=datetools.bday)
        assert_frame_equal(shifted, shifted3)

        self.assertRaises(ValueError, ps.tshift, freq='M')

        # DatetimeIndex
        shifted = self.tsframe.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_frame_equal(self.tsframe, unshifted)

        shifted2 = self.tsframe.tshift(freq=self.tsframe.index.freq)
        assert_frame_equal(shifted, shifted2)

        inferred_ts = DataFrame(self.tsframe.values,
                                Index(np.asarray(self.tsframe.index)),
                                columns=self.tsframe.columns)
        shifted = inferred_ts.tshift(1)
        unshifted = shifted.tshift(-1)
        assert_frame_equal(shifted, self.tsframe.tshift(1))
        assert_frame_equal(unshifted, inferred_ts)

        no_freq = self.tsframe.ix[[0, 5, 7], :]
        self.assertRaises(ValueError, no_freq.tshift)

    def test_apply(self):
        # ufunc
        applied = self.frame.apply(np.sqrt)
        assert_series_equal(np.sqrt(self.frame['A']), applied['A'])

        # aggregator
        applied = self.frame.apply(np.mean)
        self.assertEqual(applied['A'], np.mean(self.frame['A']))

        d = self.frame.index[0]
        applied = self.frame.apply(np.mean, axis=1)
        self.assertEqual(applied[d], np.mean(self.frame.xs(d)))
        self.assert_(applied.index is self.frame.index)  # want this

        # invalid axis
        df = DataFrame(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=['a', 'a', 'c'])
        self.assertRaises(ValueError, df.apply, lambda x: x, 2)

    def test_apply_empty(self):
        # empty
        applied = self.empty.apply(np.sqrt)
        self.assert_(applied.empty)

        applied = self.empty.apply(np.mean)
        self.assert_(applied.empty)

        no_rows = self.frame[:0]
        result = no_rows.apply(lambda x: x.mean())
        expected = Series(np.nan, index=self.frame.columns)
        assert_series_equal(result, expected)

        no_cols = self.frame.ix[:, []]
        result = no_cols.apply(lambda x: x.mean(), axis=1)
        expected = Series(np.nan, index=self.frame.index)
        assert_series_equal(result, expected)

        # 2476
        xp = DataFrame(index=['a'])
        rs = xp.apply(lambda x: x['a'], axis=1)
        assert_frame_equal(xp, rs)

    def test_apply_standard_nonunique(self):
        df = DataFrame(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]], index=['a', 'a', 'c'])
        rs = df.apply(lambda s: s[0], axis=1)
        xp = Series([1, 4, 7], ['a', 'a', 'c'])
        assert_series_equal(rs, xp)

        rs = df.T.apply(lambda s: s[0], axis=0)
        assert_series_equal(rs, xp)

    def test_apply_broadcast(self):
        broadcasted = self.frame.apply(np.mean, broadcast=True)
        agged = self.frame.apply(np.mean)

        for col, ts in broadcasted.iteritems():
            self.assert_((ts == agged[col]).all())

        broadcasted = self.frame.apply(np.mean, axis=1, broadcast=True)
        agged = self.frame.apply(np.mean, axis=1)
        for idx in broadcasted.index:
            self.assert_((broadcasted.xs(idx) == agged[idx]).all())

    def test_apply_raw(self):
        result0 = self.frame.apply(np.mean, raw=True)
        result1 = self.frame.apply(np.mean, axis=1, raw=True)

        expected0 = self.frame.apply(lambda x: x.values.mean())
        expected1 = self.frame.apply(lambda x: x.values.mean(), axis=1)

        assert_series_equal(result0, expected0)
        assert_series_equal(result1, expected1)

        # no reduction
        result = self.frame.apply(lambda x: x * 2, raw=True)
        expected = self.frame * 2
        assert_frame_equal(result, expected)

    def test_apply_axis1(self):
        d = self.frame.index[0]
        tapplied = self.frame.apply(np.mean, axis=1)
        self.assertEqual(tapplied[d], np.mean(self.frame.xs(d)))

    def test_apply_ignore_failures(self):
        result = self.mixed_frame._apply_standard(np.mean, 0,
                                                  ignore_failures=True)
        expected = self.mixed_frame._get_numeric_data().apply(np.mean)
        assert_series_equal(result, expected)

        # test with hierarchical index

    def test_apply_mixed_dtype_corner(self):
        df = DataFrame({'A': ['foo'],
                        'B': [1.]})
        result = df[:0].apply(np.mean, axis=1)
        # the result here is actually kind of ambiguous, should it be a Series
        # or a DataFrame?
        expected = Series(np.nan, index=[])
        assert_series_equal(result, expected)

    def test_apply_empty_infer_type(self):
        no_cols = DataFrame(index=['a', 'b', 'c'])
        no_index = DataFrame(columns=['a', 'b', 'c'])

        def _check(df, f):
            test_res = f(np.array([], dtype='f8'))
            is_reduction = not isinstance(test_res, np.ndarray)

            def _checkit(axis=0, raw=False):
                res = df.apply(f, axis=axis, raw=raw)
                if is_reduction:
                    agg_axis = df._get_agg_axis(axis)
                    self.assert_(isinstance(res, Series))
                    self.assert_(res.index is agg_axis)
                else:
                    self.assert_(isinstance(res, DataFrame))

            _checkit()
            _checkit(axis=1)
            _checkit(raw=True)
            _checkit(axis=0, raw=True)

        _check(no_cols, lambda x: x)
        _check(no_cols, lambda x: x.mean())
        _check(no_index, lambda x: x)
        _check(no_index, lambda x: x.mean())

        result = no_cols.apply(lambda x: x.mean(), broadcast=True)
        self.assert_(isinstance(result, DataFrame))

    def test_apply_with_args_kwds(self):
        def add_some(x, howmuch=0):
            return x + howmuch

        def agg_and_add(x, howmuch=0):
            return x.mean() + howmuch

        def subtract_and_divide(x, sub, divide=1):
            return (x - sub) / divide

        result = self.frame.apply(add_some, howmuch=2)
        exp = self.frame.apply(lambda x: x + 2)
        assert_frame_equal(result, exp)

        result = self.frame.apply(agg_and_add, howmuch=2)
        exp = self.frame.apply(lambda x: x.mean() + 2)
        assert_series_equal(result, exp)

        res = self.frame.apply(subtract_and_divide, args=(2,), divide=2)
        exp = self.frame.apply(lambda x: (x - 2.) / 2.)
        assert_frame_equal(res, exp)

    def test_apply_yield_list(self):
        result = self.frame.apply(list)
        assert_frame_equal(result, self.frame)

    def test_apply_reduce_Series(self):
        self.frame.ix[::2, 'A'] = np.nan
        result = self.frame.apply(np.mean, axis=1)
        expected = self.frame.mean(1)
        assert_series_equal(result, expected)

    def test_apply_differently_indexed(self):
        df = DataFrame(np.random.randn(20, 10))

        result0 = df.apply(Series.describe, axis=0)
        expected0 = DataFrame(dict((i, v.describe())
                                   for i, v in df.iteritems()),
                              columns=df.columns)
        assert_frame_equal(result0, expected0)

        result1 = df.apply(Series.describe, axis=1)
        expected1 = DataFrame(dict((i, v.describe())
                                   for i, v in df.T.iteritems()),
                              columns=df.index).T
        assert_frame_equal(result1, expected1)

    def test_apply_modify_traceback(self):
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

        data['C'][4] = np.nan

        def transform(row):
            if row['C'].startswith('shin') and row['A'] == 'foo':
                row['D'] = 7
            return row

        def transform2(row):
            if (notnull(row['C']) and row['C'].startswith('shin')
                    and row['A'] == 'foo'):
                row['D'] = 7
            return row

        try:
            transformed = data.apply(transform, axis=1)
        except Exception, e:
            self.assertEqual(len(e.args), 2)
            self.assertEqual(e.args[1], 'occurred at index 4')

    def test_swapaxes(self):
        df = DataFrame(np.random.randn(10, 5))
        assert_frame_equal(df.T, df.swapaxes(0, 1))
        assert_frame_equal(df.T, df.swapaxes(1, 0))
        assert_frame_equal(df, df.swapaxes(0, 0))
        self.assertRaises(ValueError, df.swapaxes, 2, 5)

    def test_apply_convert_objects(self):
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

        result = data.apply(lambda x: x, axis=1)
        assert_frame_equal(result.convert_objects(), data)

    def test_apply_attach_name(self):
        result = self.frame.apply(lambda x: x.name)
        expected = Series(self.frame.columns, index=self.frame.columns)
        assert_series_equal(result, expected)

        result = self.frame.apply(lambda x: x.name, axis=1)
        expected = Series(self.frame.index, index=self.frame.index)
        assert_series_equal(result, expected)

        # non-reductions
        result = self.frame.apply(lambda x: np.repeat(x.name, len(x)))
        expected = DataFrame(np.tile(self.frame.columns,
                                     (len(self.frame.index), 1)),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(result, expected)

        result = self.frame.apply(lambda x: np.repeat(x.name, len(x)),
                                  axis=1)
        expected = DataFrame(np.tile(self.frame.index,
                                     (len(self.frame.columns), 1)).T,
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(result, expected)

    def test_apply_multi_index(self):
        s = DataFrame([[1,2], [3,4], [5,6]])
        s.index = MultiIndex.from_arrays([['a','a','b'], ['c','d','d']])
        s.columns = ['col1','col2']
        res = s.apply(lambda x: Series({'min': min(x), 'max': max(x)}), 1)
        self.assert_(isinstance(res.index, MultiIndex))

    def test_applymap(self):
        applied = self.frame.applymap(lambda x: x * 2)
        assert_frame_equal(applied, self.frame * 2)
        result = self.frame.applymap(type)

        # GH #465, function returning tuples
        result = self.frame.applymap(lambda x: (x, x))
        self.assert_(isinstance(result['A'][0], tuple))

        # GH 2909, object conversion to float in constructor?
        df = DataFrame(data=[1,'a'])
        result = df.applymap(lambda x: x)
        self.assert_(result.dtypes[0] == object)

        df = DataFrame(data=[1.,'a'])
        result = df.applymap(lambda x: x)
        self.assert_(result.dtypes[0] == object)

        # GH2786
        df  = DataFrame(np.random.random((3,4)))
        df2 = df.copy()
        cols = ['a','a','a','a']
        df.columns = cols

        expected = df2.applymap(str)
        expected.columns = cols
        result = df.applymap(str)
        assert_frame_equal(result,expected)

    def test_filter(self):
        # items

        filtered = self.frame.filter(['A', 'B', 'E'])
        self.assertEqual(len(filtered.columns), 2)
        self.assert_('E' not in filtered)

        # like
        fcopy = self.frame.copy()
        fcopy['AA'] = 1

        filtered = fcopy.filter(like='A')
        self.assertEqual(len(filtered.columns), 2)
        self.assert_('AA' in filtered)

        # like with ints in column names
        df = DataFrame(0., index=[0, 1, 2], columns=[0, 1, '_A', '_B'])
        filtered = df.filter(like='_')
        self.assertEqual(len(filtered.columns), 2)

        # pass in None
        self.assertRaises(Exception, self.frame.filter, items=None)

        # objects
        filtered = self.mixed_frame.filter(like='foo')
        self.assert_('foo' in filtered)

        # unicode columns, won't ascii-encode
        df = self.frame.rename(columns={'B': u'\u2202'})
        filtered = df.filter(like='C')
        self.assertTrue('C' in filtered)

    def test_filter_regex_search(self):
        fcopy = self.frame.copy()
        fcopy['AA'] = 1

        # regex
        filtered = fcopy.filter(regex='[A]+')
        self.assertEqual(len(filtered.columns), 2)
        self.assert_('AA' in filtered)

        # doesn't have to be at beginning
        df = DataFrame({'aBBa': [1, 2],
                        'BBaBB': [1, 2],
                        'aCCa': [1, 2],
                        'aCCaBB': [1, 2]})

        result = df.filter(regex='BB')
        exp = df[[x for x in df.columns if 'BB' in x]]
        assert_frame_equal(result, exp)

    def test_filter_corner(self):
        empty = DataFrame()

        result = empty.filter([])
        assert_frame_equal(result, empty)

        result = empty.filter(like='foo')
        assert_frame_equal(result, empty)

    def test_select(self):
        f = lambda x: x.weekday() == 2
        result = self.tsframe.select(f, axis=0)
        expected = self.tsframe.reindex(
            index=self.tsframe.index[[f(x) for x in self.tsframe.index]])
        assert_frame_equal(result, expected)

        result = self.frame.select(lambda x: x in ('B', 'D'), axis=1)
        expected = self.frame.reindex(columns=['B', 'D'])

        assert_frame_equal(result, expected, check_names=False)  # TODO should reindex check_names?

    def test_sort_index(self):
        frame = DataFrame(np.random.randn(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # axis=0
        unordered = frame.ix[[3, 2, 4, 1]]
        sorted_df = unordered.sort_index()
        expected = frame
        assert_frame_equal(sorted_df, expected)

        sorted_df = unordered.sort_index(ascending=False)
        expected = frame[::-1]
        assert_frame_equal(sorted_df, expected)

        # axis=1
        unordered = frame.ix[:, ['D', 'B', 'C', 'A']]
        sorted_df = unordered.sort_index(axis=1)
        expected = frame
        assert_frame_equal(sorted_df, expected)

        sorted_df = unordered.sort_index(axis=1, ascending=False)
        expected = frame.ix[:, ::-1]
        assert_frame_equal(sorted_df, expected)

        # by column
        sorted_df = frame.sort_index(by='A')
        indexer = frame['A'].argsort().values
        expected = frame.ix[frame.index[indexer]]
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort_index(by='A', ascending=False)
        indexer = indexer[::-1]
        expected = frame.ix[frame.index[indexer]]
        assert_frame_equal(sorted_df, expected)

        # check for now
        sorted_df = frame.sort(columns='A')
        expected = frame.sort_index(by='A')
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort(columns='A', ascending=False)
        expected = frame.sort_index(by='A', ascending=False)
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort(columns=['A', 'B'], ascending=False)
        expected = frame.sort_index(by=['A', 'B'], ascending=False)
        assert_frame_equal(sorted_df, expected)

        self.assertRaises(ValueError, frame.sort_index, axis=2, inplace=True)

    def test_sort_index_multicolumn(self):
        import random
        A = np.arange(5).repeat(20)
        B = np.tile(np.arange(5), 20)
        random.shuffle(A)
        random.shuffle(B)
        frame = DataFrame({'A': A, 'B': B,
                           'C': np.random.randn(100)})

        result = frame.sort_index(by=['A', 'B'])
        indexer = np.lexsort((frame['B'], frame['A']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

        result = frame.sort_index(by=['A', 'B'], ascending=False)
        indexer = np.lexsort((frame['B'].rank(ascending=False),
                              frame['A'].rank(ascending=False)))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

        result = frame.sort_index(by=['B', 'A'])
        indexer = np.lexsort((frame['A'], frame['B']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

    def test_sort_index_inplace(self):
        frame = DataFrame(np.random.randn(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        # axis=0
        unordered = frame.ix[[3, 2, 4, 1]]
        a_id = id(unordered['A'])
        df = unordered.copy()
        df.sort_index(inplace=True)
        expected = frame
        assert_frame_equal(df, expected)
        self.assert_(a_id != id(df['A']))

        df = unordered.copy()
        df.sort_index(ascending=False, inplace=True)
        expected = frame[::-1]
        assert_frame_equal(df, expected)

        # axis=1
        unordered = frame.ix[:, ['D', 'B', 'C', 'A']]
        df = unordered.copy()
        df.sort_index(axis=1, inplace=True)
        expected = frame
        assert_frame_equal(df, expected)

        df = unordered.copy()
        df.sort_index(axis=1, ascending=False, inplace=True)
        expected = frame.ix[:, ::-1]
        assert_frame_equal(df, expected)

    def test_sort_index_different_sortorder(self):
        import random
        A = np.arange(20).repeat(5)
        B = np.tile(np.arange(5), 20)

        indexer = np.random.permutation(100)
        A = A.take(indexer)
        B = B.take(indexer)

        df = DataFrame({'A': A, 'B': B,
                        'C': np.random.randn(100)})

        result = df.sort_index(by=['A', 'B'], ascending=[1, 0])

        ex_indexer = np.lexsort((df.B.max() - df.B, df.A))
        expected = df.take(ex_indexer)
        assert_frame_equal(result, expected)

        # test with multiindex, too
        idf = df.set_index(['A', 'B'])

        result = idf.sort_index(ascending=[1, 0])
        expected = idf.take(ex_indexer)
        assert_frame_equal(result, expected)

        # also, Series!
        result = idf['C'].sort_index(ascending=[1, 0])
        assert_series_equal(result, expected['C'])

    def test_sort_inplace(self):
        frame = DataFrame(np.random.randn(4, 4), index=[1, 2, 3, 4],
                          columns=['A', 'B', 'C', 'D'])

        sorted_df = frame.copy()
        sorted_df.sort(columns='A', inplace=True)
        expected = frame.sort_index(by='A')
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.copy()
        sorted_df.sort(columns='A', ascending=False, inplace=True)
        expected = frame.sort_index(by='A', ascending=False)
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.copy()
        sorted_df.sort(columns=['A', 'B'], ascending=False, inplace=True)
        expected = frame.sort_index(by=['A', 'B'], ascending=False)
        assert_frame_equal(sorted_df, expected)

    def test_sort_index_duplicates(self):
        df = DataFrame([[1, 2], [3, 4]], columns=['a', 'a'])

        try:
            df.sort_index(by='a')
        except Exception, e:
            self.assertTrue('duplicate' in str(e))

        try:
            df.sort_index(by=['a'])
        except Exception, e:
            self.assertTrue('duplicate' in str(e))

    def test_sort_datetimes(self):

        # GH 3461, argsort / lexsort differences for a datetime column
        df = DataFrame(['a','a','a','b','c','d','e','f','g'],
                       columns=['A'],
                       index=date_range('20130101',periods=9))
        dts = [Timestamp(x)
               for x in  ['2004-02-11','2004-01-21','2004-01-26',
                          '2005-09-20','2010-10-04','2009-05-12',
                          '2008-11-12','2010-09-28','2010-09-28']]
        df['B'] = dts[::2] + dts[1::2]
        df['C'] = 2.
        df['A1'] = 3.

        df1 = df.sort(columns='A')
        df2 = df.sort(columns=['A'])
        assert_frame_equal(df1,df2)

        df1 = df.sort(columns='B')
        df2 = df.sort(columns=['B'])
        assert_frame_equal(df1,df2)

    def test_frame_column_inplace_sort_exception(self):
        s = self.frame['A']
        self.assertRaises(Exception, s.sort)

        cp = s.copy()
        cp.sort()  # it works!

    def test_combine_first(self):
        # disjoint
        head, tail = self.frame[:5], self.frame[5:]

        combined = head.combine_first(tail)
        reordered_frame = self.frame.reindex(combined.index)
        assert_frame_equal(combined, reordered_frame)
        self.assert_(tm.equalContents(combined.columns, self.frame.columns))
        assert_series_equal(combined['A'], reordered_frame['A'])

        # same index
        fcopy = self.frame.copy()
        fcopy['A'] = 1
        del fcopy['C']

        fcopy2 = self.frame.copy()
        fcopy2['B'] = 0
        del fcopy2['D']

        combined = fcopy.combine_first(fcopy2)

        self.assert_((combined['A'] == 1).all())
        assert_series_equal(combined['B'], fcopy['B'])
        assert_series_equal(combined['C'], fcopy2['C'])
        assert_series_equal(combined['D'], fcopy['D'])

        # overlap
        head, tail = reordered_frame[:10].copy(), reordered_frame
        head['A'] = 1

        combined = head.combine_first(tail)
        self.assert_((combined['A'][:10] == 1).all())

        # reverse overlap
        tail['A'][:10] = 0
        combined = tail.combine_first(head)
        self.assert_((combined['A'][:10] == 0).all())

        # no overlap
        f = self.frame[:10]
        g = self.frame[10:]
        combined = f.combine_first(g)
        assert_series_equal(combined['A'].reindex(f.index), f['A'])
        assert_series_equal(combined['A'].reindex(g.index), g['A'])

        # corner cases
        comb = self.frame.combine_first(self.empty)
        assert_frame_equal(comb, self.frame)

        comb = self.empty.combine_first(self.frame)
        assert_frame_equal(comb, self.frame)

        comb = self.frame.combine_first(DataFrame(index=["faz", "boo"]))
        self.assertTrue("faz" in comb.index)

        # #2525
        df = DataFrame({'a': [1]}, index=[datetime(2012, 1, 1)])
        df2 = DataFrame({}, columns=['b'])
        result = df.combine_first(df2)
        self.assertTrue('b' in result)

    def test_combine_first_mixed_bug(self):
        idx = Index(['a', 'b', 'c', 'e'])
        ser1 = Series([5.0, -9.0, 4.0, 100.], index=idx)
        ser2 = Series(['a', 'b', 'c', 'e'], index=idx)
        ser3 = Series([12, 4, 5, 97], index=idx)

        frame1 = DataFrame({"col0": ser1,
                            "col2": ser2,
                            "col3": ser3})

        idx = Index(['a', 'b', 'c', 'f'])
        ser1 = Series([5.0, -9.0, 4.0, 100.], index=idx)
        ser2 = Series(['a', 'b', 'c', 'f'], index=idx)
        ser3 = Series([12, 4, 5, 97], index=idx)

        frame2 = DataFrame({"col1": ser1,
                            "col2": ser2,
                            "col5": ser3})

        combined = frame1.combine_first(frame2)
        self.assertEqual(len(combined.columns), 5)

        # gh 3016 (same as in update)
        df = DataFrame([[1.,2.,False, True],[4.,5.,True,False]],
                       columns=['A','B','bool1','bool2'])

        other = DataFrame([[45,45]],index=[0],columns=['A','B'])
        result = df.combine_first(other)
        assert_frame_equal(result, df)

        df.ix[0,'A'] = np.nan
        result = df.combine_first(other)
        df.ix[0,'A'] = 45
        assert_frame_equal(result, df)

        # doc example
        df1 = DataFrame({'A' : [1., np.nan, 3., 5., np.nan],
                         'B' : [np.nan, 2., 3., np.nan, 6.]})

        df2 = DataFrame({'A' : [5., 2., 4., np.nan, 3., 7.],
                         'B' : [np.nan, np.nan, 3., 4., 6., 8.]})

        result = df1.combine_first(df2)
        expected = DataFrame({ 'A' : [1,2,3,5,3,7.], 'B' : [np.nan,2,3,4,6,8] })
        assert_frame_equal(result,expected)

        # GH3552, return object dtype with bools
        df1 = DataFrame([[np.nan, 3.,True], [-4.6, np.nan, True], [np.nan, 7., False]])
        df2 = DataFrame([[-42.6, np.nan, True], [-5., 1.6, False]], index=[1, 2])

        result = df1.combine_first(df2)[2]
        expected = Series([True,True,False])
        assert_series_equal(result,expected)

        # GH 3593, converting datetime64[ns] incorrecly
        df0 = DataFrame({"a":[datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)]})
        df1 = DataFrame({"a":[None, None, None]})
        df2 = df1.combine_first(df0)
        assert_frame_equal(df2,df0)

        df2 = df0.combine_first(df1)
        assert_frame_equal(df2,df0)

        df0 = DataFrame({"a":[datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)]})
        df1 = DataFrame({"a":[datetime(2000, 1, 2), None, None]})
        df2 = df1.combine_first(df0)
        result = df0.copy()
        result.iloc[0,:] = df1.iloc[0,:]
        assert_frame_equal(df2,result)

        df2 = df0.combine_first(df1)
        assert_frame_equal(df2,df0)

    def test_update(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other)

        expected = DataFrame([[1.5, nan, 3],
                              [3.6, 2, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 7.]])
        assert_frame_equal(df, expected)

    def test_update_dtypes(self):

        # gh 3016
        df = DataFrame([[1.,2.,False, True],[4.,5.,True,False]],
                       columns=['A','B','bool1','bool2'])

        other = DataFrame([[45,45]],index=[0],columns=['A','B'])
        df.update(other)

        expected = DataFrame([[45.,45.,False, True],[4.,5.,True,False]],
                             columns=['A','B','bool1','bool2'])
        assert_frame_equal(df, expected)

    def test_update_nooverwrite(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, overwrite=False)

        expected = DataFrame([[1.5, nan, 3],
                              [1.5, 2, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 3.]])
        assert_frame_equal(df, expected)

    def test_update_filtered(self):
        df = DataFrame([[1.5, nan, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[3.6, 2., np.nan],
                           [np.nan, np.nan, 7]], index=[1, 3])

        df.update(other, filter_func=lambda x: x > 2)

        expected = DataFrame([[1.5, nan, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 3],
                              [1.5, nan, 7.]])
        assert_frame_equal(df, expected)

    def test_update_raise(self):
        df = DataFrame([[1.5, 1, 3.],
                        [1.5, nan, 3.],
                        [1.5, nan, 3],
                        [1.5, nan, 3]])

        other = DataFrame([[2., nan],
                           [nan, 7]], index=[1, 3], columns=[1, 2])

        np.testing.assert_raises(Exception, df.update, *(other,),
                                 **{'raise_conflict': True})

    def test_update_from_non_df(self):
        d = {'a': Series([1, 2, 3, 4]), 'b': Series([5, 6, 7, 8])}
        df = DataFrame(d)

        d['a'] = Series([5, 6, 7, 8])
        df.update(d)

        expected = DataFrame(d)

        assert_frame_equal(df, expected)

        d = {'a': [1, 2, 3, 4], 'b': [5, 6, 7, 8]}
        df = DataFrame(d)

        d['a'] = [5, 6, 7, 8]
        df.update(d)

        expected = DataFrame(d)

        assert_frame_equal(df, expected)

    def test_combineAdd(self):
        # trivial
        comb = self.frame.combineAdd(self.frame)
        assert_frame_equal(comb, self.frame * 2)

        # more rigorous
        a = DataFrame([[1., nan, nan, 2., nan]],
                      columns=np.arange(5))
        b = DataFrame([[2., 3., nan, 2., 6., nan]],
                      columns=np.arange(6))
        expected = DataFrame([[3., 3., nan, 4., 6., nan]],
                             columns=np.arange(6))

        result = a.combineAdd(b)
        assert_frame_equal(result, expected)
        result2 = a.T.combineAdd(b.T)
        assert_frame_equal(result2, expected.T)

        expected2 = a.combine(b, operator.add, fill_value=0.)
        assert_frame_equal(expected, expected2)

        # corner cases
        comb = self.frame.combineAdd(self.empty)
        assert_frame_equal(comb, self.frame)

        comb = self.empty.combineAdd(self.frame)
        assert_frame_equal(comb, self.frame)

        # integer corner case
        df1 = DataFrame({'x': [5]})
        df2 = DataFrame({'x': [1]})
        df3 = DataFrame({'x': [6]})
        comb = df1.combineAdd(df2)
        assert_frame_equal(comb, df3)

        # mixed type GH2191
        df1 = DataFrame({'A': [1, 2], 'B': [3, 4]})
        df2 = DataFrame({'A': [1, 2], 'C': [5, 6]})
        rs = df1.combineAdd(df2)
        xp = DataFrame({'A': [2, 4], 'B': [3, 4.], 'C': [5, 6.]})
        assert_frame_equal(xp, rs)

        # TODO: test integer fill corner?

    def test_combineMult(self):
        # trivial
        comb = self.frame.combineMult(self.frame)

        assert_frame_equal(comb, self.frame ** 2)

        # corner cases
        comb = self.frame.combineMult(self.empty)
        assert_frame_equal(comb, self.frame)

        comb = self.empty.combineMult(self.frame)
        assert_frame_equal(comb, self.frame)

    def test_combine_generic(self):
        df1 = self.frame
        df2 = self.frame.ix[:-5, ['A', 'B', 'C']]

        combined = df1.combine(df2, np.add)
        combined2 = df2.combine(df1, np.add)
        self.assert_(combined['D'].isnull().all())
        self.assert_(combined2['D'].isnull().all())

        chunk = combined.ix[:-5, ['A', 'B', 'C']]
        chunk2 = combined2.ix[:-5, ['A', 'B', 'C']]

        exp = self.frame.ix[:-5, ['A', 'B', 'C']].reindex_like(chunk) * 2
        assert_frame_equal(chunk, exp)
        assert_frame_equal(chunk2, exp)

    def test_clip(self):
        median = self.frame.median().median()

        capped = self.frame.clip_upper(median)
        self.assert_(not (capped.values > median).any())

        floored = self.frame.clip_lower(median)
        self.assert_(not (floored.values < median).any())

        double = self.frame.clip(upper=median, lower=median)
        self.assert_(not (double.values != median).any())

    def test_dataframe_clip(self):

        # GH #2747
        df = DataFrame(np.random.randn(1000,2))

        for lb, ub in [(-1,1),(1,-1)]:
            clipped_df = df.clip(lb, ub)

            lb, ub = min(lb,ub), max(ub,lb)
            lb_mask = df.values <= lb
            ub_mask = df.values >= ub
            mask = ~lb_mask & ~ub_mask
            self.assert_((clipped_df.values[lb_mask] == lb).all() == True)
            self.assert_((clipped_df.values[ub_mask] == ub).all() == True)
            self.assert_((clipped_df.values[mask] == df.values[mask]).all() == True)

    def test_get_X_columns(self):
        # numeric and object columns

        df = DataFrame({'a': [1, 2, 3],
                        'b' : [True, False, True],
                        'c': ['foo', 'bar', 'baz'],
                        'd': [None, None, None],
                        'e': [3.14, 0.577, 2.773]})

        self.assert_(np.array_equal(df._get_numeric_data().columns,
                                    ['a', 'b', 'e']))

    def test_get_numeric_data(self):
        intname = np.dtype(np.int_).name
        floatname = np.dtype(np.float_).name
        datetime64name = np.dtype('M8[ns]').name
        objectname = np.dtype(np.object_).name

        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo', 'f' : Timestamp('20010102')},
                       index=np.arange(10))
        result = df.get_dtype_counts()
        expected = Series({'int64': 1, 'float64' : 1, datetime64name: 1, objectname : 1})
        result.sort()
        expected.sort()
        assert_series_equal(result, expected)

        df = DataFrame({'a': 1., 'b': 2, 'c': 'foo',
                        'd' : np.array([1.]*10,dtype='float32'),
                        'e' : np.array([1]*10,dtype='int32'),
                        'f' : np.array([1]*10,dtype='int16'),
                        'g' : Timestamp('20010102')},
                       index=np.arange(10))

        result = df._get_numeric_data()
        expected = df.ix[:, ['a', 'b','d','e','f']]
        assert_frame_equal(result, expected)

        only_obj = df.ix[:, ['c','g']]
        result = only_obj._get_numeric_data()
        expected = df.ix[:, []]
        assert_frame_equal(result, expected)

    def test_bool_describe_in_mixed_frame(self):
        df = DataFrame({
            'string_data': ['a', 'b', 'c', 'd', 'e'],
            'bool_data': [True, True, False, False, False],
            'int_data': [10, 20, 30, 40, 50],
        })

        # Boolean data and integer data is included in .describe() output, string data isn't
        self.assert_(np.array_equal(df.describe().columns, ['bool_data', 'int_data']))

        bool_describe = df.describe()['bool_data']

        # Both the min and the max values should stay booleans
        self.assert_(bool_describe['min'].dtype == np.bool_)
        self.assert_(bool_describe['max'].dtype == np.bool_)

        self.assert_(bool_describe['min'] == False)
        self.assert_(bool_describe['max'] == True)

        # For numeric operations, like mean or median, the values True/False are cast to
        # the integer values 1 and 0
        assert_almost_equal(bool_describe['mean'], 0.4)
        assert_almost_equal(bool_describe['50%'], 0)

    def test_count(self):
        f = lambda s: notnull(s).sum()
        self._check_stat_op('count', f,
                            has_skipna=False,
                            has_numeric_only=True,
                            check_dtypes=False)

        # corner case
        frame = DataFrame()
        ct1 = frame.count(1)
        self.assert_(isinstance(ct1, Series))

        ct2 = frame.count(0)
        self.assert_(isinstance(ct2, Series))

        # GH #423
        df = DataFrame(index=range(10))
        result = df.count(1)
        expected = Series(0, index=df.index)
        assert_series_equal(result, expected)

        df = DataFrame(columns=range(10))
        result = df.count(0)
        expected = Series(0, index=df.columns)
        assert_series_equal(result, expected)

        df = DataFrame()
        result = df.count()
        expected = Series(0, index=[])
        assert_series_equal(result, expected)

    def test_sum(self):
        self._check_stat_op('sum', np.sum, has_numeric_only=True)

    def test_sum_mixed_numeric(self):
        raise nose.SkipTest
        # mixed types
        self._check_stat_op('sum', np.sum, frame = self.mixed_float, has_numeric_only=True)

    def test_stat_operators_attempt_obj_array(self):
        data = {
            'a': [-0.00049987540199591344, -0.0016467257772919831,
                  0.00067695870775883013],
            'b': [-0, -0, 0.0],
            'c': [0.00031111847529610595, 0.0014902627951905339,
                  -0.00094099200035979691]
        }
        df1 = DataFrame(data, index=['foo', 'bar', 'baz'],
                        dtype='O')
        methods = ['sum', 'mean', 'prod', 'var', 'std', 'skew', 'min', 'max']

        # GH #676
        df2 = DataFrame({0: [np.nan, 2], 1: [np.nan, 3],
                        2: [np.nan, 4]}, dtype=object)

        for df in [df1, df2]:
            for meth in methods:
                self.assert_(df.values.dtype == np.object_)
                result = getattr(df, meth)(1)
                expected = getattr(df.astype('f8'), meth)(1)
                assert_series_equal(result, expected)

    def test_mean(self):
        self._check_stat_op('mean', np.mean)

    def test_product(self):
        self._check_stat_op('product', np.prod)

    def test_median(self):
        def wrapper(x):
            if isnull(x).any():
                return np.nan
            return np.median(x)

        self._check_stat_op('median', wrapper)

    def test_min(self):
        self._check_stat_op('min', np.min)
        self._check_stat_op('min', np.min, frame=self.intframe)

    def test_cummin(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cummin = self.tsframe.cummin()
        expected = self.tsframe.apply(Series.cummin)
        assert_frame_equal(cummin, expected)

        # axis = 1
        cummin = self.tsframe.cummin(axis=1)
        expected = self.tsframe.apply(Series.cummin, axis=1)
        assert_frame_equal(cummin, expected)

        # works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cummin()

        # fix issue
        cummin_xs = self.tsframe.cummin(axis=1)
        self.assertEqual(np.shape(cummin_xs), np.shape(self.tsframe))

    def test_cummax(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cummax = self.tsframe.cummax()
        expected = self.tsframe.apply(Series.cummax)
        assert_frame_equal(cummax, expected)

        # axis = 1
        cummax = self.tsframe.cummax(axis=1)
        expected = self.tsframe.apply(Series.cummax, axis=1)
        assert_frame_equal(cummax, expected)

        # works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cummax()

        # fix issue
        cummax_xs = self.tsframe.cummax(axis=1)
        self.assertEqual(np.shape(cummax_xs), np.shape(self.tsframe))

    def test_max(self):
        self._check_stat_op('max', np.max)
        self._check_stat_op('max', np.max, frame=self.intframe)

    def test_mad(self):
        f = lambda x: np.abs(x - x.mean()).mean()
        self._check_stat_op('mad', f)

    def test_var_std(self):
        alt = lambda x: np.var(x, ddof=1)
        self._check_stat_op('var', alt)

        alt = lambda x: np.std(x, ddof=1)
        self._check_stat_op('std', alt)

        result = self.tsframe.std(ddof=4)
        expected = self.tsframe.apply(lambda x: x.std(ddof=4))
        assert_almost_equal(result, expected)

        result = self.tsframe.var(ddof=4)
        expected = self.tsframe.apply(lambda x: x.var(ddof=4))
        assert_almost_equal(result, expected)

        arr = np.repeat(np.random.random((1, 1000)), 1000, 0)
        result = nanops.nanvar(arr, axis=0)
        self.assertFalse((result < 0).any())
        if nanops._USE_BOTTLENECK:
            nanops._USE_BOTTLENECK = False
            result = nanops.nanvar(arr, axis=0)
            self.assertFalse((result < 0).any())
            nanops._USE_BOTTLENECK = True

    def test_skew(self):
        _skip_if_no_scipy()
        from scipy.stats import skew

        def alt(x):
            if len(x) < 3:
                return np.nan
            return skew(x, bias=False)

        self._check_stat_op('skew', alt)

    def test_kurt(self):
        _skip_if_no_scipy()

        from scipy.stats import kurtosis

        def alt(x):
            if len(x) < 4:
                return np.nan
            return kurtosis(x, bias=False)

        self._check_stat_op('kurt', alt)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        df = DataFrame(np.random.randn(6, 3), index=index)
        assert_series_equal(df.kurt(), df.kurt(level=0).xs('bar'))

    def _check_stat_op(self, name, alternative, frame=None, has_skipna=True,
                       has_numeric_only=False, check_dtypes=True):
        if frame is None:
            frame = self.frame
            # set some NAs
            frame.ix[5:10] = np.nan
            frame.ix[15:20, -2:] = np.nan

        f = getattr(frame, name)

        if not ('max' in name or 'min' in name or 'count' in name):
            df = DataFrame({'b': date_range('1/1/2001', periods=2)})
            _f = getattr(df, name)
            print (df)
            self.assertFalse(len(_f()))

            df['a'] = range(len(df))
            self.assert_(len(getattr(df, name)()))

        if has_skipna:
            def skipna_wrapper(x):
                nona = x.dropna().values
                if len(nona) == 0:
                    return np.nan
                return alternative(nona)

            def wrapper(x):
                return alternative(x.values)

            result0 = f(axis=0, skipna=False)
            result1 = f(axis=1, skipna=False)
            assert_series_equal(result0, frame.apply(wrapper))
            assert_series_equal(result1, frame.apply(wrapper, axis=1),
                                check_dtype=False)  # HACK: win32
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        result0 = f(axis=0)
        result1 = f(axis=1)
        assert_series_equal(result0, frame.apply(skipna_wrapper))
        assert_series_equal(result1, frame.apply(skipna_wrapper, axis=1),
                            check_dtype=False)

        # check dtypes
        if check_dtypes:
            lcd_dtype = frame.values.dtype
            self.assert_(lcd_dtype == result0.dtype)
            self.assert_(lcd_dtype == result1.dtype)

        # result = f(axis=1)
        # comp = frame.apply(alternative, axis=1).reindex(result.index)
        # assert_series_equal(result, comp)

        self.assertRaises(Exception, f, axis=2)

        # make sure works on mixed-type frame
        getattr(self.mixed_frame, name)(axis=0)
        getattr(self.mixed_frame, name)(axis=1)

        if has_numeric_only:
            getattr(self.mixed_frame, name)(axis=0, numeric_only=True)
            getattr(self.mixed_frame, name)(axis=1, numeric_only=True)
            getattr(self.frame, name)(axis=0, numeric_only=False)
            getattr(self.frame, name)(axis=1, numeric_only=False)

        # all NA case
        if has_skipna:
            all_na = self.frame * np.NaN
            r0 = getattr(all_na, name)(axis=0)
            r1 = getattr(all_na, name)(axis=1)
            self.assert_(np.isnan(r0).all())
            self.assert_(np.isnan(r1).all())

    def test_sum_corner(self):
        axis0 = self.empty.sum(0)
        axis1 = self.empty.sum(1)
        self.assert_(isinstance(axis0, Series))
        self.assert_(isinstance(axis1, Series))
        self.assertEquals(len(axis0), 0)
        self.assertEquals(len(axis1), 0)

    def test_sum_object(self):
        values = self.frame.values.astype(int)
        frame = DataFrame(values, index=self.frame.index,
                          columns=self.frame.columns)
        deltas = frame * timedelta(1)
        deltas.sum()

    def test_sum_bool(self):
        # ensure this works, bug report
        bools = np.isnan(self.frame)
        bools.sum(1)
        bools.sum(0)

    def test_mean_corner(self):
        # unit test when have object data
        the_mean = self.mixed_frame.mean(axis=0)
        the_sum = self.mixed_frame.sum(axis=0, numeric_only=True)
        self.assert_(the_sum.index.equals(the_mean.index))
        self.assert_(len(the_mean.index) < len(self.mixed_frame.columns))

        # xs sum mixed type, just want to know it works...
        the_mean = self.mixed_frame.mean(axis=1)
        the_sum = self.mixed_frame.sum(axis=1, numeric_only=True)
        self.assert_(the_sum.index.equals(the_mean.index))

        # take mean of boolean column
        self.frame['bool'] = self.frame['A'] > 0
        means = self.frame.mean(0)
        self.assertEqual(means['bool'], self.frame['bool'].values.mean())

    def test_stats_mixed_type(self):
        # don't blow up
        self.mixed_frame.std(1)
        self.mixed_frame.var(1)
        self.mixed_frame.mean(1)
        self.mixed_frame.skew(1)

    def test_median_corner(self):
        def wrapper(x):
            if isnull(x).any():
                return np.nan
            return np.median(x)

        self._check_stat_op('median', wrapper, frame=self.intframe, check_dtypes=False)

    def test_quantile(self):
        from pandas.compat.scipy import scoreatpercentile

        q = self.tsframe.quantile(0.1, axis=0)
        self.assertEqual(q['A'], scoreatpercentile(self.tsframe['A'], 10))
        q = self.tsframe.quantile(0.9, axis=1)
        q = self.intframe.quantile(0.1)
        self.assertEqual(q['A'], scoreatpercentile(self.intframe['A'], 10))

        # test degenerate case
        q = DataFrame({'x': [], 'y': []}).quantile(0.1, axis=0)
        assert(np.isnan(q['x']) and np.isnan(q['y']))

        # non-numeric exclusion
        df = DataFrame({'col1':['A','A','B','B'], 'col2':[1,2,3,4]})
        rs = df.quantile(0.5)
        xp = df.median()
        assert_series_equal(rs, xp)

    def test_cumsum(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cumsum = self.tsframe.cumsum()
        expected = self.tsframe.apply(Series.cumsum)
        assert_frame_equal(cumsum, expected)

        # axis = 1
        cumsum = self.tsframe.cumsum(axis=1)
        expected = self.tsframe.apply(Series.cumsum, axis=1)
        assert_frame_equal(cumsum, expected)

        # works
        df = DataFrame({'A': np.arange(20)}, index=np.arange(20))
        result = df.cumsum()

        # fix issue
        cumsum_xs = self.tsframe.cumsum(axis=1)
        self.assertEqual(np.shape(cumsum_xs), np.shape(self.tsframe))

    def test_cumprod(self):
        self.tsframe.ix[5:10, 0] = nan
        self.tsframe.ix[10:15, 1] = nan
        self.tsframe.ix[15:, 2] = nan

        # axis = 0
        cumprod = self.tsframe.cumprod()
        expected = self.tsframe.apply(Series.cumprod)
        assert_frame_equal(cumprod, expected)

        # axis = 1
        cumprod = self.tsframe.cumprod(axis=1)
        expected = self.tsframe.apply(Series.cumprod, axis=1)
        assert_frame_equal(cumprod, expected)

        # fix issue
        cumprod_xs = self.tsframe.cumprod(axis=1)
        self.assertEqual(np.shape(cumprod_xs), np.shape(self.tsframe))

        # ints
        df = self.tsframe.fillna(0).astype(int)
        df.cumprod(0)
        df.cumprod(1)

        # ints32
        df = self.tsframe.fillna(0).astype(np.int32)
        df.cumprod(0)
        df.cumprod(1)

    def test_rank(self):
        from pandas.compat.scipy import rankdata

        self.frame['A'][::2] = np.nan
        self.frame['B'][::3] = np.nan
        self.frame['C'][::4] = np.nan
        self.frame['D'][::5] = np.nan

        ranks0 = self.frame.rank()
        ranks1 = self.frame.rank(1)
        mask = np.isnan(self.frame.values)

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fvals)
        exp0[mask] = np.nan

        exp1 = np.apply_along_axis(rankdata, 1, fvals)
        exp1[mask] = np.nan

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # integers
        df = DataFrame(np.random.randint(0, 5, size=40).reshape((10, 4)))

        result = df.rank()
        exp = df.astype(float).rank()
        assert_frame_equal(result, exp)

        result = df.rank(1)
        exp = df.astype(float).rank(1)
        assert_frame_equal(result, exp)

    def test_rank2(self):
        from datetime import datetime

        df = DataFrame([['b', 'c', 'a'], ['a', 'c', 'b']])
        expected = DataFrame([[2.0, 3.0, 1.0], [1, 3, 2]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        expected = DataFrame([[2.0, 1.5, 1.0], [1, 1.5, 2]])
        result = df.rank(0, numeric_only=False)
        assert_frame_equal(result, expected)

        df = DataFrame([['b', np.nan, 'a'], ['a', 'c', 'b']])
        expected = DataFrame([[2.0, nan, 1.0], [1.0, 3.0, 2.0]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        expected = DataFrame([[2.0, nan, 1.0], [1.0, 1.0, 2.0]])
        result = df.rank(0, numeric_only=False)
        assert_frame_equal(result, expected)

        # f7u12, this does not work without extensive workaround
        data = [[datetime(2001, 1, 5), nan, datetime(2001, 1, 2)],
                [datetime(2000, 1, 2), datetime(2000, 1, 3),
                 datetime(2000, 1, 1)]]
        df = DataFrame(data)
        expected = DataFrame([[2., nan, 1.],
                              [2., 3., 1.]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        # mixed-type frames
        self.mixed_frame['foo'] = datetime.now()
        result = self.mixed_frame.rank(1)
        expected = self.mixed_frame.rank(1, numeric_only=True)
        assert_frame_equal(result, expected)

    def test_rank_na_option(self):
        from pandas.compat.scipy import rankdata

        self.frame['A'][::2] = np.nan
        self.frame['B'][::3] = np.nan
        self.frame['C'][::4] = np.nan
        self.frame['D'][::5] = np.nan

        # bottom
        ranks0 = self.frame.rank(na_option='bottom')
        ranks1 = self.frame.rank(1, na_option='bottom')

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fvals)
        exp1 = np.apply_along_axis(rankdata, 1, fvals)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # top
        ranks0 = self.frame.rank(na_option='top')
        ranks1 = self.frame.rank(1, na_option='top')

        fval0 = self.frame.fillna((self.frame.min() - 1).to_dict()).values
        fval1 = self.frame.T
        fval1 = fval1.fillna((fval1.min() - 1).to_dict()).T
        fval1 = fval1.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fval0)
        exp1 = np.apply_along_axis(rankdata, 1, fval1)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # descending

        # bottom
        ranks0 = self.frame.rank(na_option='top', ascending=False)
        ranks1 = self.frame.rank(1, na_option='top', ascending=False)

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, -fvals)
        exp1 = np.apply_along_axis(rankdata, 1, -fvals)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

        # descending

        # top
        ranks0 = self.frame.rank(na_option='bottom', ascending=False)
        ranks1 = self.frame.rank(1, na_option='bottom', ascending=False)

        fval0 = self.frame.fillna((self.frame.min() - 1).to_dict()).values
        fval1 = self.frame.T
        fval1 = fval1.fillna((fval1.min() - 1).to_dict()).T
        fval1 = fval1.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, -fval0)
        exp1 = np.apply_along_axis(rankdata, 1, -fval1)

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

    def test_describe(self):
        desc = self.tsframe.describe()
        desc = self.mixed_frame.describe()
        desc = self.frame.describe()

    def test_describe_percentiles(self):
        desc = self.frame.describe(percentile_width=50)
        assert '75%' in desc.index
        assert '25%' in desc.index

        desc = self.frame.describe(percentile_width=95)
        assert '97.5%' in desc.index
        assert '2.5%' in desc.index

    def test_describe_no_numeric(self):
        df = DataFrame({'A': ['foo', 'foo', 'bar'] * 8,
                        'B': ['a', 'b', 'c', 'd'] * 6})
        desc = df.describe()
        expected = DataFrame(dict((k, v.describe())
                                  for k, v in df.iteritems()),
                             columns=df.columns)
        assert_frame_equal(desc, expected)

        df = DataFrame({'time': self.tsframe.index})
        desc = df.describe()
        assert(desc.time['first'] == min(self.tsframe.index))

    def test_describe_empty_int_columns(self):
        df = DataFrame([[0, 1], [1, 2]])
        desc = df[df[0] < 0].describe()  # works
        assert_series_equal(desc.xs('count'),
                            Series([0, 0], dtype=float, name='count'))
        self.assert_(isnull(desc.ix[1:]).all().all())

    def test_get_axis_etc(self):
        f = self.frame

        self.assertEquals(f._get_axis_number(0), 0)
        self.assertEquals(f._get_axis_number(1), 1)
        self.assertEquals(f._get_axis_name(0), 'index')
        self.assertEquals(f._get_axis_name(1), 'columns')

        self.assert_(f._get_axis(0) is f.index)
        self.assert_(f._get_axis(1) is f.columns)
        self.assertRaises(Exception, f._get_axis_number, 2)

    def test_axis_aliases(self):

        f = self.frame

        # reg name
        expected = f.sum(axis=0)
        result = f.sum(axis='index')
        assert_series_equal(result, expected)

        expected = f.sum(axis=1)
        result = f.sum(axis='columns')
        assert_series_equal(result, expected)

    def test_combine_first_mixed(self):
        a = Series(['a', 'b'], index=range(2))
        b = Series(range(2), index=range(2))
        f = DataFrame({'A': a, 'B': b})

        a = Series(['a', 'b'], index=range(5, 7))
        b = Series(range(2), index=range(5, 7))
        g = DataFrame({'A': a, 'B': b})

        combined = f.combine_first(g)

    def test_more_asMatrix(self):
        values = self.mixed_frame.as_matrix()
        self.assertEqual(values.shape[1], len(self.mixed_frame.columns))

    def test_reindex_boolean(self):
        frame = DataFrame(np.ones((10, 2), dtype=bool),
                          index=np.arange(0, 20, 2),
                          columns=[0, 2])

        reindexed = frame.reindex(np.arange(10))
        self.assert_(reindexed.values.dtype == np.object_)
        self.assert_(isnull(reindexed[0][1]))

        reindexed = frame.reindex(columns=range(3))
        self.assert_(reindexed.values.dtype == np.object_)
        self.assert_(isnull(reindexed[1]).all())

    def test_reindex_objects(self):
        reindexed = self.mixed_frame.reindex(columns=['foo', 'A', 'B'])
        self.assert_('foo' in reindexed)

        reindexed = self.mixed_frame.reindex(columns=['A', 'B'])
        self.assert_('foo' not in reindexed)

    def test_reindex_corner(self):
        index = Index(['a', 'b', 'c'])
        dm = self.empty.reindex(index=[1, 2, 3])
        reindexed = dm.reindex(columns=index)
        self.assert_(reindexed.columns.equals(index))

        # ints are weird

        smaller = self.intframe.reindex(columns=['A', 'B', 'E'])
        self.assert_(smaller['E'].dtype == np.float64)

    def test_reindex_axis(self):
        cols = ['A', 'B', 'E']
        reindexed1 = self.intframe.reindex_axis(cols, axis=1)
        reindexed2 = self.intframe.reindex(columns=cols)
        assert_frame_equal(reindexed1, reindexed2)

        rows = self.intframe.index[0:5]
        reindexed1 = self.intframe.reindex_axis(rows, axis=0)
        reindexed2 = self.intframe.reindex(index=rows)
        assert_frame_equal(reindexed1, reindexed2)

        self.assertRaises(ValueError, self.intframe.reindex_axis, rows, axis=2)

        # no-op case
        cols = self.frame.columns.copy()
        newFrame = self.frame.reindex_axis(cols, axis=1)
        assert_frame_equal(newFrame, self.frame)

    def test_reindex_with_nans(self):
        df = DataFrame([[1, 2], [3, 4], [np.nan, np.nan], [7, 8], [9, 10]],
                       columns=['a', 'b'],
                       index=[100.0, 101.0, np.nan, 102.0, 103.0])

        result = df.reindex(index=[101.0, 102.0, 103.0])
        expected = df.ix[[1, 3, 4]]
        assert_frame_equal(result, expected)

        result = df.reindex(index=[103.0])
        expected = df.ix[[4]]
        assert_frame_equal(result, expected)

        result = df.reindex(index=[101.0])
        expected = df.ix[[1]]
        assert_frame_equal(result, expected)

    def test_reindex_multi(self):
        df = DataFrame(np.random.randn(3, 3))

        result = df.reindex(range(4), range(4))
        expected = df.reindex(range(4)).reindex(columns=range(4))

        assert_frame_equal(result, expected)

        df = DataFrame(np.random.randint(0, 10, (3, 3)))

        result = df.reindex(range(4), range(4))
        expected = df.reindex(range(4)).reindex(columns=range(4))

        assert_frame_equal(result, expected)

        df = DataFrame(np.random.randint(0, 10, (3, 3)))

        result = df.reindex(range(2), range(2))
        expected = df.reindex(range(2)).reindex(columns=range(2))

        assert_frame_equal(result, expected)

        df = DataFrame(np.random.randn(5, 3) + 1j, columns=['a', 'b', 'c'])

        result = df.reindex(index=[0, 1], columns=['a', 'b'])
        expected = df.reindex([0, 1]).reindex(columns=['a', 'b'])

        assert_frame_equal(result, expected)

    def test_rename_objects(self):
        renamed = self.mixed_frame.rename(columns=str.upper)
        self.assert_('FOO' in renamed)
        self.assert_('foo' not in renamed)

    def test_fill_corner(self):
        self.mixed_frame['foo'][5:20] = nan
        self.mixed_frame['A'][-10:] = nan

        filled = self.mixed_frame.fillna(value=0)
        self.assert_((filled['foo'][5:20] == 0).all())
        del self.mixed_frame['foo']

        empty_float = self.frame.reindex(columns=[])
        result = empty_float.fillna(value=0)

    def test_count_objects(self):
        dm = DataFrame(self.mixed_frame._series)
        df = DataFrame(self.mixed_frame._series)

        tm.assert_series_equal(dm.count(), df.count())
        tm.assert_series_equal(dm.count(1), df.count(1))

    def test_cumsum_corner(self):
        dm = DataFrame(np.arange(20).reshape(4, 5),
                       index=range(4), columns=range(5))
        result = dm.cumsum()

    #----------------------------------------------------------------------
    # Stacking / unstacking

    def test_stack_unstack(self):
        stacked = self.frame.stack()
        stacked_df = DataFrame({'foo': stacked, 'bar': stacked})

        unstacked = stacked.unstack()
        unstacked_df = stacked_df.unstack()

        assert_frame_equal(unstacked, self.frame)
        assert_frame_equal(unstacked_df['bar'], self.frame)

        unstacked_cols = stacked.unstack(0)
        unstacked_cols_df = stacked_df.unstack(0)
        assert_frame_equal(unstacked_cols.T, self.frame)
        assert_frame_equal(unstacked_cols_df['bar'].T, self.frame)

    def test_unstack_bool(self):
        df = DataFrame([False, False],
                       index=MultiIndex.from_arrays([['a', 'b'], ['c', 'l']]),
                       columns=['col'])
        rs = df.unstack()
        xp = DataFrame(np.array([[False, np.nan], [np.nan, False]],
                                dtype=object),
                       index=['a', 'b'],
                       columns=MultiIndex.from_arrays([['col', 'col'],
                                                       ['c', 'l']]))
        assert_frame_equal(rs, xp)

    def test_unstack_to_series(self):
        # check reversibility
        data = self.frame.unstack()

        self.assertTrue(isinstance(data, Series))
        undo = data.unstack().T
        assert_frame_equal(undo, self.frame)

        # check NA handling
        data = DataFrame({'x': [1, 2, np.NaN], 'y': [3.0, 4, np.NaN]})
        data.index = Index(['a', 'b', 'c'])
        result = data.unstack()

        midx = MultiIndex(levels=[['x', 'y'], ['a', 'b', 'c']],
                          labels=[[0, 0, 0, 1, 1, 1], [0, 1, 2, 0, 1, 2]])
        expected = Series([1, 2, np.NaN, 3, 4, np.NaN], index=midx)

        assert_series_equal(result, expected)

        # check composability of unstack
        old_data = data.copy()
        for _ in xrange(4):
            data = data.unstack()
        assert_frame_equal(old_data, data)

    def test_unstack_dtypes(self):

        # GH 2929
        rows = [[1, 1, 3, 4],
                [1, 2, 3, 4],
                [2, 1, 3, 4],
                [2, 2, 3, 4]]

        df = DataFrame(rows, columns=list('ABCD'))
        result = df.get_dtype_counts()
        expected = Series({'int64' : 4})
        assert_series_equal(result, expected)

        # single dtype
        df2 = df.set_index(['A','B'])
        df3 = df2.unstack('B')
        result = df3.get_dtype_counts()
        expected = Series({'int64' : 4})
        assert_series_equal(result, expected)

        # mixed
        df2 = df.set_index(['A','B'])
        df2['C'] = 3.
        df3 = df2.unstack('B')
        result = df3.get_dtype_counts()
        expected = Series({'int64' : 2, 'float64' : 2})
        assert_series_equal(result, expected)

        df2['D'] = 'foo'
        df3 = df2.unstack('B')
        result = df3.get_dtype_counts()
        expected = Series({'float64' : 2, 'object' : 2})
        assert_series_equal(result, expected)


    def test_reset_index(self):
        stacked = self.frame.stack()[::2]
        stacked = DataFrame({'foo': stacked, 'bar': stacked})

        names = ['first', 'second']
        stacked.index.names = names
        deleveled = stacked.reset_index()
        for i, (lev, lab) in enumerate(zip(stacked.index.levels,
                                           stacked.index.labels)):
            values = lev.take(lab)
            name = names[i]
            assert_almost_equal(values, deleveled[name])

        stacked.index.names = [None, None]
        deleveled2 = stacked.reset_index()
        self.assert_(np.array_equal(deleveled['first'],
                                    deleveled2['level_0']))
        self.assert_(np.array_equal(deleveled['second'],
                                    deleveled2['level_1']))

        # default name assigned
        rdf = self.frame.reset_index()
        self.assert_(np.array_equal(rdf['index'], self.frame.index.values))

        # default name assigned, corner case
        df = self.frame.copy()
        df['index'] = 'foo'
        rdf = df.reset_index()
        self.assert_(np.array_equal(rdf['level_0'], self.frame.index.values))

        # but this is ok
        self.frame.index.name = 'index'
        deleveled = self.frame.reset_index()
        self.assert_(np.array_equal(deleveled['index'],
                                    self.frame.index.values))
        self.assert_(np.array_equal(deleveled.index,
                                    np.arange(len(deleveled))))

        # preserve column names
        self.frame.columns.name = 'columns'
        resetted = self.frame.reset_index()
        self.assertEqual(resetted.columns.name, 'columns')

        # only remove certain columns
        frame = self.frame.reset_index().set_index(['index', 'A', 'B'])
        rs = frame.reset_index(['A', 'B'])

        assert_frame_equal(rs, self.frame, check_names=False)  # TODO should reset_index check_names ?

        rs = frame.reset_index(['index', 'A', 'B'])
        assert_frame_equal(rs, self.frame.reset_index(), check_names=False)

        rs = frame.reset_index(['index', 'A', 'B'])
        assert_frame_equal(rs, self.frame.reset_index(), check_names=False)

        rs = frame.reset_index('A')
        xp = self.frame.reset_index().set_index(['index', 'B'])
        assert_frame_equal(rs, xp, check_names=False)

        # test resetting in place
        df = self.frame.copy()
        resetted = self.frame.reset_index()
        df.reset_index(inplace=True)
        assert_frame_equal(df, resetted, check_names=False)

        frame = self.frame.reset_index().set_index(['index', 'A', 'B'])
        rs = frame.reset_index('A', drop=True)
        xp = self.frame.copy()
        del xp['A']
        xp = xp.set_index(['B'], append=True)
        assert_frame_equal(rs, xp, check_names=False)

    def test_reset_index_right_dtype(self):
        time = np.arange(0.0, 10, np.sqrt(2) / 2)
        s1 = Series((9.81 * time ** 2) / 2,
                    index=Index(time, name='time'),
                    name='speed')
        df = DataFrame(s1)

        resetted = s1.reset_index()
        self.assert_(resetted['time'].dtype == np.float64)

        resetted = df.reset_index()
        self.assert_(resetted['time'].dtype == np.float64)

    def test_reset_index_multiindex_col(self):
        vals = np.random.randn(3, 3).astype(object)
        idx = ['x', 'y', 'z']
        full = np.hstack(([[x] for x in idx], vals))
        df = DataFrame(vals, Index(idx, name='a'),
                       columns=[['b', 'b', 'c'], ['mean', 'median', 'mean']])
        rs = df.reset_index()
        xp = DataFrame(full, columns=[['a', 'b', 'b', 'c'],
                                      ['', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index(col_fill=None)
        xp = DataFrame(full, columns=[['a', 'b', 'b', 'c'],
                                      ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index(col_level=1, col_fill='blah')
        xp = DataFrame(full, columns=[['blah', 'b', 'b', 'c'],
                                      ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        df = DataFrame(vals,
                       MultiIndex.from_arrays([[0, 1, 2], ['x', 'y', 'z']],
                                              names=['d', 'a']),
                       columns=[['b', 'b', 'c'], ['mean', 'median', 'mean']])
        rs = df.reset_index('a', )
        xp = DataFrame(full, Index([0, 1, 2], name='d'),
                       columns=[['a', 'b', 'b', 'c'],
                                ['', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index('a', col_fill=None)
        xp = DataFrame(full, Index(range(3), name='d'),
                       columns=[['a', 'b', 'b', 'c'],
                                ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

        rs = df.reset_index('a', col_fill='blah', col_level=1)
        xp = DataFrame(full, Index(range(3), name='d'),
                       columns=[['blah', 'b', 'b', 'c'],
                                ['a', 'mean', 'median', 'mean']])
        assert_frame_equal(rs, xp)

    #----------------------------------------------------------------------
    # Tests to cope with refactored internals
    def test_as_matrix_numeric_cols(self):
        self.frame['foo'] = 'bar'

        values = self.frame.as_matrix(['A', 'B', 'C', 'D'])
        self.assert_(values.dtype == np.float64)

    def test_as_matrix_lcd(self):

        # mixed lcd
        values = self.mixed_float.as_matrix(['A', 'B', 'C', 'D'])
        self.assert_(values.dtype == np.float64)

        values = self.mixed_float.as_matrix(['A', 'B', 'C' ])
        self.assert_(values.dtype == np.float32)

        values = self.mixed_float.as_matrix(['C'])
        self.assert_(values.dtype == np.float16)

        values = self.mixed_int.as_matrix(['A','B','C','D'])
        self.assert_(values.dtype == np.int64)

        values = self.mixed_int.as_matrix(['A','D'])
        self.assert_(values.dtype == np.int64)

        # guess all ints are cast to uints....
        values = self.mixed_int.as_matrix(['A','B','C'])
        self.assert_(values.dtype == np.int64)

        values = self.mixed_int.as_matrix(['A','C'])
        self.assert_(values.dtype == np.int32)

        values = self.mixed_int.as_matrix(['C','D'])
        self.assert_(values.dtype == np.int64)

        values = self.mixed_int.as_matrix(['A'])
        self.assert_(values.dtype == np.int32)

        values = self.mixed_int.as_matrix(['C'])
        self.assert_(values.dtype == np.uint8)

    def test_constructor_with_convert(self):
        # this is actually mostly a test of lib.maybe_convert_objects
        # #2845
        df = DataFrame({'A' : [2**63-1] })
        result = df['A']
        expected = Series(np.asarray([2**63-1], np.int64))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [2**63] })
        result = df['A']
        expected = Series(np.asarray([2**63], np.object_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [datetime(2005, 1, 1), True] })
        result = df['A']
        expected = Series(np.asarray([datetime(2005, 1, 1), True], np.object_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [None, 1] })
        result = df['A']
        expected = Series(np.asarray([np.nan, 1], np.float_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0, 2] })
        result = df['A']
        expected = Series(np.asarray([1.0, 2], np.float_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, 3] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, 3], np.complex_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, 3.0] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, 3.0], np.complex_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, True] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, True], np.object_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0, None] })
        result = df['A']
        expected = Series(np.asarray([1.0, np.nan], np.float_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [1.0+2.0j, None] })
        result = df['A']
        expected = Series(np.asarray([1.0+2.0j, np.nan], np.complex_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [2.0, 1, True, None] })
        result = df['A']
        expected = Series(np.asarray([2.0, 1, True, None], np.object_))
        assert_series_equal(result, expected)

        df = DataFrame({'A' : [2.0, 1, datetime(2006, 1, 1), None] })
        result = df['A']
        expected = Series(np.asarray([2.0, 1, datetime(2006, 1, 1),
                                      None], np.object_))
        assert_series_equal(result, expected)

    def test_constructor_frame_copy(self):
        cop = DataFrame(self.frame, copy=True)
        cop['A'] = 5
        self.assert_((cop['A'] == 5).all())
        self.assert_(not (self.frame['A'] == 5).all())

    def test_constructor_ndarray_copy(self):
        df = DataFrame(self.frame.values)

        self.frame.values[5] = 5
        self.assert_((df.values[5] == 5).all())

        df = DataFrame(self.frame.values, copy=True)
        self.frame.values[6] = 6
        self.assert_(not (df.values[6] == 6).all())

    def test_constructor_series_copy(self):
        series = self.frame._series

        df = DataFrame({'A': series['A']})
        df['A'][:] = 5

        self.assert_(not (series['A'] == 5).all())

    def test_assign_columns(self):
        self.frame['hi'] = 'there'

        frame = self.frame.copy()
        frame.columns = ['foo', 'bar', 'baz', 'quux', 'foo2']
        assert_series_equal(self.frame['C'], frame['baz'])
        assert_series_equal(self.frame['hi'], frame['foo2'])

    def test_columns_with_dups(self):

        # GH 3468 related

        # basic
        df = DataFrame([[1,2]], columns=['a','a'])
        df.columns = ['a','a.1']
        str(df)
        expected = DataFrame([[1,2]], columns=['a','a.1'])
        assert_frame_equal(df, expected)

        df = DataFrame([[1,2,3]], columns=['b','a','a'])
        df.columns = ['b','a','a.1']
        str(df)
        expected = DataFrame([[1,2,3]], columns=['b','a','a.1'])
        assert_frame_equal(df, expected)

        # with a dup index
        df = DataFrame([[1,2]], columns=['a','a'])
        df.columns = ['b','b']
        str(df)
        expected = DataFrame([[1,2]], columns=['b','b'])
        assert_frame_equal(df, expected)

        # multi-dtype
        df = DataFrame([[1,2,1.,2.,3.,'foo','bar']], columns=['a','a','b','b','d','c','c'])
        df.columns = list('ABCDEFG')
        str(df)
        expected = DataFrame([[1,2,1.,2.,3.,'foo','bar']], columns=list('ABCDEFG'))
        assert_frame_equal(df, expected)

        # this is an error because we cannot disambiguate the dup columns
        self.assertRaises(Exception, lambda x: DataFrame([[1,2,'foo','bar']], columns=['a','a','a','a']))

        # dups across blocks
        df_float  = DataFrame(np.random.randn(10, 3),dtype='float64')
        df_int    = DataFrame(np.random.randn(10, 3),dtype='int64')
        df_bool   = DataFrame(True,index=df_float.index,columns=df_float.columns)
        df_object = DataFrame('foo',index=df_float.index,columns=df_float.columns)
        df_dt     = DataFrame(Timestamp('20010101'),index=df_float.index,columns=df_float.columns)
        df        = pan.concat([ df_float, df_int, df_bool, df_object, df_dt ], axis=1)

        result = df._data._set_ref_locs()
        self.assert_(len(result) == len(df.columns))

        # testing iget
        for i in range(len(df.columns)):
             df.iloc[:,i]

        # dup columns across dtype GH 2079/2194
        vals = [[1, -1, 2.], [2, -2, 3.]]
        rs = DataFrame(vals, columns=['A', 'A', 'B'])
        xp = DataFrame(vals)
        xp.columns = ['A', 'A', 'B']
        assert_frame_equal(rs, xp)

    def test_insert_column_bug_4032(self):

        # GH4032, inserting a column and renaming causing errors
        df = DataFrame({'b': [1.1, 2.2]})
        df = df.rename(columns={})
        df.insert(0, 'a', [1, 2])

        result = df.rename(columns={})
        str(result)
        expected = DataFrame([[1,1.1],[2, 2.2]],columns=['a','b'])
        assert_frame_equal(result,expected)
        df.insert(0, 'c', [1.3, 2.3])

        result = df.rename(columns={})
        str(result)

        expected = DataFrame([[1.3,1,1.1],[2.3,2, 2.2]],columns=['c','a','b'])
        assert_frame_equal(result,expected)

    def test_cast_internals(self):
        casted = DataFrame(self.frame._data, dtype=int)
        expected = DataFrame(self.frame._series, dtype=int)
        assert_frame_equal(casted, expected)

        casted = DataFrame(self.frame._data, dtype=np.int32)
        expected = DataFrame(self.frame._series, dtype=np.int32)
        assert_frame_equal(casted, expected)

    def test_consolidate(self):
        self.frame['E'] = 7.
        consolidated = self.frame.consolidate()
        self.assert_(len(consolidated._data.blocks) == 1)

        # Ensure copy, do I want this?
        recons = consolidated.consolidate()
        self.assert_(recons is not consolidated)
        assert_frame_equal(recons, consolidated)

        self.frame['F'] = 8.
        self.assert_(len(self.frame._data.blocks) == 3)
        self.frame.consolidate(inplace=True)
        self.assert_(len(self.frame._data.blocks) == 1)

    def test_consolidate_inplace(self):
        frame = self.frame.copy()

        # triggers in-place consolidation
        for letter in range(ord('A'), ord('Z')):
            self.frame[chr(letter)] = chr(letter)

    def test_as_matrix_consolidate(self):
        self.frame['E'] = 7.
        self.assert_(not self.frame._data.is_consolidated())
        _ = self.frame.as_matrix()
        self.assert_(self.frame._data.is_consolidated())

    def test_modify_values(self):
        self.frame.values[5] = 5
        self.assert_((self.frame.values[5] == 5).all())

        # unconsolidated
        self.frame['E'] = 7.
        self.frame.values[6] = 6
        self.assert_((self.frame.values[6] == 6).all())

    def test_boolean_set_uncons(self):
        self.frame['E'] = 7.

        expected = self.frame.values.copy()
        expected[expected > 1] = 2

        self.frame[self.frame > 1] = 2
        assert_almost_equal(expected, self.frame.values)

    def test_xs_view(self):
        dm = DataFrame(np.arange(20.).reshape(4, 5),
                       index=range(4), columns=range(5))

        dm.xs(2, copy=False)[:] = 5
        self.assert_((dm.xs(2) == 5).all())

        dm.xs(2)[:] = 10
        self.assert_((dm.xs(2) == 5).all())

        # TODO (?): deal with mixed-type fiasco?
        self.assertRaises(Exception, self.mixed_frame.xs,
                          self.mixed_frame.index[2], copy=False)

        # unconsolidated
        dm['foo'] = 6.
        dm.xs(3, copy=False)[:] = 10
        self.assert_((dm.xs(3) == 10).all())

    def test_boolean_indexing(self):
        idx = range(3)
        cols = ['A','B','C']
        df1 = DataFrame(index=idx, columns=cols,
                        data=np.array([[0.0, 0.5, 1.0],
                                       [1.5, 2.0, 2.5],
                                       [3.0, 3.5, 4.0]],
                                      dtype=float))
        df2 = DataFrame(index=idx, columns=cols,
                        data=np.ones((len(idx), len(cols))))

        expected = DataFrame(index=idx, columns=cols,
                             data=np.array([[0.0, 0.5, 1.0],
                                            [1.5, 2.0, -1],
                                            [-1, -1, -1]], dtype=float))

        df1[df1 > 2.0 * df2] = -1
        assert_frame_equal(df1, expected)

    def test_boolean_indexing_mixed(self):
        df = DataFrame(
            {0L: {35: np.nan, 40: np.nan, 43: np.nan, 49: np.nan, 50: np.nan},
             1L: {35: np.nan,
                  40: 0.32632316859446198,
                  43: np.nan,
                  49: 0.32632316859446198,
                  50: 0.39114724480578139},
             2L: {35: np.nan, 40: np.nan, 43: 0.29012581014105987, 49: np.nan, 50: np.nan},
             3L: {35: np.nan, 40: np.nan, 43: np.nan, 49: np.nan, 50: np.nan},
             4L: {35: 0.34215328467153283, 40: np.nan, 43: np.nan, 49: np.nan, 50: np.nan},
             'y': {35: 0, 40: 0, 43: 0, 49: 0, 50: 1}})

        # mixed int/float ok
        df2 = df.copy()
        df2[df2>0.3] = 1
        expected = df.copy()
        expected.loc[40,1] = 1
        expected.loc[49,1] = 1
        expected.loc[50,1] = 1
        expected.loc[35,4] = 1
        assert_frame_equal(df2,expected)

        # add object, should this raise?
        df['foo'] = 'test'
        self.assertRaises(ValueError, df.__setitem__, df>0.3, 1)

    def test_sum_bools(self):
        df = DataFrame(index=range(1), columns=range(10))
        bools = isnull(df)
        self.assert_(bools.sum(axis=1)[0] == 10)

    def test_fillna_col_reordering(self):
        idx = range(20)
        cols = ["COL." + str(i) for i in range(5, 0, -1)]
        data = np.random.rand(20, 5)
        df = DataFrame(index=range(20), columns=cols, data=data)
        filled = df.fillna(method='ffill')
        self.assert_(df.columns.tolist() == filled.columns.tolist())

    def test_take(self):

        # homogeneous
        #----------------------------------------
        order  = [3, 1, 2, 0]
        for df in [self.frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['D', 'B', 'C', 'A']]
            assert_frame_equal(result, expected, check_names=False)

        # neg indicies
        order = [2,1,-1]
        for df in [self.frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['C', 'B', 'D']]
            assert_frame_equal(result, expected, check_names=False)

        # illegal indices
        self.assertRaises(IndexError, df.take, [3,1,2,30], axis=0)
        self.assertRaises(IndexError, df.take, [3,1,2,-31], axis=0)
        self.assertRaises(IndexError, df.take, [3,1,2,5], axis=1)
        self.assertRaises(IndexError, df.take, [3,1,2,-5], axis=1)

        # mixed-dtype
        #----------------------------------------
        order  = [4, 1, 2, 0, 3]
        for df in [self.mixed_frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['foo', 'B', 'C', 'A', 'D']]
            assert_frame_equal(result, expected)

        # neg indicies
        order = [4,1,-2]
        for df in [self.mixed_frame]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['foo', 'B', 'D']]
            assert_frame_equal(result, expected)

        # by dtype
        order = [1, 2, 0, 3]
        for df in [self.mixed_float,self.mixed_int]:

            result = df.take(order, axis=0)
            expected = df.reindex(df.index.take(order))
            assert_frame_equal(result, expected)

            # axis = 1
            result = df.take(order, axis=1)
            expected = df.ix[:, ['B', 'C', 'A', 'D']]
            assert_frame_equal(result, expected)

    def test_iterkv_names(self):
        for k, v in self.mixed_frame.iterkv():
            self.assertEqual(v.name, k)

    def test_series_put_names(self):
        series = self.mixed_frame._series
        for k, v in series.iteritems():
            self.assertEqual(v.name, k)

    def test_dot(self):
        a = DataFrame(np.random.randn(3, 4), index=['a', 'b', 'c'],
                      columns=['p', 'q', 'r', 's'])
        b = DataFrame(np.random.randn(4, 2), index=['p', 'q', 'r', 's'],
                      columns=['one', 'two'])

        result = a.dot(b)
        expected = DataFrame(np.dot(a.values, b.values),
                             index=['a', 'b', 'c'],
                             columns=['one', 'two'])
        # Check alignment
        b1 = b.reindex(index=reversed(b.index))
        result = a.dot(b)
        assert_frame_equal(result, expected)

        # Check series argument
        result = a.dot(b['one'])
        assert_series_equal(result, expected['one'])
        result = a.dot(b1['one'])
        assert_series_equal(result, expected['one'])

        # can pass correct-length arrays
        row = a.ix[0].values

        result = a.dot(row)
        exp = a.dot(a.ix[0])
        assert_series_equal(result, exp)

        self.assertRaises(Exception, a.dot, row[:-1])

        a = np.random.rand(1, 5)
        b = np.random.rand(5, 1)
        A = DataFrame(a)
        B = DataFrame(b)

        # it works
        result = A.dot(b)

        # unaligned
        df = DataFrame(randn(3, 4), index=[1, 2, 3], columns=range(4))
        df2 = DataFrame(randn(5, 3), index=range(5), columns=[1, 2, 3])

        self.assertRaises(ValueError, df.dot, df2)

    def test_idxmin(self):
        frame = self.frame
        frame.ix[5:10] = np.nan
        frame.ix[15:20, -2:] = np.nan
        for skipna in [True, False]:
            for axis in [0, 1]:
                for df in [frame, self.intframe]:
                    result = df.idxmin(axis=axis, skipna=skipna)
                    expected = df.apply(
                        Series.idxmin, axis=axis, skipna=skipna)
                    assert_series_equal(result, expected)

        self.assertRaises(Exception, frame.idxmin, axis=2)

    def test_idxmax(self):
        frame = self.frame
        frame.ix[5:10] = np.nan
        frame.ix[15:20, -2:] = np.nan
        for skipna in [True, False]:
            for axis in [0, 1]:
                for df in [frame, self.intframe]:
                    result = df.idxmax(axis=axis, skipna=skipna)
                    expected = df.apply(
                        Series.idxmax, axis=axis, skipna=skipna)
                    assert_series_equal(result, expected)

        self.assertRaises(Exception, frame.idxmax, axis=2)

    def test_stale_cached_series_bug_473(self):
        Y = DataFrame(np.random.random((4, 4)), index=('a', 'b', 'c', 'd'),
                      columns=('e', 'f', 'g', 'h'))
        repr(Y)
        Y['e'] = Y['e'].astype('object')
        Y['g']['c'] = np.NaN
        repr(Y)
        result = Y.sum()
        exp = Y['g'].sum()
        self.assert_(isnull(Y['g']['c']))

    def test_index_namedtuple(self):
        try:
            from collections import namedtuple
        except ImportError:
            raise nose.SkipTest
        IndexType = namedtuple("IndexType", ["a", "b"])
        idx1 = IndexType("foo", "bar")
        idx2 = IndexType("baz", "bof")
        index = Index([idx1, idx2], name="composite_index")
        df = DataFrame([(1, 2), (3, 4)], index=index, columns=["A", "B"])
        self.assertEqual(df.ix[IndexType("foo", "bar")]["A"], 1)

    def test_bool_empty_nonzero(self):
        df = DataFrame([1, 2, 3])
        self.assertTrue(bool(df))
        self.assertFalse(df.empty)
        df = DataFrame(index=['a', 'b'], columns=['c', 'd']).dropna()
        self.assertFalse(bool(df))
        self.assertFalse(bool(df.T))
        self.assertTrue(df.empty)
        self.assertTrue(df.T.empty)

    def test_any_all(self):
        self._check_bool_op('any', np.any, has_skipna=True, has_bool_only=True)
        self._check_bool_op('all', np.all, has_skipna=True, has_bool_only=True)

        df = DataFrame(randn(10, 4)) > 0
        df.any(1)
        df.all(1)
        df.any(1, bool_only=True)
        df.all(1, bool_only=True)

        # skip pathological failure cases
        # class CantNonzero(object):

        #     def __nonzero__(self):
        #         raise ValueError

        # df[4] = CantNonzero()

        # it works!
        # df.any(1)
        # df.all(1)
        # df.any(1, bool_only=True)
        # df.all(1, bool_only=True)

        # df[4][4] = np.nan
        # df.any(1)
        # df.all(1)
        # df.any(1, bool_only=True)
        # df.all(1, bool_only=True)

    def test_consolidate_datetime64(self):
        # numpy vstack bug

        data = """\
starting,ending,measure
2012-06-21 00:00,2012-06-23 07:00,77
2012-06-23 07:00,2012-06-23 16:30,65
2012-06-23 16:30,2012-06-25 08:00,77
2012-06-25 08:00,2012-06-26 12:00,0
2012-06-26 12:00,2012-06-27 08:00,77
"""
        df = read_csv(StringIO(data), parse_dates=[0, 1])

        ser_starting = df.starting
        ser_starting.index = ser_starting.values
        ser_starting = ser_starting.tz_localize('US/Eastern')
        ser_starting = ser_starting.tz_convert('UTC')

        ser_ending = df.ending
        ser_ending.index = ser_ending.values
        ser_ending = ser_ending.tz_localize('US/Eastern')
        ser_ending = ser_ending.tz_convert('UTC')

        df.starting = ser_starting.index
        df.ending = ser_ending.index

        assert_array_equal(df.starting.values, ser_starting.index.values)
        assert_array_equal(df.ending.values, ser_ending.index.values)

    def _check_bool_op(self, name, alternative, frame=None, has_skipna=True,
                       has_bool_only=False):
        if frame is None:
            frame = self.frame > 0
            # set some NAs
            frame = DataFrame(frame.values.astype(object), frame.index,
                              frame.columns)
            frame.ix[5:10] = np.nan
            frame.ix[15:20, -2:] = np.nan

        f = getattr(frame, name)

        if has_skipna:
            def skipna_wrapper(x):
                nona = x.dropna().values
                return alternative(nona)

            def wrapper(x):
                return alternative(x.values)

            result0 = f(axis=0, skipna=False)
            result1 = f(axis=1, skipna=False)
            assert_series_equal(result0, frame.apply(wrapper))
            assert_series_equal(result1, frame.apply(wrapper, axis=1),
                                check_dtype=False)  # HACK: win32
        else:
            skipna_wrapper = alternative
            wrapper = alternative

        result0 = f(axis=0)
        result1 = f(axis=1)
        assert_series_equal(result0, frame.apply(skipna_wrapper))
        assert_series_equal(result1, frame.apply(skipna_wrapper, axis=1),
                            check_dtype=False)

        # result = f(axis=1)
        # comp = frame.apply(alternative, axis=1).reindex(result.index)
        # assert_series_equal(result, comp)

        self.assertRaises(Exception, f, axis=2)

        # make sure works on mixed-type frame
        mixed = self.mixed_frame
        mixed['_bool_'] = np.random.randn(len(mixed)) > 0
        getattr(mixed, name)(axis=0)
        getattr(mixed, name)(axis=1)

        class NonzeroFail:

            def __nonzero__(self):
                raise ValueError

        mixed['_nonzero_fail_'] = NonzeroFail()

        if has_bool_only:
            getattr(mixed, name)(axis=0, bool_only=True)
            getattr(mixed, name)(axis=1, bool_only=True)
            getattr(frame, name)(axis=0, bool_only=False)
            getattr(frame, name)(axis=1, bool_only=False)

        # all NA case
        if has_skipna:
            all_na = frame * np.NaN
            r0 = getattr(all_na, name)(axis=0)
            r1 = getattr(all_na, name)(axis=1)
            if name == 'any':
                self.assert_(not r0.any())
                self.assert_(not r1.any())
            else:
                self.assert_(r0.all())
                self.assert_(r1.all())

    def test_strange_column_corruption_issue(self):
        df = DataFrame(index=[0, 1])
        df[0] = nan
        wasCol = {}
        # uncommenting these makes the results match
        # for col in xrange(100, 200):
        #    wasCol[col] = 1
        #    df[col] = nan

        for i, dt in enumerate(df.index):
            for col in xrange(100, 200):
                if not col in wasCol:
                    wasCol[col] = 1
                    df[col] = nan
                df[col][dt] = i

        myid = 100

        first = len(df.ix[isnull(df[myid]), [myid]])
        second = len(df.ix[isnull(df[myid]), [myid]])
        self.assertTrue(first == second == 0)

    def test_inplace_return_self(self):
        # re #1893

        data = DataFrame({'a': ['foo', 'bar', 'baz', 'qux'],
                          'b': [0, 0, 1, 1],
                          'c': [1, 2, 3, 4]})

        def _check_f(base, f):
            result = f(base)
            self.assertTrue(result is None)

        # -----DataFrame-----

        # set_index
        f = lambda x: x.set_index('a', inplace=True)
        _check_f(data.copy(), f)

        # reset_index
        f = lambda x: x.reset_index(inplace=True)
        _check_f(data.set_index('a'), f)

        # drop_duplicates
        f = lambda x: x.drop_duplicates(inplace=True)
        _check_f(data.copy(), f)

        # sort
        f = lambda x: x.sort('b', inplace=True)
        _check_f(data.copy(), f)

        # sort_index
        f = lambda x: x.sort_index(inplace=True)
        _check_f(data.copy(), f)

        # sortlevel
        f = lambda x: x.sortlevel(0, inplace=True)
        _check_f(data.set_index(['a', 'b']), f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(data.copy(), f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(data.copy(), f)

        # rename
        f = lambda x: x.rename({1: 'foo'}, inplace=True)
        _check_f(data.copy(), f)

        # -----Series-----

        # reset_index
        f = lambda x: x.reset_index(inplace=True, drop=True)
        _check_f(data.set_index('a')['c'], f)

        # fillna
        f = lambda x: x.fillna(0, inplace=True)
        _check_f(data.copy()['c'], f)

        # replace
        f = lambda x: x.replace(1, 0, inplace=True)
        _check_f(data.copy()['c'], f)

        # rename
        f = lambda x: x.rename({1: 'foo'}, inplace=True)
        _check_f(data.copy()['c'], f)


if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--ipdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
