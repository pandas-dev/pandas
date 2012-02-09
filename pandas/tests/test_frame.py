# pylint: disable-msg=W0612,E1101
from copy import deepcopy
from datetime import datetime, timedelta
from StringIO import StringIO
import cPickle as pickle
import operator
import os
import sys
import unittest

import nose

from numpy import random, nan
from numpy.random import randn
import numpy as np
import numpy.ma as ma

import pandas as pan
import pandas.core.common as com
import pandas.core.format as fmt
import pandas.core.datetools as datetools
from pandas.core.index import NULL_INDEX
from pandas.core.api import (DataFrame, Index, Series, notnull, isnull,
                             MultiIndex)
from pandas.io.parsers import (ExcelFile, ExcelWriter)

from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal)

import pandas.util.testing as tm
import pandas._tseries as lib

#-------------------------------------------------------------------------------
# DataFrame test cases

JOIN_TYPES = ['inner', 'outer', 'left', 'right']

class CheckIndexing(object):

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

    def test_getitem_boolean_list(self):
        df = DataFrame(np.arange(12).reshape(3,4))
        def _checkit(lst):
            result = df[lst]
            expected = df.ix[df.index[lst]]
            assert_frame_equal(result, expected)

        _checkit([True, False, True])
        _checkit([True, True, True])
        _checkit([False, False, False])

    def test_getattr(self):
        tm.assert_series_equal(self.frame.A, self.frame['A'])
        self.assertRaises(AttributeError, getattr, self.frame,
                          'NONEXISTENT_NAME')

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

        df[df > 0] = 5
        values[values > 0] = 5
        assert_almost_equal(df.values, values)

        df[df == 5] = 0
        values[values == 5] = 0
        assert_almost_equal(df.values, values)

        self.assertRaises(Exception, df.__setitem__, df[:-1] > 0, 2)
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
        self.frame['B'] = 0
        self.assert_(self.frame['B'].dtype == np.float64)

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
        df = DataFrame({'B' : [1., 2., 3.],
                         'C' : ['a', 'b', 'c']},
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

        dm['C'] = 1
        self.assertEqual(dm['C'].dtype, np.int64)

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
        data = {"title" : ['foobar','bar','foobar'] + ['foobar'] * 17 ,
                "cruft" : np.random.random(20)}

        df = DataFrame(data)
        ix = df[df['title'] == 'bar'].index

        df.ix[ix, ['title']] = 'foobar'
        df.ix[ix, ['cruft']] = 0

        assert( df.ix[1, 'title'] == 'foobar' )
        assert( df.ix[1, 'cruft'] == 0 )

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
                       index=[0,1,2,3])
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

        self.assertRaises(Exception, self.mixed_frame.ix.__setitem__,
                          5, self.mixed_frame.ix[6])

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
                expected.values[i,j] = val
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
        result = self.frame.ix[[1,4,7]]
        expected = self.frame.ix[self.frame.index[[1,4,7]]]
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

    def test_setitem_single_column_mixed(self):
        df = DataFrame(randn(5, 3), index=['a', 'b', 'c', 'd', 'e'],
                       columns=['foo', 'bar', 'baz'])
        df['str'] = 'qux'
        df.ix[::2, 'str'] = nan
        expected = [nan, 'qux', nan, 'qux', nan]
        assert_almost_equal(df['str'].values, expected)

    def test_setitem_fancy_exceptions(self):
        pass

    def test_getitem_boolean_missing(self):
        pass

    def test_setitem_boolean_missing(self):
        pass

    def test_get_value(self):
        for idx in self.frame.index:
            for col in self.frame.columns:
                result = self.frame.get_value(idx, col)
                expected = self.frame[col][idx]
                assert_almost_equal(result, expected)

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

        df = DataFrame({'label' : ['a', 'b', 'a', 'c'],
                        'mask_a' : [True, True, False, True],
                        'mask_b' : [True, False, False, False],
                        'mask_c' : [False, True, False, True]})
        df['mask'] = df.lookup(df.index, 'mask_' + df['label'])
        exp_mask = alt(df, df.index, 'mask_' + df['label'])
        assert_almost_equal(df['mask'], exp_mask)
        self.assert_(df['mask'].dtype == np.bool_)

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

    def test_iget_value(self):
        for i, row in enumerate(self.frame.index):
            for j, col in enumerate(self.frame.columns):
                result = self.frame.iget_value(i, j)
                expected = self.frame.get_value(row, col)
                assert_almost_equal(result, expected)

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
        assert_frame_equal(joined, self.frame)

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

    def setUp(self):
        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.intframe = _intframe.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()

        self.ts1 = tm.makeTimeSeries()
        self.ts2 = tm.makeTimeSeries()[5:]
        self.ts3 = tm.makeTimeSeries()[-5:]
        self.ts4 = tm.makeTimeSeries()[1:-1]

        self.ts_dict = {
            'col1' : self.ts1,
            'col2' : self.ts2,
            'col3' : self.ts3,
            'col4' : self.ts4,
        }
        self.empty = DataFrame({})

        arr = np.array([[1., 2., 3.],
                        [4., 5., 6.],
                        [7., 8., 9.]])

        self.simple = DataFrame(arr, columns=['one', 'two', 'three'],
                                 index=['a', 'b', 'c'])

    def test_get_axis(self):
        self.assert_(DataFrame._get_axis_name(0) == 'index')
        self.assert_(DataFrame._get_axis_name(1) == 'columns')
        self.assert_(DataFrame._get_axis_name('index') == 'index')
        self.assert_(DataFrame._get_axis_name('columns') == 'columns')
        self.assertRaises(Exception, DataFrame._get_axis_name, 'foo')
        self.assertRaises(Exception, DataFrame._get_axis_name, None)

        self.assert_(DataFrame._get_axis_number(0) == 0)
        self.assert_(DataFrame._get_axis_number(1) == 1)
        self.assert_(DataFrame._get_axis_number('index') == 0)
        self.assert_(DataFrame._get_axis_number('columns') == 1)
        self.assertRaises(Exception, DataFrame._get_axis_number, 2)
        self.assertRaises(Exception, DataFrame._get_axis_number, None)

        self.assert_(self.frame._get_axis(0) is self.frame.index)
        self.assert_(self.frame._get_axis(1) is self.frame.columns)

    def test_set_index(self):
        idx = Index(np.arange(len(self.mixed_frame)))

        # cache it
        _ = self.mixed_frame['foo']
        self.mixed_frame.index = idx
        self.assert_(self.mixed_frame['foo'].index  is idx)
        self.assertRaises(Exception, setattr, self.mixed_frame, 'index',
                          idx[::2])

    def test_set_columns(self):
        cols = Index(np.arange(len(self.mixed_frame.columns)))
        self.mixed_frame.columns = cols
        self.assertRaises(Exception, setattr, self.mixed_frame, 'columns',
                          cols[::2])

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
        df = DataFrame({0 : np.ones(10, dtype=bool),
                        1 : np.zeros(10, dtype=bool)})
        self.assertEqual(df.values.dtype, np.bool_)

    def test_is_mixed_type(self):
        self.assert_(not self.frame._is_mixed_type)
        self.assert_(self.mixed_frame._is_mixed_type)

    def test_constructor_dict(self):
        frame = DataFrame({'col1' : self.ts1,
                            'col2' : self.ts2})

        tm.assert_dict_equal(self.ts1, frame['col1'], compare_keys=False)
        tm.assert_dict_equal(self.ts2, frame['col2'], compare_keys=False)

        frame = DataFrame({'col1' : self.ts1,
                            'col2' : self.ts2},
                           columns=['col2', 'col3', 'col4'])

        self.assertEqual(len(frame), len(self.ts2))
        self.assert_('col1' not in frame)
        self.assert_(np.isnan(frame['col3']).all())

        # Corner cases
        self.assertEqual(len(DataFrame({})), 0)
        self.assertRaises(Exception, lambda x: DataFrame([self.ts1, self.ts2]))

        # mix dict and array, wrong size
        self.assertRaises(Exception, DataFrame,
                          {'A' : {'a' : 'a', 'b' : 'b'},
                           'B' : ['a', 'b', 'c']})


        # Length-one dict micro-optimization
        frame = DataFrame({'A' : {'1' : 1, '2' : 2}})
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
        frame = DataFrame({'A' : [], 'B' : []}, columns=['A', 'B'])
        self.assert_(frame.index is NULL_INDEX)

    def test_constructor_subclass_dict(self):
        # Test for passing dict subclass to constructor
        data = {'col1': tm.TestSubDict((x, 10.0 * x) for x in xrange(10)),
                'col2': tm.TestSubDict((x, 20.0 * x) for x in xrange(10))}
        df = DataFrame(data)
        refdf = DataFrame(dict((col, dict(val.iteritems())) for col, val in data.iteritems()))
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
        df = DataFrame({'d' : [4.],'c' : [3.],'b' : [2.],'a' : [1.]},
                       columns=['d', 'c', 'b', 'a'])
        assert_almost_equal(df.values, expected)

    def test_constructor_dict_cast(self):
        # cast float tests
        test_data = {
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
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
                'A' : dict(zip(range(20), tm.makeDateIndex(20))),
                'B' : dict(zip(range(15), randn(15)))
        }
        frame = DataFrame(test_data, dtype=float)
        self.assertEqual(len(frame), 20)
        self.assert_(frame['A'].dtype == np.object_)
        self.assert_(frame['B'].dtype == np.float64)

    def test_constructor_dict_dont_upcast(self):
        d = {'Col1': {'Row1': 'A String', 'Row2': np.nan}}
        df = DataFrame(d)
        self.assert_(isinstance(df['Col1']['Row2'], float))

        dm = DataFrame([[1,2],['a','b']], index=[1,2], columns=[1,2])
        self.assert_(isinstance(dm[1][1], int))

    def test_constructor_ndarray(self):
        mat = np.zeros((2, 3), dtype=float)

        # 2-D input
        frame = DataFrame(mat, columns=['A', 'B', 'C'], index=[1, 2])

        self.assertEqual(len(frame.index), 2)
        self.assertEqual(len(frame.columns), 3)

        # cast type
        frame = DataFrame(mat, columns=['A', 'B', 'C'],
                           index=[1, 2], dtype=int)
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
        self.assert_(frame.index is NULL_INDEX)

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
                           index=[1, 2], dtype=int)
        self.assert_(frame.values.dtype == np.int64)

        # Check non-masked values
        mat2 = ma.copy(mat)
        mat2[0,0] = 1.0
        mat2[1,2] = 2.0
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
        self.assert_(frame.index is NULL_INDEX)

        frame = DataFrame(ma.masked_all((3, 0)))
        self.assert_(len(frame.columns) == 0)

    def test_constructor_corner(self):
        df = DataFrame(index=[])
        self.assertEqual(df.values.shape, (0, 0))

        # empty but with specified dtype
        df = DataFrame(index=range(10), columns=['a','b'], dtype=object)
        self.assert_(df.values.dtype == np.object_)

        # does not error but ends up float
        df = DataFrame(index=range(10), columns=['a','b'], dtype=int)
        self.assert_(df.values.dtype == np.float64)

    def test_constructor_scalar_inference(self):
        data = {'int' : 1, 'bool' : True,
                'float' : 3., 'object' : 'foo'}
        df = DataFrame(data, index=np.arange(10))

        self.assert_(df['int'].dtype == np.int64)
        self.assert_(df['bool'].dtype == np.bool_)
        self.assert_(df['float'].dtype == np.float64)
        self.assert_(df['object'].dtype == np.object_)

    def test_constructor_DataFrame(self):
        df = DataFrame(self.frame)
        assert_frame_equal(df, self.frame)

        df_casted = DataFrame(self.frame, dtype=int)
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
        dm = DataFrame({'A' : np.ones(10, dtype=int),
                         'B' : np.ones(10, dtype=float)},
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
        data = [{'a': 1.5, 'b': 3, 'c':4, 'd':6},
                {'a': 1.5, 'b': 3, 'd':6},
                {'a': 1.5, 'd':6},
                {},
                {'a': 1.5, 'b': 3, 'c':4},
                {'b': 3, 'c':4, 'd':6}]

        result = DataFrame(data)
        expected = DataFrame.from_dict(dict(zip(range(len(data)), data)),
                                       orient='index')
        assert_frame_equal(result, expected.reindex(result.index))

        result = DataFrame([{}])
        expected = DataFrame([])
        assert_frame_equal(result, expected)

    def test_constructor_ragged(self):
        data = {'A' : randn(10),
                'B' : randn(8)}
        self.assertRaises(Exception, DataFrame, data)

    def test_constructor_scalar(self):
        idx = Index(range(3))
        df = DataFrame({"a" : 0}, index=idx)
        expected = DataFrame({"a" : [0, 0, 0]}, index=idx)
        assert_frame_equal(df, expected)

    def test_constructor_Series_copy_bug(self):
        df = DataFrame(self.frame['A'], index=self.frame.index, columns=['A'])
        df.copy()

    def test_constructor_mixed_dict_and_Series(self):
        data = {}
        data['A'] = {'foo' : 1, 'bar' : 2, 'baz' : 3}
        data['B'] = Series([4, 3, 2, 1], index=['bar', 'qux', 'baz', 'foo'])

        result = DataFrame(data)
        self.assert_(result.index.is_monotonic)

        # ordering ambiguous, raise exception
        self.assertRaises(Exception, DataFrame,
                          {'A' : ['a', 'b'], 'B' : {'a' : 'a', 'b' : 'b'}})

        # this is OK though
        result = DataFrame({'A' : ['a', 'b'],
                            'B' : Series(['a', 'b'], index=['a', 'b'])})
        expected = DataFrame({'A' : ['a', 'b'], 'B' : ['a', 'b']},
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

    def test_constructor_Series_named(self):
        a = Series([1,2,3], index=['a','b','c'], name='x')
        df = DataFrame(a)
        self.assert_(df.columns[0] == 'x')
        self.assert_(df.index.equals(a.index))

    def test_constructor_Series_differently_indexed(self):
        # name
        s1 = Series([1, 2, 3], index=['a','b','c'], name='x')

        # no name
        s2 = Series([1, 2, 3], index=['a','b','c'])

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
        arr = lib.list_to_object_array([('bar', 'baz')] * len(self.mixed_frame))
        self.mixed_frame['foo'] = arr
        row_items = [(idx, list(self.mixed_frame.xs(idx)))
                     for idx in self.mixed_frame.index]
        recons = DataFrame.from_items(row_items,
                                      columns=self.mixed_frame.columns,
                                      orient='index')
        assert_frame_equal(recons, self.mixed_frame)
        self.assert_(isinstance(recons['foo'][0], tuple))

    def test_constructor_mix_series_nonseries(self):
        df = DataFrame({'A' : self.frame['A'],
                        'B' : list(self.frame['B'])}, columns=['A', 'B'])
        assert_frame_equal(df, self.frame.ix[:, ['A', 'B']])

        self.assertRaises(Exception, DataFrame,
                          {'A' : self.frame['A'],
                           'B' : list(self.frame['B'])[:-2]})

    def test_astype(self):
        casted = self.frame.astype(int)
        expected = DataFrame(self.frame.values.astype(int),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

        self.frame['foo'] = '5'
        casted = self.frame.astype(int)
        expected = DataFrame(self.frame.values.astype(int),
                             index=self.frame.index,
                             columns=self.frame.columns)
        assert_frame_equal(casted, expected)

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
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
        }
        recons_data = DataFrame(test_data).to_dict()

        for k, v in test_data.iteritems():
            for k2, v2 in v.iteritems():
                self.assertEqual(v2, recons_data[k][k2])

    def test_from_records_to_records(self):
        # from numpy documentation
        arr = np.zeros((2,),dtype=('i4,f4,a10'))
        arr[:] = [(1,2.,'Hello'),(2,3.,"World")]

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

    def test_from_records_sequencelike(self):
        df = DataFrame({'A' : np.random.randn(6),
                        'B' : np.arange(6),
                        'C' : ['foo'] * 6,
                        'D' : np.array([True, False] * 3, dtype=bool)})

        tuples = [tuple(x) for x in df.values]
        lists = [list(x) for x in tuples]
        asdict = dict((x,y) for x, y in df.iteritems())

        result = DataFrame.from_records(tuples, columns=df.columns)
        result2 = DataFrame.from_records(lists, columns=df.columns)
        result3 = DataFrame.from_records(asdict, columns=df.columns)

        assert_frame_equal(result, df)
        assert_frame_equal(result2, df)
        assert_frame_equal(result3, df)

        result = DataFrame.from_records(tuples)
        self.assert_(np.array_equal(result.columns, range(4)))

        # test exclude parameter
        result = DataFrame.from_records(tuples, exclude=[0,1,3])
        result.columns = ['C']
        assert_frame_equal(result, df[['C']])

        # empty case
        result = DataFrame.from_records([], columns=['foo', 'bar', 'baz'])
        self.assertEqual(len(result), 0)
        self.assert_(np.array_equal(result.columns, ['foo', 'bar', 'baz']))

        result = DataFrame.from_records([])
        self.assertEqual(len(result), 0)
        self.assertEqual(len(result.columns), 0)

    def test_from_records_with_index_data(self):
        df = DataFrame(np.random.randn(10,3), columns=['A', 'B', 'C'])

        data = np.random.randn(10)
        df1 = DataFrame.from_records(df, index=data)
        assert(df1.index.equals(Index(data)))

    def test_from_records_bad_index_column(self):
        df = DataFrame(np.random.randn(10,3), columns=['A', 'B', 'C'])

        # should pass
        df1 = DataFrame.from_records(df, index=['C'])
        assert(df1.index.equals(Index(df.C)))

        df1 = DataFrame.from_records(df, index='C')
        assert(df1.index.equals(Index(df.C)))

        # should fail
        self.assertRaises(Exception, DataFrame.from_records, df, index=[2])
        self.assertRaises(KeyError, DataFrame.from_records, df, index=2)

    def test_get_agg_axis(self):
        cols = self.frame._get_agg_axis(0)
        self.assert_(cols is self.frame.columns)

        idx = self.frame._get_agg_axis(1)
        self.assert_(idx is self.frame.index)

        self.assertRaises(Exception, self.frame._get_agg_axis, 2)

    def test_nonzero(self):
        self.assertFalse(self.empty)

        self.assert_(self.frame)
        self.assert_(self.mixed_frame)

        # corner case
        df = DataFrame({'A' : [1., 2., 3.],
                         'B' : ['a', 'b', 'c']},
                        index=np.arange(3))
        del df['A']
        self.assert_(df)

    def test_repr(self):
        buf = StringIO()

        # empty
        foo = repr(self.empty)

        # empty with index
        frame = DataFrame(index=np.arange(1000))
        foo = repr(frame)

        # small one
        foo = repr(self.frame)
        self.frame.info(verbose=False, buf=buf)

        # even smaller
        self.frame.reindex(columns=['A']).info(verbose=False, buf=buf)
        self.frame.reindex(columns=['A', 'B']).info(verbose=False, buf=buf)

        # big one
        biggie = DataFrame(np.zeros((200, 4)), columns=range(4),
                            index=range(200))
        foo = repr(biggie)

        # mixed
        foo = repr(self.mixed_frame)
        self.mixed_frame.info(verbose=False, buf=buf)

        # big mixed
        biggie = DataFrame({'A' : randn(200),
                             'B' : tm.makeStringIndex(200)},
                            index=range(200))
        biggie['A'][:20] = nan
        biggie['B'][:20] = nan

        foo = repr(biggie)

        # exhausting cases in DataFrame.info

        # columns but no index
        no_index = DataFrame(columns=[0, 1, 3])
        foo = repr(no_index)

        # no columns or index
        self.empty.info(buf=buf)

        # columns are not sortable

        unsortable = DataFrame({'foo' : [1] * 50,
                                datetime.today() : [1] * 50,
                                'bar' : ['bar'] * 50,
                                datetime.today() + timedelta(1) : ['bar'] * 50},
                               index=np.arange(50))
        foo = repr(unsortable)

        fmt.set_printoptions(precision=3, column_space=10)
        repr(self.frame)

        fmt.set_printoptions(max_rows=10, max_columns=2)
        repr(self.frame)
        fmt.reset_printoptions()

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
        A = self.frame.pop('A')
        self.assert_('A' not in self.frame)

        self.frame['foo'] = 'bar'
        foo = self.frame.pop('foo')
        self.assert_('foo' not in self.frame)

    def test_iter(self):
        self.assert_(tm.equalContents(list(self.frame), self.frame.columns))

    def test_iterrows(self):
        for i, (k, v) in enumerate(self.frame.iterrows()):
            exp = self.frame.xs(self.frame.index[i])
            assert_series_equal(v, exp)

        for i, (k, v) in enumerate(self.mixed_frame.iterrows()):
            exp = self.mixed_frame.xs(self.mixed_frame.index[i])
            assert_series_equal(v, exp)

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
               'b': {'a': False, 'b': True, 'c': False, 'd': False, 'e': False},
               'c': {'a': False, 'b': False, 'c': True, 'd': False, 'e': False},
               'd': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True},
               'e': {'a': True, 'b': False, 'c': False, 'd': True, 'e': True}}

        df2 = {'a': {'a': True, 'b': False, 'c': True, 'd': False, 'e': False},
               'b': {'a': False, 'b': True, 'c': False, 'd': False, 'e': False},
               'c': {'a': True, 'b': False, 'c': True, 'd': False, 'e': False},
               'd': {'a': False, 'b': False, 'c': False, 'd': True, 'e': False},
               'e': {'a': False, 'b': False, 'c': False, 'd': False, 'e': True}}

        df1 = DataFrame(df1)
        df2 = DataFrame(df2)

        _check_bin_op(operator.and_)
        _check_bin_op(operator.or_)
        _check_bin_op(operator.xor)

        _check_unary_op(operator.neg)

    def test_neg(self):
        # what to do?
        assert_frame_equal(-self.frame, -1 * self.frame)

    def test_first_last_valid(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan
        mat[-5:] = nan

        frame = DataFrame({'foo' : mat}, index=self.frame.index)
        index = frame.first_valid_index()

        self.assert_(index == frame.index[5])

        index = frame.last_valid_index()
        self.assert_(index == frame.index[-6])

    def test_arith_flex_frame(self):
        res_add = self.frame.add(self.frame)
        res_sub = self.frame.sub(self.frame)
        res_mul = self.frame.mul(self.frame)
        res_div = self.frame.div(2 * self.frame)

        assert_frame_equal(res_add, self.frame + self.frame)
        assert_frame_equal(res_sub, self.frame - self.frame)
        assert_frame_equal(res_mul, self.frame * self.frame)
        assert_frame_equal(res_div, self.frame / (2 * self.frame))

        const_add = self.frame.add(1)
        assert_frame_equal(const_add, self.frame + 1)

        # corner cases
        result = self.frame.add(self.frame[:0])
        assert_frame_equal(result, self.frame * np.nan)

        result = self.frame[:0].add(self.frame)
        assert_frame_equal(result, self.frame * np.nan)

    def test_arith_flex_series(self):
        df = self.simple

        row = df.xs('a')
        col = df['two']

        assert_frame_equal(df.add(row), df + row)
        assert_frame_equal(df.add(row, axis=None), df + row)
        assert_frame_equal(df.sub(row), df - row)
        assert_frame_equal(df.div(row), df / row)
        assert_frame_equal(df.mul(row), df * row)

        assert_frame_equal(df.add(col, axis=0), (df.T + col).T)
        assert_frame_equal(df.sub(col, axis=0), (df.T - col).T)
        assert_frame_equal(df.div(col, axis=0), (df.T / col).T)
        assert_frame_equal(df.mul(col, axis=0), (df.T * col).T)

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
        self.assert_(not empty_empty)

        # out of order
        reverse = self.frame.reindex(columns=self.frame.columns[::-1])

        assert_frame_equal(reverse + self.frame, self.frame * 2)

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

        # TimeSeries
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

    def test_combineFunc(self):
        result = self.frame * 2
        self.assert_(np.array_equal(result.values, self.frame.values * 2))

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

    def test_to_csv_from_csv(self):
        path = '__tmp__'

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
        dm = DataFrame({'s1' : Series(range(3),range(3)),
                        's2' : Series(range(2),range(2))})
        dm.to_csv(path)
        recons = DataFrame.from_csv(path)
        assert_frame_equal(dm, recons)

        os.remove(path)

    def test_to_csv_multiindex(self):
        path = '__tmp__'

        frame = self.frame
        old_index = frame.index
        arrays = np.arange(len(old_index)*2).reshape(2,-1)
        new_index = MultiIndex.from_arrays(arrays, names=['first', 'second'])
        frame.index = new_index
        frame.to_csv(path, header=False)
        frame.to_csv(path, cols=['A', 'B'])

        # round trip
        frame.to_csv(path)
        df = DataFrame.from_csv(path, index_col=[0,1], parse_dates=False)

        assert_frame_equal(frame, df)
        self.assertEqual(frame.index.names, df.index.names)
        self.frame.index = old_index # needed if setUP becomes a classmethod

        # try multiindex with dates
        tsframe = self.tsframe
        old_index = tsframe.index
        new_index = [old_index, np.arange(len(old_index))]
        tsframe.index = MultiIndex.from_arrays(new_index)

        tsframe.to_csv(path, index_label = ['time','foo'])
        recons = DataFrame.from_csv(path, index_col=[0,1])
        assert_frame_equal(tsframe, recons)

        # do not load index
        tsframe.to_csv(path)
        recons = DataFrame.from_csv(path, index_col=None)
        np.testing.assert_equal(len(recons.columns), len(tsframe.columns) + 2)

        # no index
        tsframe.to_csv(path, index=False)
        recons = DataFrame.from_csv(path, index_col=None)
        assert_almost_equal(recons.values, self.tsframe.values)
        self.tsframe.index = old_index # needed if setUP becomes classmethod

        os.remove(path)

        # empty
        tsframe[:0].to_csv(path)
        recons = DataFrame.from_csv(path)
        assert_frame_equal(recons, tsframe[:0])

    def test_to_csv_float32_nanrep(self):
        df = DataFrame(np.random.randn(1, 4).astype(np.float32))
        df[1] = np.nan

        pth = '__tmp__.csv'
        df.to_csv(pth, na_rep=999)

        lines = open(pth).readlines()
        self.assert_(lines[1].split(',')[2] == '999')
        os.remove(pth)

    def test_to_csv_withcommas(self):

        path = '__tmp__'
        # Commas inside fields should be correctly escaped when saving as CSV.

        df = DataFrame({'A':[1,2,3], 'B':['5,6','7,8','9,0']})
        df.to_csv(path)
        df2 = DataFrame.from_csv(path)
        assert_frame_equal(df2, df)

        os.remove(path)

    def test_to_csv_bug(self):
        path = '__tmp__.csv'
        f1 = StringIO('a,1.0\nb,2.0')
        df = DataFrame.from_csv(f1,header=None)
        newdf = DataFrame({'t': df[df.columns[0]]})
        newdf.to_csv(path)

        recons = pan.read_csv(path, index_col=0)
        assert_frame_equal(recons, newdf)

        os.remove(path)

    def test_to_csv_unicode(self):
        path = '__tmp__.csv'
        df = DataFrame({u'c/\u03c3':[1,2,3]})
        df.to_csv(path, encoding='UTF-8')
        df2 = pan.read_csv(path, index_col=0, encoding='UTF-8')
        assert_frame_equal(df, df2)

        df.to_csv(path, encoding='UTF-8', index=False)
        df2 = pan.read_csv(path, index_col=None, encoding='UTF-8')
        assert_frame_equal(df, df2)

        os.remove(path)

    def test_to_csv_stringio(self):
        buf = StringIO()
        self.frame.to_csv(buf)
        buf.seek(0)
        recons = pan.read_csv(buf, index_col=0)
        assert_frame_equal(recons, self.frame)

    def test_to_excel_from_excel(self):
        try:
            import xlwt
            import xlrd
            import openpyxl
        except ImportError:
            raise nose.SkipTest

        for ext in ['xls', 'xlsx']:
            path = '__tmp__.' + ext

            self.frame['A'][:5] = nan

            self.frame.to_excel(path,'test1')
            self.frame.to_excel(path,'test1', cols=['A', 'B'])
            self.frame.to_excel(path,'test1', header=False)
            self.frame.to_excel(path,'test1', index=False)

            # test roundtrip
            self.frame.to_excel(path,'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0)
            assert_frame_equal(self.frame, recons)

            self.frame.to_excel(path,'test1', index=False)
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=None)
            recons.index = self.frame.index
            assert_frame_equal(self.frame, recons)

            self.frame.to_excel(path,'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0, skiprows=[1])
            assert_frame_equal(self.frame.ix[1:], recons)

            self.frame.to_excel(path,'test1',na_rep='NA')
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0, na_values=['NA'])
            assert_frame_equal(self.frame, recons)

            self.mixed_frame.to_excel(path,'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0)
            assert_frame_equal(self.mixed_frame, recons)

            self.tsframe.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1')
            assert_frame_equal(self.tsframe, recons)

            #Test np.int64
            frame = DataFrame(np.random.randn(10,2))
            frame.to_excel(path,'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1')
            assert_frame_equal(frame, recons)

            # Test writing to separate sheets
            writer = ExcelWriter(path)
            self.frame.to_excel(writer,'test1')
            self.tsframe.to_excel(writer,'test2')
            writer.save()
            reader = ExcelFile(path)
            recons = reader.parse('test1',index_col=0)
            assert_frame_equal(self.frame, recons)
            recons = reader.parse('test2',index_col=0)
            assert_frame_equal(self.tsframe, recons)
            np.testing.assert_equal(2, len(reader.sheet_names))
            np.testing.assert_equal('test1', reader.sheet_names[0])
            np.testing.assert_equal('test2', reader.sheet_names[1])


            os.remove(path)

        # datetime.date, not sure what to test here exactly
        path = '__tmp__.xls'
        tsf = self.tsframe.copy()
        tsf.index = [x.date() for x in self.tsframe.index]
        tsf.to_excel(path, 'test1')
        reader = ExcelFile(path)
        recons = reader.parse('test1')
        assert_frame_equal(self.tsframe, recons)
        os.remove(path)

    def test_to_excel_multiindex(self):
        try:
            import xlwt
            import xlrd
            import openpyxl
        except ImportError:
            raise nose.SkipTest

        for ext in ['xls', 'xlsx']:
            path = '__tmp__.' + ext

            frame = self.frame
            old_index = frame.index
            arrays = np.arange(len(old_index)*2).reshape(2,-1)
            new_index = MultiIndex.from_arrays(arrays,
                                               names=['first', 'second'])
            frame.index = new_index
            frame.to_excel(path, 'test1', header=False)
            frame.to_excel(path, 'test1', cols=['A', 'B'])

            # round trip
            frame.to_excel(path, 'test1')
            reader = ExcelFile(path)
            df = reader.parse('test1', index_col=[0,1], parse_dates=False)
            assert_frame_equal(frame, df)
            self.assertEqual(frame.index.names, df.index.names)
            self.frame.index = old_index # needed if setUP becomes a classmethod

            # try multiindex with dates
            tsframe = self.tsframe
            old_index = tsframe.index
            new_index = [old_index, np.arange(len(old_index))]
            tsframe.index = MultiIndex.from_arrays(new_index)

            tsframe.to_excel(path, 'test1', index_label = ['time','foo'])
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=[0,1])
            assert_frame_equal(tsframe, recons)

            # infer index
            tsframe.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1')
            assert_frame_equal(tsframe, recons)

            # no index
            tsframe.index.names = ['first', 'second']
            tsframe.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1')
            assert_almost_equal(tsframe.values,
                                recons.ix[:, tsframe.columns].values)
            self.assertEqual(len(tsframe.columns) + 2, len(recons.columns))

            tsframe.index.names = [None, None]

            # no index
            tsframe.to_excel(path, 'test1', index=False)
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=None)
            assert_almost_equal(recons.values, self.tsframe.values)
            self.tsframe.index = old_index # needed if setUP becomes classmethod

            # write a big DataFrame
            df = DataFrame(np.random.randn(1005, 1))
            df.to_excel(path, 'test1')

            os.remove(path)

    def test_info(self):
        io = StringIO()
        self.frame.info(buf=io)
        self.tsframe.info(buf=io)

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

    def test_convert_objects_no_conversion(self):
        mixed1 = DataFrame({'a': [1,2,3], 'b': [4.0, 5, 6], 'c': ['x','y','z']})
        mixed2 = mixed1.convert_objects()
        assert_frame_equal(mixed1, mixed2)

    def test_append_series_dict(self):
        df = DataFrame(np.random.randn(5, 4),
                       columns=['foo', 'bar', 'baz', 'qux'])

        series = df.ix[4]
        self.assertRaises(Exception, df.append, series)
        series.name = None
        self.assertRaises(Exception, df.append, series)

        result = df.append(series[::-1], ignore_index=True)
        expected = df.append(DataFrame({0 : series[::-1]}, index=df.columns).T, ignore_index=True)
        assert_frame_equal(result, expected)

        # dict
        result = df.append(series.to_dict(), ignore_index=True)
        assert_frame_equal(result, expected)

        result = df.append(series[::-1][:3], ignore_index=True)
        expected = df.append(DataFrame({0 : series[::-1][:3]}).T,
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

    def test_asfreq(self):
        offset_monthly = self.tsframe.asfreq(datetools.bmonthEnd)
        rule_monthly = self.tsframe.asfreq('EOM')

        assert_almost_equal(offset_monthly['A'], rule_monthly['A'])

        filled = rule_monthly.asfreq('WEEKDAY', method='pad')
        # TODO: actually check that this worked.

        # don't forget!
        filled_dep = rule_monthly.asfreq('WEEKDAY', method='pad')

        # test does not blow up on length-0 DataFrame
        zero_length = self.tsframe.reindex([])
        result = zero_length.asfreq('EOM')
        self.assert_(result is not zero_length)

    def test_asfreq_DateRange(self):
        from pandas.core.daterange import DateRange
        df = DataFrame({'A': [1,2,3]},
                       index=[datetime(2011,11,01), datetime(2011,11,2),
                              datetime(2011,11,3)])
        df = df.asfreq('WEEKDAY')
        self.assert_(isinstance(df.index, DateRange))

        ts = df['A'].asfreq('WEEKDAY')
        self.assert_(isinstance(ts.index, DateRange))

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

        # single block corner case
        mat = self.frame.as_matrix(['A', 'B'])
        expected = self.frame.reindex(columns=['A', 'B']).values
        assert_almost_equal(mat, expected)

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

    def test_corr(self):
        self.frame['A'][:5] = nan
        self.frame['B'][:10] = nan

        def _check_method(method='pearson'):
            correls = self.frame.corr(method=method)
            exp = self.frame['A'].corr(self.frame['C'], method=method)
            assert_almost_equal(correls['A']['C'], exp)

        _check_method('pearson')
        _check_method('kendall')
        _check_method('spearman')

        # exclude non-numeric types
        result = self.mixed_frame.corr()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].corr()
        assert_frame_equal(result, expected)

    def test_cov(self):
        self.frame['A'][:5] = nan
        self.frame['B'][:10] = nan
        cov = self.frame.cov()

        assert_almost_equal(cov['A']['C'],
                            self.frame['A'].cov(self.frame['C']))

        # exclude non-numeric types
        result = self.mixed_frame.cov()
        expected = self.mixed_frame.ix[:, ['A', 'B', 'C', 'D']].cov()
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

    def test_dropEmptyRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo' : mat}, index=self.frame.index)

        smaller_frame = frame.dropna(how='all')
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

        smaller_frame = frame.dropna(how='all', subset=['foo'])
        self.assert_(np.array_equal(smaller_frame['foo'], mat[5:]))

    def test_dropIncompleteRows(self):
        N = len(self.frame.index)
        mat = randn(N)
        mat[:5] = nan

        frame = DataFrame({'foo' : mat}, index=self.frame.index)
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

    def test_drop_duplicates(self):
        df = DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                               'foo', 'bar', 'bar', 'foo'],
                        'B' : ['one', 'one', 'two', 'two',
                               'two', 'two', 'one', 'two'],
                        'C' : [1, 1, 2, 2, 2, 2, 1, 2],
                        'D' : range(8)})

        # single column
        result = df.drop_duplicates('A')
        expected = df[:2]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates('A', take_last=True)
        expected = df.ix[[6, 7]]
        assert_frame_equal(result, expected)

        # multi column
        result = df.drop_duplicates(['A', 'B'])
        expected = df.ix[[0, 1, 2, 3]]
        assert_frame_equal(result, expected)

        result = df.drop_duplicates(['A', 'B'], take_last=True)
        expected = df.ix[[0, 5, 6, 7]]
        assert_frame_equal(result, expected)

        # consider everything
        df2 = df.ix[:, ['A', 'B', 'C']]

        result = df2.drop_duplicates()
        # in this case only
        expected = df2.drop_duplicates(['A', 'B'])
        assert_frame_equal(result, expected)

        result = df2.drop_duplicates(take_last=True)
        expected = df2.drop_duplicates(['A', 'B'], take_last=True)
        assert_frame_equal(result, expected)

    def test_drop_col_still_multiindex(self):
        arrays = [[  'a',   'b',   'c',    'top'],
                  [  '',    '',    '',     'OD' ],
                  [  '',    '',    '',     'wx' ]]

        tuples = zip(*arrays)
        tuples.sort()
        index = MultiIndex.from_tuples(tuples)

        df = DataFrame(randn(3,4), columns=index)
        del df[('a','','')]
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
                'A' : {'1' : 1, '2' : 2},
                'B' : {'1' : '1', '2' : '2', '3' : '3'},
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

    def test_pivot(self):
        data = {
            'index' : ['A', 'B', 'C', 'C', 'B', 'A'],
            'columns' : ['One', 'One', 'One', 'Two', 'Two', 'Two'],
            'values' : [1., 2., 3., 3., 2., 1.]
        }

        frame = DataFrame(data)
        pivoted = frame.pivot(index='index', columns='columns', values='values')

        expected = DataFrame({
            'One' : {'A' : 1., 'B' : 2., 'C' : 3.},
            'Two' : {'A' : 1., 'B' : 2., 'C' : 3.}
        })

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
        data = DataFrame({'a' : ['bar', 'bar', 'foo', 'foo', 'foo'],
                          'b' : ['one', 'two', 'one', 'one', 'two'],
                          'c' : [1., 2., 3., 3., 4.]})
        self.assertRaises(Exception, data.pivot, 'a', 'b', 'c')

    def test_pivot_empty(self):
        df = DataFrame({}, columns=['a', 'b', 'c'])
        result = df.pivot('a', 'b', 'c')
        expected = DataFrame({})
        assert_frame_equal(result, expected)

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
        self.assert_(not newFrame)
        self.assertEqual(len(newFrame.columns), len(self.frame.columns))

        # length zero with columns reindexed with non-empty index
        newFrame = self.frame.reindex([])
        newFrame = newFrame.reindex(self.frame.index)
        self.assertEqual(len(newFrame.index), len(self.frame.index))
        self.assertEqual(len(newFrame.columns), len(self.frame.columns))

        # pass non-Index
        newFrame = self.frame.reindex(list(self.ts1.index))
        self.assert_(newFrame.index.equals(self.ts1.index))

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
        self.assert_(not newFrame)

    def test_add_index(self):
        df = DataFrame({'A' : ['foo', 'foo', 'foo', 'bar', 'bar'],
                        'B' : ['one', 'two', 'three', 'one', 'two'],
                        'C' : ['a', 'b', 'c', 'd', 'e'],
                        'D' : np.random.randn(5),
                        'E' : np.random.randn(5)})

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
        self.assertRaises(Exception, df.set_index, 'A')

    def test_align(self):
        af, bf = self.frame.align(self.frame)
        self.assert_(af._data is not self.frame._data)

        af, bf = self.frame.align(self.frame, copy=False)
        self.assert_(af._data is self.frame._data)

        # axis = 0
        other = self.frame.ix[:-5, :3]
        af, bf = self.frame.align(other, axis=0)
        self.assert_(bf.columns.equals(other.columns))

        af, bf = self.frame.align(other, join='right', axis=0)
        self.assert_(bf.columns.equals(other.columns))
        self.assert_(bf.index.equals(other.index))
        self.assert_(af.index.equals(other.index))

        # axis = 1
        other = self.frame.ix[:-5, :3].copy()
        af, bf = self.frame.align(other, axis=1)
        self.assert_(bf.columns.equals(self.frame.columns))
        self.assert_(bf.index.equals(other.index))

        af, bf = self.frame.align(other, join='inner', axis=1)
        self.assert_(bf.columns.equals(other.columns))

        # try to align dataframe to series along bad axis
        self.assertRaises(ValueError, self.frame.align, af.ix[0,:3],
                          join='inner', axis=2)

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
            'A' : 'a',
            'B' : 'b',
            'C' : 'c',
            'D' : 'd'
        }
        bad_mapping = {
            'A' : 'a',
            'B' : 'b',
            'C' : 'b',
            'D' : 'd'
        }

        renamed = self.frame.rename(columns=mapping)
        renamed2 = self.frame.rename(columns=str.lower)

        assert_frame_equal(renamed, renamed2)
        assert_frame_equal(renamed2.rename(columns=str.upper),
                           self.frame)

        self.assertRaises(Exception, self.frame.rename,
                          columns=bad_mapping)

        # index

        data = {
            'A' : {'foo' : 0, 'bar' : 1}
        }

        # gets sorted alphabetical
        df = DataFrame(data)
        renamed = df.rename(index={'foo' : 'bar', 'bar' : 'foo'})
        self.assert_(np.array_equal(renamed.index, ['foo', 'bar']))

        renamed = df.rename(index=str.upper)
        self.assert_(np.array_equal(renamed.index, ['BAR', 'FOO']))

        # have to pass something
        self.assertRaises(Exception, self.frame.rename)

        # partial columns
        renamed = self.frame.rename(columns={'C' : 'foo', 'D' : 'bar'})
        self.assert_(np.array_equal(renamed.columns, ['A', 'B', 'foo', 'bar']))

        # other axis
        renamed = self.frame.T.rename(index={'C' : 'foo', 'D' : 'bar'})
        self.assert_(np.array_equal(renamed.index, ['A', 'B', 'foo', 'bar']))

    def test_rename_nocopy(self):
        renamed = self.frame.rename(columns={'C' : 'foo'}, copy=False)
        renamed['foo'] = 1.
        self.assert_((self.frame['C'] == 1.).all())

    #----------------------------------------------------------------------
    # Time series related

    def test_diff(self):
        the_diff = self.tsframe.diff(1)

        assert_series_equal(the_diff['A'],
                            self.tsframe['A'] - self.tsframe['A'].shift(1))

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
        shiftedFrame = self.tsframe.shift(5, offset=datetools.BDay())
        self.assert_(len(shiftedFrame) == len(self.tsframe))

        shiftedFrame2 = self.tsframe.shift(5, timeRule='WEEKDAY')
        assert_frame_equal(shiftedFrame, shiftedFrame2)

        d = self.tsframe.index[0]
        shifted_d = d + datetools.BDay(5)
        assert_series_equal(self.tsframe.xs(d),
                            shiftedFrame.xs(shifted_d))

        # shift int frame
        int_shifted = self.intframe.shift(1)

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
        self.assert_(applied.index is self.frame.index) # want this

        # empty
        applied = self.empty.apply(np.sqrt)
        self.assert_(not applied)

        applied = self.empty.apply(np.mean)
        self.assert_(not applied)

        no_rows = self.frame[:0]
        result = no_rows.apply(lambda x: x.mean())
        expected = Series(np.nan, index=self.frame.columns)
        assert_series_equal(result, expected)

        no_cols = self.frame.ix[:, []]
        result = no_cols.apply(lambda x: x.mean(), axis=1)
        expected = Series(np.nan, index=self.frame.index)
        assert_series_equal(result, expected)

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
        df = DataFrame({'A' : ['foo'],
                        'B' : [1.]})
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

        data['C'][4] = np.nan

        def transform(row):
            if row['C'].startswith('shin') and row['A'] == 'foo':
                row['D'] = 7
            return row

        def transform2(row):
            if (notnull(row['C']) and  row['C'].startswith('shin')
                and row['A'] == 'foo'):
                row['D'] = 7
            return row

        try:
            transformed = data.apply(transform, axis=1)
        except Exception, e:
            self.assertEqual(len(e.args), 2)
            self.assertEqual(e.args[1], 'occurred at index 4')


    def test_apply_convert_objects(self):
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

        result = data.apply(lambda x: x, axis=1)
        assert_frame_equal(result, data)

    def test_applymap(self):
        applied = self.frame.applymap(lambda x: x * 2)
        assert_frame_equal(applied, self.frame * 2)
        result = self.frame.applymap(type)

        # GH #465, function returning tuples
        result = self.frame.applymap(lambda x: (x, x))
        self.assert_(isinstance(result['A'][0], tuple))

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

        # regex
        filtered = fcopy.filter(regex='[A]+')
        self.assertEqual(len(filtered.columns), 2)
        self.assert_('AA' in filtered)

        # pass in None
        self.assertRaises(Exception, self.frame.filter, items=None)

        # objects
        filtered = self.mixed_frame.filter(like='foo')
        self.assert_('foo' in filtered)

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
        assert_frame_equal(result, expected)

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
        sorted_df = frame.sort(column='A')
        expected = frame.sort_index(by='A')
        assert_frame_equal(sorted_df, expected)

        sorted_df = frame.sort(column='A', ascending=False)
        expected = frame.sort_index(by='A', ascending=False)
        assert_frame_equal(sorted_df, expected)

    def test_sort_index_multicolumn(self):
        import random
        A = np.arange(5).repeat(20)
        B = np.tile(np.arange(5), 20)
        random.shuffle(A)
        random.shuffle(B)
        frame = DataFrame({'A' : A, 'B' : B,
                           'C' : np.random.randn(100)})

        result = frame.sort_index(by=['A', 'B'])
        indexer = np.lexsort((frame['B'], frame['A']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

        result = frame.sort_index(by=['A', 'B'], ascending=False)
        expected = frame.take(indexer[::-1])
        assert_frame_equal(result, expected)

        result = frame.sort_index(by=['B', 'A'])
        indexer = np.lexsort((frame['A'], frame['B']))
        expected = frame.take(indexer)
        assert_frame_equal(result, expected)

    def test_frame_column_inplace_sort_exception(self):
        s = self.frame['A']
        self.assertRaises(Exception, s.sort)

        cp = s.copy()
        cp.sort() # it works!

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

    def test_combine_first_mixed_bug(self):
        idx = Index(['a','b','c','e'])
        ser1 = Series([5.0,-9.0,4.0,100.],index=idx)
        ser2 = Series(['a', 'b', 'c', 'e'], index=idx)
        ser3 = Series([12,4,5,97], index=idx)

        frame1 = DataFrame({"col0" : ser1,
                                "col2" : ser2,
                                "col3" : ser3})

        idx = Index(['a','b','c','f'])
        ser1 = Series([5.0,-9.0,4.0,100.], index=idx)
        ser2 = Series(['a','b','c','f'], index=idx)
        ser3 = Series([12,4,5,97],index=idx)

        frame2 = DataFrame({"col1" : ser1,
                                 "col2" : ser2,
                                 "col5" : ser3})


        combined = frame1.combine_first(frame2)
        self.assertEqual(len(combined.columns), 5)

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
        df1 = DataFrame({'x':[5]})
        df2 = DataFrame({'x':[1]})
        df3 = DataFrame({'x':[6]})
        comb = df1.combineAdd(df2)
        assert_frame_equal(comb, df3)

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

    def test_get_X_columns(self):
        # numeric and object columns

        # Booleans get casted to float in DataFrame, so skip for now
        df = DataFrame({'a' : [1, 2, 3],
                         # 'b' : [True, False, True],
                         'c' : ['foo', 'bar', 'baz'],
                         'd' : [None, None, None],
                         'e' : [3.14, 0.577, 2.773]})

        self.assert_(np.array_equal(df._get_numeric_data().columns,
                                    ['a', 'e']))

    def test_get_numeric_data(self):
        df = DataFrame({'a' : 1., 'b' : 2, 'c' : 'foo'},
                       index=np.arange(10))

        result = df._get_numeric_data()
        expected = df.ix[:, ['a', 'b']]
        assert_frame_equal(result, expected)

        only_obj = df.ix[:, ['c']]
        result = only_obj._get_numeric_data()
        expected = df.ix[:, []]
        assert_frame_equal(result, expected)

    def test_count(self):
        f = lambda s: notnull(s).sum()
        self._check_stat_op('count', f,
                            has_skipna=False,
                            has_numeric_only=True)

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
        df = DataFrame({'A' : np.arange(20)}, index=np.arange(20))
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
        df = DataFrame({'A' : np.arange(20)}, index=np.arange(20))
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

    def test_var(self):
        alt = lambda x: np.var(x, ddof=1)
        self._check_stat_op('var', alt)

    def test_std(self):
        alt = lambda x: np.std(x, ddof=1)
        self._check_stat_op('std', alt)

    def test_skew(self):
        from scipy.stats import skew

        def alt(x):
            if len(x) < 3:
                return np.nan
            return skew(x, bias=False)

        self._check_stat_op('skew', alt)

    def _check_stat_op(self, name, alternative, frame=None, has_skipna=True,
                       has_numeric_only=False):
        if frame is None:
            frame = self.frame
            # set some NAs
            frame.ix[5:10] = np.nan
            frame.ix[15:20, -2:] = np.nan

        f = getattr(frame, name)

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
                                check_dtype=False) # HACK: win32
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

        self._check_stat_op('median', wrapper, frame=self.intframe)

    def test_quantile(self):
        try:
            from scipy.stats import scoreatpercentile
        except ImportError:
            return

        q = self.tsframe.quantile(0.1, axis=0)
        self.assertEqual(q['A'], scoreatpercentile(self.tsframe['A'], 10))
        q = self.tsframe.quantile(0.9, axis=1)
        q = self.intframe.quantile(0.1)
        self.assertEqual(q['A'], scoreatpercentile(self.intframe['A'], 10))

        # test degenerate case
        q = DataFrame({'x':[],'y':[]}).quantile(0.1, axis=0)
        assert(np.isnan(q['x']) and np.isnan(q['y']))

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
        df = DataFrame({'A' : np.arange(20)}, index=np.arange(20))
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
        df = self.tsframe.astype(int)
        df.cumprod(0)
        df.cumprod(1)

    def test_rank(self):
        from scipy.stats import rankdata

        self.frame['A'][::2] = np.nan
        self.frame['B'][::3] = np.nan
        self.frame['C'][::4] = np.nan
        self.frame['D'][::5] = np.nan


        ranks0 = self.frame.rank()
        ranks1 = self.frame.rank(1)
        mask =  np.isnan(self.frame.values)

        fvals = self.frame.fillna(np.inf).values

        exp0 = np.apply_along_axis(rankdata, 0, fvals)
        exp0[mask] = np.nan

        exp1 = np.apply_along_axis(rankdata, 1, fvals)
        exp1[mask] = np.nan

        assert_almost_equal(ranks0.values, exp0)
        assert_almost_equal(ranks1.values, exp1)

    def test_rank2(self):
        from datetime import datetime

        df = DataFrame([['b','c','a'],['a','c','b']])
        expected = DataFrame([[2.0, 3.0, 1.0], [1, 3, 2]])
        result = df.rank(1, numeric_only=False)
        assert_frame_equal(result, expected)

        expected = DataFrame([[2.0, 1.5, 1.0], [1, 1.5, 2]])
        result = df.rank(0, numeric_only=False)
        assert_frame_equal(result, expected)

        df = DataFrame([['b',np.nan,'a'],['a','c','b']])
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

    def test_describe(self):
        desc = self.tsframe.describe()
        desc = self.mixed_frame.describe()
        desc = self.frame.describe()

    def test_describe_no_numeric(self):
        df = DataFrame({'A' : ['foo', 'foo', 'bar'] * 8,
                        'B' : ['a', 'b', 'c', 'd'] * 6})
        desc = df.describe()
        expected = DataFrame(dict((k, v.describe())
                                  for k, v in df.iteritems()),
                             columns=df.columns)
        assert_frame_equal(desc, expected)

    def test_get_axis_etc(self):
        f = self.frame

        self.assertEquals(f._get_axis_number(0), 0)
        self.assertEquals(f._get_axis_number(1), 1)
        self.assertEquals(f._get_axis_name(0), 'index')
        self.assertEquals(f._get_axis_name(1), 'columns')

        self.assert_(f._get_axis(0) is f.index)
        self.assert_(f._get_axis(1) is f.columns)
        self.assertRaises(Exception, f._get_axis_number, 2)

    def test_combine_first_mixed(self):
        a = Series(['a','b'], index=range(2))
        b = Series(range(2), index=range(2))
        f = DataFrame({'A' : a, 'B' : b})

        a = Series(['a','b'], index=range(5, 7))
        b = Series(range(2), index=range(5, 7))
        g = DataFrame({'A' : a, 'B' : b})

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
        stacked_df = DataFrame({'foo' : stacked, 'bar' : stacked})

        unstacked = stacked.unstack()
        unstacked_df = stacked_df.unstack()

        assert_frame_equal(unstacked, self.frame)
        assert_frame_equal(unstacked_df['bar'], self.frame)

        unstacked_cols = stacked.unstack(0)
        unstacked_cols_df = stacked_df.unstack(0)
        assert_frame_equal(unstacked_cols.T, self.frame)
        assert_frame_equal(unstacked_cols_df['bar'].T, self.frame)

    def test_unstack_to_series(self):
        # check reversibility
        data = self.frame.unstack()

        self.assertTrue(isinstance(data, Series))
        undo = data.unstack().T
        assert_frame_equal(undo, self.frame)

        # check NA handling
        data = DataFrame({'x': [1, 2, np.NaN], 'y': [3.0, 4, np.NaN]})
        data.index = Index(['a','b','c'])
        result = data.unstack()

        midx = MultiIndex(levels=[['x','y'],['a','b','c']],
                          labels=[[0,0,0,1,1,1],[0,1,2,0,1,2]])
        expected = Series([1,2,np.NaN,3,4,np.NaN], index=midx)

        assert_series_equal(result, expected)

        # check composability of unstack
        old_data = data.copy()
        for _ in xrange(4):
            data = data.unstack()
        assert_frame_equal(old_data, data)

    def test_reset_index(self):
        stacked = self.frame.stack()[::2]
        stacked = DataFrame({'foo' : stacked, 'bar' : stacked})

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

    def test_reset_index_right_dtype(self):
        time = np.arange(0.0, 10, np.sqrt(2)/2)
        s1 = Series((9.81 * time ** 2) /2,
                    index=Index(time, name='time'),
                    name='speed')
        df = DataFrame(s1)

        resetted = s1.reset_index()
        self.assert_(resetted['time'].dtype == np.float64)

        resetted = df.reset_index()
        self.assert_(resetted['time'].dtype == np.float64)

    #----------------------------------------------------------------------
    # Tests to cope with refactored internals

    def test_as_matrix_numeric_cols(self):
        self.frame['foo'] = 'bar'

        values = self.frame.as_matrix(['A', 'B', 'C', 'D'])
        self.assert_(values.dtype == np.float64)

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

        df = DataFrame({'A' : series['A']})
        df['A'][:] = 5

        self.assert_(not (series['A'] == 5).all())

    def test_assign_columns(self):
        self.frame['hi'] = 'there'

        frame = self.frame.copy()
        frame.columns = ['foo', 'bar', 'baz', 'quux', 'foo2']
        assert_series_equal(self.frame['C'], frame['baz'])
        assert_series_equal(self.frame['hi'], frame['foo2'])

    def test_cast_internals(self):
        casted = DataFrame(self.frame._data, dtype=int)
        expected = DataFrame(self.frame._series, dtype=int)
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

    def test_boolean_set_mixed_type(self):
        bools = self.mixed_frame.applymap(lambda x: x != 2).astype(bool)
        self.assertRaises(Exception, self.mixed_frame.__setitem__, bools, 2)

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
        cols = range(3)
        df1 = DataFrame(index=idx, columns=cols, \
                           data=np.array([[0.0, 0.5, 1.0],
                                          [1.5, 2.0, 2.5],
                                          [3.0, 3.5, 4.0]], dtype=float))
        df2 = DataFrame(index=idx, columns=cols, data=np.ones((len(idx), len(cols))))

        expected = DataFrame(index=idx, columns=cols, \
                           data=np.array([[0.0, 0.5, 1.0],
                                          [1.5, 2.0, -1],
                                          [-1,  -1,  -1]], dtype=float))

        df1[df1 > 2.0 * df2] = -1
        assert_frame_equal(df1, expected)

    def test_sum_bools(self):
        df = DataFrame(index=range(1), columns=range(10))
        bools = np.isnan(df)
        self.assert_(bools.sum(axis=1)[0] == 10)

    def test_fillna_col_reordering(self):
        idx = range(20)
        cols = ["COL." + str(i) for i in range(5, 0, -1)]
        data = np.random.rand(20, 5)
        df = DataFrame(index=range(20), columns=cols, data=data)
        self.assert_(df.columns.tolist() == df.fillna().columns.tolist())

    def test_take(self):
        # homogeneous
        #----------------------------------------

        # mixed-dtype
        #----------------------------------------
        order = [4, 1, 2, 0, 3]

        result = self.mixed_frame.take(order, axis=0)
        expected = self.mixed_frame.reindex(self.mixed_frame.index.take(order))
        assert_frame_equal(result, expected)

        # axis = 1
        result = self.mixed_frame.take(order, axis=1)
        expected = self.mixed_frame.ix[:, ['foo', 'B', 'C', 'A', 'D']]
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
        assert_frame_equal(result, expected)

    def test_idxmin(self):
        frame = self.frame
        frame.ix[5:10] = np.nan
        frame.ix[15:20, -2:] = np.nan
        for skipna in [True, False]:
            for axis in [0, 1]:
                for df in [frame, self.intframe]:
                    result = df.idxmin(axis=axis, skipna=skipna)
                    expected = df.apply(Series.idxmin, axis=axis, skipna=skipna)
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
                    expected = df.apply(Series.idxmax, axis=axis, skipna=skipna)
                    assert_series_equal(result, expected)

        self.assertRaises(Exception, frame.idxmax, axis=2)

    def test_stale_cached_series_bug_473(self):
        Y = DataFrame(np.random.random((4, 4)), index=('a', 'b','c','d'),
                      columns=('e','f','g','h'))
        repr(Y)
        Y['e'] = Y['e'].astype('object')
        Y['g']['c'] = np.NaN
        repr(Y)
        result = Y.sum()
        exp = Y['g'].sum()
        self.assert_(isnull(Y['g']['c']))

if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--ipdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
