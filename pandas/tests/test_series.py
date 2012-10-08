# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta, date
import os
import operator
import unittest

import nose

from numpy import nan
import numpy as np
import numpy.ma as ma

from pandas import (Index, Series, TimeSeries, DataFrame, isnull, notnull,
                    bdate_range, date_range)
from pandas.core.index import MultiIndex
from pandas.tseries.index import Timestamp, DatetimeIndex
import pandas.core.series as smod
import pandas.lib as lib

import pandas.core.datetools as datetools
import pandas.core.nanops as nanops

from pandas.util.py3compat import StringIO
from pandas.util import py3compat
from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

def _skip_if_no_scipy():
    try:
        import scipy.stats
    except ImportError:
        raise nose.SkipTest


#-------------------------------------------------------------------------------
# Series test cases

JOIN_TYPES = ['inner', 'outer', 'left', 'right']

class CheckNameIntegration(object):

    def test_scalarop_preserve_name(self):
        result = self.ts * 2
        self.assertEquals(result.name, self.ts.name)

    def test_copy_name(self):
        result = self.ts.copy()
        self.assertEquals(result.name, self.ts.name)

    # def test_copy_index_name_checking(self):
    #     # don't want to be able to modify the index stored elsewhere after
    #     # making a copy

    #     self.ts.index.name = None
    #     cp = self.ts.copy()
    #     cp.index.name = 'foo'
    #     self.assert_(self.ts.index.name is None)

    def test_append_preserve_name(self):
        result = self.ts[:5].append(self.ts[5:])
        self.assertEquals(result.name, self.ts.name)

    def test_binop_maybe_preserve_name(self):
        # names match, preserve
        result = self.ts * self.ts
        self.assertEquals(result.name, self.ts.name)

        result = self.ts * self.ts[:-2]
        self.assertEquals(result.name, self.ts.name)

        # names don't match, don't preserve
        cp = self.ts.copy()
        cp.name = 'something else'
        result = self.ts + cp
        self.assert_(result.name is None)

    def test_combine_first_name(self):
        result = self.ts.combine_first(self.ts[:5])
        self.assertEquals(result.name, self.ts.name)

    def test_getitem_preserve_name(self):
        result = self.ts[self.ts > 0]
        self.assertEquals(result.name, self.ts.name)

        result = self.ts[[0, 2, 4]]
        self.assertEquals(result.name, self.ts.name)

        result = self.ts[5:10]
        self.assertEquals(result.name, self.ts.name)

    def test_getitem_setitem_ellipsis(self):
        s = Series(np.random.randn(10))

        np.fix(s)

        result = s[...]
        assert_series_equal(result, s)

        s[...] = 5
        self.assert_((result == 5).all())

    def test_multilevel_name_print(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        s = Series(range(0,len(index)), index=index, name='sth')
        expected = ["first  second",
                    "foo    one       0",
                    "       two       1",
                    "       three     2",
                    "bar    one       3",
                    "       two       4",
                    "baz    two       5",
                    "       three     6",
                    "qux    one       7",
                    "       two       8",
                    "       three     9",
                    "Name: sth"]
        expected = "\n".join(expected)
        self.assertEquals(repr(s), expected)

    def test_multilevel_preserve_name(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        s = Series(np.random.randn(len(index)), index=index, name='sth')

        result = s['foo']
        result2 = s.ix['foo']
        self.assertEquals(result.name, s.name)
        self.assertEquals(result2.name, s.name)

    def test_name_printing(self):
        # test small series
        s = Series([0, 1, 2])
        s.name = "test"
        self.assert_("Name: test" in repr(s))
        s.name = None
        self.assert_(not "Name:" in repr(s))
        # test big series (diff code path)
        s = Series(range(0,1000))
        s.name = "test"
        self.assert_("Name: test" in repr(s))
        s.name = None
        self.assert_(not "Name:" in repr(s))

    def test_pickle_preserve_name(self):
        unpickled = self._pickle_roundtrip(self.ts)
        self.assertEquals(unpickled.name, self.ts.name)

    def _pickle_roundtrip(self, obj):
        obj.save('__tmp__')
        unpickled = Series.load('__tmp__')
        os.remove('__tmp__')
        return unpickled

    def test_argsort_preserve_name(self):
        result = self.ts.argsort()
        self.assertEquals(result.name, self.ts.name)

    def test_sort_index_name(self):
        result = self.ts.sort_index(ascending=False)
        self.assertEquals(result.name, self.ts.name)

    def test_to_sparse_pass_name(self):
        result = self.ts.to_sparse()
        self.assertEquals(result.name, self.ts.name)

class TestNanops(unittest.TestCase):

    def test_comparisons(self):
        left = np.random.randn(10)
        right = np.random.randn(10)
        left[:3] = np.nan

        result = nanops.nangt(left, right)
        expected = (left > right).astype('O')
        expected[:3] = np.nan

        assert_almost_equal(result, expected)

    def test_none_comparison(self):
        # bug brought up by #1079
        s = Series(np.random.randn(10), index=range(0, 20, 2))
        self.assertRaises(TypeError, s.__eq__, None)

    def test_sum_zero(self):
        arr = np.array([])
        self.assert_(nanops.nansum(arr) == 0)

        arr = np.empty((10, 0))
        self.assert_((nanops.nansum(arr, axis=1) == 0).all())

        # GH #844
        s = Series([], index=[])
        self.assert_(s.sum() == 0)

        df = DataFrame(np.empty((10, 0)))
        self.assert_((df.sum(1) == 0).all())

    def test_nansum_buglet(self):
        s = Series([1.0, np.nan], index=[0,1])
        result = np.nansum(s)
        assert_almost_equal(result, 1)

class SafeForSparse(object):
    pass

class TestSeries(unittest.TestCase, CheckNameIntegration):

    def setUp(self):
        self.ts = tm.makeTimeSeries()
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

        self.objSeries = tm.makeObjectSeries()
        self.objSeries.name = 'objects'

        self.empty = Series([], index=[])

    def test_constructor(self):
        # Recognize TimeSeries
        self.assert_(isinstance(self.ts, TimeSeries))

        # Pass in Series
        derived = Series(self.ts)
        self.assert_(isinstance(derived, TimeSeries))

        self.assert_(tm.equalContents(derived.index, self.ts.index))
        # Ensure new index is not created
        self.assertEquals(id(self.ts.index), id(derived.index))

        # Pass in scalar
        scalar = Series(0.5)
        self.assert_(isinstance(scalar, float))

        # Mixed type Series
        mixed = Series(['hello', np.NaN], index=[0, 1])
        self.assert_(mixed.dtype == np.object_)
        self.assert_(mixed[1] is np.NaN)

        self.assert_(not isinstance(self.empty, TimeSeries))
        self.assert_(not isinstance(Series({}), TimeSeries))

        self.assertRaises(Exception, Series, np.random.randn(3, 3),
                          index=np.arange(3))

        mixed.name = 'Series'
        rs = Series(mixed).name
        xp = 'Series'
        self.assertEqual(rs, xp)

    def test_constructor_empty(self):
        empty = Series()
        empty2 = Series([])
        assert_series_equal(empty, empty2)

        empty = Series(index=range(10))
        empty2 = Series(np.nan, index=range(10))
        assert_series_equal(empty, empty2)

    def test_constructor_series(self):
        index1 = ['d', 'b', 'a', 'c']
        index2 = sorted(index1)
        s1 = Series([4, 7, -5, 3], index=index1)
        s2 = Series(s1, index=index2)

        assert_series_equal(s2, s1.sort_index())

    def test_constructor_generator(self):
        gen = (i for i in range(10))

        result = Series(gen)
        exp = Series(range(10))
        assert_series_equal(result, exp)

        gen = (i for i in range(10))
        result = Series(gen, index=range(10, 20))
        exp.index = range(10, 20)
        assert_series_equal(result, exp)

    def test_constructor_maskedarray(self):
        data = ma.masked_all((3,), dtype=float)
        result = Series(data)
        expected = Series([nan, nan, nan])
        assert_series_equal(result, expected)

        data[0] = 0.0
        data[2] = 2.0
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([0.0, nan, 2.0], index=index)
        assert_series_equal(result, expected)

    def test_constructor_default_index(self):
        s = Series([0, 1, 2])
        assert_almost_equal(s.index, np.arange(3))

    def test_constructor_corner(self):
        df = tm.makeTimeDataFrame()
        objs = [df, df]
        s = Series(objs, index=[0, 1])
        self.assert_(isinstance(s, Series))

    def test_constructor_sanitize(self):
        s = Series(np.array([1., 1., 8.]), dtype='i8')
        self.assertEquals(s.dtype, np.dtype('i8'))

        s = Series(np.array([1., 1., np.nan]), copy=True, dtype='i8')
        self.assertEquals(s.dtype, np.dtype('f8'))

    def test_constructor_pass_none(self):
        s = Series(None, index=range(5))
        self.assert_(s.dtype == np.float64)

        s = Series(None, index=range(5), dtype=object)
        self.assert_(s.dtype == np.object_)

    def test_constructor_cast(self):
        self.assertRaises(ValueError, Series, ['a', 'b', 'c'], dtype=float)

    def test_constructor_dtype_nocast(self):
        # #1572
        s = Series([1, 2, 3])

        s2 = Series(s, dtype=np.int64)

        s2[1] = 5
        self.assertEquals(s[1], 5)

    def test_constructor_dict(self):
        d = {'a' : 0., 'b' : 1., 'c' : 2.}
        result = Series(d, index=['b', 'c', 'd', 'a'])
        expected = Series([1, 2, nan, 0], index=['b', 'c', 'd', 'a'])
        assert_series_equal(result, expected)

        pidx = tm.makePeriodIndex(100)
        d = {pidx[0] : 0, pidx[1] : 1}
        result = Series(d, index=pidx)
        expected = Series(np.nan, pidx)
        expected.ix[0] = 0
        expected.ix[1] = 1
        assert_series_equal(result, expected)

    def test_constructor_subclass_dict(self):
        data = tm.TestSubDict((x, 10.0 * x) for x in xrange(10))
        series = Series(data)
        refseries = Series(dict(data.iteritems()))
        assert_series_equal(refseries, series)

    def test_constructor_list_of_tuples(self):
        data = [(1, 1), (2, 2), (2, 3)]
        s = Series(data)
        self.assertEqual(list(s), data)

    def test_constructor_tuple_of_tuples(self):
        data = ((1, 1), (2, 2), (2, 3))
        s = Series(data)
        self.assertEqual(tuple(s), data)

    def test_constructor_set(self):
        values = set([1, 2, 3, 4, 5])

        self.assertRaises(TypeError, Series, values)

    def test_fromDict(self):
        data = {'a' : 0, 'b' : 1, 'c' : 2, 'd' : 3}

        series = Series(data)
        self.assert_(tm.is_sorted(series.index))

        data = {'a' : 0, 'b' : '1', 'c' : '2', 'd' : datetime.now()}
        series = Series(data)
        self.assert_(series.dtype == np.object_)

        data = {'a' : 0, 'b' : '1', 'c' : '2', 'd' : '3'}
        series = Series(data)
        self.assert_(series.dtype == np.object_)

        data = {'a' : '0', 'b' : '1'}
        series = Series(data, dtype=float)
        self.assert_(series.dtype == np.float64)

    def test_from_json_to_json(self):
        raise nose.SkipTest

        def _check_orient(series, orient, dtype=None, numpy=True):
            series = series.sort_index()
            unser = Series.from_json(series.to_json(orient=orient),
                                     orient=orient, numpy=numpy, dtype=dtype)
            unser = unser.sort_index()
            if series.index.dtype.type == np.datetime64:
                unser.index = DatetimeIndex(unser.index.values.astype('i8'))
            if orient == "records" or orient == "values":
                assert_almost_equal(series.values, unser.values)
            else:
                try:
                    assert_series_equal(series, unser)
                except:
                    raise
                if orient == "split":
                    self.assert_(series.name == unser.name)

        def _check_all_orients(series, dtype=None):
            _check_orient(series, "columns", dtype=dtype)
            _check_orient(series, "records", dtype=dtype)
            _check_orient(series, "split", dtype=dtype)
            _check_orient(series, "index", dtype=dtype)
            _check_orient(series, "values", dtype=dtype)

            _check_orient(series, "columns", dtype=dtype, numpy=False)
            _check_orient(series, "records", dtype=dtype, numpy=False)
            _check_orient(series, "split", dtype=dtype, numpy=False)
            _check_orient(series, "index", dtype=dtype, numpy=False)
            _check_orient(series, "values", dtype=dtype, numpy=False)

        # basic
        _check_all_orients(self.series)
        self.assertEqual(self.series.to_json(),
                         self.series.to_json(orient="index"))

        objSeries = Series([str(d) for d in self.objSeries],
                           index=self.objSeries.index,
                           name=self.objSeries.name)
        _check_all_orients(objSeries)
        _check_all_orients(self.empty)
        _check_all_orients(self.ts)

        # dtype
        s = Series(range(6), index=['a','b','c','d','e','f'])
        _check_all_orients(Series(s, dtype=np.float64), dtype=np.float64)
        _check_all_orients(Series(s, dtype=np.int), dtype=np.int)

    def test_to_json_except(self):
        raise nose.SkipTest
        s = Series([1, 2, 3])
        self.assertRaises(ValueError, s.to_json, orient="garbage")

    def test_setindex(self):
        # wrong type
        series = self.series.copy()
        self.assertRaises(TypeError, setattr, series, 'index', None)

        # wrong length
        series = self.series.copy()
        self.assertRaises(AssertionError, setattr, series, 'index',
                          np.arange(len(series) - 1))

        # works
        series = self.series.copy()
        series.index = np.arange(len(series))
        self.assert_(isinstance(series.index, Index))

    def test_array_finalize(self):
        pass

    def test_fromValue(self):
        nans = Series(np.NaN, index=self.ts.index)
        self.assert_(nans.dtype == np.float_)
        self.assertEqual(len(nans), len(self.ts))

        strings = Series('foo', index=self.ts.index)
        self.assert_(strings.dtype == np.object_)
        self.assertEqual(len(strings), len(self.ts))

        d = datetime.now()
        dates = Series(d, index=self.ts.index)
        self.assert_(dates.dtype == np.object_)
        self.assertEqual(len(dates), len(self.ts))

    def test_contains(self):
        tm.assert_contains_all(self.ts.index, self.ts)

    def test_pickle(self):
        unp_series = self._pickle_roundtrip(self.series)
        unp_ts = self._pickle_roundtrip(self.ts)
        assert_series_equal(unp_series, self.series)
        assert_series_equal(unp_ts, self.ts)

    def _pickle_roundtrip(self, obj):
        obj.save('__tmp__')
        unpickled = Series.load('__tmp__')
        os.remove('__tmp__')
        return unpickled

    def test_getitem_get(self):
        idx1 = self.series.index[5]
        idx2 = self.objSeries.index[5]

        self.assertEqual(self.series[idx1], self.series.get(idx1))
        self.assertEqual(self.objSeries[idx2], self.objSeries.get(idx2))

        self.assertEqual(self.series[idx1], self.series[5])
        self.assertEqual(self.objSeries[idx2], self.objSeries[5])

        self.assert_(self.series.get(-1) is None)
        self.assertEqual(self.series[5], self.series.get(self.series.index[5]))

        # missing
        d = self.ts.index[0] - datetools.bday
        self.assertRaises(KeyError, self.ts.__getitem__, d)

    def test_iget(self):
        s = Series(np.random.randn(10), index=range(0, 20, 2))
        for i in range(len(s)):
            result = s.iget(i)
            exp = s[s.index[i]]
            assert_almost_equal(result, exp)

        # pass a slice
        result = s.iget(slice(1, 3))
        expected = s.ix[2:4]
        assert_series_equal(result, expected)

        # test slice is a view
        result[:] = 0
        self.assert_((s[1:3] == 0).all())

        # list of integers
        result = s.iget([0, 2, 3, 4, 5])
        expected = s.reindex(s.index[[0, 2, 3, 4, 5]])
        assert_series_equal(result, expected)

    def test_iget_nonunique(self):
        s = Series([0, 1, 2], index=[0, 1, 0])
        self.assertEqual(s.iget(2), 2)

    def test_getitem_regression(self):
        s = Series(range(5), index=range(5))
        result = s[range(5)]
        assert_series_equal(result, s)

    def test_getitem_setitem_slice_bug(self):
        s = Series(range(10), range(10))
        result = s[-12:]
        assert_series_equal(result, s)

        result = s[-7:]
        assert_series_equal(result, s[3:])

        result = s[:-12]
        assert_series_equal(result, s[:0])

        s = Series(range(10), range(10))
        s[-12:] = 0
        self.assert_((s == 0).all())

        s[:-12] = 5
        self.assert_((s == 0).all())

    def test_getitem_int64(self):
        idx = np.int64(5)
        self.assertEqual(self.ts[idx], self.ts[5])

    def test_getitem_fancy(self):
        slice1 = self.series[[1,2,3]]
        slice2 = self.objSeries[[1,2,3]]
        self.assertEqual(self.series.index[2], slice1.index[1])
        self.assertEqual(self.objSeries.index[2], slice2.index[1])
        self.assertEqual(self.series[2], slice1[1])
        self.assertEqual(self.objSeries[2], slice2[1])

    def test_getitem_boolean(self):
        s = self.series
        mask = s > s.median()

        # passing list is OK
        result = s[list(mask)]
        expected = s[mask]
        assert_series_equal(result, expected)
        self.assert_(np.array_equal(result.index, s.index[mask]))

    def test_getitem_generator(self):
        gen = (x > 0 for x in self.series)
        result = self.series[gen]
        result2 = self.series[iter(self.series > 0)]
        expected = self.series[self.series > 0]
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

    def test_getitem_boolean_object(self):
        # using column from DataFrame
        s = self.series
        mask = s > s.median()
        omask = mask.astype(object)

        # getitem
        result = s[omask]
        expected = s[mask]
        assert_series_equal(result, expected)

        # setitem
        cop = s.copy()
        cop[omask] = 5
        s[mask] = 5
        assert_series_equal(cop, s)

        # nans raise exception
        omask[5:10] = np.nan
        self.assertRaises(Exception, s.__getitem__, omask)
        self.assertRaises(Exception, s.__setitem__, omask, 5)

    def test_getitem_setitem_boolean_corner(self):
        ts = self.ts
        mask_shifted = ts.shift(1, freq=datetools.bday) > ts.median()
        self.assertRaises(Exception, ts.__getitem__, mask_shifted)
        self.assertRaises(Exception, ts.__setitem__, mask_shifted, 1)

        self.assertRaises(Exception, ts.ix.__getitem__, mask_shifted)
        self.assertRaises(Exception, ts.ix.__setitem__, mask_shifted, 1)

    def test_getitem_setitem_slice_integers(self):
        s = Series(np.random.randn(8), index=[2, 4, 6, 8, 10, 12, 14, 16])

        result = s[:4]
        expected = s.reindex([2, 4, 6, 8])
        assert_series_equal(result, expected)

        s[:4] = 0
        self.assert_((s[:4] == 0).all())
        self.assert_(not (s[4:] == 0).any())

    def test_getitem_out_of_bounds(self):
        # don't segfault, GH #495
        self.assertRaises(IndexError, self.ts.__getitem__, len(self.ts))

        # GH #917
        s = Series([])
        self.assertRaises(IndexError, s.__getitem__, -1)

    def test_getitem_setitem_integers(self):
        # caused bug without test
        s = Series([1,2,3], ['a','b','c'])

        self.assertEqual(s.ix[0], s['a'])
        s.ix[0] = 5
        self.assertAlmostEqual(s['a'], 5)

    def test_getitem_box_float64(self):
        value = self.ts[5]
        self.assert_(isinstance(value, np.float64))

    def test_getitem_ambiguous_keyerror(self):
        s = Series(range(10), index=range(0, 20, 2))
        self.assertRaises(KeyError, s.__getitem__, 1)
        self.assertRaises(KeyError, s.ix.__getitem__, 1)

    def test_getitem_unordered_dup(self):
        obj = Series(range(5), index=['c', 'a', 'a', 'b', 'b'])
        self.assert_(np.isscalar(obj['c']))
        self.assert_(obj['c'] == 0)

    def test_setitem_ambiguous_keyerror(self):
        s = Series(range(10), index=range(0, 20, 2))
        self.assertRaises(KeyError, s.__setitem__, 1, 5)
        self.assertRaises(KeyError, s.ix.__setitem__, 1, 5)

    def test_setitem_float_labels(self):
        # note labels are floats
        s = Series(['a','b','c'],index=[0,0.5,1])
        tmp = s.copy()

        s.ix[1] = 'zoo'
        tmp.values[1] = 'zoo'

        assert_series_equal(s, tmp)

    def test_slice(self):
        numSlice = self.series[10:20]
        numSliceEnd = self.series[-10:]
        objSlice = self.objSeries[10:20]

        self.assert_(self.series.index[9] not in numSlice.index)
        self.assert_(self.objSeries.index[9] not in objSlice.index)

        self.assertEqual(len(numSlice), len(numSlice.index))
        self.assertEqual(self.series[numSlice.index[0]],
                         numSlice[numSlice.index[0]])

        self.assertEqual(numSlice.index[1], self.series.index[11])

        self.assert_(tm.equalContents(numSliceEnd,
                                          np.array(self.series)[-10:]))

        # test return view
        sl = self.series[10:20]
        sl[:] = 0
        self.assert_((self.series[10:20] == 0).all())

    def test_slice_can_reorder_not_uniquely_indexed(self):
        s = Series(1, index=['a', 'a', 'b', 'b', 'c'])
        result = s[::-1] # it works!

    def test_slice_float_get_set(self):
        result = self.ts[4.0:10.0]
        expected = self.ts[4:10]
        assert_series_equal(result, expected)

        self.ts[4.0:10.0] = 0
        self.assert_((self.ts[4:10] == 0).all())

        self.assertRaises(TypeError, self.ts.__getitem__, slice(4.5, 10.0))
        self.assertRaises(TypeError, self.ts.__setitem__, slice(4.5, 10.0), 0)

    def test_slice_floats2(self):
        s = Series(np.random.rand(10), index=np.arange(10,20,dtype=float))

        self.assert_(len(s.ix[12.0:]) == 8)
        self.assert_(len(s.ix[12.5:]) == 7)

        i = np.arange(10,20,dtype=float)
        i[2] = 12.2
        s.index = i
        self.assert_(len(s.ix[12.0:]) == 8)
        self.assert_(len(s.ix[12.5:]) == 7)

    def test_slice_float64(self):
        values = np.arange(10., 50., 2)
        index = Index(values)

        start, end = values[[5, 15]]

        s = Series(np.random.randn(20), index=index)

        result = s[start:end]
        expected = s.ix[5:16]
        assert_series_equal(result, expected)

        result = s.ix[start:end]
        assert_series_equal(result, expected)

        df = DataFrame(np.random.randn(20, 3), index=index)

        result = df[start:end]
        expected = df.ix[5:16]
        tm.assert_frame_equal(result, expected)

        result = df.ix[start:end]
        tm.assert_frame_equal(result, expected)

    def test_setitem(self):
        self.ts[self.ts.index[5]] = np.NaN
        self.ts[[1,2,17]] = np.NaN
        self.ts[6] = np.NaN
        self.assert_(np.isnan(self.ts[6]))
        self.assert_(np.isnan(self.ts[2]))
        self.ts[np.isnan(self.ts)] = 5
        self.assert_(not np.isnan(self.ts[2]))

        # caught this bug when writing tests
        series = Series(tm.makeIntIndex(20).astype(float),
                        index=tm.makeIntIndex(20))

        series[::2] = 0
        self.assert_((series[::2] == 0).all())

        # set item that's not contained
        self.assertRaises(Exception, self.series.__setitem__,
                          'foobar', 1)

    def test_set_value(self):
        idx = self.ts.index[10]
        res = self.ts.set_value(idx, 0)
        self.assert_(res is self.ts)
        self.assertEqual(self.ts[idx], 0)

        res = self.series.set_value('foobar', 0)
        self.assert_(res is not self.series)
        self.assert_(res.index[-1] == 'foobar')
        self.assertEqual(res['foobar'], 0)

    def test_setslice(self):
        sl = self.ts[5:20]
        self.assertEqual(len(sl), len(sl.index))
        self.assert_(sl.index.is_unique)

    def test_basic_getitem_setitem_corner(self):
        # invalid tuples, e.g. self.ts[:, None] vs. self.ts[:, 2]
        self.assertRaises(Exception, self.ts.__getitem__,
                          (slice(None, None), 2))
        self.assertRaises(Exception, self.ts.__setitem__,
                          (slice(None, None), 2), 2)

        # weird lists. [slice(0, 5)] will work but not two slices
        result = self.ts[[slice(None, 5)]]
        expected = self.ts[:5]
        assert_series_equal(result, expected)

        # OK
        self.assertRaises(Exception, self.ts.__getitem__,
                          [5, slice(None, None)])
        self.assertRaises(Exception, self.ts.__setitem__,
                          [5, slice(None, None)], 2)

    def test_reshape_non_2d(self):
        x = Series(np.random.random(201), name='x')
        self.assertRaises(TypeError, x.reshape, (len(x),))

    def test_reshape_2d_return_array(self):
        x = Series(np.random.random(201), name='x')
        result = x.reshape((-1,1))
        self.assert_(not isinstance(result, Series))

        result2 = np.reshape(x, (-1,1))
        self.assert_(not isinstance(result, Series))

        result = x[:, None]
        expected = x.reshape((-1,1))
        assert_almost_equal(result, expected)

    def test_basic_getitem_with_labels(self):
        indices = self.ts.index[[5, 10, 15]]

        result = self.ts[indices]
        expected = self.ts.reindex(indices)
        assert_series_equal(result, expected)

        result = self.ts[indices[0]:indices[2]]
        expected = self.ts.ix[indices[0]:indices[2]]
        assert_series_equal(result, expected)

        # integer indexes, be careful
        s = Series(np.random.randn(10), index=range(0, 20, 2))
        inds = [0, 2, 5, 7, 8]
        arr_inds = np.array([0, 2, 5, 7, 8])
        result = s[inds]
        expected = s.reindex(inds)
        assert_series_equal(result, expected)

        result = s[arr_inds]
        expected = s.reindex(arr_inds)
        assert_series_equal(result, expected)

    def test_basic_setitem_with_labels(self):
        indices = self.ts.index[[5, 10, 15]]

        cp = self.ts.copy()
        exp = self.ts.copy()
        cp[indices] = 0
        exp.ix[indices] = 0
        assert_series_equal(cp, exp)

        cp = self.ts.copy()
        exp = self.ts.copy()
        cp[indices[0]:indices[2]] = 0
        exp.ix[indices[0]:indices[2]] = 0
        assert_series_equal(cp, exp)

        # integer indexes, be careful
        s = Series(np.random.randn(10), index=range(0, 20, 2))
        inds = [0, 4, 6]
        arr_inds = np.array([0, 4, 6])

        cp = s.copy()
        exp = s.copy()
        s[inds] = 0
        s.ix[inds] = 0
        assert_series_equal(cp, exp)

        cp = s.copy()
        exp = s.copy()
        s[arr_inds] = 0
        s.ix[arr_inds] = 0
        assert_series_equal(cp, exp)

        inds_notfound = [0, 4, 5, 6]
        arr_inds_notfound = np.array([0, 4, 5, 6])
        self.assertRaises(Exception, s.__setitem__, inds_notfound, 0)
        self.assertRaises(Exception, s.__setitem__, arr_inds_notfound, 0)

    def test_ix_getitem(self):
        inds = self.series.index[[3,4,7]]
        assert_series_equal(self.series.ix[inds], self.series.reindex(inds))
        assert_series_equal(self.series.ix[5::2], self.series[5::2])

        # slice with indices
        d1, d2 = self.ts.index[[5, 15]]
        result = self.ts.ix[d1:d2]
        expected = self.ts.truncate(d1, d2)
        assert_series_equal(result, expected)

        # boolean
        mask = self.series > self.series.median()
        assert_series_equal(self.series.ix[mask], self.series[mask])

        # ask for index value
        self.assertEquals(self.ts.ix[d1], self.ts[d1])
        self.assertEquals(self.ts.ix[d2], self.ts[d2])

    def test_ix_getitem_not_monotonic(self):
        d1, d2 = self.ts.index[[5, 15]]

        ts2 = self.ts[::2][::-1]

        self.assertRaises(KeyError, ts2.ix.__getitem__, slice(d1, d2))
        self.assertRaises(KeyError, ts2.ix.__setitem__, slice(d1, d2), 0)

    def test_ix_getitem_setitem_integer_slice_keyerrors(self):
        s = Series(np.random.randn(10), index=range(0, 20, 2))

        # this is OK
        cp = s.copy()
        cp.ix[4:10] = 0
        self.assert_((cp.ix[4:10] == 0).all())

        # so is this
        cp = s.copy()
        cp.ix[3:11] = 0
        self.assert_((cp.ix[3:11] == 0).values.all())

        result = s.ix[4:10]
        result2 = s.ix[3:11]
        expected = s.reindex([4, 6, 8, 10])

        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        # non-monotonic, raise KeyError
        s2 = s[::-1]
        self.assertRaises(KeyError, s2.ix.__getitem__, slice(3, 11))
        self.assertRaises(KeyError, s2.ix.__setitem__, slice(3, 11), 0)

    def test_ix_getitem_iterator(self):
        idx = iter(self.series.index[:10])
        result = self.series.ix[idx]
        assert_series_equal(result, self.series[:10])

    def test_ix_setitem(self):
        inds = self.series.index[[3,4,7]]

        result = self.series.copy()
        result.ix[inds] = 5

        expected = self.series.copy()
        expected[[3,4,7]] = 5
        assert_series_equal(result, expected)

        result.ix[5:10] = 10
        expected[5:10] = 10
        assert_series_equal(result, expected)

        # set slice with indices
        d1, d2 = self.series.index[[5, 15]]
        result.ix[d1:d2] = 6
        expected[5:16] = 6 # because it's inclusive
        assert_series_equal(result, expected)

        # set index value
        self.series.ix[d1] = 4
        self.series.ix[d2] = 6
        self.assertEquals(self.series[d1], 4)
        self.assertEquals(self.series[d2], 6)

    def test_ix_setitem_boolean(self):
        mask = self.series > self.series.median()

        result = self.series.copy()
        result.ix[mask] = 0
        expected = self.series
        expected[mask] = 0
        assert_series_equal(result, expected)

    def test_ix_setitem_corner(self):
        inds = list(self.series.index[[5, 8, 12]])
        self.series.ix[inds] = 5
        self.assertRaises(Exception, self.series.ix.__setitem__,
                          inds + ['foo'], 5)

    def test_get_set_boolean_different_order(self):
        ordered = self.series.order()

        # setting
        copy = self.series.copy()
        copy[ordered > 0] = 0

        expected = self.series.copy()
        expected[expected > 0] = 0

        assert_series_equal(copy, expected)

        # getting
        sel = self.series[ordered > 0]
        exp = self.series[self.series > 0]
        assert_series_equal(sel, exp)

    def test_repr(self):
        str(self.ts)
        str(self.series)
        str(self.series.astype(int))
        str(self.objSeries)

        str(Series(tm.randn(1000), index=np.arange(1000)))
        str(Series(tm.randn(1000), index=np.arange(1000, 0, step=-1)))

        # empty
        str(self.empty)

        # with NaNs
        self.series[5:7] = np.NaN
        str(self.series)

        # with Nones
        ots = self.ts.astype('O')
        ots[::2] = None
        repr(ots)

        # tuple name, e.g. from hierarchical index
        self.series.name = ('foo', 'bar', 'baz')
        repr(self.series)

        biggie = Series(tm.randn(1000), index=np.arange(1000),
                        name=('foo', 'bar', 'baz'))
        repr(biggie)

        # 0 as name
        ser = Series(np.random.randn(100), name=0)
        rep_str = repr(ser)
        self.assert_("Name: 0" in rep_str)

        # tidy repr
        ser = Series(np.random.randn(1001), name=0)
        rep_str = repr(ser)
        self.assert_("Name: 0" in rep_str)

    def test_repr_bool_fails(self):
        s = Series([DataFrame(np.random.randn(2,2)) for i in range(5)])

        import sys

        buf = StringIO()
        sys.stderr = buf
        # it works (with no Cython exception barf)!
        repr(s)
        sys.stderr = sys.__stderr__
        self.assertEquals(buf.getvalue(), '')

    def test_repr_name_iterable_indexable(self):
        s = Series([1, 2, 3], name=np.int64(3))

        # it works!
        repr(s)

        s.name = (u"\u05d0",) * 2
        repr(s)

    def test_timeseries_repr_object_dtype(self):
        index = Index([datetime(2000, 1, 1) + timedelta(i)
                       for i in range(1000)], dtype=object)
        ts = Series(np.random.randn(len(index)), index)
        repr(ts)

        ts = tm.makeTimeSeries(1000)
        self.assert_(repr(ts).splitlines()[-1].startswith('Freq:'))

        ts2 = ts.ix[np.random.randint(0, len(ts)-1, 400)]
        repr(ts).splitlines()[-1]

    def test_iter(self):
        for i, val in enumerate(self.series):
            self.assertEqual(val, self.series[i])

        for i, val in enumerate(self.ts):
            self.assertEqual(val, self.ts[i])

    def test_keys(self):
        # HACK: By doing this in two stages, we avoid 2to3 wrapping the call
        # to .keys() in a list()
        getkeys = self.ts.keys
        self.assert_(getkeys() is self.ts.index)

    def test_values(self):
        self.assert_(np.array_equal(self.ts, self.ts.values))

    def test_iteritems(self):
        for idx, val in self.series.iteritems():
            self.assertEqual(val, self.series[idx])

        for idx, val in self.ts.iteritems():
            self.assertEqual(val, self.ts[idx])

    def test_sum(self):
        self._check_stat_op('sum', np.sum)

    def test_sum_inf(self):
        s = Series(np.random.randn(10))
        s2 = s.copy()
        s[5:8] = np.inf
        s2[5:8] = np.nan
        assert_almost_equal(s.sum(), s2.sum())

        import pandas.core.nanops as nanops
        arr = np.random.randn(100, 100).astype('f4')
        arr[:, 2] = np.inf
        res = nanops.nansum(arr, axis=1)
        expected = nanops._nansum(arr, axis=1)
        assert_almost_equal(res, expected)

    def test_mean(self):
        self._check_stat_op('mean', np.mean)

    def test_median(self):
        self._check_stat_op('median', np.median)

        # test with integers, test failure
        int_ts = TimeSeries(np.ones(10, dtype=int), index=range(10))
        self.assertAlmostEqual(np.median(int_ts), int_ts.median())

    def test_prod(self):
        self._check_stat_op('prod', np.prod)

    def test_min(self):
        self._check_stat_op('min', np.min, check_objects=True)

    def test_max(self):
        self._check_stat_op('max', np.max, check_objects=True)

    def test_var_std(self):
        alt = lambda x: np.std(x, ddof=1)
        self._check_stat_op('std', alt)

        alt = lambda x: np.var(x, ddof=1)
        self._check_stat_op('var', alt)

        result = self.ts.std(ddof=4)
        expected = np.std(self.ts.values, ddof=4)
        assert_almost_equal(result, expected)

        result = self.ts.var(ddof=4)
        expected = np.var(self.ts.values, ddof=4)
        assert_almost_equal(result, expected)

    def test_skew(self):
        _skip_if_no_scipy()

        from scipy.stats import skew
        alt =lambda x: skew(x, bias=False)
        self._check_stat_op('skew', alt)

    def test_kurt(self):
        _skip_if_no_scipy()

        from scipy.stats import kurtosis
        alt = lambda x: kurtosis(x, bias=False)
        self._check_stat_op('kurt', alt)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        self.assertAlmostEqual(s.kurt(), s.kurt(level=0)['bar'])

    def test_argsort(self):
        self._check_accum_op('argsort')
        argsorted = self.ts.argsort()
        self.assert_(issubclass(argsorted.dtype.type, np.integer))

    def test_argsort_stable(self):
        s = Series(np.random.randint(0, 100, size=10000))
        mindexer = s.argsort(kind='mergesort')
        qindexer = s.argsort()

        mexpected = np.argsort(s.values, kind='mergesort')
        qexpected = np.argsort(s.values, kind='quicksort')

        self.assert_(np.array_equal(mindexer, mexpected))
        self.assert_(np.array_equal(qindexer, qexpected))
        self.assert_(not np.array_equal(qindexer, mindexer))

    def test_cumsum(self):
        self._check_accum_op('cumsum')

    def test_cumprod(self):
        self._check_accum_op('cumprod')

    def test_cummin(self):
        self.assert_(np.array_equal(self.ts.cummin(),
                                    np.minimum.accumulate(np.array(self.ts))))
        ts = self.ts.copy()
        ts[::2]  = np.NaN
        result   = ts.cummin()[1::2]
        expected = np.minimum.accumulate(ts.valid())

        self.assert_(np.array_equal(result, expected))

    def test_cummax(self):
        self.assert_(np.array_equal(self.ts.cummax(),
                                    np.maximum.accumulate(np.array(self.ts))))
        ts = self.ts.copy()
        ts[::2]  = np.NaN
        result   = ts.cummax()[1::2]
        expected = np.maximum.accumulate(ts.valid())

        self.assert_(np.array_equal(result, expected))

    def test_npdiff(self):
        s = Series(np.arange(5))
        r = np.diff(s)
        assert_series_equal(Series([nan, 0, 0, 0, nan]), r)

    def _check_stat_op(self, name, alternate, check_objects=False):
        import pandas.core.nanops as nanops

        def testit():
            f = getattr(Series, name)

            # add some NaNs
            self.series[5:15] = np.NaN

            # skipna or no
            self.assert_(notnull(f(self.series)))
            self.assert_(isnull(f(self.series, skipna=False)))

            # check the result is correct
            nona = self.series.dropna()
            assert_almost_equal(f(nona), alternate(nona))

            allna = self.series * nan
            self.assert_(np.isnan(f(allna)))

            # dtype=object with None, it works!
            s = Series([1, 2, 3, None, 5])
            f(s)

            # check date range
            if check_objects:
                s = Series(bdate_range('1/1/2000', periods=10))
                res = f(s)
                exp = alternate(s)
                self.assertEqual(res, exp)

        testit()

        try:
            import bottleneck as bn
            nanops._USE_BOTTLENECK = False
            testit()
            nanops._USE_BOTTLENECK = True
        except ImportError:
            pass

    def _check_accum_op(self, name):
        func = getattr(np, name)
        self.assert_(np.array_equal(func(self.ts), func(np.array(self.ts))))

        # with missing values
        ts = self.ts.copy()
        ts[::2] = np.NaN

        result = func(ts)[1::2]
        expected = func(np.array(ts.valid()))

        self.assert_(np.array_equal(result, expected))

    def test_round(self):
        # numpy.round doesn't preserve metadata, probably a numpy bug,
        # re: GH #314
        result = np.round(self.ts, 2)
        expected = Series(np.round(self.ts.values, 2), index=self.ts.index)
        assert_series_equal(result, expected)
        self.assertEqual(result.name, self.ts.name)

    def test_prod_numpy16_bug(self):
        s = Series([1., 1., 1.] , index=range(3))
        result = s.prod()
        self.assert_(not isinstance(result, Series))

    def test_quantile(self):
        from pandas.compat.scipy import scoreatpercentile

        q = self.ts.quantile(0.1)
        self.assertEqual(q, scoreatpercentile(self.ts.valid(), 10))

        q = self.ts.quantile(0.9)
        self.assertEqual(q, scoreatpercentile(self.ts.valid(), 90))

    def test_describe(self):
        _ = self.series.describe()
        _ = self.ts.describe()

    def test_describe_percentiles(self):
        desc = self.series.describe(percentile_width=50)
        assert '75%' in desc.index
        assert '25%' in desc.index

        desc = self.series.describe(percentile_width=95)
        assert '97.5%' in desc.index
        assert '2.5%' in desc.index

    def test_describe_objects(self):
        s = Series(['a', 'b', 'b', np.nan, np.nan, np.nan, 'c', 'd', 'a', 'a'])
        result = s.describe()
        expected = Series({'count' : 7, 'unique' : 4,
                           'top' : 'a', 'freq' : 3}, index=result.index)
        assert_series_equal(result, expected)

        dt = list(self.ts.index)
        dt.append(dt[0])
        ser = Series(dt)
        rs = ser.describe()
        min_date = min(dt)
        max_date = max(dt)
        xp = Series({'count' : len(dt),
                     'unique' : len(self.ts.index),
                     'first' : min_date, 'last' : max_date, 'freq' : 2,
                     'top' : min_date}, index=rs.index)
        assert_series_equal(rs, xp)

    def test_describe_empty(self):
        result = self.empty.describe()

        self.assert_(result['count'] == 0)
        self.assert_(result.drop('count').isnull().all())

        nanSeries = Series([np.nan])
        nanSeries.name = 'NaN'
        result = nanSeries.describe()
        self.assert_(result['count'] == 0)
        self.assert_(result.drop('count').isnull().all())

    def test_describe_none(self):
        noneSeries = Series([None])
        noneSeries.name = 'None'
        assert_series_equal(noneSeries.describe(),
                            Series([0, 0], index=['count', 'unique']))

    def test_append(self):
        appendedSeries = self.series.append(self.objSeries)
        for idx, value in appendedSeries.iteritems():
            if idx in self.series.index:
                self.assertEqual(value, self.series[idx])
            elif idx in self.objSeries.index:
                self.assertEqual(value, self.objSeries[idx])
            else:
                self.fail("orphaned index!")

        self.assertRaises(Exception, self.ts.append, self.ts,
                          verify_integrity=True)

    def test_append_many(self):
        pieces = [self.ts[:5], self.ts[5:10], self.ts[10:]]

        result = pieces[0].append(pieces[1:])
        assert_series_equal(result, self.ts)

    def test_all_any(self):
        np.random.seed(12345)
        ts = tm.makeTimeSeries()
        bool_series = ts > 0
        self.assert_(not bool_series.all())
        self.assert_(bool_series.any())

    def test_operators(self):

        def _check_op(series, other, op, pos_only=False):
            left = np.abs(series) if pos_only else series
            right = np.abs(other) if pos_only else other

            cython_or_numpy = op(left, right)
            python = left.combine(right, op)
            tm.assert_almost_equal(cython_or_numpy, python)

        def check(series, other):
            simple_ops = ['add', 'sub', 'mul', 'truediv', 'floordiv']

            for opname in simple_ops:
                _check_op(series, other, getattr(operator, opname))

            _check_op(series, other, operator.pow, pos_only=True)

            _check_op(series, other, lambda x, y: operator.add(y, x))
            _check_op(series, other, lambda x, y: operator.sub(y, x))
            _check_op(series, other, lambda x, y: operator.truediv(y, x))
            _check_op(series, other, lambda x, y: operator.floordiv(y, x))
            _check_op(series, other, lambda x, y: operator.mul(y, x))
            _check_op(series, other, lambda x, y: operator.pow(y, x),
                      pos_only=True)

        check(self.ts, self.ts * 2)
        check(self.ts, self.ts * 0)
        check(self.ts, self.ts[::2])
        check(self.ts, 5)

        def check_comparators(series, other):
            _check_op(series, other, operator.gt)
            _check_op(series, other, operator.ge)
            _check_op(series, other, operator.eq)
            _check_op(series, other, operator.lt)
            _check_op(series, other, operator.le)

        check_comparators(self.ts, 5)
        check_comparators(self.ts, self.ts + 1)

    def test_operators_empty_int_corner(self):
        s1 = Series([], [], dtype=np.int32)
        s2 = Series({'x' : 0.})

        # it works!
        _ = s1 * s2

    # NumPy limitiation =(

    # def test_logical_range_select(self):
    #     np.random.seed(12345)
    #     selector = -0.5 <= self.ts <= 0.5
    #     expected = (self.ts >= -0.5) & (self.ts <= 0.5)
    #     assert_series_equal(selector, expected)

    def test_operators_na_handling(self):
        from decimal import Decimal
        from datetime import date
        s = Series([Decimal('1.3'), Decimal('2.3')],
                   index=[date(2012,1,1), date(2012,1,2)])

        result = s + s.shift(1)
        self.assert_(isnull(result[0]))

        s = Series(['foo', 'bar', 'baz', np.nan])
        result = 'prefix_' + s
        expected = Series(['prefix_foo', 'prefix_bar', 'prefix_baz', np.nan])
        assert_series_equal(result, expected)

        result = s + '_suffix'
        expected = Series(['foo_suffix', 'bar_suffix', 'baz_suffix', np.nan])
        assert_series_equal(result, expected)

    def test_object_comparisons(self):
        s = Series(['a', 'b', np.nan, 'c', 'a'])

        result = s == 'a'
        expected = Series([True, False, False, False, True])
        assert_series_equal(result, expected)

        result = s < 'a'
        expected = Series([False, False, False, False, False])
        assert_series_equal(result, expected)

        result = s != 'a'
        expected = -(s == 'a')
        assert_series_equal(result, expected)

    def test_comparison_operators_with_nas(self):
        s = Series(bdate_range('1/1/2000', periods=10), dtype=object)
        s[::2] = np.nan

        # test that comparions work
        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            val = s[5]

            f = getattr(operator, op)
            result = f(s, val)

            expected = f(s.dropna(), val).reindex(s.index)

            if op == 'ne':
                expected = expected.fillna(True).astype(bool)
            else:
                expected = expected.fillna(False).astype(bool)

            assert_series_equal(result, expected)

            # fffffffuuuuuuuuuuuu
            # result = f(val, s)
            # expected = f(val, s.dropna()).reindex(s.index)
            # assert_series_equal(result, expected)

        # boolean &, |, ^ should work with object arrays and propagate NAs

        ops = ['and_', 'or_', 'xor']
        mask = s.isnull()
        for bool_op in ops:
            f = getattr(operator, bool_op)

            filled = s.fillna(s[0])

            result = f(s < s[9], s > s[3])

            expected = f(filled < filled[9], filled > filled[3])
            expected[mask] = False
            assert_series_equal(result, expected)

    def test_comparison_object_numeric_nas(self):
        s = Series(np.random.randn(10), dtype=object)
        shifted = s.shift(2)

        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            f = getattr(operator, op)

            result = f(s, shifted)
            expected = f(s.astype(float), shifted.astype(float))
            assert_series_equal(result, expected)

    def test_more_na_comparisons(self):
        left = Series(['a', np.nan, 'c'])
        right = Series(['a', np.nan, 'd'])

        result = left == right
        expected = Series([True, False, False])
        assert_series_equal(result, expected)

        result = left != right
        expected = Series([False, True, True])
        assert_series_equal(result, expected)

        result = left == np.nan
        expected = Series([False, False, False])
        assert_series_equal(result, expected)

        result = left != np.nan
        expected = Series([True, True, True])
        assert_series_equal(result, expected)

    def test_comparison_different_length(self):
        a = Series(['a', 'b', 'c'])
        b = Series(['b', 'a'])
        self.assertRaises(ValueError, a.__lt__, b)

    def test_between(self):
        s = Series(bdate_range('1/1/2000', periods=20).asobject)
        s[::2] = np.nan

        result = s[s.between(s[3], s[17])]
        expected = s[3:18].dropna()
        assert_series_equal(result, expected)

        result = s[s.between(s[3], s[17], inclusive=False)]
        expected = s[5:16].dropna()
        assert_series_equal(result, expected)

    def test_setitem_na_exception(self):
        def testme1():
            s = Series([2,3,4,5,6,7,8,9,10])
            s[::2] = np.nan

        def testme2():
            s = Series([True, True, False, False])
            s[::2] = np.nan

        def testme3():
            s = Series(np.arange(10))
            s[:5] = np.nan

        self.assertRaises(Exception, testme1)
        self.assertRaises(Exception, testme2)
        self.assertRaises(Exception, testme3)

    def test_scalar_na_cmp_corners(self):
        s = Series([2,3,4,5,6,7,8,9,10])

        def tester(a, b):
            return a & b

        self.assertRaises(ValueError, tester, s, datetime(2005,1,1))

        s = Series([2,3,4,5,6,7,8,9,datetime(2005,1,1)])
        s[::2] = np.nan

        assert_series_equal(tester(s, list(s)), s)

        d = DataFrame({'A':s})
        self.assertRaises(TypeError, tester, s, d)

    def test_idxmin(self):
        # test idxmin
        # _check_stat_op approach can not be used here because of isnull check.

        # add some NaNs
        self.series[5:15] = np.NaN

        # skipna or no
        self.assertEqual(self.series[self.series.idxmin()], self.series.min())
        self.assert_(isnull(self.series.idxmin(skipna=False)))

        # no NaNs
        nona = self.series.dropna()
        self.assertEqual(nona[nona.idxmin()], nona.min())
        self.assertEqual(nona.index.values.tolist().index(nona.idxmin()),
                         nona.values.argmin())

        # all NaNs
        allna = self.series * nan
        self.assert_(isnull(allna.idxmin()))

    def test_idxmax(self):
        # test idxmax
        # _check_stat_op approach can not be used here because of isnull check.

        # add some NaNs
        self.series[5:15] = np.NaN

        # skipna or no
        self.assertEqual(self.series[self.series.idxmax()], self.series.max())
        self.assert_(isnull(self.series.idxmax(skipna=False)))

        # no NaNs
        nona = self.series.dropna()
        self.assertEqual(nona[nona.idxmax()], nona.max())
        self.assertEqual(nona.index.values.tolist().index(nona.idxmax()),
                         nona.values.argmax())

        # all NaNs
        allna = self.series * nan
        self.assert_(isnull(allna.idxmax()))

    def test_operators_date(self):
        result = self.objSeries + timedelta(1)
        result = self.objSeries - timedelta(1)

    def test_operators_corner(self):
        series = self.ts

        empty = Series([], index=Index([]))

        result = series + empty
        self.assert_(np.isnan(result).all())

        result = empty + Series([], index=Index([]))
        self.assert_(len(result) == 0)

        # TODO: this returned NotImplemented earlier, what to do?
        # deltas = Series([timedelta(1)] * 5, index=np.arange(5))
        # sub_deltas = deltas[::2]
        # deltas5 = deltas * 5
        # deltas = deltas + sub_deltas

        # float + int
        int_ts = self.ts.astype(int)[:-5]
        added = self.ts + int_ts
        expected = self.ts.values[:-5] + int_ts.values
        self.assert_(np.array_equal(added[:-5], expected))

    def test_operators_reverse_object(self):
        # GH 56
        arr = Series(np.random.randn(10), index=np.arange(10),
                     dtype=object)

        def _check_op(arr, op):
            result = op(1., arr)
            expected = op(1., arr.astype(float))
            assert_series_equal(result.astype(float), expected)

        _check_op(arr, operator.add)
        _check_op(arr, operator.sub)
        _check_op(arr, operator.mul)
        _check_op(arr, operator.truediv)
        _check_op(arr, operator.floordiv)

    def test_series_frame_radd_bug(self):
        from pandas.util.testing import rands
        import operator

        # GH 353
        vals = Series([rands(5) for _ in xrange(10)])
        result = 'foo_' + vals
        expected = vals.map(lambda x: 'foo_' + x)
        assert_series_equal(result, expected)

        frame = DataFrame({'vals' : vals})
        result = 'foo_' + frame
        expected = DataFrame({'vals' : vals.map(lambda x: 'foo_' + x)})
        tm.assert_frame_equal(result, expected)

        # really raise this time
        self.assertRaises(TypeError, operator.add, datetime.now(), self.ts)

    def test_operators_frame(self):
        # rpow does not work with DataFrame
        df = DataFrame({'A' : self.ts})

        tm.assert_almost_equal(self.ts + self.ts, (self.ts + df)['A'])
        tm.assert_almost_equal(self.ts ** self.ts, (self.ts ** df)['A'])
        tm.assert_almost_equal(self.ts < self.ts, (self.ts < df)['A'])

    def test_operators_combine(self):
        def _check_fill(meth, op, a, b, fill_value=0):
            exp_index = a.index.union(b.index)
            a = a.reindex(exp_index)
            b = b.reindex(exp_index)

            amask = isnull(a)
            bmask = isnull(b)

            exp_values = []
            for i in range(len(exp_index)):
                if amask[i]:
                    if bmask[i]:
                        exp_values.append(nan)
                        continue
                    exp_values.append(op(fill_value, b[i]))
                elif bmask[i]:
                    if amask[i]:
                        exp_values.append(nan)
                        continue
                    exp_values.append(op(a[i], fill_value))
                else:
                    exp_values.append(op(a[i], b[i]))

            result = meth(a, b, fill_value=fill_value)
            expected = Series(exp_values, exp_index)
            assert_series_equal(result, expected)

        a = Series([nan, 1., 2., 3., nan], index=np.arange(5))
        b = Series([nan, 1, nan, 3, nan, 4.], index=np.arange(6))

        ops = [Series.add, Series.sub, Series.mul, Series.div]
        equivs = [operator.add, operator.sub, operator.mul]
        if py3compat.PY3:
            equivs.append(operator.truediv)
        else:
            equivs.append(operator.div)
        fillvals = [0, 0, 1, 1]

        for op, equiv_op, fv in zip(ops, equivs, fillvals):
            result = op(a, b)
            exp = equiv_op(a, b)
            assert_series_equal(result, exp)
            _check_fill(op, equiv_op, a, b, fill_value=fv)

    def test_combine_first(self):
        values = tm.makeIntIndex(20).values.astype(float)
        series = Series(values, index=tm.makeIntIndex(20))

        series_copy = series * 2
        series_copy[::2] = np.NaN

        # nothing used from the input
        combined = series.combine_first(series_copy)

        self.assert_(np.array_equal(combined, series))

        # Holes filled from input
        combined = series_copy.combine_first(series)
        self.assert_(np.isfinite(combined).all())

        self.assert_(np.array_equal(combined[::2], series[::2]))
        self.assert_(np.array_equal(combined[1::2], series_copy[1::2]))

        # mixed types
        index = tm.makeStringIndex(20)
        floats = Series(tm.randn(20), index=index)
        strings = Series(tm.makeStringIndex(10), index=index[::2])

        combined = strings.combine_first(floats)

        tm.assert_dict_equal(strings, combined, compare_keys=False)
        tm.assert_dict_equal(floats[1::2], combined, compare_keys=False)

        # corner case
        s = Series([1., 2, 3], index=[0, 1, 2])
        result = s.combine_first(Series([], index=[]))
        assert_series_equal(s, result)

    def test_update(self):
        s = Series([1.5, nan, 3., 4., nan])
        s2 = Series([nan, 3.5, nan, 5.])
        s.update(s2)

        expected = Series([1.5, 3.5, 3., 5., np.nan])
        assert_series_equal(s, expected)

    def test_corr(self):
        _skip_if_no_scipy()

        import scipy.stats as stats

        # full overlap
        self.assertAlmostEqual(self.ts.corr(self.ts), 1)

        # partial overlap
        self.assertAlmostEqual(self.ts[:15].corr(self.ts[5:]), 1)

        # No overlap
        self.assert_(np.isnan(self.ts[::2].corr(self.ts[1::2])))

        # all NA
        cp = self.ts[:10].copy()
        cp[:] = np.nan
        self.assert_(isnull(cp.corr(cp)))

        A = tm.makeTimeSeries()
        B = tm.makeTimeSeries()
        result = A.corr(B)
        expected, _ = stats.pearsonr(A, B)
        self.assertAlmostEqual(result, expected)

    def test_corr_rank(self):
        _skip_if_no_scipy()

        import scipy
        import scipy.stats as stats

        # kendall and spearman
        A = tm.makeTimeSeries()
        B = tm.makeTimeSeries()
        A[-5:] = A[:5]
        result = A.corr(B, method='kendall')
        expected = stats.kendalltau(A, B)[0]
        self.assertAlmostEqual(result, expected)

        result = A.corr(B, method='spearman')
        expected = stats.spearmanr(A, B)[0]
        self.assertAlmostEqual(result, expected)

        # these methods got rewritten in 0.8
        if int(scipy.__version__.split('.')[1]) < 9:
            raise nose.SkipTest

        # results from R
        A = Series([-0.89926396,  0.94209606, -1.03289164, -0.95445587,
                    0.76910310, -0.06430576, -2.09704447, 0.40660407,
                    -0.89926396,  0.94209606])
        B = Series([-1.01270225, -0.62210117, -1.56895827,  0.59592943,
                    -0.01680292,  1.17258718, -1.06009347, -0.10222060,
                    -0.89076239,  0.89372375])
        kexp = 0.4319297
        sexp = 0.5853767
        self.assertAlmostEqual(A.corr(B, method='kendall'), kexp)
        self.assertAlmostEqual(A.corr(B, method='spearman'), sexp)

    def test_cov(self):
        # full overlap
        self.assertAlmostEqual(self.ts.cov(self.ts), self.ts.std()**2)

        # partial overlap
        self.assertAlmostEqual(self.ts[:15].cov(self.ts[5:]), self.ts[5:15].std()**2)

        # No overlap
        self.assert_(np.isnan(self.ts[::2].cov(self.ts[1::2])))

        # all NA
        cp = self.ts[:10].copy()
        cp[:] = np.nan
        self.assert_(isnull(cp.cov(cp)))

    def test_copy(self):
        ts = self.ts.copy()

        ts[::2] = np.NaN

        # Did not modify original Series
        self.assertFalse(np.isnan(self.ts[0]))

    def test_count(self):
        self.assertEqual(self.ts.count(), len(self.ts))

        self.ts[::2] = np.NaN

        self.assertEqual(self.ts.count(), np.isfinite(self.ts).sum())

    def test_value_counts_nunique(self):
        s = Series(['a', 'b', 'b', 'b', 'b', 'a', 'c', 'd', 'd', 'a'])
        hist = s.value_counts()
        expected = Series([4, 3, 2, 1], index=['b', 'a', 'd', 'c'])
        assert_series_equal(hist, expected)

        self.assertEquals(s.nunique(), 4)

        # handle NA's properly
        s[5:7] = np.nan
        hist = s.value_counts()
        expected = s.dropna().value_counts()
        assert_series_equal(hist, expected)

        s = Series({})
        hist = s.value_counts()
        expected = Series([])
        assert_series_equal(hist, expected)

    def test_unique(self):
        # 714 also, dtype=float
        s = Series([1.2345] * 100)
        s[::2] = np.nan
        result = s.unique()
        self.assert_(len(result) == 2)

        s = Series([1.2345] * 100, dtype='f4')
        s[::2] = np.nan
        result = s.unique()
        self.assert_(len(result) == 2)

        # NAs in object arrays #714
        s = Series(['foo'] * 100, dtype='O')
        s[::2] = np.nan
        result = s.unique()
        self.assert_(len(result) == 2)

        # integers
        s = Series(np.random.randint(0, 100, size=100))
        result = np.sort(s.unique())
        expected = np.unique(s.values)
        self.assert_(np.array_equal(result, expected))

        s = Series(np.random.randint(0, 100, size=100).astype(np.int32))
        result = np.sort(s.unique())
        expected = np.unique(s.values)
        self.assert_(np.array_equal(result, expected))

        # test string arrays for coverage
        strings = np.tile(np.array([tm.rands(10) for _ in xrange(10)]), 10)
        result = np.sort(nanops.unique1d(strings))
        expected = np.unique(strings)
        self.assert_(np.array_equal(result, expected))

        # decision about None

        s = Series([1, 2, 3, None, None, None], dtype=object)
        result = s.unique()
        expected = np.array([1, 2, 3, None], dtype=object)
        self.assert_(np.array_equal(result, expected))

    def test_sort(self):
        ts = self.ts.copy()
        ts.sort()

        self.assert_(np.array_equal(ts, self.ts.order()))
        self.assert_(np.array_equal(ts.index, self.ts.order().index))

    def test_sort_index(self):
        import random

        rindex = list(self.ts.index)
        random.shuffle(rindex)

        random_order = self.ts.reindex(rindex)
        sorted_series = random_order.sort_index()
        assert_series_equal(sorted_series, self.ts)


        # descending
        sorted_series = random_order.sort_index(ascending=False)
        assert_series_equal(sorted_series,
                            self.ts.reindex(self.ts.index[::-1]))

    def test_order(self):
        ts = self.ts.copy()
        ts[:5] = np.NaN
        vals = ts.values

        result = ts.order()
        self.assert_(np.isnan(result[-5:]).all())
        self.assert_(np.array_equal(result[:-5], np.sort(vals[5:])))

        result = ts.order(na_last=False)
        self.assert_(np.isnan(result[:5]).all())
        self.assert_(np.array_equal(result[5:], np.sort(vals[5:])))

        # something object-type
        ser = Series(['A', 'B'], [1, 2])
        # no failure
        ser.order()

        # ascending=False
        ordered = ts.order(ascending=False)
        expected = np.sort(ts.valid().values)[::-1]
        assert_almost_equal(expected, ordered.valid().values)
        ordered = ts.order(ascending=False, na_last=False)
        assert_almost_equal(expected, ordered.valid().values)

    def test_rank(self):
        from pandas.compat.scipy import rankdata

        self.ts[::2] = np.nan
        self.ts[:10][::3] = 4.

        ranks = self.ts.rank()
        oranks = self.ts.astype('O').rank()

        assert_series_equal(ranks, oranks)

        mask =  np.isnan(self.ts)
        filled = self.ts.fillna(np.inf)

        exp = rankdata(filled)
        exp[mask] = np.nan

        assert_almost_equal(ranks, exp)

        iseries = Series(np.arange(5).repeat(2))

        iranks = iseries.rank()
        exp = iseries.astype(float).rank()
        assert_series_equal(iranks, exp)

    def test_from_csv(self):
        self.ts.to_csv('_foo')
        ts = Series.from_csv('_foo')
        assert_series_equal(self.ts, ts)

        self.series.to_csv('_foo')
        series = Series.from_csv('_foo')
        self.assert_(series.name is None)
        self.assert_(series.index.name is None)
        assert_series_equal(self.series, series)

        outfile = open('_foo', 'w')
        outfile.write('1998-01-01|1.0\n1999-01-01|2.0')
        outfile.close()
        series = Series.from_csv('_foo',sep='|')
        checkseries = Series({datetime(1998,1,1): 1.0, datetime(1999,1,1): 2.0})
        assert_series_equal(checkseries, series)

        series = Series.from_csv('_foo',sep='|',parse_dates=False)
        checkseries = Series({'1998-01-01': 1.0, '1999-01-01': 2.0})
        assert_series_equal(checkseries, series)

        os.remove('_foo')

    def test_to_csv(self):
        self.ts.to_csv('_foo')

        lines = open('_foo', 'U').readlines()
        assert(lines[1] != '\n')

        self.ts.to_csv('_foo', index=False)
        arr = np.loadtxt('_foo')
        assert_almost_equal(arr, self.ts.values)

        os.remove('_foo')

    def test_to_csv_stringio(self):
        buf = StringIO()
        self.ts.to_csv(buf, index=False)
        buf.seek(0)
        arr = np.loadtxt(buf)
        assert_almost_equal(arr, self.ts.values)

    def test_to_dict(self):
        self.assert_(np.array_equal(Series(self.ts.to_dict()), self.ts))

    def test_to_csv_float_format(self):
        filename = '__tmp__.csv'
        ser = Series([0.123456, 0.234567, 0.567567])
        ser.to_csv(filename, float_format='%.2f')

        rs = Series.from_csv(filename)
        xp = Series([0.12, 0.23, 0.57])
        assert_series_equal(rs, xp)
        os.remove(filename)

    def test_to_csv_list_entries(self):
        s = Series(['jack and jill','jesse and frank'])

        split = s.str.split(r'\s+and\s+')

        buf = StringIO()
        split.to_csv(buf)

    def test_clip(self):
        val = self.ts.median()

        self.assertEqual(self.ts.clip_lower(val).min(), val)
        self.assertEqual(self.ts.clip_upper(val).max(), val)

        self.assertEqual(self.ts.clip(lower=val).min(), val)
        self.assertEqual(self.ts.clip(upper=val).max(), val)

        result = self.ts.clip(-0.5, 0.5)
        expected = np.clip(self.ts, -0.5, 0.5)
        assert_series_equal(result, expected)
        self.assert_(isinstance(expected, Series))

    def test_valid(self):
        ts = self.ts.copy()
        ts[::2] = np.NaN

        result = ts.valid()
        self.assertEqual(len(result), ts.count())

        tm.assert_dict_equal(result, ts, compare_keys=False)

    def test_isnull(self):
        ser = Series([0,5.4,3,nan,-0.001])
        assert_series_equal(ser.isnull(), Series([False,False,False,True,False]))
        ser = Series(["hi","",nan])
        assert_series_equal(ser.isnull(), Series([False,False,True]))

    def test_notnull(self):
        ser = Series([0,5.4,3,nan,-0.001])
        assert_series_equal(ser.notnull(), Series([True,True,True,False,True]))
        ser = Series(["hi","",nan])
        assert_series_equal(ser.notnull(), Series([True,True,False]))

    def test_shift(self):
        shifted = self.ts.shift(1)
        unshifted = shifted.shift(-1)

        tm.assert_dict_equal(unshifted.valid(), self.ts, compare_keys=False)

        offset = datetools.bday
        shifted = self.ts.shift(1, freq=offset)
        unshifted = shifted.shift(-1, freq=offset)

        assert_series_equal(unshifted, self.ts)

        unshifted = self.ts.shift(0, freq=offset)
        assert_series_equal(unshifted, self.ts)

        shifted = self.ts.shift(1, freq='B')
        unshifted = shifted.shift(-1, freq='B')

        assert_series_equal(unshifted, self.ts)

        # corner case
        unshifted = self.ts.shift(0)
        assert_series_equal(unshifted, self.ts)

        # Shifting with PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        tm.assert_dict_equal(unshifted.valid(), ps, compare_keys=False)

        shifted2 = ps.shift(1, 'B')
        shifted3 = ps.shift(1, datetools.bday)
        assert_series_equal(shifted2, shifted3)
        assert_series_equal(ps, shifted2.shift(-1, 'B'))

        self.assertRaises(ValueError, ps.shift, freq='D')

        # legacy support
        smod._SHOW_WARNINGS = False
        shifted4 = ps.shift(1, timeRule='B')
        assert_series_equal(shifted2, shifted4)

        shifted5 = ps.shift(1, offset=datetools.bday)
        assert_series_equal(shifted5, shifted4)
        smod._SHOW_WARNINGS = True

    def test_tshift(self):
        # PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_series_equal(unshifted, ps)

        shifted2 = ps.tshift(freq='B')
        assert_series_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=datetools.bday)
        assert_series_equal(shifted, shifted3)

        self.assertRaises(ValueError, ps.tshift, freq='M')

        # DatetimeIndex
        shifted = self.ts.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_series_equal(self.ts, unshifted)

        shifted2 = self.ts.tshift(freq=self.ts.index.freq)
        assert_series_equal(shifted, shifted2)

        inferred_ts = Series(self.ts.values, Index(np.asarray(self.ts.index)))
        shifted = inferred_ts.tshift(1)
        unshifted = shifted.tshift(-1)
        assert_series_equal(shifted, self.ts.tshift(1))
        assert_series_equal(unshifted, inferred_ts)

        no_freq = self.ts[[0, 5, 7]]
        self.assertRaises(ValueError, no_freq.tshift)

    def test_shift_int(self):
        ts = self.ts.astype(int)
        shifted = ts.shift(1)
        expected = ts.astype(float).shift(1)
        assert_series_equal(shifted, expected)

    def test_truncate(self):
        offset = datetools.bday

        ts = self.ts[::3]

        start, end = self.ts.index[3], self.ts.index[6]
        start_missing, end_missing = self.ts.index[2], self.ts.index[7]

        # neither specified
        truncated = ts.truncate()
        assert_series_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        assert_series_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        assert_series_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        assert_series_equal(truncated, expected)

        # corner case, empty series returned
        truncated = ts.truncate(after=self.ts.index[0] - offset)
        assert(len(truncated) == 0)

        truncated = ts.truncate(before=self.ts.index[-1] + offset)
        assert(len(truncated) == 0)

        self.assertRaises(Exception, ts.truncate,
                          before=self.ts.index[-1] + offset,
                          after=self.ts.index[0] - offset)

    def test_ptp(self):
        N = 1000
        arr = np.random.randn(N)
        ser = Series(arr)
        self.assertEqual(np.ptp(ser), np.ptp(arr))

    def test_asof(self):
        # array or list or dates
        N = 50
        rng = date_range('1/1/1990', periods=N, freq='53s')
        ts = Series(np.random.randn(N), index=rng)
        ts[15:30] = np.nan
        dates = date_range('1/1/1990', periods=N * 3, freq='25s')

        result = ts.asof(dates)
        self.assert_(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        result = ts.asof(list(dates))
        self.assert_(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        mask = (result.index >= lb) & (result.index < ub)
        rs = result[mask]
        self.assert_((rs == ts[lb]).all())

        val = result[result.index[result.index >= ub][0]]
        self.assertEqual(ts[ub], val)

        self.ts[5:10] = np.NaN
        self.ts[15:20] = np.NaN

        val1 = self.ts.asof(self.ts.index[7])
        val2 = self.ts.asof(self.ts.index[19])

        self.assertEqual(val1, self.ts[4])
        self.assertEqual(val2, self.ts[14])

        # accepts strings
        val1 = self.ts.asof(str(self.ts.index[7]))
        self.assertEqual(val1, self.ts[4])

        # in there
        self.assertEqual(self.ts.asof(self.ts.index[3]), self.ts[3])

        # no as of value
        d = self.ts.index[0] - datetools.bday
        self.assert_(np.isnan(self.ts.asof(d)))

    def test_asof_periodindex(self):
        from pandas import period_range, PeriodIndex
        # array or list or dates
        N = 50
        rng = period_range('1/1/1990', periods=N, freq='H')
        ts = Series(np.random.randn(N), index=rng)
        ts[15:30] = np.nan
        dates = date_range('1/1/1990', periods=N * 3, freq='37min')

        result = ts.asof(dates)
        self.assert_(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        result = ts.asof(list(dates))
        self.assert_(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        pix = PeriodIndex(result.index.values, freq='H')
        mask = (pix >= lb) & (pix < ub)
        rs = result[mask]
        self.assert_((rs == ts[lb]).all())

        ts[5:10] = np.NaN
        ts[15:20] = np.NaN

        val1 = ts.asof(ts.index[7])
        val2 = ts.asof(ts.index[19])

        self.assertEqual(val1, ts[4])
        self.assertEqual(val2, ts[14])

        # accepts strings
        val1 = ts.asof(str(ts.index[7]))
        self.assertEqual(val1, ts[4])

        # in there
        self.assertEqual(ts.asof(ts.index[3]), ts[3])

        # no as of value
        d = ts.index[0].to_timestamp() - datetools.bday
        self.assert_(np.isnan(ts.asof(d)))

    def test_asof_more(self):
        from pandas import date_range
        s = Series([nan, nan, 1, 2, nan, nan, 3, 4, 5],
                   index=date_range('1/1/2000', periods=9))

        dates = s.index[[4, 5, 6, 2, 1]]

        result = s.asof(dates)
        expected = Series([2, 2, 3, 1, np.nan], index=dates)

        assert_series_equal(result, expected)

        s = Series([1.5, 2.5, 1, 2, nan, nan, 3, 4, 5],
                   index=date_range('1/1/2000', periods=9))
        result = s.asof(s.index[0])
        self.assertEqual(result, s[0])

    def test_astype_cast_nan_int(self):
        df = Series([1.0, 2.0, 3.0, np.nan])
        self.assertRaises(ValueError, df.astype, np.int64)

    def test_astype_cast_object_int(self):
        arr = Series(["car", "house", "tree","1"])

        self.assertRaises(ValueError, arr.astype, int)
        self.assertRaises(ValueError, arr.astype, np.int64)
        self.assertRaises(ValueError, arr.astype, np.int8)

        arr = Series(['1', '2', '3', '4'], dtype=object)
        result = arr.astype(int)
        self.assert_(np.array_equal(result, np.arange(1, 5)))

    def test_map(self):
        index, data = tm.getMixedTypeDict()

        source = Series(data['B'], index=data['C'])
        target = Series(data['C'][:4], index=data['D'][:4])

        merged = target.map(source)

        for k, v in merged.iteritems():
            self.assertEqual(v, source[target[k]])

        # input could be a dict
        merged = target.map(source.to_dict())

        for k, v in merged.iteritems():
            self.assertEqual(v, source[target[k]])

        # function
        result = self.ts.map(lambda x: x * 2)
        self.assert_(np.array_equal(result, self.ts * 2))

    def test_map_int(self):
        left = Series({'a' : 1., 'b' : 2., 'c' : 3., 'd' : 4})
        right = Series({1 : 11, 2 : 22, 3 : 33})

        self.assert_(left.dtype == np.float_)
        self.assert_(issubclass(right.dtype.type, np.integer))

        merged = left.map(right)
        self.assert_(merged.dtype == np.float_)
        self.assert_(isnull(merged['d']))
        self.assert_(not isnull(merged['c']))

    def test_map_type_inference(self):
        s = Series(range(3))
        s2 = s.map(lambda x: np.where(x == 0, 0, 1))
        self.assert_(issubclass(s2.dtype.type, np.integer))

    def test_map_decimal(self):
        from decimal import Decimal

        result = self.series.map(lambda x: Decimal(str(x)))
        self.assert_(result.dtype == np.object_)
        self.assert_(isinstance(result[0], Decimal))

    def test_map_na_exclusion(self):
        s = Series([1.5, np.nan, 3, np.nan, 5])

        result = s.map(lambda x: x * 2, na_action='ignore')
        exp = s * 2
        assert_series_equal(result, exp)

    def test_apply(self):
        assert_series_equal(self.ts.apply(np.sqrt), np.sqrt(self.ts))

        # elementwise-apply
        import math
        assert_series_equal(self.ts.apply(math.exp), np.exp(self.ts))

        # does not return Series
        result = self.ts.apply(lambda x: x.values * 2)
        assert_series_equal(result, self.ts * 2)

    def test_apply_same_length_inference_bug(self):
        s = Series([1, 2])
        f = lambda x: (x, x + 1)

        result = s.apply(f)
        expected = s.map(f)
        assert_series_equal(result, expected)

        s = Series([1, 2, 3])
        result = s.apply(f)
        expected = s.map(f)
        assert_series_equal(result, expected)

    def test_apply_dont_convert_dtype(self):
        s = Series(np.random.randn(10))

        f = lambda x: x if x > 0 else np.nan
        result = s.apply(f, convert_dtype=False)
        self.assert_(result.dtype == object)

    def test_apply_args(self):
        s = Series(['foo,bar'])

        result = s.apply(str.split, args=(',',))
        self.assert_(result[0] == ['foo', 'bar'])

    def test_align(self):
        def _check_align(a, b, how='left', fill=None):
            aa, ab = a.align(b, join=how, fill_value=fill)

            join_index = a.index.join(b.index, how=how)
            if fill is not None:
                diff_a = aa.index.diff(join_index)
                diff_b = ab.index.diff(join_index)
                if len(diff_a) > 0:
                    self.assert_((aa.reindex(diff_a) == fill).all())
                if len(diff_b) > 0:
                    self.assert_((ab.reindex(diff_b) == fill).all())

            ea = a.reindex(join_index)
            eb = b.reindex(join_index)

            if fill is not None:
                ea = ea.fillna(fill)
                eb = eb.fillna(fill)

            assert_series_equal(aa, ea)
            assert_series_equal(ab, eb)

        for kind in JOIN_TYPES:
            _check_align(self.ts[2:], self.ts[:-5], how=kind)
            _check_align(self.ts[2:], self.ts[:-5], how=kind, fill=-1)

            # empty left
            _check_align(self.ts[:0], self.ts[:-5], how=kind)

            # empty right
            _check_align(self.ts[:-5], self.ts[:0], how=kind)

            # both empty
            _check_align(self.ts[:0], self.ts[:0], how=kind)

    def test_align_fill_method(self):
        def _check_align(a, b, how='left', method='pad', limit=None):
            aa, ab = a.align(b, join=how, method=method, limit=limit)

            join_index = a.index.join(b.index, how=how)
            ea = a.reindex(join_index)
            eb = b.reindex(join_index)

            ea = ea.fillna(method=method, limit=limit)
            eb = eb.fillna(method=method, limit=limit)

            assert_series_equal(aa, ea)
            assert_series_equal(ab, eb)

        for kind in JOIN_TYPES:
            for meth in ['pad', 'bfill']:
                _check_align(self.ts[2:], self.ts[:-5], how=kind, method=meth)
                _check_align(self.ts[2:], self.ts[:-5], how=kind,
                             method=meth, limit=1)

                # empty left
                _check_align(self.ts[:0], self.ts[:-5], how=kind, method=meth)
                _check_align(self.ts[:0], self.ts[:-5], how=kind, method=meth,
                             limit=1)

                # empty right
                _check_align(self.ts[:-5], self.ts[:0], how=kind, method=meth)
                _check_align(self.ts[:-5], self.ts[:0], how=kind, method=meth,
                             limit=1)

                # both empty
                _check_align(self.ts[:0], self.ts[:0], how=kind, method=meth)
                _check_align(self.ts[:0], self.ts[:0], how=kind, method=meth,
                             limit=1)

    def test_align_nocopy(self):
        b = self.ts[:5].copy()

        # do copy
        a = self.ts.copy()
        ra, _ = a.align(b, join='left')
        ra[:5] = 5
        self.assert_(not (a[:5] == 5).any())

        # do not copy
        a = self.ts.copy()
        ra, _ = a.align(b, join='left', copy=False)
        ra[:5] = 5
        self.assert_((a[:5] == 5).all())

        # do copy
        a = self.ts.copy()
        b = self.ts[:5].copy()
        _, rb = a.align(b, join='right')
        rb[:3] = 5
        self.assert_(not (b[:3] == 5).any())

        # do not copy
        a = self.ts.copy()
        b = self.ts[:5].copy()
        _, rb = a.align(b, join='right', copy=False)
        rb[:2] = 5
        self.assert_((b[:2] == 5).all())

    def test_align_sameindex(self):
        a, b = self.ts.align(self.ts, copy=False)
        self.assert_(a.index is self.ts.index)
        self.assert_(b.index is self.ts.index)

        # a, b = self.ts.align(self.ts, copy=True)
        # self.assert_(a.index is not self.ts.index)
        # self.assert_(b.index is not self.ts.index)

    def test_reindex(self):
        identity = self.series.reindex(self.series.index)
        self.assertEqual(id(self.series.index), id(identity.index))

        subIndex = self.series.index[10:20]
        subSeries = self.series.reindex(subIndex)

        for idx, val in subSeries.iteritems():
            self.assertEqual(val, self.series[idx])

        subIndex2 = self.ts.index[10:20]
        subTS = self.ts.reindex(subIndex2)

        for idx, val in subTS.iteritems():
            self.assertEqual(val, self.ts[idx])
        stuffSeries = self.ts.reindex(subIndex)

        self.assert_(np.isnan(stuffSeries).all())

        # This is extremely important for the Cython code to not screw up
        nonContigIndex = self.ts.index[::2]
        subNonContig = self.ts.reindex(nonContigIndex)
        for idx, val in subNonContig.iteritems():
            self.assertEqual(val, self.ts[idx])

        self.assertRaises(ValueError, self.ts.reindex)

    def test_reindex_corner(self):
        # (don't forget to fix this) I think it's fixed
        reindexed_dep = self.empty.reindex(self.ts.index, method='pad')

        # corner case: pad empty series
        reindexed = self.empty.reindex(self.ts.index, method='pad')

        # pass non-Index
        reindexed = self.ts.reindex(list(self.ts.index))
        assert_series_equal(self.ts, reindexed)

        # bad fill method
        ts = self.ts[::2]
        self.assertRaises(Exception, ts.reindex, self.ts.index, method='foo')

    def test_reindex_pad(self):
        s = Series(np.arange(10), np.arange(10))

        s2 = s[::2]

        reindexed = s2.reindex(s.index, method='pad')
        reindexed2 = s2.reindex(s.index, method='ffill')
        assert_series_equal(reindexed, reindexed2)

        # used platform int above, need to pass int explicitly here per #1219
        expected = Series([0, 0, 2, 2, 4, 4, 6, 6, 8, 8], dtype=int,
                          index=np.arange(10))
        assert_series_equal(reindexed, expected)

    def test_reindex_backfill(self):
        pass

    def test_reindex_int(self):
        ts = self.ts[::2]
        int_ts = Series(np.zeros(len(ts), dtype=int), index=ts.index)

        # this should work fine
        reindexed_int = int_ts.reindex(self.ts.index)

        # if NaNs introduced
        self.assert_(reindexed_int.dtype == np.float_)

        # NO NaNs introduced
        reindexed_int = int_ts.reindex(int_ts.index[::2])
        self.assert_(reindexed_int.dtype == np.int_)

    def test_reindex_bool(self):

        # A series other than float, int, string, or object
        ts = self.ts[::2]
        bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)

        # this should work fine
        reindexed_bool = bool_ts.reindex(self.ts.index)

        # if NaNs introduced
        self.assert_(reindexed_bool.dtype == np.object_)

        # NO NaNs introduced
        reindexed_bool = bool_ts.reindex(bool_ts.index[::2])
        self.assert_(reindexed_bool.dtype == np.bool_)

    def test_reindex_bool_pad(self):
        # fail
        ts = self.ts[5:]
        bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)
        filled_bool = bool_ts.reindex(self.ts.index, method='pad')
        self.assert_(isnull(filled_bool[:5]).all())

    def test_reindex_like(self):
        other = self.ts[::2]
        assert_series_equal(self.ts.reindex(other.index),
                            self.ts.reindex_like(other))

    def test_reindex_fill_value(self):
        #------------------------------------------------------------
        # floats
        floats = Series([1., 2., 3.])
        result = floats.reindex([1, 2, 3])
        expected = Series([2., 3., np.nan], index=[1, 2, 3])
        assert_series_equal(result, expected)

        result = floats.reindex([1, 2, 3], fill_value=0)
        expected = Series([2., 3., 0], index=[1, 2, 3])
        assert_series_equal(result, expected)

        #------------------------------------------------------------
        # ints
        ints = Series([1, 2, 3])

        result = ints.reindex([1, 2, 3])
        expected = Series([2., 3., np.nan], index=[1, 2, 3])
        assert_series_equal(result, expected)

        # don't upcast
        result = ints.reindex([1, 2, 3], fill_value=0)
        expected = Series([2, 3, 0], index=[1, 2, 3])
        self.assert_(issubclass(result.dtype.type, np.integer))
        assert_series_equal(result, expected)

        #------------------------------------------------------------
        # objects
        objects = Series([1, 2, 3], dtype=object)

        result = objects.reindex([1, 2, 3])
        expected = Series([2, 3, np.nan], index=[1, 2, 3], dtype=object)
        assert_series_equal(result, expected)

        result = objects.reindex([1, 2, 3], fill_value='foo')
        expected = Series([2, 3, 'foo'], index=[1, 2, 3], dtype=object)
        assert_series_equal(result, expected)

        #------------------------------------------------------------
        # bools
        bools = Series([True, False, True])

        result = bools.reindex([1, 2, 3])
        expected = Series([False, True, np.nan], index=[1, 2, 3], dtype=object)
        assert_series_equal(result, expected)

        result = bools.reindex([1, 2, 3], fill_value=False)
        expected = Series([False, True, False], index=[1, 2, 3])
        assert_series_equal(result, expected)

    def test_rename(self):
        renamer = lambda x: x.strftime('%Y%m%d')
        renamed = self.ts.rename(renamer)
        self.assertEqual(renamed.index[0], renamer(self.ts.index[0]))

        # dict
        rename_dict = dict(zip(self.ts.index, renamed.index))
        renamed2 = self.ts.rename(rename_dict)
        assert_series_equal(renamed, renamed2)

        # partial dict
        s = Series(np.arange(4), index=['a', 'b', 'c', 'd'])
        renamed = s.rename({'b' : 'foo', 'd' : 'bar'})
        self.assert_(np.array_equal(renamed.index, ['a', 'foo', 'c', 'bar']))

    def test_rename_inplace(self):
        renamer = lambda x: x.strftime('%Y%m%d')
        expected = renamer(self.ts.index[0])
        self.ts.rename(renamer, inplace=True)
        self.assertEqual(self.ts.index[0], expected)

    def test_preserveRefs(self):
        seq = self.ts[[5,10,15]]
        seq[1] = np.NaN
        self.assertFalse(np.isnan(self.ts[10]))

    def test_ne(self):
        ts = TimeSeries([3, 4, 5, 6, 7], [3, 4, 5, 6, 7], dtype=float)
        expected = [True, True, False, True, True]
        self.assert_(tm.equalContents(ts.index != 5, expected))
        self.assert_(tm.equalContents(~(ts.index == 5), expected))

    def test_pad_nan(self):
        x = TimeSeries([np.nan, 1., np.nan, 3., np.nan],
                       ['z', 'a', 'b', 'c', 'd'], dtype=float)
        x.fillna(method='pad', inplace=True)
        expected = TimeSeries([np.nan, 1.0, 1.0, 3.0, 3.0],
                                ['z', 'a', 'b', 'c', 'd'], dtype=float)
        assert_series_equal(x[1:], expected[1:])
        self.assert_(np.isnan(x[0]), np.isnan(expected[0]))

    def test_unstack(self):
        from numpy import nan
        from pandas.util.testing import assert_frame_equal

        index = MultiIndex(levels=[['bar', 'foo'], ['one', 'three', 'two']],
                           labels=[[1, 1, 0, 0], [0, 1, 0, 2]])

        s = Series(np.arange(4.), index=index)
        unstacked = s.unstack()

        expected = DataFrame([[2., nan, 3.], [0., 1., nan]],
                             index=['bar', 'foo'],
                             columns=['one', 'three', 'two'])

        assert_frame_equal(unstacked, expected)

        unstacked = s.unstack(level=0)
        assert_frame_equal(unstacked, expected.T)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        exp_index = MultiIndex(levels=[['one', 'two', 'three'], [0, 1]],
                               labels=[[0, 1, 2, 0, 1, 2],
                                       [0, 1, 0, 1, 0, 1]])
        expected = DataFrame({'bar' : s.values}, index=exp_index).sortlevel(0)
        unstacked = s.unstack(0)
        assert_frame_equal(unstacked, expected)

    def test_head_tail(self):
        assert_series_equal(self.series.head(), self.series[:5])
        assert_series_equal(self.series.tail(), self.series[-5:])

    def test_isin(self):
        s = Series(['A', 'B', 'C', 'a', 'B', 'B', 'A', 'C'])

        result = s.isin(['A', 'C'])
        expected = Series([True, False, True, False, False, False, True, True])
        assert_series_equal(result, expected)

    def test_fillna_int(self):
        s = Series(np.random.randint(-100, 100, 50))
        self.assert_(s.fillna(inplace=True) is s)
        assert_series_equal(s.fillna(inplace=False), s)

#-------------------------------------------------------------------------------
# TimeSeries-specific

    def test_fillna(self):
        ts = Series([0., 1., 2., 3., 4.], index=tm.makeDateIndex(5))

        self.assert_(np.array_equal(ts, ts.fillna()))

        ts[2] = np.NaN

        self.assert_(np.array_equal(ts.fillna(), [0., 1., 1., 3., 4.]))
        self.assert_(np.array_equal(ts.fillna(method='backfill'),
                                    [0., 1., 3., 3., 4.]))

        self.assert_(np.array_equal(ts.fillna(value=5), [0., 1., 5., 3., 4.]))

    def test_fillna_bug(self):
        x = Series([nan, 1., nan, 3., nan],['z','a','b','c','d'])
        filled = x.fillna(method='ffill')
        expected = Series([nan, 1., 1., 3., 3.], x.index)
        assert_series_equal(filled, expected)

        filled = x.fillna(method='bfill')
        expected = Series([1., 1., 3., 3., nan], x.index)
        assert_series_equal(filled, expected)

    def test_fillna_inplace(self):
        x = Series([nan, 1., nan, 3., nan],['z','a','b','c','d'])
        y = x.copy()

        y2 = y.fillna(value=0, inplace=True)
        self.assert_(y is y2)

        expected = x.fillna(value=0)
        assert_series_equal(y2, expected)

    def test_fillna_invalid_method(self):
        try:
            self.ts.fillna(method='ffil')
        except ValueError, inst:
            self.assert_('ffil' in str(inst))

    def test_replace(self):
        N = 100
        ser = Series(np.random.randn(N))
        ser[0:4] = np.nan
        ser[6:10] = 0

        # replace list with a single value
        rs = ser.replace([np.nan], -1, inplace=True)
        exp = ser.fillna(-1)
        assert_series_equal(rs, exp)

        rs = ser.replace(0., np.nan)
        ser[ser == 0.] = np.nan
        assert_series_equal(rs, ser)

        ser = Series(np.fabs(np.random.randn(N)), tm.makeDateIndex(N),
                     dtype=object)
        ser[:5] = np.nan
        ser[6:10] = 'foo'
        ser[20:30] = 'bar'

        # replace list with a single value
        rs = ser.replace([np.nan, 'foo', 'bar'], -1)

        self.assert_((rs[:5] == -1).all())
        self.assert_((rs[6:10] == -1).all())
        self.assert_((rs[20:30] == -1).all())
        self.assert_((isnull(ser[:5])).all())

        # replace with different values
        rs = ser.replace({np.nan : -1, 'foo' : -2, 'bar' : -3})

        self.assert_((rs[:5] == -1).all())
        self.assert_((rs[6:10] == -2).all())
        self.assert_((rs[20:30] == -3).all())
        self.assert_((isnull(ser[:5])).all())

        # replace with different values with 2 lists
        rs2 = ser.replace([np.nan, 'foo', 'bar'], [-1, -2, -3])
        assert_series_equal(rs, rs2)

        # replace with forward fill not considering np.nan missing
        s2 = ser.copy()
        s2[5] = np.nan
        rs3 = s2.replace(['foo', 'bar'])
        self.assert_(isnull(rs3[6]))

        # replace with back fill considering np.nan as missing
        rs4 = ser.replace([np.nan, 'foo', 'bar'], method='bfill')
        assert_almost_equal(rs4[4], ser[5])

        # replace inplace
        ser.replace([np.nan, 'foo', 'bar'], -1, inplace=True)
        self.assert_((ser[:5] == -1).all())
        self.assert_((ser[6:10] == -1).all())
        self.assert_((ser[20:30] == -1).all())

        ser = Series([np.nan, 0, np.inf])
        assert_series_equal(ser.replace(np.nan, 0), ser.fillna(0))

        ser = Series([np.nan, 0, 'foo', 'bar', np.inf, None, lib.NaT])
        assert_series_equal(ser.replace(np.nan, 0), ser.fillna(0))
        filled = ser.copy()
        filled[4] = 0
        assert_series_equal(ser.replace(np.inf, 0), filled)

        ser = Series(self.ts.index)
        assert_series_equal(ser.replace(np.nan, 0), ser.fillna(0))

        # malformed
        self.assertRaises(ValueError, ser.replace, [1,2,3], [np.nan, 0])
        self.assertRaises(ValueError, ser.replace, xrange(1,3), [np.nan, 0])

        ser = Series([0, 1, 2, 3, 4])
        result = ser.replace([0, 1, 2, 3, 4], [4, 3, 2, 1, 0])
        assert_series_equal(result, Series([4, 3, 2, 1, 0]))

    def test_asfreq(self):
        ts = Series([0., 1., 2.], index=[datetime(2009, 10, 30),
                                         datetime(2009, 11, 30),
                                         datetime(2009, 12, 31)])

        daily_ts = ts.asfreq('B')
        monthly_ts = daily_ts.asfreq('BM')
        self.assert_(np.array_equal(monthly_ts, ts))

        daily_ts = ts.asfreq('B', method='pad')
        monthly_ts = daily_ts.asfreq('BM')
        self.assert_(np.array_equal(monthly_ts, ts))

        daily_ts = ts.asfreq(datetools.bday)
        monthly_ts = daily_ts.asfreq(datetools.bmonthEnd)
        self.assert_(np.array_equal(monthly_ts, ts))

        result = ts[:0].asfreq('M')
        self.assert_(len(result) == 0)
        self.assert_(result is not ts)

    def test_interpolate(self):
        ts = Series(np.arange(len(self.ts), dtype=float), self.ts.index)

        ts_copy = ts.copy()
        ts_copy[5:10] = np.NaN

        linear_interp = ts_copy.interpolate(method='linear')
        self.assert_(np.array_equal(linear_interp, ts))

        ord_ts = Series([d.toordinal() for d in self.ts.index],
                        index=self.ts.index).astype(float)

        ord_ts_copy = ord_ts.copy()
        ord_ts_copy[5:10] = np.NaN

        time_interp = ord_ts_copy.interpolate(method='time')
        self.assert_(np.array_equal(time_interp, ord_ts))

        # try time interpolation on a non-TimeSeries
        self.assertRaises(Exception, self.series.interpolate, method='time')

    def test_interpolate_index_values(self):
        s = Series(np.nan, index=np.sort(np.random.rand(30)))
        s[::3] = np.random.randn(10)

        vals = s.index.values.astype(float)

        result = s.interpolate(method='values')

        expected = s.copy()
        bad = isnull(expected.values)
        good = -bad
        expected[bad] = np.interp(vals[bad], vals[good], s.values[good])

        assert_series_equal(result, expected)

    def test_weekday(self):
        # Just run the function
        weekdays = self.ts.weekday

    def test_diff(self):
        # Just run the function
        self.ts.diff()

    def test_pct_change(self):
        rs = self.ts.pct_change(fill_method=None)
        assert_series_equal(rs, self.ts / self.ts.shift(1) - 1)

        rs = self.ts.pct_change(2)
        filled = self.ts.fillna(method='pad')
        assert_series_equal(rs, filled / filled.shift(2) - 1)

        rs = self.ts.pct_change(fill_method='bfill', limit=1)
        filled = self.ts.fillna(method='bfill', limit=1)
        assert_series_equal(rs, filled / filled.shift(1) - 1)

        rs = self.ts.pct_change(freq='5D')
        filled = self.ts.fillna(method='pad')
        assert_series_equal(rs, filled / filled.shift(freq='5D') - 1)

    def test_pct_change_shift_over_nas(self):
        s = Series([1., 1.5, np.nan, 2.5, 3.])

        chg = s.pct_change()
        expected = Series([np.nan, 0.5, np.nan, 2.5/1.5 -1, .2])
        assert_series_equal(chg, expected)

    def test_autocorr(self):
        # Just run the function
        self.ts.autocorr()

    def test_first_last_valid(self):
        ts = self.ts.copy()
        ts[:5] = np.NaN

        index = ts.first_valid_index()
        self.assertEqual(index, ts.index[5])

        ts[-5:] = np.NaN
        index = ts.last_valid_index()
        self.assertEqual(index, ts.index[-6])

        ts[:] = np.nan
        self.assert_(ts.last_valid_index() is None)
        self.assert_(ts.first_valid_index() is None)

        ser = Series([], index=[])
        self.assert_(ser.last_valid_index() is None)
        self.assert_(ser.first_valid_index() is None)

    def test_mpl_compat_hack(self):
        result = self.ts[:, np.newaxis]
        expected = self.ts.values[:, np.newaxis]
        assert_almost_equal(result, expected)

#-------------------------------------------------------------------------------
# GroupBy

    def test_select(self):
        n = len(self.ts)
        result = self.ts.select(lambda x: x >= self.ts.index[n // 2])
        expected = self.ts.reindex(self.ts.index[n//2:])
        assert_series_equal(result, expected)

        result = self.ts.select(lambda x: x.weekday() == 2)
        expected = self.ts[self.ts.weekday == 2]
        assert_series_equal(result, expected)

#----------------------------------------------------------------------
# Misc not safe for sparse

    def test_dropna_preserve_name(self):
        self.ts[:5] = np.nan
        result = self.ts.dropna()
        self.assertEquals(result.name, self.ts.name)

    def test_numpy_unique(self):
        # it works!
        result = np.unique(self.ts)

class TestSeriesNonUnique(unittest.TestCase):

    def setUp(self):
        pass

    def test_basic_indexing(self):
        s = Series(np.random.randn(5), index=['a', 'b', 'a', 'a', 'b'])

        self.assertRaises(IndexError, s.__getitem__, 5)
        self.assertRaises(IndexError, s.__setitem__, 5, 0)

        self.assertRaises(KeyError, s.__getitem__, 'c')
        self.assertRaises(KeyError, s.__setitem__, 'c', 0)

        s = s.sort_index()

        self.assertRaises(IndexError, s.__getitem__, 5)
        self.assertRaises(IndexError, s.__setitem__, 5, 0)

        self.assertRaises(KeyError, s.__getitem__, 'c')
        self.assertRaises(KeyError, s.__setitem__, 'c', 0)

    def test_int_indexing(self):
        s = Series(np.random.randn(6), index=[0, 0, 1, 1, 2, 2])

        self.assertRaises(KeyError, s.__getitem__, 5)
        self.assertRaises(KeyError, s.__setitem__, 5, 0)

        self.assertRaises(KeyError, s.__getitem__, 'c')
        self.assertRaises(KeyError, s.__setitem__, 'c', 0)

        # not monotonic
        s = Series(np.random.randn(6), index=[2, 2, 0, 0, 1, 1])

        self.assertRaises(KeyError, s.__getitem__, 5)
        self.assertRaises(KeyError, s.__setitem__, 5, 0)

        self.assertRaises(KeyError, s.__getitem__, 'c')
        self.assertRaises(KeyError, s.__setitem__, 'c', 0)

    def test_datetime_indexing(self):
        from pandas import date_range

        index = date_range('1/1/2000', '1/7/2000')
        index = index.repeat(3)

        s = Series(len(index), index=index)
        stamp = Timestamp('1/8/2000')

        self.assertRaises(KeyError, s.__getitem__, stamp)
        self.assertRaises(KeyError, s.__setitem__, stamp, 0)

        # not monotonic
        s = s[::-1]

        self.assertRaises(KeyError, s.__getitem__, stamp)
        self.assertRaises(KeyError, s.__setitem__, stamp, 0)

    def test_reset_index(self):
        df = tm.makeDataFrame()[:5]
        ser = df.stack()
        ser.index.names = ['hash', 'category']

        ser.name = 'value'
        df = ser.reset_index()
        self.assert_('value' in df)

        df = ser.reset_index(name='value2')
        self.assert_('value2' in df)

        #check inplace
        s = ser.reset_index(drop=True)
        s2 = ser
        s2.reset_index(drop=True, inplace=True)
        assert_series_equal(s, s2)

        #level
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        rs = s.reset_index(level=1)
        self.assert_(len(rs.columns) == 2)

        rs = s.reset_index(level=[0, 2], drop=True)
        self.assert_(rs.index.equals(Index(index.get_level_values(1))))
        self.assert_(isinstance(rs, Series))

    def test_timeseries_coercion(self):
        idx = tm.makeDateIndex(10000)
        ser = Series(np.random.randn(len(idx)), idx.astype(object))
        self.assert_(isinstance(ser, TimeSeries))
        self.assert_(isinstance(ser.index, DatetimeIndex))

    def test_replace(self):
        N = 100
        ser = Series(np.fabs(np.random.randn(N)), tm.makeDateIndex(N),
                     dtype=object)
        ser[:5] = np.nan
        ser[6:10] = 'foo'
        ser[20:30] = 'bar'

        # replace list with a single value
        rs = ser.replace([np.nan, 'foo', 'bar'], -1)

        self.assert_((rs[:5] == -1).all())
        self.assert_((rs[6:10] == -1).all())
        self.assert_((rs[20:30] == -1).all())
        self.assert_((isnull(ser[:5])).all())

        # replace with different values
        rs = ser.replace({np.nan : -1, 'foo' : -2, 'bar' : -3})

        self.assert_((rs[:5] == -1).all())
        self.assert_((rs[6:10] == -2).all())
        self.assert_((rs[20:30] == -3).all())
        self.assert_((isnull(ser[:5])).all())

        # replace with different values with 2 lists
        rs2 = ser.replace([np.nan, 'foo', 'bar'], [-1, -2, -3])
        assert_series_equal(rs, rs2)

        # replace with forward fill not considering np.nan missing
        s2 = ser.copy()
        s2[5] = np.nan
        rs3 = s2.replace(['foo', 'bar'])
        self.assert_(isnull(rs3[6]))

        # replace with back fill considering np.nan as missing
        rs4 = ser.replace([np.nan, 'foo', 'bar'], method='bfill')
        assert_almost_equal(rs4[4], ser[5])

        # replace inplace
        ser.replace([np.nan, 'foo', 'bar'], -1, inplace=True)
        self.assert_((ser[:5] == -1).all())
        self.assert_((ser[6:10] == -1).all())
        self.assert_((ser[20:30] == -1).all())

    def test_repeat(self):
        s = Series(np.random.randn(3), index=['a', 'b', 'c'])

        reps = s.repeat(5)
        exp = Series(s.values.repeat(5), index=s.index.values.repeat(5))
        assert_series_equal(reps, exp)

        to_rep = [2, 3, 4]
        reps = s.repeat(to_rep)
        exp = Series(s.values.repeat(to_rep),
                     index=s.index.values.repeat(to_rep))
        assert_series_equal(reps, exp)

    def test_unique_data_ownership(self):
        # it works! #1807
        Series(Series(["a","c","b"]).unique()).sort()


if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
