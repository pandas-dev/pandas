# -*- coding: utf-8 -*-

from datetime import datetime, timedelta

# TODO(wesm): fix long line flake8 issues
# flake8: noqa

import pandas.util.testing as tm
from pandas.indexes.api import Index, MultiIndex
from .common import Base

from pandas.compat import (is_platform_windows, range, lrange, lzip, u,
                           zip, PY3)
import operator
import os

import numpy as np

from pandas import (period_range, date_range, Series,
                    Float64Index, Int64Index,
                    CategoricalIndex, DatetimeIndex, TimedeltaIndex,
                    PeriodIndex)
from pandas.util.testing import assert_almost_equal
from pandas.compat.numpy import np_datetime64_compat

import pandas.core.config as cf

from pandas.tseries.index import _to_m8

import pandas as pd
from pandas.lib import Timestamp


class TestIndex(Base, tm.TestCase):
    _holder = Index
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(unicodeIndex=tm.makeUnicodeIndex(100),
                            strIndex=tm.makeStringIndex(100),
                            dateIndex=tm.makeDateIndex(100),
                            periodIndex=tm.makePeriodIndex(100),
                            tdIndex=tm.makeTimedeltaIndex(100),
                            intIndex=tm.makeIntIndex(100),
                            rangeIndex=tm.makeIntIndex(100),
                            floatIndex=tm.makeFloatIndex(100),
                            boolIndex=Index([True, False]),
                            catIndex=tm.makeCategoricalIndex(100),
                            empty=Index([]),
                            tuples=MultiIndex.from_tuples(lzip(
                                ['foo', 'bar', 'baz'], [1, 2, 3])))
        self.setup_indices()

    def create_index(self):
        return Index(list('abcde'))

    def test_new_axis(self):
        new_index = self.dateIndex[None, :]
        self.assertEqual(new_index.ndim, 2)
        tm.assertIsInstance(new_index, np.ndarray)

    def test_copy_and_deepcopy(self):
        super(TestIndex, self).test_copy_and_deepcopy()

        new_copy2 = self.intIndex.copy(dtype=int)
        self.assertEqual(new_copy2.dtype.kind, 'i')

    def test_constructor(self):
        # regular instance creation
        tm.assert_contains_all(self.strIndex, self.strIndex)
        tm.assert_contains_all(self.dateIndex, self.dateIndex)

        # casting
        arr = np.array(self.strIndex)
        index = Index(arr)
        tm.assert_contains_all(arr, index)
        tm.assert_numpy_array_equal(self.strIndex, index)

        # copy
        arr = np.array(self.strIndex)
        index = Index(arr, copy=True, name='name')
        tm.assertIsInstance(index, Index)
        self.assertEqual(index.name, 'name')
        tm.assert_numpy_array_equal(arr, index)
        arr[0] = "SOMEBIGLONGSTRING"
        self.assertNotEqual(index[0], "SOMEBIGLONGSTRING")

        # what to do here?
        # arr = np.array(5.)
        # self.assertRaises(Exception, arr.view, Index)

    def test_constructor_corner(self):
        # corner case
        self.assertRaises(TypeError, Index, 0)

    def test_construction_list_mixed_tuples(self):
        # 10697
        # if we are constructing from a mixed list of tuples, make sure that we
        # are independent of the sorting order
        idx1 = Index([('A', 1), 'B'])
        self.assertIsInstance(idx1, Index) and self.assertNotInstance(
            idx1, MultiIndex)
        idx2 = Index(['B', ('A', 1)])
        self.assertIsInstance(idx2, Index) and self.assertNotInstance(
            idx2, MultiIndex)

    def test_constructor_from_index_datetimetz(self):
        idx = pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                            tz='US/Eastern')
        result = pd.Index(idx)
        tm.assert_index_equal(result, idx)
        self.assertEqual(result.tz, idx.tz)

        result = pd.Index(idx.asobject)
        tm.assert_index_equal(result, idx)
        self.assertEqual(result.tz, idx.tz)

    def test_constructor_from_index_timedelta(self):
        idx = pd.timedelta_range('1 days', freq='D', periods=3)
        result = pd.Index(idx)
        tm.assert_index_equal(result, idx)

        result = pd.Index(idx.asobject)
        tm.assert_index_equal(result, idx)

    def test_constructor_from_index_period(self):
        idx = pd.period_range('2015-01-01', freq='D', periods=3)
        result = pd.Index(idx)
        tm.assert_index_equal(result, idx)

        result = pd.Index(idx.asobject)
        tm.assert_index_equal(result, idx)

    def test_constructor_from_series_datetimetz(self):
        idx = pd.date_range('2015-01-01 10:00', freq='D', periods=3,
                            tz='US/Eastern')
        result = pd.Index(pd.Series(idx))
        tm.assert_index_equal(result, idx)
        self.assertEqual(result.tz, idx.tz)

    def test_constructor_from_series_timedelta(self):
        idx = pd.timedelta_range('1 days', freq='D', periods=3)
        result = pd.Index(pd.Series(idx))
        tm.assert_index_equal(result, idx)

    def test_constructor_from_series_period(self):
        idx = pd.period_range('2015-01-01', freq='D', periods=3)
        result = pd.Index(pd.Series(idx))
        tm.assert_index_equal(result, idx)

    def test_constructor_from_series(self):

        expected = DatetimeIndex([Timestamp('20110101'), Timestamp('20120101'),
                                  Timestamp('20130101')])
        s = Series([Timestamp('20110101'), Timestamp('20120101'), Timestamp(
            '20130101')])
        result = Index(s)
        self.assertTrue(result.equals(expected))
        result = DatetimeIndex(s)
        self.assertTrue(result.equals(expected))

        # GH 6273
        # create from a series, passing a freq
        s = Series(pd.to_datetime(['1-1-1990', '2-1-1990', '3-1-1990',
                                   '4-1-1990', '5-1-1990']))
        result = DatetimeIndex(s, freq='MS')
        expected = DatetimeIndex(
            ['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990'
             ], freq='MS')
        self.assertTrue(result.equals(expected))

        df = pd.DataFrame(np.random.rand(5, 3))
        df['date'] = ['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990',
                      '5-1-1990']
        result = DatetimeIndex(df['date'], freq='MS')
        self.assertTrue(result.equals(expected))
        self.assertEqual(df['date'].dtype, object)

        exp = pd.Series(
            ['1-1-1990', '2-1-1990', '3-1-1990', '4-1-1990', '5-1-1990'
             ], name='date')
        self.assert_series_equal(df['date'], exp)

        # GH 6274
        # infer freq of same
        result = pd.infer_freq(df['date'])
        self.assertEqual(result, 'MS')

    def test_constructor_ndarray_like(self):
        # GH 5460#issuecomment-44474502
        # it should be possible to convert any object that satisfies the numpy
        # ndarray interface directly into an Index
        class ArrayLike(object):

            def __init__(self, array):
                self.array = array

            def __array__(self, dtype=None):
                return self.array

        for array in [np.arange(5), np.array(['a', 'b', 'c']),
                      date_range('2000-01-01', periods=3).values]:
            expected = pd.Index(array)
            result = pd.Index(ArrayLike(array))
            self.assertTrue(result.equals(expected))

    def test_index_ctor_infer_periodindex(self):
        xp = period_range('2012-1-1', freq='M', periods=3)
        rs = Index(xp)
        tm.assert_numpy_array_equal(rs, xp)
        tm.assertIsInstance(rs, PeriodIndex)

    def test_constructor_simple_new(self):
        idx = Index([1, 2, 3, 4, 5], name='int')
        result = idx._simple_new(idx, 'int')
        self.assertTrue(result.equals(idx))

        idx = Index([1.1, np.nan, 2.2, 3.0], name='float')
        result = idx._simple_new(idx, 'float')
        self.assertTrue(result.equals(idx))

        idx = Index(['A', 'B', 'C', np.nan], name='obj')
        result = idx._simple_new(idx, 'obj')
        self.assertTrue(result.equals(idx))

    def test_constructor_dtypes(self):

        for idx in [Index(np.array([1, 2, 3], dtype=int)),
                    Index(np.array([1, 2, 3], dtype=int), dtype=int),
                    Index([1, 2, 3], dtype=int)]:
            self.assertIsInstance(idx, Int64Index)

        # these should coerce
        for idx in [Index(np.array([1., 2., 3.], dtype=float), dtype=int),
                    Index([1., 2., 3.], dtype=int)]:
            self.assertIsInstance(idx, Int64Index)

        for idx in [Index(np.array([1., 2., 3.], dtype=float)),
                    Index(np.array([1, 2, 3], dtype=int), dtype=float),
                    Index(np.array([1., 2., 3.], dtype=float), dtype=float),
                    Index([1, 2, 3], dtype=float),
                    Index([1., 2., 3.], dtype=float)]:
            self.assertIsInstance(idx, Float64Index)

        for idx in [Index(np.array([True, False, True], dtype=bool)),
                    Index([True, False, True]),
                    Index(np.array([True, False, True], dtype=bool), dtype=bool),
                    Index([True, False, True], dtype=bool)]:
            self.assertIsInstance(idx, Index)
            self.assertEqual(idx.dtype, object)

        for idx in [Index(np.array([1, 2, 3], dtype=int), dtype='category'),
                    Index([1, 2, 3], dtype='category'),
                    Index(np.array([np_datetime64_compat('2011-01-01'),
                                    np_datetime64_compat('2011-01-02')]), dtype='category'),
                    Index([datetime(2011, 1, 1), datetime(2011, 1, 2)], dtype='category')]:
            self.assertIsInstance(idx, CategoricalIndex)

        for idx in [Index(np.array([np_datetime64_compat('2011-01-01'),
                                    np_datetime64_compat('2011-01-02')])),
                    Index([datetime(2011, 1, 1), datetime(2011, 1, 2)])]:
            self.assertIsInstance(idx, DatetimeIndex)

        for idx in [Index(np.array([np_datetime64_compat('2011-01-01'),
                                    np_datetime64_compat('2011-01-02')]), dtype=object),
                    Index([datetime(2011, 1, 1),
                           datetime(2011, 1, 2)], dtype=object)]:
            self.assertNotIsInstance(idx, DatetimeIndex)
            self.assertIsInstance(idx, Index)
            self.assertEqual(idx.dtype, object)

        for idx in [Index(np.array([np.timedelta64(1, 'D'), np.timedelta64(
                1, 'D')])), Index([timedelta(1), timedelta(1)])]:
            self.assertIsInstance(idx, TimedeltaIndex)

        for idx in [Index(np.array([np.timedelta64(1, 'D'),
                                    np.timedelta64(1, 'D')]), dtype=object),
                    Index([timedelta(1), timedelta(1)], dtype=object)]:
            self.assertNotIsInstance(idx, TimedeltaIndex)
            self.assertIsInstance(idx, Index)
            self.assertEqual(idx.dtype, object)

    def test_view_with_args(self):

        restricted = ['unicodeIndex', 'strIndex', 'catIndex', 'boolIndex',
                      'empty']

        for i in restricted:
            ind = self.indices[i]

            # with arguments
            self.assertRaises(TypeError, lambda: ind.view('i8'))

        # these are ok
        for i in list(set(self.indices.keys()) - set(restricted)):
            ind = self.indices[i]

            # with arguments
            ind.view('i8')

    def test_legacy_pickle_identity(self):

        # GH 8431
        pth = tm.get_data_path()
        s1 = pd.read_pickle(os.path.join(pth, 's1-0.12.0.pickle'))
        s2 = pd.read_pickle(os.path.join(pth, 's2-0.12.0.pickle'))
        self.assertFalse(s1.index.identical(s2.index))
        self.assertFalse(s1.index.equals(s2.index))

    def test_astype(self):
        casted = self.intIndex.astype('i8')

        # it works!
        casted.get_loc(5)

        # pass on name
        self.intIndex.name = 'foobar'
        casted = self.intIndex.astype('i8')
        self.assertEqual(casted.name, 'foobar')

    def test_equals(self):
        # same
        self.assertTrue(Index(['a', 'b', 'c']).equals(Index(['a', 'b', 'c'])))

        # different length
        self.assertFalse(Index(['a', 'b', 'c']).equals(Index(['a', 'b'])))

        # same length, different values
        self.assertFalse(Index(['a', 'b', 'c']).equals(Index(['a', 'b', 'd'])))

        # Must also be an Index
        self.assertFalse(Index(['a', 'b', 'c']).equals(['a', 'b', 'c']))

    def test_insert(self):

        # GH 7256
        # validate neg/pos inserts
        result = Index(['b', 'c', 'd'])

        # test 0th element
        self.assertTrue(Index(['a', 'b', 'c', 'd']).equals(result.insert(0,
                                                                         'a')))

        # test Nth element that follows Python list behavior
        self.assertTrue(Index(['b', 'c', 'e', 'd']).equals(result.insert(-1,
                                                                         'e')))

        # test loc +/- neq (0, -1)
        self.assertTrue(result.insert(1, 'z').equals(result.insert(-2, 'z')))

        # test empty
        null_index = Index([])
        self.assertTrue(Index(['a']).equals(null_index.insert(0, 'a')))

    def test_delete(self):
        idx = Index(['a', 'b', 'c', 'd'], name='idx')

        expected = Index(['b', 'c', 'd'], name='idx')
        result = idx.delete(0)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)

        expected = Index(['a', 'b', 'c'], name='idx')
        result = idx.delete(-1)
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)

        with tm.assertRaises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = idx.delete(5)

    def test_identical(self):

        # index
        i1 = Index(['a', 'b', 'c'])
        i2 = Index(['a', 'b', 'c'])

        self.assertTrue(i1.identical(i2))

        i1 = i1.rename('foo')
        self.assertTrue(i1.equals(i2))
        self.assertFalse(i1.identical(i2))

        i2 = i2.rename('foo')
        self.assertTrue(i1.identical(i2))

        i3 = Index([('a', 'a'), ('a', 'b'), ('b', 'a')])
        i4 = Index([('a', 'a'), ('a', 'b'), ('b', 'a')], tupleize_cols=False)
        self.assertFalse(i3.identical(i4))

    def test_is_(self):
        ind = Index(range(10))
        self.assertTrue(ind.is_(ind))
        self.assertTrue(ind.is_(ind.view().view().view().view()))
        self.assertFalse(ind.is_(Index(range(10))))
        self.assertFalse(ind.is_(ind.copy()))
        self.assertFalse(ind.is_(ind.copy(deep=False)))
        self.assertFalse(ind.is_(ind[:]))
        self.assertFalse(ind.is_(ind.view(np.ndarray).view(Index)))
        self.assertFalse(ind.is_(np.array(range(10))))

        # quasi-implementation dependent
        self.assertTrue(ind.is_(ind.view()))
        ind2 = ind.view()
        ind2.name = 'bob'
        self.assertTrue(ind.is_(ind2))
        self.assertTrue(ind2.is_(ind))
        # doesn't matter if Indices are *actually* views of underlying data,
        self.assertFalse(ind.is_(Index(ind.values)))
        arr = np.array(range(1, 11))
        ind1 = Index(arr, copy=False)
        ind2 = Index(arr, copy=False)
        self.assertFalse(ind1.is_(ind2))

    def test_asof(self):
        d = self.dateIndex[0]
        self.assertEqual(self.dateIndex.asof(d), d)
        self.assertTrue(np.isnan(self.dateIndex.asof(d - timedelta(1))))

        d = self.dateIndex[-1]
        self.assertEqual(self.dateIndex.asof(d + timedelta(1)), d)

        d = self.dateIndex[0].to_datetime()
        tm.assertIsInstance(self.dateIndex.asof(d), Timestamp)

    def test_asof_datetime_partial(self):
        idx = pd.date_range('2010-01-01', periods=2, freq='m')
        expected = Timestamp('2010-02-28')
        result = idx.asof('2010-02')
        self.assertEqual(result, expected)
        self.assertFalse(isinstance(result, Index))

    def test_nanosecond_index_access(self):
        s = Series([Timestamp('20130101')]).values.view('i8')[0]
        r = DatetimeIndex([s + 50 + i for i in range(100)])
        x = Series(np.random.randn(100), index=r)

        first_value = x.asof(x.index[0])

        # this does not yet work, as parsing strings is done via dateutil
        # self.assertEqual(first_value,
        #                  x['2013-01-01 00:00:00.000000050+0000'])

        self.assertEqual(
            first_value,
            x[Timestamp(np_datetime64_compat('2013-01-01 00:00:00.000000050+0000',
                                             'ns'))])

    def test_comparators(self):
        index = self.dateIndex
        element = index[len(index) // 2]
        element = _to_m8(element)

        arr = np.array(index)

        def _check(op):
            arr_result = op(arr, element)
            index_result = op(index, element)

            self.assertIsInstance(index_result, np.ndarray)
            tm.assert_numpy_array_equal(arr_result, index_result)

        _check(operator.eq)
        _check(operator.ne)
        _check(operator.gt)
        _check(operator.lt)
        _check(operator.ge)
        _check(operator.le)

    def test_booleanindex(self):
        boolIdx = np.repeat(True, len(self.strIndex)).astype(bool)
        boolIdx[5:30:2] = False

        subIndex = self.strIndex[boolIdx]

        for i, val in enumerate(subIndex):
            self.assertEqual(subIndex.get_loc(val), i)

        subIndex = self.strIndex[list(boolIdx)]
        for i, val in enumerate(subIndex):
            self.assertEqual(subIndex.get_loc(val), i)

    def test_fancy(self):
        sl = self.strIndex[[1, 2, 3]]
        for i in sl:
            self.assertEqual(i, sl[sl.get_loc(i)])

    def test_empty_fancy(self):
        empty_farr = np.array([], dtype=np.float_)
        empty_iarr = np.array([], dtype=np.int_)
        empty_barr = np.array([], dtype=np.bool_)

        # pd.DatetimeIndex is excluded, because it overrides getitem and should
        # be tested separately.
        for idx in [self.strIndex, self.intIndex, self.floatIndex]:
            empty_idx = idx.__class__([])

            self.assertTrue(idx[[]].identical(empty_idx))
            self.assertTrue(idx[empty_iarr].identical(empty_idx))
            self.assertTrue(idx[empty_barr].identical(empty_idx))

            # np.ndarray only accepts ndarray of int & bool dtypes, so should
            # Index.
            self.assertRaises(IndexError, idx.__getitem__, empty_farr)

    def test_getitem(self):
        arr = np.array(self.dateIndex)
        exp = self.dateIndex[5]
        exp = _to_m8(exp)

        self.assertEqual(exp, arr[5])

    def test_intersection(self):
        first = self.strIndex[:20]
        second = self.strIndex[:10]
        intersect = first.intersection(second)
        self.assertTrue(tm.equalContents(intersect, second))

        # Corner cases
        inter = first.intersection(first)
        self.assertIs(inter, first)

        idx1 = Index([1, 2, 3, 4, 5], name='idx')
        # if target has the same name, it is preserved
        idx2 = Index([3, 4, 5, 6, 7], name='idx')
        expected2 = Index([3, 4, 5], name='idx')
        result2 = idx1.intersection(idx2)
        self.assertTrue(result2.equals(expected2))
        self.assertEqual(result2.name, expected2.name)

        # if target name is different, it will be reset
        idx3 = Index([3, 4, 5, 6, 7], name='other')
        expected3 = Index([3, 4, 5], name=None)
        result3 = idx1.intersection(idx3)
        self.assertTrue(result3.equals(expected3))
        self.assertEqual(result3.name, expected3.name)

        # non monotonic
        idx1 = Index([5, 3, 2, 4, 1], name='idx')
        idx2 = Index([4, 7, 6, 5, 3], name='idx')
        result2 = idx1.intersection(idx2)
        self.assertTrue(tm.equalContents(result2, expected2))
        self.assertEqual(result2.name, expected2.name)

        idx3 = Index([4, 7, 6, 5, 3], name='other')
        result3 = idx1.intersection(idx3)
        self.assertTrue(tm.equalContents(result3, expected3))
        self.assertEqual(result3.name, expected3.name)

        # non-monotonic non-unique
        idx1 = Index(['A', 'B', 'A', 'C'])
        idx2 = Index(['B', 'D'])
        expected = Index(['B'], dtype='object')
        result = idx1.intersection(idx2)
        self.assertTrue(result.equals(expected))

        # preserve names
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        first.name = 'A'
        second.name = 'A'
        intersect = first.intersection(second)
        self.assertEqual(intersect.name, 'A')

        second.name = 'B'
        intersect = first.intersection(second)
        self.assertIsNone(intersect.name)

        first.name = None
        second.name = 'B'
        intersect = first.intersection(second)
        self.assertIsNone(intersect.name)

    def test_union(self):
        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        everything = self.strIndex[:20]
        union = first.union(second)
        self.assertTrue(tm.equalContents(union, everything))

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.union(case)
            self.assertTrue(tm.equalContents(result, everything))

        # Corner cases
        union = first.union(first)
        self.assertIs(union, first)

        union = first.union([])
        self.assertIs(union, first)

        union = Index([]).union(first)
        self.assertIs(union, first)

        # preserve names
        first = Index(list('ab'), name='A')
        second = Index(list('ab'), name='B')
        union = first.union(second)
        self.assertIsNone(union.name)

        first = Index(list('ab'), name='A')
        second = Index([], name='B')
        union = first.union(second)
        self.assertIsNone(union.name)

        first = Index([], name='A')
        second = Index(list('ab'), name='B')
        union = first.union(second)
        self.assertIsNone(union.name)

        first = Index(list('ab'))
        second = Index(list('ab'), name='B')
        union = first.union(second)
        self.assertEqual(union.name, 'B')

        first = Index([])
        second = Index(list('ab'), name='B')
        union = first.union(second)
        self.assertEqual(union.name, 'B')

        first = Index(list('ab'))
        second = Index([], name='B')
        union = first.union(second)
        self.assertEqual(union.name, 'B')

        first = Index(list('ab'), name='A')
        second = Index(list('ab'))
        union = first.union(second)
        self.assertEqual(union.name, 'A')

        first = Index(list('ab'), name='A')
        second = Index([])
        union = first.union(second)
        self.assertEqual(union.name, 'A')

        first = Index([], name='A')
        second = Index(list('ab'))
        union = first.union(second)
        self.assertEqual(union.name, 'A')

    def test_add(self):

        # - API change GH 8226
        with tm.assert_produces_warning():
            self.strIndex + self.strIndex
        with tm.assert_produces_warning():
            self.strIndex + self.strIndex.tolist()
        with tm.assert_produces_warning():
            self.strIndex.tolist() + self.strIndex

        with tm.assert_produces_warning(RuntimeWarning):
            firstCat = self.strIndex.union(self.dateIndex)
        secondCat = self.strIndex.union(self.strIndex)

        if self.dateIndex.dtype == np.object_:
            appended = np.append(self.strIndex, self.dateIndex)
        else:
            appended = np.append(self.strIndex, self.dateIndex.astype('O'))

        self.assertTrue(tm.equalContents(firstCat, appended))
        self.assertTrue(tm.equalContents(secondCat, self.strIndex))
        tm.assert_contains_all(self.strIndex, firstCat)
        tm.assert_contains_all(self.strIndex, secondCat)
        tm.assert_contains_all(self.dateIndex, firstCat)

        # test add and radd
        idx = Index(list('abc'))
        expected = Index(['a1', 'b1', 'c1'])
        self.assert_index_equal(idx + '1', expected)
        expected = Index(['1a', '1b', '1c'])
        self.assert_index_equal('1' + idx, expected)

    def test_append_multiple(self):
        index = Index(['a', 'b', 'c', 'd', 'e', 'f'])

        foos = [index[:2], index[2:4], index[4:]]
        result = foos[0].append(foos[1:])
        self.assertTrue(result.equals(index))

        # empty
        result = index.append([])
        self.assertTrue(result.equals(index))

    def test_append_empty_preserve_name(self):
        left = Index([], name='foo')
        right = Index([1, 2, 3], name='foo')

        result = left.append(right)
        self.assertEqual(result.name, 'foo')

        left = Index([], name='foo')
        right = Index([1, 2, 3], name='bar')

        result = left.append(right)
        self.assertIsNone(result.name)

    def test_add_string(self):
        # from bug report
        index = Index(['a', 'b', 'c'])
        index2 = index + 'foo'

        self.assertNotIn('a', index2)
        self.assertIn('afoo', index2)

    def test_iadd_string(self):
        index = pd.Index(['a', 'b', 'c'])
        # doesn't fail test unless there is a check before `+=`
        self.assertIn('a', index)

        index += '_x'
        self.assertIn('a_x', index)

    def test_difference(self):

        first = self.strIndex[5:20]
        second = self.strIndex[:10]
        answer = self.strIndex[10:20]
        first.name = 'name'
        # different names
        result = first.difference(second)

        self.assertTrue(tm.equalContents(result, answer))
        self.assertEqual(result.name, None)

        # same names
        second.name = 'name'
        result = first.difference(second)
        self.assertEqual(result.name, 'name')

        # with empty
        result = first.difference([])
        self.assertTrue(tm.equalContents(result, first))
        self.assertEqual(result.name, first.name)

        # with everythin
        result = first.difference(first)
        self.assertEqual(len(result), 0)
        self.assertEqual(result.name, first.name)

    def test_symmetric_difference(self):
        # smoke
        idx1 = Index([1, 2, 3, 4], name='idx1')
        idx2 = Index([2, 3, 4, 5])
        result = idx1.symmetric_difference(idx2)
        expected = Index([1, 5])
        self.assertTrue(tm.equalContents(result, expected))
        self.assertIsNone(result.name)

        # __xor__ syntax
        expected = idx1 ^ idx2
        self.assertTrue(tm.equalContents(result, expected))
        self.assertIsNone(result.name)

        # multiIndex
        idx1 = MultiIndex.from_tuples(self.tuples)
        idx2 = MultiIndex.from_tuples([('foo', 1), ('bar', 3)])
        result = idx1.symmetric_difference(idx2)
        expected = MultiIndex.from_tuples([('bar', 2), ('baz', 3), ('bar', 3)])
        self.assertTrue(tm.equalContents(result, expected))

        # nans:
        # GH #6444, sorting of nans. Make sure the number of nans is right
        # and the correct non-nan values are there. punt on sorting.
        idx1 = Index([1, 2, 3, np.nan])
        idx2 = Index([0, 1, np.nan])
        result = idx1.symmetric_difference(idx2)
        # expected = Index([0.0, np.nan, 2.0, 3.0, np.nan])

        nans = pd.isnull(result)
        self.assertEqual(nans.sum(), 1)
        self.assertEqual((~nans).sum(), 3)
        [self.assertIn(x, result) for x in [0.0, 2.0, 3.0]]

        # other not an Index:
        idx1 = Index([1, 2, 3, 4], name='idx1')
        idx2 = np.array([2, 3, 4, 5])
        expected = Index([1, 5])
        result = idx1.symmetric_difference(idx2)
        self.assertTrue(tm.equalContents(result, expected))
        self.assertEqual(result.name, 'idx1')

        result = idx1.symmetric_difference(idx2, result_name='new_name')
        self.assertTrue(tm.equalContents(result, expected))
        self.assertEqual(result.name, 'new_name')

    def test_is_numeric(self):
        self.assertFalse(self.dateIndex.is_numeric())
        self.assertFalse(self.strIndex.is_numeric())
        self.assertTrue(self.intIndex.is_numeric())
        self.assertTrue(self.floatIndex.is_numeric())
        self.assertFalse(self.catIndex.is_numeric())

    def test_is_object(self):
        self.assertTrue(self.strIndex.is_object())
        self.assertTrue(self.boolIndex.is_object())
        self.assertFalse(self.catIndex.is_object())
        self.assertFalse(self.intIndex.is_object())
        self.assertFalse(self.dateIndex.is_object())
        self.assertFalse(self.floatIndex.is_object())

    def test_is_all_dates(self):
        self.assertTrue(self.dateIndex.is_all_dates)
        self.assertFalse(self.strIndex.is_all_dates)
        self.assertFalse(self.intIndex.is_all_dates)

    def test_summary(self):
        self._check_method_works(Index.summary)
        # GH3869
        ind = Index(['{other}%s', "~:{range}:0"], name='A')
        result = ind.summary()
        # shouldn't be formatted accidentally.
        self.assertIn('~:{range}:0', result)
        self.assertIn('{other}%s', result)

    def test_format(self):
        self._check_method_works(Index.format)

        index = Index([datetime.now()])

        # windows has different precision on datetime.datetime.now (it doesn't
        # include us since the default for Timestamp shows these but Index
        # formating does not we are skipping
        if not is_platform_windows():
            formatted = index.format()
            expected = [str(index[0])]
            self.assertEqual(formatted, expected)

        # 2845
        index = Index([1, 2.0 + 3.0j, np.nan])
        formatted = index.format()
        expected = [str(index[0]), str(index[1]), u('NaN')]
        self.assertEqual(formatted, expected)

        # is this really allowed?
        index = Index([1, 2.0 + 3.0j, None])
        formatted = index.format()
        expected = [str(index[0]), str(index[1]), u('NaN')]
        self.assertEqual(formatted, expected)

        self.strIndex[:0].format()

    def test_format_with_name_time_info(self):
        # bug I fixed 12/20/2011
        inc = timedelta(hours=4)
        dates = Index([dt + inc for dt in self.dateIndex], name='something')

        formatted = dates.format(name=True)
        self.assertEqual(formatted[0], 'something')

    def test_format_datetime_with_time(self):
        t = Index([datetime(2012, 2, 7), datetime(2012, 2, 7, 23)])

        result = t.format()
        expected = ['2012-02-07 00:00:00', '2012-02-07 23:00:00']
        self.assertEqual(len(result), 2)
        self.assertEqual(result, expected)

    def test_format_none(self):
        values = ['a', 'b', 'c', None]

        idx = Index(values)
        idx.format()
        self.assertIsNone(idx[3])

    def test_logical_compat(self):
        idx = self.create_index()
        self.assertEqual(idx.all(), idx.values.all())
        self.assertEqual(idx.any(), idx.values.any())

    def _check_method_works(self, method):
        method(self.empty)
        method(self.dateIndex)
        method(self.unicodeIndex)
        method(self.strIndex)
        method(self.intIndex)
        method(self.tuples)
        method(self.catIndex)

    def test_get_indexer(self):
        idx1 = Index([1, 2, 3, 4, 5])
        idx2 = Index([2, 4, 6])

        r1 = idx1.get_indexer(idx2)
        assert_almost_equal(r1, [1, 3, -1])

        r1 = idx2.get_indexer(idx1, method='pad')
        e1 = [-1, 0, 0, 1, 1]
        assert_almost_equal(r1, e1)

        r2 = idx2.get_indexer(idx1[::-1], method='pad')
        assert_almost_equal(r2, e1[::-1])

        rffill1 = idx2.get_indexer(idx1, method='ffill')
        assert_almost_equal(r1, rffill1)

        r1 = idx2.get_indexer(idx1, method='backfill')
        e1 = [0, 0, 1, 1, 2]
        assert_almost_equal(r1, e1)

        rbfill1 = idx2.get_indexer(idx1, method='bfill')
        assert_almost_equal(r1, rbfill1)

        r2 = idx2.get_indexer(idx1[::-1], method='backfill')
        assert_almost_equal(r2, e1[::-1])

    def test_get_indexer_invalid(self):
        # GH10411
        idx = Index(np.arange(10))

        with tm.assertRaisesRegexp(ValueError, 'tolerance argument'):
            idx.get_indexer([1, 0], tolerance=1)

        with tm.assertRaisesRegexp(ValueError, 'limit argument'):
            idx.get_indexer([1, 0], limit=1)

    def test_get_indexer_nearest(self):
        idx = Index(np.arange(10))

        all_methods = ['pad', 'backfill', 'nearest']
        for method in all_methods:
            actual = idx.get_indexer([0, 5, 9], method=method)
            tm.assert_numpy_array_equal(actual, [0, 5, 9])

            actual = idx.get_indexer([0, 5, 9], method=method, tolerance=0)
            tm.assert_numpy_array_equal(actual, [0, 5, 9])

        for method, expected in zip(all_methods, [[0, 1, 8], [1, 2, 9], [0, 2,
                                                                         9]]):
            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method)
            tm.assert_numpy_array_equal(actual, expected)

            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method,
                                     tolerance=1)
            tm.assert_numpy_array_equal(actual, expected)

        for method, expected in zip(all_methods, [[0, -1, -1], [-1, 2, -1],
                                                  [0, 2, -1]]):
            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method,
                                     tolerance=0.2)
            tm.assert_numpy_array_equal(actual, expected)

        with tm.assertRaisesRegexp(ValueError, 'limit argument'):
            idx.get_indexer([1, 0], method='nearest', limit=1)

    def test_get_indexer_nearest_decreasing(self):
        idx = Index(np.arange(10))[::-1]

        all_methods = ['pad', 'backfill', 'nearest']
        for method in all_methods:
            actual = idx.get_indexer([0, 5, 9], method=method)
            tm.assert_numpy_array_equal(actual, [9, 4, 0])

        for method, expected in zip(all_methods, [[8, 7, 0], [9, 8, 1], [9, 7,
                                                                         0]]):
            actual = idx.get_indexer([0.2, 1.8, 8.5], method=method)
            tm.assert_numpy_array_equal(actual, expected)

    def test_get_indexer_strings(self):
        idx = pd.Index(['b', 'c'])

        actual = idx.get_indexer(['a', 'b', 'c', 'd'], method='pad')
        expected = [-1, 0, 1, 1]
        tm.assert_numpy_array_equal(actual, expected)

        actual = idx.get_indexer(['a', 'b', 'c', 'd'], method='backfill')
        expected = [0, 0, 1, -1]
        tm.assert_numpy_array_equal(actual, expected)

        with tm.assertRaises(TypeError):
            idx.get_indexer(['a', 'b', 'c', 'd'], method='nearest')

        with tm.assertRaises(TypeError):
            idx.get_indexer(['a', 'b', 'c', 'd'], method='pad', tolerance=2)

    def test_get_loc(self):
        idx = pd.Index([0, 1, 2])
        all_methods = [None, 'pad', 'backfill', 'nearest']
        for method in all_methods:
            self.assertEqual(idx.get_loc(1, method=method), 1)
            if method is not None:
                self.assertEqual(idx.get_loc(1, method=method, tolerance=0), 1)
            with tm.assertRaises(TypeError):
                idx.get_loc([1, 2], method=method)

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc(1.1, method), loc)

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc(1.1, method, tolerance=1), loc)

        for method in ['pad', 'backfill', 'nearest']:
            with tm.assertRaises(KeyError):
                idx.get_loc(1.1, method, tolerance=0.05)

        with tm.assertRaisesRegexp(ValueError, 'must be numeric'):
            idx.get_loc(1.1, 'nearest', tolerance='invalid')
        with tm.assertRaisesRegexp(ValueError, 'tolerance .* valid if'):
            idx.get_loc(1.1, tolerance=1)

        idx = pd.Index(['a', 'c'])
        with tm.assertRaises(TypeError):
            idx.get_loc('a', method='nearest')
        with tm.assertRaises(TypeError):
            idx.get_loc('a', method='pad', tolerance='invalid')

    def test_slice_locs(self):
        for dtype in [int, float]:
            idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=dtype))
            n = len(idx)

            self.assertEqual(idx.slice_locs(start=2), (2, n))
            self.assertEqual(idx.slice_locs(start=3), (3, n))
            self.assertEqual(idx.slice_locs(3, 8), (3, 6))
            self.assertEqual(idx.slice_locs(5, 10), (3, n))
            self.assertEqual(idx.slice_locs(end=8), (0, 6))
            self.assertEqual(idx.slice_locs(end=9), (0, 7))

            # reversed
            idx2 = idx[::-1]
            self.assertEqual(idx2.slice_locs(8, 2), (2, 6))
            self.assertEqual(idx2.slice_locs(7, 3), (2, 5))

        # float slicing
        idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=float))
        n = len(idx)
        self.assertEqual(idx.slice_locs(5.0, 10.0), (3, n))
        self.assertEqual(idx.slice_locs(4.5, 10.5), (3, 8))
        idx2 = idx[::-1]
        self.assertEqual(idx2.slice_locs(8.5, 1.5), (2, 6))
        self.assertEqual(idx2.slice_locs(10.5, -1), (0, n))

        # int slicing with floats
        # GH 4892, these are all TypeErrors
        idx = Index(np.array([0, 1, 2, 5, 6, 7, 9, 10], dtype=int))
        self.assertRaises(TypeError,
                          lambda: idx.slice_locs(5.0, 10.0), (3, n))
        self.assertRaises(TypeError,
                          lambda: idx.slice_locs(4.5, 10.5), (3, 8))
        idx2 = idx[::-1]
        self.assertRaises(TypeError,
                          lambda: idx2.slice_locs(8.5, 1.5), (2, 6))
        self.assertRaises(TypeError,
                          lambda: idx2.slice_locs(10.5, -1), (0, n))

    def test_slice_locs_dup(self):
        idx = Index(['a', 'a', 'b', 'c', 'd', 'd'])
        self.assertEqual(idx.slice_locs('a', 'd'), (0, 6))
        self.assertEqual(idx.slice_locs(end='d'), (0, 6))
        self.assertEqual(idx.slice_locs('a', 'c'), (0, 4))
        self.assertEqual(idx.slice_locs('b', 'd'), (2, 6))

        idx2 = idx[::-1]
        self.assertEqual(idx2.slice_locs('d', 'a'), (0, 6))
        self.assertEqual(idx2.slice_locs(end='a'), (0, 6))
        self.assertEqual(idx2.slice_locs('d', 'b'), (0, 4))
        self.assertEqual(idx2.slice_locs('c', 'a'), (2, 6))

        for dtype in [int, float]:
            idx = Index(np.array([10, 12, 12, 14], dtype=dtype))
            self.assertEqual(idx.slice_locs(12, 12), (1, 3))
            self.assertEqual(idx.slice_locs(11, 13), (1, 3))

            idx2 = idx[::-1]
            self.assertEqual(idx2.slice_locs(12, 12), (1, 3))
            self.assertEqual(idx2.slice_locs(13, 11), (1, 3))

    def test_slice_locs_na(self):
        idx = Index([np.nan, 1, 2])
        self.assertRaises(KeyError, idx.slice_locs, start=1.5)
        self.assertRaises(KeyError, idx.slice_locs, end=1.5)
        self.assertEqual(idx.slice_locs(1), (1, 3))
        self.assertEqual(idx.slice_locs(np.nan), (0, 3))

        idx = Index([0, np.nan, np.nan, 1, 2])
        self.assertEqual(idx.slice_locs(np.nan), (1, 5))

    def test_slice_locs_negative_step(self):
        idx = Index(list('bcdxy'))

        SLC = pd.IndexSlice

        def check_slice(in_slice, expected):
            s_start, s_stop = idx.slice_locs(in_slice.start, in_slice.stop,
                                             in_slice.step)
            result = idx[s_start:s_stop:in_slice.step]
            expected = pd.Index(list(expected))
            self.assertTrue(result.equals(expected))

        for in_slice, expected in [
                (SLC[::-1], 'yxdcb'), (SLC['b':'y':-1], ''),
                (SLC['b'::-1], 'b'), (SLC[:'b':-1], 'yxdcb'),
                (SLC[:'y':-1], 'y'), (SLC['y'::-1], 'yxdcb'),
                (SLC['y'::-4], 'yb'),
                # absent labels
                (SLC[:'a':-1], 'yxdcb'), (SLC[:'a':-2], 'ydb'),
                (SLC['z'::-1], 'yxdcb'), (SLC['z'::-3], 'yc'),
                (SLC['m'::-1], 'dcb'), (SLC[:'m':-1], 'yx'),
                (SLC['a':'a':-1], ''), (SLC['z':'z':-1], ''),
                (SLC['m':'m':-1], '')
        ]:
            check_slice(in_slice, expected)

    def test_drop(self):
        n = len(self.strIndex)

        drop = self.strIndex[lrange(5, 10)]
        dropped = self.strIndex.drop(drop)
        expected = self.strIndex[lrange(5) + lrange(10, n)]
        self.assertTrue(dropped.equals(expected))

        self.assertRaises(ValueError, self.strIndex.drop, ['foo', 'bar'])
        self.assertRaises(ValueError, self.strIndex.drop, ['1', 'bar'])

        # errors='ignore'
        mixed = drop.tolist() + ['foo']
        dropped = self.strIndex.drop(mixed, errors='ignore')
        expected = self.strIndex[lrange(5) + lrange(10, n)]
        self.assert_index_equal(dropped, expected)

        dropped = self.strIndex.drop(['foo', 'bar'], errors='ignore')
        expected = self.strIndex[lrange(n)]
        self.assert_index_equal(dropped, expected)

        dropped = self.strIndex.drop(self.strIndex[0])
        expected = self.strIndex[1:]
        self.assert_index_equal(dropped, expected)

        ser = Index([1, 2, 3])
        dropped = ser.drop(1)
        expected = Index([2, 3])
        self.assert_index_equal(dropped, expected)

        # errors='ignore'
        self.assertRaises(ValueError, ser.drop, [3, 4])

        dropped = ser.drop(4, errors='ignore')
        expected = Index([1, 2, 3])
        self.assert_index_equal(dropped, expected)

        dropped = ser.drop([3, 4, 5], errors='ignore')
        expected = Index([1, 2])
        self.assert_index_equal(dropped, expected)

    def test_tuple_union_bug(self):
        import pandas
        import numpy as np

        aidx1 = np.array([(1, 'A'), (2, 'A'), (1, 'B'), (2, 'B')],
                         dtype=[('num', int), ('let', 'a1')])
        aidx2 = np.array([(1, 'A'), (2, 'A'), (1, 'B'),
                          (2, 'B'), (1, 'C'), (2, 'C')],
                         dtype=[('num', int), ('let', 'a1')])

        idx1 = pandas.Index(aidx1)
        idx2 = pandas.Index(aidx2)

        # intersection broken?
        int_idx = idx1.intersection(idx2)
        # needs to be 1d like idx1 and idx2
        expected = idx1[:4]  # pandas.Index(sorted(set(idx1) & set(idx2)))
        self.assertEqual(int_idx.ndim, 1)
        self.assertTrue(int_idx.equals(expected))

        # union broken
        union_idx = idx1.union(idx2)
        expected = idx2
        self.assertEqual(union_idx.ndim, 1)
        self.assertTrue(union_idx.equals(expected))

    def test_is_monotonic_incomparable(self):
        index = Index([5, datetime.now(), 7])
        self.assertFalse(index.is_monotonic)
        self.assertFalse(index.is_monotonic_decreasing)

    def test_get_set_value(self):
        values = np.random.randn(100)
        date = self.dateIndex[67]

        assert_almost_equal(self.dateIndex.get_value(values, date), values[67])

        self.dateIndex.set_value(values, date, 10)
        self.assertEqual(values[67], 10)

    def test_isin(self):
        values = ['foo', 'bar', 'quux']

        idx = Index(['qux', 'baz', 'foo', 'bar'])
        result = idx.isin(values)
        expected = np.array([False, False, True, True])
        tm.assert_numpy_array_equal(result, expected)

        # set
        result = idx.isin(set(values))
        tm.assert_numpy_array_equal(result, expected)

        # empty, return dtype bool
        idx = Index([])
        result = idx.isin(values)
        self.assertEqual(len(result), 0)
        self.assertEqual(result.dtype, np.bool_)

    def test_isin_nan(self):
        tm.assert_numpy_array_equal(
            Index(['a', np.nan]).isin([np.nan]), [False, True])
        tm.assert_numpy_array_equal(
            Index(['a', pd.NaT]).isin([pd.NaT]), [False, True])
        tm.assert_numpy_array_equal(
            Index(['a', np.nan]).isin([float('nan')]), [False, False])
        tm.assert_numpy_array_equal(
            Index(['a', np.nan]).isin([pd.NaT]), [False, False])
        # Float64Index overrides isin, so must be checked separately
        tm.assert_numpy_array_equal(
            Float64Index([1.0, np.nan]).isin([np.nan]), [False, True])
        tm.assert_numpy_array_equal(
            Float64Index([1.0, np.nan]).isin([float('nan')]), [False, True])
        tm.assert_numpy_array_equal(
            Float64Index([1.0, np.nan]).isin([pd.NaT]), [False, True])

    def test_isin_level_kwarg(self):
        def check_idx(idx):
            values = idx.tolist()[-2:] + ['nonexisting']

            expected = np.array([False, False, True, True])
            tm.assert_numpy_array_equal(expected, idx.isin(values, level=0))
            tm.assert_numpy_array_equal(expected, idx.isin(values, level=-1))

            self.assertRaises(IndexError, idx.isin, values, level=1)
            self.assertRaises(IndexError, idx.isin, values, level=10)
            self.assertRaises(IndexError, idx.isin, values, level=-2)

            self.assertRaises(KeyError, idx.isin, values, level=1.0)
            self.assertRaises(KeyError, idx.isin, values, level='foobar')

            idx.name = 'foobar'
            tm.assert_numpy_array_equal(expected,
                                        idx.isin(values, level='foobar'))

            self.assertRaises(KeyError, idx.isin, values, level='xyzzy')
            self.assertRaises(KeyError, idx.isin, values, level=np.nan)

        check_idx(Index(['qux', 'baz', 'foo', 'bar']))
        # Float64Index overrides isin, so must be checked separately
        check_idx(Float64Index([1.0, 2.0, 3.0, 4.0]))

    def test_boolean_cmp(self):
        values = [1, 2, 3, 4]

        idx = Index(values)
        res = (idx == values)

        tm.assert_numpy_array_equal(res, np.array(
            [True, True, True, True], dtype=bool))

    def test_get_level_values(self):
        result = self.strIndex.get_level_values(0)
        self.assertTrue(result.equals(self.strIndex))

    def test_slice_keep_name(self):
        idx = Index(['a', 'b'], name='asdf')
        self.assertEqual(idx.name, idx[1:].name)

    def test_join_self(self):
        # instance attributes of the form self.<name>Index
        indices = 'unicode', 'str', 'date', 'int', 'float'
        kinds = 'outer', 'inner', 'left', 'right'
        for index_kind in indices:
            res = getattr(self, '{0}Index'.format(index_kind))

            for kind in kinds:
                joined = res.join(res, how=kind)
                self.assertIs(res, joined)

    def test_str_attribute(self):
        # GH9068
        methods = ['strip', 'rstrip', 'lstrip']
        idx = Index([' jack', 'jill ', ' jesse ', 'frank'])
        for method in methods:
            expected = Index([getattr(str, method)(x) for x in idx.values])
            tm.assert_index_equal(
                getattr(Index.str, method)(idx.str), expected)

        # create a few instances that are not able to use .str accessor
        indices = [Index(range(5)), tm.makeDateIndex(10),
                   MultiIndex.from_tuples([('foo', '1'), ('bar', '3')]),
                   PeriodIndex(start='2000', end='2010', freq='A')]
        for idx in indices:
            with self.assertRaisesRegexp(AttributeError,
                                         'only use .str accessor'):
                idx.str.repeat(2)

        idx = Index(['a b c', 'd e', 'f'])
        expected = Index([['a', 'b', 'c'], ['d', 'e'], ['f']])
        tm.assert_index_equal(idx.str.split(), expected)
        tm.assert_index_equal(idx.str.split(expand=False), expected)

        expected = MultiIndex.from_tuples([('a', 'b', 'c'), ('d', 'e', np.nan),
                                           ('f', np.nan, np.nan)])
        tm.assert_index_equal(idx.str.split(expand=True), expected)

        # test boolean case, should return np.array instead of boolean Index
        idx = Index(['a1', 'a2', 'b1', 'b2'])
        expected = np.array([True, True, False, False])
        tm.assert_numpy_array_equal(idx.str.startswith('a'), expected)
        self.assertIsInstance(idx.str.startswith('a'), np.ndarray)
        s = Series(range(4), index=idx)
        expected = Series(range(2), index=['a1', 'a2'])
        tm.assert_series_equal(s[s.index.str.startswith('a')], expected)

    def test_tab_completion(self):
        # GH 9910
        idx = Index(list('abcd'))
        self.assertTrue('str' in dir(idx))

        idx = Index(range(4))
        self.assertTrue('str' not in dir(idx))

    def test_indexing_doesnt_change_class(self):
        idx = Index([1, 2, 3, 'a', 'b', 'c'])

        self.assertTrue(idx[1:3].identical(pd.Index([2, 3], dtype=np.object_)))
        self.assertTrue(idx[[0, 1]].identical(pd.Index(
            [1, 2], dtype=np.object_)))

    def test_outer_join_sort(self):
        left_idx = Index(np.random.permutation(15))
        right_idx = tm.makeDateIndex(10)

        with tm.assert_produces_warning(RuntimeWarning):
            joined = left_idx.join(right_idx, how='outer')

        # right_idx in this case because DatetimeIndex has join precedence over
        # Int64Index
        with tm.assert_produces_warning(RuntimeWarning):
            expected = right_idx.astype(object).union(left_idx.astype(object))
        tm.assert_index_equal(joined, expected)

    def test_nan_first_take_datetime(self):
        idx = Index([pd.NaT, Timestamp('20130101'), Timestamp('20130102')])
        res = idx.take([-1, 0, 1])
        exp = Index([idx[-1], idx[0], idx[1]])
        tm.assert_index_equal(res, exp)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.Index(list('ABC'), name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.Index(list('BAC'), name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.Index(['B', 'A', np.nan], name='xxx')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.Index(['B', 'A', 'C'], name='xxx')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_reindex_preserves_name_if_target_is_list_or_ndarray(self):
        # GH6552
        idx = pd.Index([0, 1, 2])

        dt_idx = pd.date_range('20130101', periods=3)

        idx.name = None
        self.assertEqual(idx.reindex([])[0].name, None)
        self.assertEqual(idx.reindex(np.array([]))[0].name, None)
        self.assertEqual(idx.reindex(idx.tolist())[0].name, None)
        self.assertEqual(idx.reindex(idx.tolist()[:-1])[0].name, None)
        self.assertEqual(idx.reindex(idx.values)[0].name, None)
        self.assertEqual(idx.reindex(idx.values[:-1])[0].name, None)

        # Must preserve name even if dtype changes.
        self.assertEqual(idx.reindex(dt_idx.values)[0].name, None)
        self.assertEqual(idx.reindex(dt_idx.tolist())[0].name, None)

        idx.name = 'foobar'
        self.assertEqual(idx.reindex([])[0].name, 'foobar')
        self.assertEqual(idx.reindex(np.array([]))[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.tolist())[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.tolist()[:-1])[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.values)[0].name, 'foobar')
        self.assertEqual(idx.reindex(idx.values[:-1])[0].name, 'foobar')

        # Must preserve name even if dtype changes.
        self.assertEqual(idx.reindex(dt_idx.values)[0].name, 'foobar')
        self.assertEqual(idx.reindex(dt_idx.tolist())[0].name, 'foobar')

    def test_reindex_preserves_type_if_target_is_empty_list_or_array(self):
        # GH7774
        idx = pd.Index(list('abc'))

        def get_reindex_type(target):
            return idx.reindex(target)[0].dtype.type

        self.assertEqual(get_reindex_type([]), np.object_)
        self.assertEqual(get_reindex_type(np.array([])), np.object_)
        self.assertEqual(get_reindex_type(np.array([], dtype=np.int64)),
                         np.object_)

    def test_reindex_doesnt_preserve_type_if_target_is_empty_index(self):
        # GH7774
        idx = pd.Index(list('abc'))

        def get_reindex_type(target):
            return idx.reindex(target)[0].dtype.type

        self.assertEqual(get_reindex_type(pd.Int64Index([])), np.int64)
        self.assertEqual(get_reindex_type(pd.Float64Index([])), np.float64)
        self.assertEqual(get_reindex_type(pd.DatetimeIndex([])), np.datetime64)

        reindexed = idx.reindex(pd.MultiIndex(
            [pd.Int64Index([]), pd.Float64Index([])], [[], []]))[0]
        self.assertEqual(reindexed.levels[0].dtype.type, np.int64)
        self.assertEqual(reindexed.levels[1].dtype.type, np.float64)

    def test_groupby(self):
        idx = Index(range(5))
        groups = idx.groupby(np.array([1, 1, 2, 2, 2]))
        exp = {1: [0, 1], 2: [2, 3, 4]}
        tm.assert_dict_equal(groups, exp)

    def test_equals_op_multiindex(self):
        # GH9785
        # test comparisons of multiindex
        from pandas.compat import StringIO
        df = pd.read_csv(StringIO('a,b,c\n1,2,3\n4,5,6'), index_col=[0, 1])
        tm.assert_numpy_array_equal(df.index == df.index,
                                    np.array([True, True]))

        mi1 = MultiIndex.from_tuples([(1, 2), (4, 5)])
        tm.assert_numpy_array_equal(df.index == mi1, np.array([True, True]))
        mi2 = MultiIndex.from_tuples([(1, 2), (4, 6)])
        tm.assert_numpy_array_equal(df.index == mi2, np.array([True, False]))
        mi3 = MultiIndex.from_tuples([(1, 2), (4, 5), (8, 9)])
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            df.index == mi3

        index_a = Index(['foo', 'bar', 'baz'])
        with tm.assertRaisesRegexp(ValueError, "Lengths must match"):
            df.index == index_a
        tm.assert_numpy_array_equal(index_a == mi3,
                                    np.array([False, False, False]))

    def test_conversion_preserves_name(self):
        # GH 10875
        i = pd.Index(['01:02:03', '01:02:04'], name='label')
        self.assertEqual(i.name, pd.to_datetime(i).name)
        self.assertEqual(i.name, pd.to_timedelta(i).name)

    def test_string_index_repr(self):
        # py3/py2 repr can differ because of "u" prefix
        # which also affects to displayed element size

        # suppress flake8 warnings
        if PY3:
            coerce = lambda x: x
        else:
            coerce = unicode

        # short
        idx = pd.Index(['a', 'bb', 'ccc'])
        if PY3:
            expected = u"""Index(['a', 'bb', 'ccc'], dtype='object')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'a', u'bb', u'ccc'], dtype='object')"""
            self.assertEqual(coerce(idx), expected)

        # multiple lines
        idx = pd.Index(['a', 'bb', 'ccc'] * 10)
        if PY3:
            expected = u"""\
Index(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc',
       'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc',
       'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
      dtype='object')"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""\
Index([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
       u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb',
       u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
      dtype='object')"""

            self.assertEqual(coerce(idx), expected)

        # truncated
        idx = pd.Index(['a', 'bb', 'ccc'] * 100)
        if PY3:
            expected = u"""\
Index(['a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a',
       ...
       'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc', 'a', 'bb', 'ccc'],
      dtype='object', length=300)"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""\
Index([u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a',
       ...
       u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc', u'a', u'bb', u'ccc'],
      dtype='object', length=300)"""

            self.assertEqual(coerce(idx), expected)

        # short
        idx = pd.Index([u'あ', u'いい', u'ううう'])
        if PY3:
            expected = u"""Index(['あ', 'いい', 'ううう'], dtype='object')"""
            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""\
Index([u'あ', u'いい', u'ううう'], dtype='object')"""
            self.assertEqual(coerce(idx), expected)

        # multiple lines
        idx = pd.Index([u'あ', u'いい', u'ううう'] * 10)
        if PY3:
            expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
      dtype='object')"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
      dtype='object')"""

            self.assertEqual(coerce(idx), expected)

        # truncated
        idx = pd.Index([u'あ', u'いい', u'ううう'] * 100)
        if PY3:
            expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ',
       ...
       'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう'],
      dtype='object', length=300)"""

            self.assertEqual(repr(idx), expected)
        else:
            expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       ...
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
      dtype='object', length=300)"""

            self.assertEqual(coerce(idx), expected)

        # Emable Unicode option -----------------------------------------
        with cf.option_context('display.unicode.east_asian_width', True):

            # short
            idx = pd.Index([u'あ', u'いい', u'ううう'])
            if PY3:
                expected = u"""Index(['あ', 'いい', 'ううう'], dtype='object')"""
                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""Index([u'あ', u'いい', u'ううう'], dtype='object')"""
                self.assertEqual(coerce(idx), expected)

            # multiple lines
            idx = pd.Index([u'あ', u'いい', u'ううう'] * 10)
            if PY3:
                expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ', 'いい', 'ううう'],
      dtype='object')"""

                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう'],
      dtype='object')"""

                self.assertEqual(coerce(idx), expected)

            # truncated
            idx = pd.Index([u'あ', u'いい', u'ううう'] * 100)
            if PY3:
                expected = u"""Index(['あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう',
       'あ',
       ...
       'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい', 'ううう', 'あ', 'いい',
       'ううう'],
      dtype='object', length=300)"""

                self.assertEqual(repr(idx), expected)
            else:
                expected = u"""Index([u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい',
       u'ううう', u'あ',
       ...
       u'ううう', u'あ', u'いい', u'ううう', u'あ', u'いい', u'ううう', u'あ',
       u'いい', u'ううう'],
      dtype='object', length=300)"""

                self.assertEqual(coerce(idx), expected)


def test_get_combined_index():
    from pandas.core.index import _get_combined_index
    result = _get_combined_index([])
    assert (result.equals(Index([])))
