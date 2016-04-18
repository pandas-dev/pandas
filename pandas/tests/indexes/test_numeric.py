# -*- coding: utf-8 -*-

from datetime import datetime
from pandas import compat
from pandas.compat import range, lrange, u, PY3

import numpy as np

from pandas import (date_range, Series, DataFrame,
                    Index, Float64Index, Int64Index, RangeIndex)
from pandas.util.testing import assertRaisesRegexp

import pandas.util.testing as tm
import pandas.core.config as cf

import pandas as pd
from pandas.lib import Timestamp

from .common import Base


class Numeric(Base):

    def test_numeric_compat(self):

        idx = self.create_index()
        didx = idx * idx

        result = idx * 1
        tm.assert_index_equal(result, idx)

        result = 1 * idx
        tm.assert_index_equal(result, idx)

        # in general not true for RangeIndex
        if not isinstance(idx, RangeIndex):
            result = idx * idx
            tm.assert_index_equal(result, idx ** 2)

        # truediv under PY3
        result = idx / 1
        expected = idx
        if PY3:
            expected = expected.astype('float64')
        tm.assert_index_equal(result, expected)

        result = idx / 2
        if PY3:
            expected = expected.astype('float64')
        expected = Index(idx.values / 2)
        tm.assert_index_equal(result, expected)

        result = idx // 1
        tm.assert_index_equal(result, idx)

        result = idx * np.array(5, dtype='int64')
        tm.assert_index_equal(result, idx * 5)

        result = idx * np.arange(5, dtype='int64')
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5, dtype='int64'))
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5, dtype='float64') + 0.1)
        expected = Float64Index(np.arange(5, dtype='float64') *
                                (np.arange(5, dtype='float64') + 0.1))
        tm.assert_index_equal(result, expected)

        # invalid
        self.assertRaises(TypeError,
                          lambda: idx * date_range('20130101', periods=5))
        self.assertRaises(ValueError, lambda: idx * idx[0:3])
        self.assertRaises(ValueError, lambda: idx * np.array([1, 2]))

    def test_explicit_conversions(self):

        # GH 8608
        # add/sub are overriden explicity for Float/Int Index
        idx = self._holder(np.arange(5, dtype='int64'))

        # float conversions
        arr = np.arange(5, dtype='int64') * 3.2
        expected = Float64Index(arr)
        fidx = idx * 3.2
        tm.assert_index_equal(fidx, expected)
        fidx = 3.2 * idx
        tm.assert_index_equal(fidx, expected)

        # interops with numpy arrays
        expected = Float64Index(arr)
        a = np.zeros(5, dtype='float64')
        result = fidx - a
        tm.assert_index_equal(result, expected)

        expected = Float64Index(-arr)
        a = np.zeros(5, dtype='float64')
        result = a - fidx
        tm.assert_index_equal(result, expected)

    def test_ufunc_compat(self):
        idx = self._holder(np.arange(5, dtype='int64'))
        result = np.sin(idx)
        expected = Float64Index(np.sin(np.arange(5, dtype='int64')))
        tm.assert_index_equal(result, expected)

    def test_index_groupby(self):
        int_idx = Index(range(6))
        float_idx = Index(np.arange(0, 0.6, 0.1))
        obj_idx = Index('A B C D E F'.split())
        dt_idx = pd.date_range('2013-01-01', freq='M', periods=6)

        for idx in [int_idx, float_idx, obj_idx, dt_idx]:
            to_groupby = np.array([1, 2, np.nan, np.nan, 2, 1])
            self.assertEqual(idx.groupby(to_groupby),
                             {1.0: [idx[0], idx[5]], 2.0: [idx[1], idx[4]]})

            to_groupby = Index([datetime(2011, 11, 1),
                                datetime(2011, 12, 1),
                                pd.NaT,
                                pd.NaT,
                                datetime(2011, 12, 1),
                                datetime(2011, 11, 1)],
                               tz='UTC').values

            ex_keys = pd.tslib.datetime_to_datetime64(np.array([Timestamp(
                '2011-11-01'), Timestamp('2011-12-01')]))
            expected = {ex_keys[0][0]: [idx[0], idx[5]],
                        ex_keys[0][1]: [idx[1], idx[4]]}
            self.assertEqual(idx.groupby(to_groupby), expected)

    def test_modulo(self):
        # GH 9244
        index = self.create_index()
        expected = Index(index.values % 2)
        self.assert_index_equal(index % 2, expected)


class TestFloat64Index(Numeric, tm.TestCase):
    _holder = Float64Index
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(mixed=Float64Index([1.5, 2, 3, 4, 5]),
                            float=Float64Index(np.arange(5) * 2.5))
        self.setup_indices()

    def create_index(self):
        return Float64Index(np.arange(5, dtype='float64'))

    def test_repr_roundtrip(self):
        for ind in (self.mixed, self.float):
            tm.assert_index_equal(eval(repr(ind)), ind)

    def check_is_index(self, i):
        self.assertIsInstance(i, Index)
        self.assertNotIsInstance(i, Float64Index)

    def check_coerce(self, a, b, is_float_index=True):
        self.assertTrue(a.equals(b))
        if is_float_index:
            self.assertIsInstance(b, Float64Index)
        else:
            self.check_is_index(b)

    def test_constructor(self):

        # explicit construction
        index = Float64Index([1, 2, 3, 4, 5])
        self.assertIsInstance(index, Float64Index)
        self.assertTrue((index.values == np.array(
            [1, 2, 3, 4, 5], dtype='float64')).all())
        index = Float64Index(np.array([1, 2, 3, 4, 5]))
        self.assertIsInstance(index, Float64Index)
        index = Float64Index([1., 2, 3, 4, 5])
        self.assertIsInstance(index, Float64Index)
        index = Float64Index(np.array([1., 2, 3, 4, 5]))
        self.assertIsInstance(index, Float64Index)
        self.assertEqual(index.dtype, float)

        index = Float64Index(np.array([1., 2, 3, 4, 5]), dtype=np.float32)
        self.assertIsInstance(index, Float64Index)
        self.assertEqual(index.dtype, np.float64)

        index = Float64Index(np.array([1, 2, 3, 4, 5]), dtype=np.float32)
        self.assertIsInstance(index, Float64Index)
        self.assertEqual(index.dtype, np.float64)

        # nan handling
        result = Float64Index([np.nan, np.nan])
        self.assertTrue(pd.isnull(result.values).all())
        result = Float64Index(np.array([np.nan]))
        self.assertTrue(pd.isnull(result.values).all())
        result = Index(np.array([np.nan]))
        self.assertTrue(pd.isnull(result.values).all())

    def test_constructor_invalid(self):

        # invalid
        self.assertRaises(TypeError, Float64Index, 0.)
        self.assertRaises(TypeError, Float64Index, ['a', 'b', 0.])
        self.assertRaises(TypeError, Float64Index, [Timestamp('20130101')])

    def test_constructor_coerce(self):

        self.check_coerce(self.mixed, Index([1.5, 2, 3, 4, 5]))
        self.check_coerce(self.float, Index(np.arange(5) * 2.5))
        self.check_coerce(self.float, Index(np.array(
            np.arange(5) * 2.5, dtype=object)))

    def test_constructor_explicit(self):

        # these don't auto convert
        self.check_coerce(self.float,
                          Index((np.arange(5) * 2.5), dtype=object),
                          is_float_index=False)
        self.check_coerce(self.mixed, Index(
            [1.5, 2, 3, 4, 5], dtype=object), is_float_index=False)

    def test_astype(self):

        result = self.float.astype(object)
        self.assertTrue(result.equals(self.float))
        self.assertTrue(self.float.equals(result))
        self.check_is_index(result)

        i = self.mixed.copy()
        i.name = 'foo'
        result = i.astype(object)
        self.assertTrue(result.equals(i))
        self.assertTrue(i.equals(result))
        self.check_is_index(result)

        # GH 12881
        # a float astype int
        for dtype in ['int16', 'int32', 'int64']:
            i = Float64Index([0, 1, 2])
            result = i.astype(dtype)
            expected = Int64Index([0, 1, 2])
            tm.assert_index_equal(result, expected)

            i = Float64Index([0, 1.1, 2])
            result = i.astype(dtype)
            expected = Int64Index([0, 1, 2])
            tm.assert_index_equal(result, expected)

        for dtype in ['float32', 'float64']:
            i = Float64Index([0, 1, 2])
            result = i.astype(dtype)
            expected = i
            tm.assert_index_equal(result, expected)

            i = Float64Index([0, 1.1, 2])
            result = i.astype(dtype)
            expected = Index(i.values.astype(dtype))
            tm.assert_index_equal(result, expected)

        # invalid
        for dtype in ['M8[ns]', 'm8[ns]']:
            self.assertRaises(TypeError, lambda: i.astype(dtype))

    def test_equals(self):

        i = Float64Index([1.0, 2.0])
        self.assertTrue(i.equals(i))
        self.assertTrue(i.identical(i))

        i2 = Float64Index([1.0, 2.0])
        self.assertTrue(i.equals(i2))

        i = Float64Index([1.0, np.nan])
        self.assertTrue(i.equals(i))
        self.assertTrue(i.identical(i))

        i2 = Float64Index([1.0, np.nan])
        self.assertTrue(i.equals(i2))

    def test_get_indexer(self):
        idx = Float64Index([0.0, 1.0, 2.0])
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = [-0.1, 0.5, 1.1]
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest'), [0, 1, 1])

    def test_get_loc(self):
        idx = Float64Index([0.0, 1.0, 2.0])
        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(1, method), 1)
            if method is not None:
                self.assertEqual(idx.get_loc(1, method, tolerance=0), 1)

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc(1.1, method), loc)
            self.assertEqual(idx.get_loc(1.1, method, tolerance=0.9), loc)

        self.assertRaises(KeyError, idx.get_loc, 'foo')
        self.assertRaises(KeyError, idx.get_loc, 1.5)
        self.assertRaises(KeyError, idx.get_loc, 1.5, method='pad',
                          tolerance=0.1)

        with tm.assertRaisesRegexp(ValueError, 'must be numeric'):
            idx.get_loc(1.4, method='nearest', tolerance='foo')

    def test_get_loc_na(self):
        idx = Float64Index([np.nan, 1, 2])
        self.assertEqual(idx.get_loc(1), 1)
        self.assertEqual(idx.get_loc(np.nan), 0)

        idx = Float64Index([np.nan, 1, np.nan])
        self.assertEqual(idx.get_loc(1), 1)

        # representable by slice [0:2:2]
        # self.assertRaises(KeyError, idx.slice_locs, np.nan)
        sliced = idx.slice_locs(np.nan)
        self.assertTrue(isinstance(sliced, tuple))
        self.assertEqual(sliced, (0, 3))

        # not representable by slice
        idx = Float64Index([np.nan, 1, np.nan, np.nan])
        self.assertEqual(idx.get_loc(1), 1)
        self.assertRaises(KeyError, idx.slice_locs, np.nan)

    def test_contains_nans(self):
        i = Float64Index([1.0, 2.0, np.nan])
        self.assertTrue(np.nan in i)

    def test_contains_not_nans(self):
        i = Float64Index([1.0, 2.0, np.nan])
        self.assertTrue(1.0 in i)

    def test_doesnt_contain_all_the_things(self):
        i = Float64Index([np.nan])
        self.assertFalse(i.isin([0]).item())
        self.assertFalse(i.isin([1]).item())
        self.assertTrue(i.isin([np.nan]).item())

    def test_nan_multiple_containment(self):
        i = Float64Index([1.0, np.nan])
        tm.assert_numpy_array_equal(i.isin([1.0]), np.array([True, False]))
        tm.assert_numpy_array_equal(i.isin([2.0, np.pi]),
                                    np.array([False, False]))
        tm.assert_numpy_array_equal(i.isin([np.nan]), np.array([False, True]))
        tm.assert_numpy_array_equal(i.isin([1.0, np.nan]),
                                    np.array([True, True]))
        i = Float64Index([1.0, 2.0])
        tm.assert_numpy_array_equal(i.isin([np.nan]), np.array([False, False]))

    def test_astype_from_object(self):
        index = Index([1.0, np.nan, 0.2], dtype='object')
        result = index.astype(float)
        expected = Float64Index([1.0, np.nan, 0.2])
        tm.assert_equal(result.dtype, expected.dtype)
        tm.assert_index_equal(result, expected)

    def test_fillna_float64(self):
        # GH 11343
        idx = Index([1.0, np.nan, 3.0], dtype=float, name='x')
        # can't downcast
        exp = Index([1.0, 0.1, 3.0], name='x')
        self.assert_index_equal(idx.fillna(0.1), exp)

        # downcast
        exp = Float64Index([1.0, 2.0, 3.0], name='x')
        self.assert_index_equal(idx.fillna(2), exp)

        # object
        exp = Index([1.0, 'obj', 3.0], name='x')
        self.assert_index_equal(idx.fillna('obj'), exp)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.Float64Index([1., 2., 3.], name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.Float64Index([2., 1., 3.], name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.Float64Index([2., 1., np.nan], name='xxx')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.Float64Index([2., 1., 3.], name='xxx')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))


class TestInt64Index(Numeric, tm.TestCase):
    _holder = Int64Index
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index=Int64Index(np.arange(0, 20, 2)))
        self.setup_indices()

    def create_index(self):
        return Int64Index(np.arange(5, dtype='int64'))

    def test_too_many_names(self):
        def testit():
            self.index.names = ["roger", "harold"]

        assertRaisesRegexp(ValueError, "^Length", testit)

    def test_constructor(self):
        # pass list, coerce fine
        index = Int64Index([-5, 0, 1, 2])
        expected = np.array([-5, 0, 1, 2], dtype=np.int64)
        tm.assert_numpy_array_equal(index, expected)

        # from iterable
        index = Int64Index(iter([-5, 0, 1, 2]))
        tm.assert_numpy_array_equal(index, expected)

        # scalar raise Exception
        self.assertRaises(TypeError, Int64Index, 5)

        # copy
        arr = self.index.values
        new_index = Int64Index(arr, copy=True)
        tm.assert_numpy_array_equal(new_index, self.index)
        val = arr[0] + 3000

        # this should not change index
        arr[0] = val
        self.assertNotEqual(new_index[0], val)

        # interpret list-like
        expected = Int64Index([5, 0])
        for cls in [Index, Int64Index]:
            for idx in [cls([5, 0], dtype='int64'),
                        cls(np.array([5, 0]), dtype='int64'),
                        cls(Series([5, 0]), dtype='int64')]:
                tm.assert_index_equal(idx, expected)

    def test_constructor_corner(self):
        arr = np.array([1, 2, 3, 4], dtype=object)
        index = Int64Index(arr)
        self.assertEqual(index.values.dtype, np.int64)
        self.assertTrue(index.equals(arr))

        # preventing casting
        arr = np.array([1, '2', 3, '4'], dtype=object)
        with tm.assertRaisesRegexp(TypeError, 'casting'):
            Int64Index(arr)

        arr_with_floats = [0, 2, 3, 4, 5, 1.25, 3, -1]
        with tm.assertRaisesRegexp(TypeError, 'casting'):
            Int64Index(arr_with_floats)

    def test_copy(self):
        i = Int64Index([], name='Foo')
        i_copy = i.copy()
        self.assertEqual(i_copy.name, 'Foo')

    def test_view(self):
        super(TestInt64Index, self).test_view()

        i = Int64Index([], name='Foo')
        i_view = i.view()
        self.assertEqual(i_view.name, 'Foo')

        i_view = i.view('i8')
        tm.assert_index_equal(i, Int64Index(i_view, name='Foo'))

        i_view = i.view(Int64Index)
        tm.assert_index_equal(i, Int64Index(i_view, name='Foo'))

    def test_coerce_list(self):
        # coerce things
        arr = Index([1, 2, 3, 4])
        tm.assertIsInstance(arr, Int64Index)

        # but not if explicit dtype passed
        arr = Index([1, 2, 3, 4], dtype=object)
        tm.assertIsInstance(arr, Index)

    def test_dtype(self):
        self.assertEqual(self.index.dtype, np.int64)

    def test_is_monotonic(self):
        self.assertTrue(self.index.is_monotonic)
        self.assertTrue(self.index.is_monotonic_increasing)
        self.assertFalse(self.index.is_monotonic_decreasing)

        index = Int64Index([4, 3, 2, 1])
        self.assertFalse(index.is_monotonic)
        self.assertTrue(index.is_monotonic_decreasing)

        index = Int64Index([1])
        self.assertTrue(index.is_monotonic)
        self.assertTrue(index.is_monotonic_increasing)
        self.assertTrue(index.is_monotonic_decreasing)

    def test_is_monotonic_na(self):
        examples = [Index([np.nan]),
                    Index([np.nan, 1]),
                    Index([1, 2, np.nan]),
                    Index(['a', 'b', np.nan]),
                    pd.to_datetime(['NaT']),
                    pd.to_datetime(['NaT', '2000-01-01']),
                    pd.to_datetime(['2000-01-01', 'NaT', '2000-01-02']),
                    pd.to_timedelta(['1 day', 'NaT']), ]
        for index in examples:
            self.assertFalse(index.is_monotonic_increasing)
            self.assertFalse(index.is_monotonic_decreasing)

    def test_equals(self):
        same_values = Index(self.index, dtype=object)
        self.assertTrue(self.index.equals(same_values))
        self.assertTrue(same_values.equals(self.index))

    def test_logical_compat(self):
        idx = self.create_index()
        self.assertEqual(idx.all(), idx.values.all())
        self.assertEqual(idx.any(), idx.values.any())

    def test_identical(self):
        i = Index(self.index.copy())
        self.assertTrue(i.identical(self.index))

        same_values_different_type = Index(i, dtype=object)
        self.assertFalse(i.identical(same_values_different_type))

        i = self.index.copy(dtype=object)
        i = i.rename('foo')
        same_values = Index(i, dtype=object)
        self.assertTrue(same_values.identical(i))

        self.assertFalse(i.identical(self.index))
        self.assertTrue(Index(same_values, name='foo', dtype=object).identical(
            i))

        self.assertFalse(self.index.copy(dtype=object)
                         .identical(self.index.copy(dtype='int64')))

    def test_get_indexer(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target)
        expected = np.array([0, -1, 1, -1, 2, -1, 3, -1, 4, -1])
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_pad(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target, method='pad')
        expected = np.array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
        tm.assert_numpy_array_equal(indexer, expected)

    def test_get_indexer_backfill(self):
        target = Int64Index(np.arange(10))
        indexer = self.index.get_indexer(target, method='backfill')
        expected = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5])
        tm.assert_numpy_array_equal(indexer, expected)

    def test_join_outer(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        # guarantee of sortedness
        res, lidx, ridx = self.index.join(other, how='outer',
                                          return_indexers=True)
        noidx_res = self.index.join(other, how='outer')
        self.assertTrue(res.equals(noidx_res))

        eres = Int64Index([0, 1, 2, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 25])
        elidx = np.array([0, -1, 1, 2, -1, 3, -1, 4, 5, 6, 7, 8, 9, -1],
                         dtype=np.int64)
        eridx = np.array([-1, 3, 4, -1, 5, -1, 0, -1, -1, 1, -1, -1, -1, 2],
                         dtype=np.int64)

        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='outer',
                                          return_indexers=True)
        noidx_res = self.index.join(other_mono, how='outer')
        self.assertTrue(res.equals(noidx_res))

        eridx = np.array([-1, 0, 1, -1, 2, -1, 3, -1, -1, 4, -1, -1, -1, 5],
                         dtype=np.int64)
        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_inner(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='inner',
                                          return_indexers=True)

        # no guarantee of sortedness, so sort for comparison purposes
        ind = res.argsort()
        res = res.take(ind)
        lidx = lidx.take(ind)
        ridx = ridx.take(ind)

        eres = Int64Index([2, 12])
        elidx = np.array([1, 6])
        eridx = np.array([4, 1])

        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='inner',
                                          return_indexers=True)

        res2 = self.index.intersection(other_mono)
        self.assertTrue(res.equals(res2))

        eridx = np.array([1, 4])
        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_left(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='left',
                                          return_indexers=True)
        eres = self.index
        eridx = np.array([-1, 4, -1, -1, -1, -1, 1, -1, -1, -1],
                         dtype=np.int64)

        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        self.assertIsNone(lidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='left',
                                          return_indexers=True)
        eridx = np.array([-1, 1, -1, -1, -1, -1, 4, -1, -1, -1],
                         dtype=np.int64)
        tm.assertIsInstance(res, Int64Index)
        self.assertTrue(res.equals(eres))
        self.assertIsNone(lidx)
        tm.assert_numpy_array_equal(ridx, eridx)

        # non-unique
        idx = Index([1, 1, 2, 5])
        idx2 = Index([1, 2, 5, 7, 9])
        res, lidx, ridx = idx2.join(idx, how='left', return_indexers=True)
        eres = Index([1, 1, 2, 5, 7, 9])  # 1 is in idx2, so it should be x2
        eridx = np.array([0, 1, 2, 3, -1, -1])
        elidx = np.array([0, 0, 1, 2, 3, 4])
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_right(self):
        other = Int64Index([7, 12, 25, 1, 2, 5])
        other_mono = Int64Index([1, 2, 5, 7, 12, 25])

        # not monotonic
        res, lidx, ridx = self.index.join(other, how='right',
                                          return_indexers=True)
        eres = other
        elidx = np.array([-1, 6, -1, -1, 1, -1], dtype=np.int64)

        tm.assertIsInstance(other, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        self.assertIsNone(ridx)

        # monotonic
        res, lidx, ridx = self.index.join(other_mono, how='right',
                                          return_indexers=True)
        eres = other_mono
        elidx = np.array([-1, 1, -1, -1, 6, -1], dtype=np.int64)
        tm.assertIsInstance(other, Int64Index)
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        self.assertIsNone(ridx)

        # non-unique
        idx = Index([1, 1, 2, 5])
        idx2 = Index([1, 2, 5, 7, 9])
        res, lidx, ridx = idx.join(idx2, how='right', return_indexers=True)
        eres = Index([1, 1, 2, 5, 7, 9])  # 1 is in idx2, so it should be x2
        elidx = np.array([0, 1, 2, 3, -1, -1])
        eridx = np.array([0, 0, 1, 2, 3, 4])
        self.assertTrue(res.equals(eres))
        tm.assert_numpy_array_equal(lidx, elidx)
        tm.assert_numpy_array_equal(ridx, eridx)

    def test_join_non_int_index(self):
        other = Index([3, 6, 7, 8, 10], dtype=object)

        outer = self.index.join(other, how='outer')
        outer2 = other.join(self.index, how='outer')
        expected = Index([0, 2, 3, 4, 6, 7, 8, 10, 12, 14,
                          16, 18], dtype=object)
        self.assertTrue(outer.equals(outer2))
        self.assertTrue(outer.equals(expected))

        inner = self.index.join(other, how='inner')
        inner2 = other.join(self.index, how='inner')
        expected = Index([6, 8, 10], dtype=object)
        self.assertTrue(inner.equals(inner2))
        self.assertTrue(inner.equals(expected))

        left = self.index.join(other, how='left')
        self.assertTrue(left.equals(self.index))

        left2 = other.join(self.index, how='left')
        self.assertTrue(left2.equals(other))

        right = self.index.join(other, how='right')
        self.assertTrue(right.equals(other))

        right2 = other.join(self.index, how='right')
        self.assertTrue(right2.equals(self.index))

    def test_join_non_unique(self):
        left = Index([4, 4, 3, 3])

        joined, lidx, ridx = left.join(left, return_indexers=True)

        exp_joined = Index([3, 3, 3, 3, 4, 4, 4, 4])
        self.assertTrue(joined.equals(exp_joined))

        exp_lidx = np.array([2, 2, 3, 3, 0, 0, 1, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(lidx, exp_lidx)

        exp_ridx = np.array([2, 3, 2, 3, 0, 1, 0, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(ridx, exp_ridx)

    def test_join_self(self):
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            joined = self.index.join(self.index, how=kind)
            self.assertIs(self.index, joined)

    def test_intersection(self):
        other = Index([1, 2, 3, 4, 5])
        result = self.index.intersection(other)
        expected = np.sort(np.intersect1d(self.index.values, other.values))
        tm.assert_numpy_array_equal(result, expected)

        result = other.intersection(self.index)
        expected = np.sort(np.asarray(np.intersect1d(self.index.values,
                                                     other.values)))
        tm.assert_numpy_array_equal(result, expected)

    def test_intersect_str_dates(self):
        dt_dates = [datetime(2012, 2, 9), datetime(2012, 2, 22)]

        i1 = Index(dt_dates, dtype=object)
        i2 = Index(['aa'], dtype=object)
        res = i2.intersection(i1)

        self.assertEqual(len(res), 0)

    def test_union_noncomparable(self):
        from datetime import datetime, timedelta
        # corner case, non-Int64Index
        now = datetime.now()
        other = Index([now + timedelta(i) for i in range(4)], dtype=object)
        result = self.index.union(other)
        expected = np.concatenate((self.index, other))
        tm.assert_numpy_array_equal(result, expected)

        result = other.union(self.index)
        expected = np.concatenate((other, self.index))
        tm.assert_numpy_array_equal(result, expected)

    def test_cant_or_shouldnt_cast(self):
        # can't
        data = ['foo', 'bar', 'baz']
        self.assertRaises(TypeError, Int64Index, data)

        # shouldn't
        data = ['0', '1', '2']
        self.assertRaises(TypeError, Int64Index, data)

    def test_view_Index(self):
        self.index.view(Index)

    def test_prevent_casting(self):
        result = self.index.astype('O')
        self.assertEqual(result.dtype, np.object_)

    def test_take_preserve_name(self):
        index = Int64Index([1, 2, 3, 4], name='foo')
        taken = index.take([3, 0, 1])
        self.assertEqual(index.name, taken.name)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.Int64Index([1, 2, 3], name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.Int64Index([2, 1, 3], name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        msg = "Unable to fill values because Int64Index cannot contain NA"
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -1]), fill_value=True)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.Int64Index([2, 1, 3], name='xxx')
        tm.assert_index_equal(result, expected)

        msg = "Unable to fill values because Int64Index cannot contain NA"
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_int_name_format(self):
        index = Index(['a', 'b', 'c'], name=0)
        s = Series(lrange(3), index)
        df = DataFrame(lrange(3), index=index)
        repr(s)
        repr(df)

    def test_print_unicode_columns(self):
        df = pd.DataFrame({u("\u05d0"): [1, 2, 3],
                           "\u05d1": [4, 5, 6],
                           "c": [7, 8, 9]})
        repr(df.columns)  # should not raise UnicodeDecodeError

    def test_repr_summary(self):
        with cf.option_context('display.max_seq_items', 10):
            r = repr(pd.Index(np.arange(1000)))
            self.assertTrue(len(r) < 200)
            self.assertTrue("..." in r)

    def test_repr_roundtrip(self):
        tm.assert_index_equal(eval(repr(self.index)), self.index)

    def test_unicode_string_with_unicode(self):
        idx = Index(lrange(1000))

        if PY3:
            str(idx)
        else:
            compat.text_type(idx)

    def test_bytestring_with_unicode(self):
        idx = Index(lrange(1000))
        if PY3:
            bytes(idx)
        else:
            str(idx)

    def test_slice_keep_name(self):
        idx = Int64Index([1, 2], name='asdf')
        self.assertEqual(idx.name, idx[1:].name)

    def test_ufunc_coercions(self):
        idx = Int64Index([1, 2, 3, 4, 5], name='x')

        result = np.sqrt(idx)
        tm.assertIsInstance(result, Float64Index)
        exp = Float64Index(np.sqrt(np.array([1, 2, 3, 4, 5])), name='x')
        tm.assert_index_equal(result, exp)

        result = np.divide(idx, 2.)
        tm.assertIsInstance(result, Float64Index)
        exp = Float64Index([0.5, 1., 1.5, 2., 2.5], name='x')
        tm.assert_index_equal(result, exp)

        # _evaluate_numeric_binop
        result = idx + 2.
        tm.assertIsInstance(result, Float64Index)
        exp = Float64Index([3., 4., 5., 6., 7.], name='x')
        tm.assert_index_equal(result, exp)

        result = idx - 2.
        tm.assertIsInstance(result, Float64Index)
        exp = Float64Index([-1., 0., 1., 2., 3.], name='x')
        tm.assert_index_equal(result, exp)

        result = idx * 1.
        tm.assertIsInstance(result, Float64Index)
        exp = Float64Index([1., 2., 3., 4., 5.], name='x')
        tm.assert_index_equal(result, exp)

        result = idx / 2.
        tm.assertIsInstance(result, Float64Index)
        exp = Float64Index([0.5, 1., 1.5, 2., 2.5], name='x')
        tm.assert_index_equal(result, exp)
