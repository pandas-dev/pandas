from pandas.compat import range
import re
from numpy import nan, ndarray
import numpy as np

import operator

from pandas.core.series import Series
from pandas.core.common import notnull
from pandas.sparse.api import SparseArray
from pandas.util.testing import assert_almost_equal, assertRaisesRegexp
import pandas.util.testing as tm


def assert_sp_array_equal(left, right):
    assert_almost_equal(left.sp_values, right.sp_values)
    assert(left.sp_index.equals(right.sp_index))
    if np.isnan(left.fill_value):
        assert(np.isnan(right.fill_value))
    else:
        assert(left.fill_value == right.fill_value)


class TestSparseArray(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.arr_data = np.array([nan, nan, 1, 2, 3, nan, 4, 5, nan, 6])
        self.arr = SparseArray(self.arr_data)
        self.zarr = SparseArray([0, 0, 1, 2, 3, 0, 4, 5, 0, 6], fill_value=0)

    def test_get_item(self):
        errmsg = re.compile("bounds")
        assertRaisesRegexp(IndexError, errmsg, lambda: self.arr[11])
        assertRaisesRegexp(IndexError, errmsg, lambda: self.arr[-11])
        self.assertEqual(self.arr[-1], self.arr[len(self.arr) - 1])

    def test_bad_take(self):
        assertRaisesRegexp(IndexError, "bounds", lambda: self.arr.take(11))
        self.assertRaises(IndexError, lambda: self.arr.take(-11))

    def test_set_item(self):
        def setitem():
            self.arr[5] = 3

        def setslice():
            self.arr[1:5] = 2
        assertRaisesRegexp(TypeError, "item assignment", setitem)
        assertRaisesRegexp(TypeError, "item assignment", setslice)

    def test_constructor_from_sparse(self):
        res = SparseArray(self.zarr)
        self.assertEqual(res.fill_value, 0)
        assert_almost_equal(res.sp_values, self.zarr.sp_values)

    def test_constructor_copy(self):
        cp = SparseArray(self.arr, copy=True)
        cp.sp_values[:3] = 0
        self.assertFalse((self.arr.sp_values[:3] == 0).any())

        not_copy = SparseArray(self.arr)
        not_copy.sp_values[:3] = 0
        self.assertTrue((self.arr.sp_values[:3] == 0).all())

    def test_astype(self):
        res = self.arr.astype('f8')
        res.sp_values[:3] = 27
        self.assertFalse((self.arr.sp_values[:3] == 27).any())

        assertRaisesRegexp(TypeError, "floating point", self.arr.astype, 'i8')

    def test_copy_shallow(self):
        arr2 = self.arr.copy(deep=False)

        def _get_base(values):
            base = values.base
            while base.base is not None:
                base = base.base
            return base

        assert(_get_base(arr2) is _get_base(self.arr))

    def test_values_asarray(self):
        assert_almost_equal(self.arr.values, self.arr_data)
        assert_almost_equal(self.arr.to_dense(), self.arr_data)
        assert_almost_equal(self.arr.sp_values, np.asarray(self.arr))

    def test_getitem(self):
        def _checkit(i):
            assert_almost_equal(self.arr[i], self.arr.values[i])

        for i in range(len(self.arr)):
            _checkit(i)
            _checkit(-i)

    def test_getslice(self):
        result = self.arr[:-3]
        exp = SparseArray(self.arr.values[:-3])
        assert_sp_array_equal(result, exp)

        result = self.arr[-4:]
        exp = SparseArray(self.arr.values[-4:])
        assert_sp_array_equal(result, exp)

        # two corner cases from Series
        result = self.arr[-12:]
        exp = SparseArray(self.arr)
        assert_sp_array_equal(result, exp)

        result = self.arr[:-12]
        exp = SparseArray(self.arr.values[:0])
        assert_sp_array_equal(result, exp)

    def test_binary_operators(self):
        data1 = np.random.randn(20)
        data2 = np.random.randn(20)
        data1[::2] = np.nan
        data2[::3] = np.nan

        arr1 = SparseArray(data1)
        arr2 = SparseArray(data2)

        data1[::2] = 3
        data2[::3] = 3
        farr1 = SparseArray(data1, fill_value=3)
        farr2 = SparseArray(data2, fill_value=3)

        def _check_op(op, first, second):
            res = op(first, second)
            exp = SparseArray(op(first.values, second.values),
                              fill_value=first.fill_value)
            tm.assertIsInstance(res, SparseArray)
            assert_almost_equal(res.values, exp.values)

            res2 = op(first, second.values)
            tm.assertIsInstance(res2, SparseArray)
            assert_sp_array_equal(res, res2)

            res3 = op(first.values, second)
            tm.assertIsInstance(res3, SparseArray)
            assert_sp_array_equal(res, res3)

            res4 = op(first, 4)
            tm.assertIsInstance(res4, SparseArray)

            # ignore this if the actual op raises (e.g. pow)
            try:
                exp = op(first.values, 4)
                exp_fv = op(first.fill_value, 4)
                assert_almost_equal(res4.fill_value, exp_fv)
                assert_almost_equal(res4.values, exp)
            except (ValueError) :
                pass

        def _check_inplace_op(op):
            tmp = arr1.copy()
            self.assertRaises(NotImplementedError, op, tmp, arr2)

        bin_ops = [operator.add, operator.sub, operator.mul, operator.truediv,
                   operator.floordiv, operator.pow]
        for op in bin_ops:
            _check_op(op, arr1, arr2)
            _check_op(op, farr1, farr2)

        inplace_ops = ['iadd', 'isub', 'imul', 'itruediv', 'ifloordiv', 'ipow']
        for op in inplace_ops:
            _check_inplace_op(getattr(operator, op))

    def test_pickle(self):
        def _check_roundtrip(obj):
            unpickled = self.round_trip_pickle(obj)
            assert_sp_array_equal(unpickled, obj)

        _check_roundtrip(self.arr)
        _check_roundtrip(self.zarr)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
