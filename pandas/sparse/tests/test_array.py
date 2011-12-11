from numpy import nan, ndarray
import numpy as np

import operator
import unittest

from pandas.sparse.api import SparseArray
from pandas.util.testing import assert_almost_equal

def assert_sp_array_equal(left, right):
    assert_almost_equal(left.sp_values, right.sp_values)
    assert(left.sp_index.equals(right.sp_index))
    if np.isnan(left.fill_value):
        assert(np.isnan(right.fill_value))
    else:
        assert(left.fill_value == right.fill_value)


class TestSparseArray(unittest.TestCase):

    def setUp(self):
        self.arr_data = np.array([nan, nan, 1, 2, 3, nan, 4, 5, nan, 6])
        self.arr = SparseArray(self.arr_data)
        self.zarr = SparseArray([0, 0, 1, 2, 3, 0, 4, 5, 0, 6], fill_value=0)

    def test_values_asarray(self):
        assert_almost_equal(self.arr.values, self.arr_data)
        assert_almost_equal(self.arr.sp_values, np.asarray(self.arr))

    def test_getslice(self):
        result = self.arr[:-3]
        exp = SparseArray(self.arr.values[:-3])
        assert_sp_array_equal(result, exp)

    def test_binary_operators(self):
        data1 = np.random.randn(20)
        data2 = np.random.randn(20)
        data1[::2] = np.nan
        data2[::3] = np.nan

        arr1 = SparseArray(data1)
        arr2 = SparseArray(data2)

        def _check_op(op):
            res = op(arr1, arr2)
            exp = SparseArray(op(arr1.values, arr2.values))
            self.assert_(isinstance(res, SparseArray))
            assert_almost_equal(res.values, exp.values)

        def _check_inplace_op(op):
            tmp = arr1.copy()
            self.assertRaises(NotImplementedError, op, tmp, arr2)

        bin_ops = [operator.add, operator.sub, operator.mul, operator.truediv,
                   operator.floordiv, operator.pow]
        for op in bin_ops:
            _check_op(op)

        inplace_ops = ['iadd', 'isub', 'imul', 'itruediv', 'ifloordiv', 'ipow']
        for op in inplace_ops:
            _check_inplace_op(getattr(operator, op))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
