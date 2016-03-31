from pandas.compat import range
import re
import operator
import warnings

from numpy import nan
import numpy as np

from pandas import _np_version_under1p8
from pandas.sparse.api import SparseArray
from pandas._sparse import IntIndex
from pandas.util.testing import assert_almost_equal, assertRaisesRegexp
import pandas.util.testing as tm


class TestSparseArray(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.arr_data = np.array([nan, nan, 1, 2, 3, nan, 4, 5, nan, 6])
        self.arr = SparseArray(self.arr_data)
        self.zarr = SparseArray([0, 0, 1, 2, 3, 0, 4, 5, 0, 6], fill_value=0)

    def test_constructor_dtype(self):
        arr = SparseArray([np.nan, 1, 2, np.nan])
        self.assertEqual(arr.dtype, np.float64)
        self.assertTrue(np.isnan(arr.fill_value))

        arr = SparseArray([np.nan, 1, 2, np.nan], fill_value=0)
        self.assertEqual(arr.dtype, np.float64)
        self.assertEqual(arr.fill_value, 0)

        arr = SparseArray([0, 1, 2, 4], dtype=np.int64)
        self.assertEqual(arr.dtype, np.int64)
        self.assertTrue(np.isnan(arr.fill_value))

        arr = SparseArray([0, 1, 2, 4], fill_value=0, dtype=np.int64)
        self.assertEqual(arr.dtype, np.int64)
        self.assertEqual(arr.fill_value, 0)

        arr = SparseArray([0, 1, 2, 4], dtype=None)
        self.assertEqual(arr.dtype, np.int64)
        self.assertTrue(np.isnan(arr.fill_value))

        arr = SparseArray([0, 1, 2, 4], fill_value=0, dtype=None)
        self.assertEqual(arr.dtype, np.int64)
        self.assertEqual(arr.fill_value, 0)

    def test_constructor_spindex_dtype(self):
        arr = SparseArray(data=[1, 2], sparse_index=IntIndex(4, [1, 2]))
        tm.assert_sp_array_equal(arr, SparseArray([np.nan, 1, 2, np.nan]))
        self.assertEqual(arr.dtype, np.float64)
        self.assertTrue(np.isnan(arr.fill_value))

        arr = SparseArray(data=[0, 1, 2, 3],
                          sparse_index=IntIndex(4, [0, 1, 2, 3]),
                          dtype=np.int64)
        exp = SparseArray([0, 1, 2, 3], dtype=np.int64)
        tm.assert_sp_array_equal(arr, exp)
        self.assertEqual(arr.dtype, np.int64)
        self.assertTrue(np.isnan(arr.fill_value))

        arr = SparseArray(data=[1, 2], sparse_index=IntIndex(4, [1, 2]),
                          fill_value=0, dtype=np.int64)
        exp = SparseArray([0, 1, 2, 0], fill_value=0, dtype=np.int64)
        tm.assert_sp_array_equal(arr, exp)
        self.assertEqual(arr.dtype, np.int64)
        self.assertEqual(arr.fill_value, 0)

        arr = SparseArray(data=[0, 1, 2, 3],
                          sparse_index=IntIndex(4, [0, 1, 2, 3]),
                          dtype=None)
        exp = SparseArray([0, 1, 2, 3], dtype=None)
        tm.assert_sp_array_equal(arr, exp)
        self.assertEqual(arr.dtype, np.int64)
        self.assertTrue(np.isnan(arr.fill_value))

        # scalar input
        arr = SparseArray(data=1,
                          sparse_index=IntIndex(1, [0]),
                          dtype=None)
        exp = SparseArray([1], dtype=None)
        tm.assert_sp_array_equal(arr, exp)
        self.assertEqual(arr.dtype, np.int64)
        self.assertTrue(np.isnan(arr.fill_value))

        arr = SparseArray(data=[1, 2], sparse_index=IntIndex(4, [1, 2]),
                          fill_value=0, dtype=None)
        exp = SparseArray([0, 1, 2, 0], fill_value=0, dtype=None)
        tm.assert_sp_array_equal(arr, exp)
        self.assertEqual(arr.dtype, np.int64)
        self.assertEqual(arr.fill_value, 0)

    def test_get_item(self):

        self.assertTrue(np.isnan(self.arr[1]))
        self.assertEqual(self.arr[2], 1)
        self.assertEqual(self.arr[7], 5)

        self.assertEqual(self.zarr[0], 0)
        self.assertEqual(self.zarr[2], 1)
        self.assertEqual(self.zarr[7], 5)

        errmsg = re.compile("bounds")
        assertRaisesRegexp(IndexError, errmsg, lambda: self.arr[11])
        assertRaisesRegexp(IndexError, errmsg, lambda: self.arr[-11])
        self.assertEqual(self.arr[-1], self.arr[len(self.arr) - 1])

    def test_take(self):
        self.assertTrue(np.isnan(self.arr.take(0)))
        self.assertTrue(np.isscalar(self.arr.take(2)))

        # np.take in < 1.8 doesn't support scalar indexing
        if not _np_version_under1p8:
            self.assertEqual(self.arr.take(2), np.take(self.arr_data, 2))
            self.assertEqual(self.arr.take(6), np.take(self.arr_data, 6))

        exp = SparseArray(np.take(self.arr_data, [2, 3]))
        tm.assert_sp_array_equal(self.arr.take([2, 3]), exp)

        exp = SparseArray(np.take(self.arr_data, [0, 1, 2]))
        tm.assert_sp_array_equal(self.arr.take([0, 1, 2]), exp)

    def test_take_fill_value(self):
        data = np.array([1, np.nan, 0, 3, 0])
        sparse = SparseArray(data, fill_value=0)

        exp = SparseArray(np.take(data, [0]), fill_value=0)
        tm.assert_sp_array_equal(sparse.take([0]), exp)

        exp = SparseArray(np.take(data, [1, 3, 4]), fill_value=0)
        tm.assert_sp_array_equal(sparse.take([1, 3, 4]), exp)

    def test_take_negative(self):
        exp = SparseArray(np.take(self.arr_data, [-1]))
        tm.assert_sp_array_equal(self.arr.take([-1]), exp)

        exp = SparseArray(np.take(self.arr_data, [-4, -3, -2]))
        tm.assert_sp_array_equal(self.arr.take([-4, -3, -2]), exp)

    def test_bad_take(self):
        assertRaisesRegexp(IndexError, "bounds", lambda: self.arr.take(11))
        self.assertRaises(IndexError, lambda: self.arr.take(-11))

    def test_take_invalid_kwargs(self):
        msg = "take\(\) got an unexpected keyword argument 'foo'"
        tm.assertRaisesRegexp(TypeError, msg, self.arr.take,
                              [2, 3], foo=2)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, self.arr.take,
                              [2, 3], out=self.arr)

        msg = "the 'mode' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, self.arr.take,
                              [2, 3], mode='clip')

    def test_take_filling(self):
        # similar tests as GH 12631
        sparse = SparseArray([np.nan, np.nan, 1, np.nan, 4])
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([np.nan, np.nan, 4])
        tm.assert_sp_array_equal(result, expected)

        # fill_value
        result = sparse.take(np.array([1, 0, -1]), fill_value=True)
        expected = SparseArray([np.nan, np.nan, np.nan])
        tm.assert_sp_array_equal(result, expected)

        # allow_fill=False
        result = sparse.take(np.array([1, 0, -1]),
                             allow_fill=False, fill_value=True)
        expected = SparseArray([np.nan, np.nan, 4])
        tm.assert_sp_array_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            sparse.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            sparse.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, -6]))
        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, 5]))
        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, 5]), fill_value=True)

    def test_take_filling_fill_value(self):
        # same tests as GH 12631
        sparse = SparseArray([np.nan, 0, 1, 0, 4], fill_value=0)
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([0, np.nan, 4], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        # fill_value
        result = sparse.take(np.array([1, 0, -1]), fill_value=True)
        expected = SparseArray([0, np.nan, 0], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        # allow_fill=False
        result = sparse.take(np.array([1, 0, -1]),
                             allow_fill=False, fill_value=True)
        expected = SparseArray([0, np.nan, 4], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            sparse.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            sparse.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, -6]))
        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, 5]))
        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, 5]), fill_value=True)

    def test_take_filling_all_nan(self):
        sparse = SparseArray([np.nan, np.nan, np.nan, np.nan, np.nan])
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([np.nan, np.nan, np.nan])
        tm.assert_sp_array_equal(result, expected)

        result = sparse.take(np.array([1, 0, -1]), fill_value=True)
        expected = SparseArray([np.nan, np.nan, np.nan])
        tm.assert_sp_array_equal(result, expected)

        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, -6]))
        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, 5]))
        with tm.assertRaises(IndexError):
            sparse.take(np.array([1, 5]), fill_value=True)

    def test_set_item(self):
        def setitem():
            self.arr[5] = 3

        def setslice():
            self.arr[1:5] = 2

        assertRaisesRegexp(TypeError, "item assignment", setitem)
        assertRaisesRegexp(TypeError, "item assignment", setslice)

    def test_constructor_from_too_large_array(self):
        assertRaisesRegexp(TypeError, "expected dimension <= 1 data",
                           SparseArray, np.arange(10).reshape((2, 5)))

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

    def test_constructor_bool(self):
        # GH 10648
        data = np.array([False, False, True, True, False, False])
        arr = SparseArray(data, fill_value=False, dtype=bool)

        self.assertEqual(arr.dtype, bool)
        tm.assert_numpy_array_equal(arr.sp_values, np.array([True, True]))
        tm.assert_numpy_array_equal(arr.sp_values, np.asarray(arr))
        tm.assert_numpy_array_equal(arr.sp_index.indices, np.array([2, 3]))

        for dense in [arr.to_dense(), arr.values]:
            self.assertEqual(dense.dtype, bool)
            tm.assert_numpy_array_equal(dense, data)

    def test_constructor_bool_fill_value(self):
        arr = SparseArray([True, False, True], dtype=None)
        self.assertEqual(arr.dtype, np.bool)
        self.assertFalse(arr.fill_value)

        arr = SparseArray([True, False, True], dtype=np.bool)
        self.assertEqual(arr.dtype, np.bool)
        self.assertFalse(arr.fill_value)

        arr = SparseArray([True, False, True], dtype=np.bool, fill_value=True)
        self.assertEqual(arr.dtype, np.bool)
        self.assertTrue(arr.fill_value)

    def test_constructor_float32(self):
        # GH 10648
        data = np.array([1., np.nan, 3], dtype=np.float32)
        arr = SparseArray(data, dtype=np.float32)

        self.assertEqual(arr.dtype, np.float32)
        tm.assert_numpy_array_equal(arr.sp_values, np.array([1, 3]))
        tm.assert_numpy_array_equal(arr.sp_values, np.asarray(arr))
        tm.assert_numpy_array_equal(arr.sp_index.indices, np.array([0, 2]))

        for dense in [arr.to_dense(), arr.values]:
            self.assertEqual(dense.dtype, np.float32)
            self.assert_numpy_array_equal(dense, data)

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

        assert (_get_base(arr2) is _get_base(self.arr))

    def test_values_asarray(self):
        assert_almost_equal(self.arr.values, self.arr_data)
        assert_almost_equal(self.arr.to_dense(), self.arr_data)
        assert_almost_equal(self.arr.sp_values, np.asarray(self.arr))

    def test_to_dense(self):
        vals = np.array([1, np.nan, np.nan, 3, np.nan])
        res = SparseArray(vals).to_dense()
        tm.assert_numpy_array_equal(res, vals)

        res = SparseArray(vals, fill_value=0).to_dense()
        tm.assert_numpy_array_equal(res, vals)

        vals = np.array([1, np.nan, 0, 3, 0])
        res = SparseArray(vals).to_dense()
        tm.assert_numpy_array_equal(res, vals)

        res = SparseArray(vals, fill_value=0).to_dense()
        tm.assert_numpy_array_equal(res, vals)

        vals = np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
        res = SparseArray(vals).to_dense()
        tm.assert_numpy_array_equal(res, vals)

        res = SparseArray(vals, fill_value=0).to_dense()
        tm.assert_numpy_array_equal(res, vals)

    def test_getitem(self):
        def _checkit(i):
            assert_almost_equal(self.arr[i], self.arr.values[i])

        for i in range(len(self.arr)):
            _checkit(i)
            _checkit(-i)

    def test_getslice(self):
        result = self.arr[:-3]
        exp = SparseArray(self.arr.values[:-3])
        tm.assert_sp_array_equal(result, exp)

        result = self.arr[-4:]
        exp = SparseArray(self.arr.values[-4:])
        tm.assert_sp_array_equal(result, exp)

        # two corner cases from Series
        result = self.arr[-12:]
        exp = SparseArray(self.arr)
        tm.assert_sp_array_equal(result, exp)

        result = self.arr[:-12]
        exp = SparseArray(self.arr.values[:0])
        tm.assert_sp_array_equal(result, exp)

    def test_getslice_tuple(self):
        dense = np.array([np.nan, 0, 3, 4, 0, 5, np.nan, np.nan, 0])

        sparse = SparseArray(dense)
        res = sparse[4:, ]
        exp = SparseArray(dense[4:, ])
        tm.assert_sp_array_equal(res, exp)

        sparse = SparseArray(dense, fill_value=0)
        res = sparse[4:, ]
        exp = SparseArray(dense[4:, ], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

        with tm.assertRaises(IndexError):
            sparse[4:, :]

        with tm.assertRaises(IndexError):
            # check numpy compat
            dense[4:, :]

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
            tm.assert_sp_array_equal(res, res2)

            res3 = op(first.values, second)
            tm.assertIsInstance(res3, SparseArray)
            tm.assert_sp_array_equal(res, res3)

            res4 = op(first, 4)
            tm.assertIsInstance(res4, SparseArray)

            # ignore this if the actual op raises (e.g. pow)
            try:
                exp = op(first.values, 4)
                exp_fv = op(first.fill_value, 4)
                assert_almost_equal(res4.fill_value, exp_fv)
                assert_almost_equal(res4.values, exp)
            except ValueError:
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
            tm.assert_sp_array_equal(unpickled, obj)

        _check_roundtrip(self.arr)
        _check_roundtrip(self.zarr)

    def test_generator_warnings(self):
        sp_arr = SparseArray([1, 2, 3])
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings(action='always',
                                    category=DeprecationWarning)
            warnings.filterwarnings(action='always',
                                    category=PendingDeprecationWarning)
            for _ in sp_arr:
                pass
            assert len(w) == 0

    def test_fillna(self):
        s = SparseArray([1, np.nan, np.nan, 3, np.nan])
        res = s.fillna(-1)
        exp = SparseArray([1, -1, -1, 3, -1], fill_value=-1)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([1, np.nan, np.nan, 3, np.nan], fill_value=0)
        res = s.fillna(-1)
        exp = SparseArray([1, -1, -1, 3, -1], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([1, np.nan, 0, 3, 0])
        res = s.fillna(-1)
        exp = SparseArray([1, -1, 0, 3, 0], fill_value=-1)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([1, np.nan, 0, 3, 0], fill_value=0)
        res = s.fillna(-1)
        exp = SparseArray([1, -1, 0, 3, 0], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([np.nan, np.nan, np.nan, np.nan])
        res = s.fillna(-1)
        exp = SparseArray([-1, -1, -1, -1], fill_value=-1)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([np.nan, np.nan, np.nan, np.nan], fill_value=0)
        res = s.fillna(-1)
        exp = SparseArray([-1, -1, -1, -1], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([0, 0, 0, 0])
        res = s.fillna(-1)
        exp = SparseArray([0, 0, 0, 0], fill_value=-1)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([0, 0, 0, 0], fill_value=0)
        res = s.fillna(-1)
        exp = SparseArray([0, 0, 0, 0], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

    def test_fillna_overlap(self):
        s = SparseArray([1, np.nan, np.nan, 3, np.nan])
        # filling with existing value doesn't replace existing value with
        # fill_value, i.e. existing 3 remains in sp_values
        res = s.fillna(3)
        exp = np.array([1, 3, 3, 3, 3])
        tm.assert_numpy_array_equal(res.to_dense(), exp)

        s = SparseArray([1, np.nan, np.nan, 3, np.nan], fill_value=0)
        res = s.fillna(3)
        exp = SparseArray([1, 3, 3, 3, 3], fill_value=0)
        tm.assert_sp_array_equal(res, exp)


class TestSparseArrayArithmetic(tm.TestCase):

    _multiprocess_can_split_ = True

    def _check_numeric_ops(self, a, b, a_dense, b_dense):
        tm.assert_numpy_array_equal((a + b).to_dense(), a_dense + b_dense)
        tm.assert_numpy_array_equal((b + a).to_dense(), b_dense + a_dense)

        tm.assert_numpy_array_equal((a - b).to_dense(), a_dense - b_dense)
        tm.assert_numpy_array_equal((b - a).to_dense(), b_dense - a_dense)

        tm.assert_numpy_array_equal((a * b).to_dense(), a_dense * b_dense)
        tm.assert_numpy_array_equal((b * a).to_dense(), b_dense * a_dense)

        tm.assert_numpy_array_equal((a / b).to_dense(), a_dense / b_dense)
        tm.assert_numpy_array_equal((b / a).to_dense(), b_dense / a_dense)

        tm.assert_numpy_array_equal((a // b).to_dense(), a_dense // b_dense)
        tm.assert_numpy_array_equal((b // a).to_dense(), b_dense // a_dense)

        tm.assert_numpy_array_equal((a % b).to_dense(), a_dense % b_dense)
        tm.assert_numpy_array_equal((b % a).to_dense(), b_dense % a_dense)

        tm.assert_numpy_array_equal((a ** b).to_dense(), a_dense ** b_dense)
        tm.assert_numpy_array_equal((b ** a).to_dense(), b_dense ** a_dense)

    def _check_comparison_ops(self, a, b, a_dense, b_dense):

        def _check(res):
            tm.assertIsInstance(res, SparseArray)
            self.assertEqual(res.dtype, np.bool)
            self.assertIsInstance(res.fill_value, bool)

        _check(a == b)
        tm.assert_numpy_array_equal((a == b).to_dense(), a_dense == b_dense)

        _check(a != b)
        tm.assert_numpy_array_equal((a != b).to_dense(), a_dense != b_dense)

        _check(a >= b)
        tm.assert_numpy_array_equal((a >= b).to_dense(), a_dense >= b_dense)

        _check(a <= b)
        tm.assert_numpy_array_equal((a <= b).to_dense(), a_dense <= b_dense)

        _check(a > b)
        tm.assert_numpy_array_equal((a > b).to_dense(), a_dense > b_dense)

        _check(a < b)
        tm.assert_numpy_array_equal((a < b).to_dense(), a_dense < b_dense)

    def test_float_scalar(self):
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])

        for kind in ['integer', 'block']:
            a = SparseArray(values, kind=kind)
            self._check_numeric_ops(a, 1, values, 1)
            self._check_numeric_ops(a, 0, values, 0)
            self._check_numeric_ops(a, 3, values, 3)

            a = SparseArray(values, kind=kind, fill_value=0)
            self._check_numeric_ops(a, 1, values, 1)
            self._check_numeric_ops(a, 0, values, 0)
            self._check_numeric_ops(a, 3, values, 3)

            a = SparseArray(values, kind=kind, fill_value=2)
            self._check_numeric_ops(a, 1, values, 1)
            self._check_numeric_ops(a, 0, values, 0)
            self._check_numeric_ops(a, 3, values, 3)

    def test_float_scalar_comparison(self):
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])

        for kind in ['integer', 'block']:
            a = SparseArray(values, kind=kind)
            self._check_comparison_ops(a, 1, values, 1)
            self._check_comparison_ops(a, 0, values, 0)
            self._check_comparison_ops(a, 3, values, 3)

            a = SparseArray(values, kind=kind, fill_value=0)
            self._check_comparison_ops(a, 1, values, 1)
            self._check_comparison_ops(a, 0, values, 0)
            self._check_comparison_ops(a, 3, values, 3)

            a = SparseArray(values, kind=kind, fill_value=2)
            self._check_comparison_ops(a, 1, values, 1)
            self._check_comparison_ops(a, 0, values, 0)
            self._check_comparison_ops(a, 3, values, 3)

    def test_float_same_index(self):
        # when sp_index are the same
        for kind in ['integer', 'block']:
            values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
            rvalues = np.array([np.nan, 2, 3, 4, np.nan, 0, 1, 3, 2, np.nan])

            a = SparseArray(values, kind=kind)
            b = SparseArray(rvalues, kind=kind)
            self._check_numeric_ops(a, b, values, rvalues)

            values = np.array([0., 1., 2., 6., 0., 0., 1., 2., 1., 0.])
            rvalues = np.array([0., 2., 3., 4., 0., 0., 1., 3., 2., 0.])

            a = SparseArray(values, kind=kind, fill_value=0)
            b = SparseArray(rvalues, kind=kind, fill_value=0)
            self._check_numeric_ops(a, b, values, rvalues)

    def test_float_same_index_comparison(self):
        # when sp_index are the same
        for kind in ['integer', 'block']:
            values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
            rvalues = np.array([np.nan, 2, 3, 4, np.nan, 0, 1, 3, 2, np.nan])

            a = SparseArray(values, kind=kind)
            b = SparseArray(rvalues, kind=kind)
            self._check_comparison_ops(a, b, values, rvalues)

            values = np.array([0., 1., 2., 6., 0., 0., 1., 2., 1., 0.])
            rvalues = np.array([0., 2., 3., 4., 0., 0., 1., 3., 2., 0.])

            a = SparseArray(values, kind=kind, fill_value=0)
            b = SparseArray(rvalues, kind=kind, fill_value=0)
            self._check_comparison_ops(a, b, values, rvalues)

    def test_float_array(self):
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, np.nan, 2, 3, np.nan, 0, 1, 5, 2, np.nan])

        for kind in ['integer', 'block']:
            a = SparseArray(values, kind=kind)
            b = SparseArray(rvalues, kind=kind)
            self._check_numeric_ops(a, b, values, rvalues)
            self._check_numeric_ops(a, b * 0, values, rvalues * 0)

            a = SparseArray(values, kind=kind, fill_value=0)
            b = SparseArray(rvalues, kind=kind)
            self._check_numeric_ops(a, b, values, rvalues)

            a = SparseArray(values, kind=kind, fill_value=0)
            b = SparseArray(rvalues, kind=kind, fill_value=0)
            self._check_numeric_ops(a, b, values, rvalues)

            a = SparseArray(values, kind=kind, fill_value=1)
            b = SparseArray(rvalues, kind=kind, fill_value=2)
            self._check_numeric_ops(a, b, values, rvalues)

    def test_float_array_different_kind(self):
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, np.nan, 2, 3, np.nan, 0, 1, 5, 2, np.nan])

        a = SparseArray(values, kind='integer')
        b = SparseArray(rvalues, kind='block')
        self._check_numeric_ops(a, b, values, rvalues)
        self._check_numeric_ops(a, b * 0, values, rvalues * 0)

        a = SparseArray(values, kind='integer', fill_value=0)
        b = SparseArray(rvalues, kind='block')
        self._check_numeric_ops(a, b, values, rvalues)

        a = SparseArray(values, kind='integer', fill_value=0)
        b = SparseArray(rvalues, kind='block', fill_value=0)
        self._check_numeric_ops(a, b, values, rvalues)

        a = SparseArray(values, kind='integer', fill_value=1)
        b = SparseArray(rvalues, kind='block', fill_value=2)
        self._check_numeric_ops(a, b, values, rvalues)

    def test_float_array_comparison(self):
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, np.nan, 2, 3, np.nan, 0, 1, 5, 2, np.nan])

        for kind in ['integer', 'block']:
            a = SparseArray(values, kind=kind)
            b = SparseArray(rvalues, kind=kind)
            self._check_comparison_ops(a, b, values, rvalues)
            self._check_comparison_ops(a, b * 0, values, rvalues * 0)

            a = SparseArray(values, kind=kind, fill_value=0)
            b = SparseArray(rvalues, kind=kind)
            self._check_comparison_ops(a, b, values, rvalues)

            a = SparseArray(values, kind=kind, fill_value=0)
            b = SparseArray(rvalues, kind=kind, fill_value=0)
            self._check_comparison_ops(a, b, values, rvalues)

            a = SparseArray(values, kind=kind, fill_value=1)
            b = SparseArray(rvalues, kind=kind, fill_value=2)
            self._check_comparison_ops(a, b, values, rvalues)


class TestSparseArrayAnalytics(tm.TestCase):
    def test_sum(self):
        data = np.arange(10).astype(float)
        out = SparseArray(data).sum()
        self.assertEqual(out, 45.0)

        data[5] = np.nan
        out = SparseArray(data, fill_value=2).sum()
        self.assertEqual(out, 40.0)

        out = SparseArray(data, fill_value=np.nan).sum()
        self.assertEqual(out, 40.0)

    def test_numpy_sum(self):
        data = np.arange(10).astype(float)
        out = np.sum(SparseArray(data))
        self.assertEqual(out, 45.0)

        data[5] = np.nan
        out = np.sum(SparseArray(data, fill_value=2))
        self.assertEqual(out, 40.0)

        out = np.sum(SparseArray(data, fill_value=np.nan))
        self.assertEqual(out, 40.0)

        msg = "the 'dtype' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.sum,
                              SparseArray(data), dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.sum,
                              SparseArray(data), out=out)

    def test_cumsum(self):
        data = np.arange(10).astype(float)
        out = SparseArray(data).cumsum()
        expected = SparseArray(data.cumsum())
        tm.assert_sp_array_equal(out, expected)

        # TODO: gh-12855 - return a SparseArray here
        data[5] = np.nan
        out = SparseArray(data, fill_value=2).cumsum()
        self.assertNotIsInstance(out, SparseArray)
        tm.assert_numpy_array_equal(out, data.cumsum())

        out = SparseArray(data, fill_value=np.nan).cumsum()
        expected = SparseArray(np.array([
            0, 1, 3, 6, 10, np.nan, 16, 23, 31, 40]))
        tm.assert_sp_array_equal(out, expected)

    def test_numpy_cumsum(self):
        data = np.arange(10).astype(float)
        out = np.cumsum(SparseArray(data))
        expected = SparseArray(data.cumsum())
        tm.assert_sp_array_equal(out, expected)

        # TODO: gh-12855 - return a SparseArray here
        data[5] = np.nan
        out = np.cumsum(SparseArray(data, fill_value=2))
        self.assertNotIsInstance(out, SparseArray)
        tm.assert_numpy_array_equal(out, data.cumsum())

        out = np.cumsum(SparseArray(data, fill_value=np.nan))
        expected = SparseArray(np.array([
            0, 1, 3, 6, 10, np.nan, 16, 23, 31, 40]))
        tm.assert_sp_array_equal(out, expected)

        msg = "the 'dtype' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.cumsum,
                              SparseArray(data), dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.cumsum,
                              SparseArray(data), out=out)

    def test_mean(self):
        data = np.arange(10).astype(float)
        out = SparseArray(data).mean()
        self.assertEqual(out, 4.5)

        data[5] = np.nan
        out = SparseArray(data).mean()
        self.assertEqual(out, 40.0 / 9)

    def test_numpy_mean(self):
        data = np.arange(10).astype(float)
        out = np.mean(SparseArray(data))
        self.assertEqual(out, 4.5)

        data[5] = np.nan
        out = np.mean(SparseArray(data))
        self.assertEqual(out, 40.0 / 9)

        msg = "the 'dtype' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.mean,
                              SparseArray(data), dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.mean,
                              SparseArray(data), out=out)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
