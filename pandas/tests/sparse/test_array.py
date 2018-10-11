from pandas.compat import range

import re
import operator
import pytest
import warnings

from numpy import nan
import numpy as np
import pandas as pd

from pandas.core.sparse.api import SparseArray, SparseSeries, SparseDtype
from pandas._libs.sparse import IntIndex
from pandas.util.testing import assert_almost_equal
import pandas.util.testing as tm
import pandas.util._test_decorators as td


@pytest.fixture(params=["integer", "block"])
def kind(request):
    return request.param


class TestSparseArray(object):

    def setup_method(self, method):
        self.arr_data = np.array([nan, nan, 1, 2, 3, nan, 4, 5, nan, 6])
        self.arr = SparseArray(self.arr_data)
        self.zarr = SparseArray([0, 0, 1, 2, 3, 0, 4, 5, 0, 6], fill_value=0)

    def test_constructor_dtype(self):
        arr = SparseArray([np.nan, 1, 2, np.nan])
        assert arr.dtype == SparseDtype(np.float64, np.nan)
        assert arr.dtype.subtype == np.float64
        assert np.isnan(arr.fill_value)

        arr = SparseArray([np.nan, 1, 2, np.nan], fill_value=0)
        assert arr.dtype == SparseDtype(np.float64, 0)
        assert arr.fill_value == 0

        arr = SparseArray([0, 1, 2, 4], dtype=np.float64)
        assert arr.dtype == SparseDtype(np.float64, np.nan)
        assert np.isnan(arr.fill_value)

        arr = SparseArray([0, 1, 2, 4], dtype=np.int64)
        assert arr.dtype == SparseDtype(np.int64, 0)
        assert arr.fill_value == 0

        arr = SparseArray([0, 1, 2, 4], fill_value=0, dtype=np.int64)
        assert arr.dtype == SparseDtype(np.int64, 0)
        assert arr.fill_value == 0

        arr = SparseArray([0, 1, 2, 4], dtype=None)
        assert arr.dtype == SparseDtype(np.int64, 0)
        assert arr.fill_value == 0

        arr = SparseArray([0, 1, 2, 4], fill_value=0, dtype=None)
        assert arr.dtype == SparseDtype(np.int64, 0)
        assert arr.fill_value == 0

    def test_constructor_dtype_str(self):
        result = SparseArray([1, 2, 3], dtype='int')
        expected = SparseArray([1, 2, 3], dtype=int)
        tm.assert_sp_array_equal(result, expected)

    def test_constructor_sparse_dtype(self):
        result = SparseArray([1, 0, 0, 1], dtype=SparseDtype('int64', -1))
        expected = SparseArray([1, 0, 0, 1], fill_value=-1, dtype=np.int64)
        tm.assert_sp_array_equal(result, expected)
        assert result.sp_values.dtype == np.dtype('int64')

    def test_constructor_sparse_dtype_str(self):
        result = SparseArray([1, 0, 0, 1], dtype='Sparse[int32]')
        expected = SparseArray([1, 0, 0, 1], dtype=np.int32)
        tm.assert_sp_array_equal(result, expected)
        assert result.sp_values.dtype == np.dtype('int32')

    def test_constructor_object_dtype(self):
        # GH 11856
        arr = SparseArray(['A', 'A', np.nan, 'B'], dtype=np.object)
        assert arr.dtype == SparseDtype(np.object)
        assert np.isnan(arr.fill_value)

        arr = SparseArray(['A', 'A', np.nan, 'B'], dtype=np.object,
                          fill_value='A')
        assert arr.dtype == SparseDtype(np.object, 'A')
        assert arr.fill_value == 'A'

        # GH 17574
        data = [False, 0, 100.0, 0.0]
        arr = SparseArray(data, dtype=np.object, fill_value=False)
        assert arr.dtype == SparseDtype(np.object, False)
        assert arr.fill_value is False
        arr_expected = np.array(data, dtype=np.object)
        it = (type(x) == type(y) and x == y for x, y in zip(arr, arr_expected))
        assert np.fromiter(it, dtype=np.bool).all()

    @pytest.mark.parametrize("dtype", [SparseDtype(int, 0), int])
    def test_constructor_na_dtype(self, dtype):
        with tm.assert_raises_regex(ValueError, "Cannot convert"):
            SparseArray([0, 1, np.nan], dtype=dtype)

    def test_constructor_spindex_dtype(self):
        arr = SparseArray(data=[1, 2], sparse_index=IntIndex(4, [1, 2]))
        # XXX: Behavior change: specifying SparseIndex no longer changes the
        # fill_value
        expected = SparseArray([0, 1, 2, 0], kind='integer')
        tm.assert_sp_array_equal(arr, expected)
        assert arr.dtype == SparseDtype(np.int64)
        assert arr.fill_value == 0

        arr = SparseArray(data=[1, 2, 3],
                          sparse_index=IntIndex(4, [1, 2, 3]),
                          dtype=np.int64, fill_value=0)
        exp = SparseArray([0, 1, 2, 3], dtype=np.int64, fill_value=0)
        tm.assert_sp_array_equal(arr, exp)
        assert arr.dtype == SparseDtype(np.int64)
        assert arr.fill_value == 0

        arr = SparseArray(data=[1, 2], sparse_index=IntIndex(4, [1, 2]),
                          fill_value=0, dtype=np.int64)
        exp = SparseArray([0, 1, 2, 0], fill_value=0, dtype=np.int64)
        tm.assert_sp_array_equal(arr, exp)
        assert arr.dtype == SparseDtype(np.int64)
        assert arr.fill_value == 0

        arr = SparseArray(data=[1, 2, 3],
                          sparse_index=IntIndex(4, [1, 2, 3]),
                          dtype=None, fill_value=0)
        exp = SparseArray([0, 1, 2, 3], dtype=None)
        tm.assert_sp_array_equal(arr, exp)
        assert arr.dtype == SparseDtype(np.int64)
        assert arr.fill_value == 0

    @pytest.mark.parametrize("sparse_index", [
        None, IntIndex(1, [0]),
    ])
    def test_constructor_spindex_dtype_scalar(self, sparse_index):
        # scalar input
        arr = SparseArray(data=1, sparse_index=sparse_index, dtype=None)
        exp = SparseArray([1], dtype=None)
        tm.assert_sp_array_equal(arr, exp)
        assert arr.dtype == SparseDtype(np.int64)
        assert arr.fill_value == 0

        arr = SparseArray(data=1, sparse_index=IntIndex(1, [0]), dtype=None)
        exp = SparseArray([1], dtype=None)
        tm.assert_sp_array_equal(arr, exp)
        assert arr.dtype == SparseDtype(np.int64)
        assert arr.fill_value == 0

    def test_constructor_spindex_dtype_scalar_broadcasts(self):
        arr = SparseArray(data=[1, 2], sparse_index=IntIndex(4, [1, 2]),
                          fill_value=0, dtype=None)
        exp = SparseArray([0, 1, 2, 0], fill_value=0, dtype=None)
        tm.assert_sp_array_equal(arr, exp)
        assert arr.dtype == SparseDtype(np.int64)
        assert arr.fill_value == 0

    @pytest.mark.parametrize('data, fill_value', [
        (np.array([1, 2]), 0),
        (np.array([1.0, 2.0]), np.nan),
        ([True, False], False),
        ([pd.Timestamp('2017-01-01')], pd.NaT),
    ])
    def test_constructor_inferred_fill_value(self, data, fill_value):
        result = SparseArray(data).fill_value

        if pd.isna(fill_value):
            assert pd.isna(result)
        else:
            assert result == fill_value

    @pytest.mark.parametrize('scalar,dtype', [
        (False, SparseDtype(bool, False)),
        (0.0, SparseDtype('float64', 0)),
        (1, SparseDtype('int64', 1)),
        ('z', SparseDtype('object', 'z'))])
    def test_scalar_with_index_infer_dtype(self, scalar, dtype):
        # GH 19163
        arr = SparseArray(scalar, index=[1, 2, 3], fill_value=scalar)
        exp = SparseArray([scalar, scalar, scalar], fill_value=scalar)

        tm.assert_sp_array_equal(arr, exp)

        assert arr.dtype == dtype
        assert exp.dtype == dtype

    @pytest.mark.parametrize("fill", [1, np.nan, 0])
    def test_sparse_series_round_trip(self, kind, fill):
        # see gh-13999
        arr = SparseArray([np.nan, 1, np.nan, 2, 3],
                          kind=kind, fill_value=fill)
        res = SparseArray(SparseSeries(arr))
        tm.assert_sp_array_equal(arr, res)

        arr = SparseArray([0, 0, 0, 1, 1, 2], dtype=np.int64,
                          kind=kind, fill_value=fill)
        res = SparseArray(SparseSeries(arr), dtype=np.int64)
        tm.assert_sp_array_equal(arr, res)

        res = SparseArray(SparseSeries(arr))
        tm.assert_sp_array_equal(arr, res)

    @pytest.mark.parametrize("fill", [True, False, np.nan])
    def test_sparse_series_round_trip2(self, kind, fill):
        # see gh-13999
        arr = SparseArray([True, False, True, True], dtype=np.bool,
                          kind=kind, fill_value=fill)
        res = SparseArray(SparseSeries(arr))
        tm.assert_sp_array_equal(arr, res)

        res = SparseArray(SparseSeries(arr))
        tm.assert_sp_array_equal(arr, res)

    def test_get_item(self):

        assert np.isnan(self.arr[1])
        assert self.arr[2] == 1
        assert self.arr[7] == 5

        assert self.zarr[0] == 0
        assert self.zarr[2] == 1
        assert self.zarr[7] == 5

        errmsg = re.compile("bounds")
        tm.assert_raises_regex(IndexError, errmsg, lambda: self.arr[11])
        tm.assert_raises_regex(IndexError, errmsg, lambda: self.arr[-11])
        assert self.arr[-1] == self.arr[len(self.arr) - 1]

    def test_take_scalar_raises(self):
        msg = "'indices' must be an array, not a scalar '2'."
        with tm.assert_raises_regex(ValueError, msg):
            self.arr.take(2)

    def test_take(self):
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
        tm.assert_raises_regex(
            IndexError, "bounds", lambda: self.arr.take([11]))

    def test_take_filling(self):
        # similar tests as GH 12631
        sparse = SparseArray([np.nan, np.nan, 1, np.nan, 4])
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([np.nan, np.nan, 4])
        tm.assert_sp_array_equal(result, expected)

        # XXX: test change: fill_value=True -> allow_fill=True
        result = sparse.take(np.array([1, 0, -1]), allow_fill=True)
        expected = SparseArray([np.nan, np.nan, np.nan])
        tm.assert_sp_array_equal(result, expected)

        # allow_fill=False
        result = sparse.take(np.array([1, 0, -1]),
                             allow_fill=False, fill_value=True)
        expected = SparseArray([np.nan, np.nan, 4])
        tm.assert_sp_array_equal(result, expected)

        msg = ("Invalid value in 'indices'")
        with tm.assert_raises_regex(ValueError, msg):
            sparse.take(np.array([1, 0, -2]), allow_fill=True)
        with tm.assert_raises_regex(ValueError, msg):
            sparse.take(np.array([1, 0, -5]), allow_fill=True)

        with pytest.raises(IndexError):
            sparse.take(np.array([1, -6]))
        with pytest.raises(IndexError):
            sparse.take(np.array([1, 5]))
        with pytest.raises(IndexError):
            sparse.take(np.array([1, 5]), allow_fill=True)

    def test_take_filling_fill_value(self):
        # same tests as GH 12631
        sparse = SparseArray([np.nan, 0, 1, 0, 4], fill_value=0)
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([0, np.nan, 4], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        # fill_value
        result = sparse.take(np.array([1, 0, -1]), allow_fill=True)
        # XXX: behavior change.
        # the old way of filling self.fill_value doesn't follow EA rules.
        # It's supposed to be self.dtype.na_value (nan in this case)
        expected = SparseArray([0, np.nan, np.nan], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        # allow_fill=False
        result = sparse.take(np.array([1, 0, -1]),
                             allow_fill=False, fill_value=True)
        expected = SparseArray([0, np.nan, 4], fill_value=0)
        tm.assert_sp_array_equal(result, expected)

        msg = ("Invalid value in 'indices'.")
        with tm.assert_raises_regex(ValueError, msg):
            sparse.take(np.array([1, 0, -2]), allow_fill=True)
        with tm.assert_raises_regex(ValueError, msg):
            sparse.take(np.array([1, 0, -5]), allow_fill=True)

        with pytest.raises(IndexError):
            sparse.take(np.array([1, -6]))
        with pytest.raises(IndexError):
            sparse.take(np.array([1, 5]))
        with pytest.raises(IndexError):
            sparse.take(np.array([1, 5]), fill_value=True)

    def test_take_filling_all_nan(self):
        sparse = SparseArray([np.nan, np.nan, np.nan, np.nan, np.nan])
        # XXX: did the default kind from take change?
        result = sparse.take(np.array([1, 0, -1]))
        expected = SparseArray([np.nan, np.nan, np.nan], kind='block')
        tm.assert_sp_array_equal(result, expected)

        result = sparse.take(np.array([1, 0, -1]), fill_value=True)
        expected = SparseArray([np.nan, np.nan, np.nan], kind='block')
        tm.assert_sp_array_equal(result, expected)

        with pytest.raises(IndexError):
            sparse.take(np.array([1, -6]))
        with pytest.raises(IndexError):
            sparse.take(np.array([1, 5]))
        with pytest.raises(IndexError):
            sparse.take(np.array([1, 5]), fill_value=True)

    def test_set_item(self):
        def setitem():
            self.arr[5] = 3

        def setslice():
            self.arr[1:5] = 2

        tm.assert_raises_regex(TypeError, "item assignment", setitem)
        tm.assert_raises_regex(TypeError, "item assignment", setslice)

    def test_constructor_from_too_large_array(self):
        tm.assert_raises_regex(TypeError, "expected dimension <= 1 data",
                               SparseArray, np.arange(10).reshape((2, 5)))

    def test_constructor_from_sparse(self):
        res = SparseArray(self.zarr)
        assert res.fill_value == 0
        assert_almost_equal(res.sp_values, self.zarr.sp_values)

    def test_constructor_copy(self):
        cp = SparseArray(self.arr, copy=True)
        cp.sp_values[:3] = 0
        assert not (self.arr.sp_values[:3] == 0).any()

        not_copy = SparseArray(self.arr)
        not_copy.sp_values[:3] = 0
        assert (self.arr.sp_values[:3] == 0).all()

    def test_constructor_bool(self):
        # GH 10648
        data = np.array([False, False, True, True, False, False])
        arr = SparseArray(data, fill_value=False, dtype=bool)

        assert arr.dtype == SparseDtype(bool)
        tm.assert_numpy_array_equal(arr.sp_values, np.array([True, True]))
        # Behavior change: np.asarray densifies.
        # tm.assert_numpy_array_equal(arr.sp_values, np.asarray(arr))
        tm.assert_numpy_array_equal(arr.sp_index.indices,
                                    np.array([2, 3], np.int32))

        for dense in [arr.to_dense(), arr.values]:
            assert dense.dtype == bool
            tm.assert_numpy_array_equal(dense, data)

    def test_constructor_bool_fill_value(self):
        arr = SparseArray([True, False, True], dtype=None)
        assert arr.dtype == SparseDtype(np.bool)
        assert not arr.fill_value

        arr = SparseArray([True, False, True], dtype=np.bool)
        assert arr.dtype == SparseDtype(np.bool)
        assert not arr.fill_value

        arr = SparseArray([True, False, True], dtype=np.bool, fill_value=True)
        assert arr.dtype == SparseDtype(np.bool, True)
        assert arr.fill_value

    def test_constructor_float32(self):
        # GH 10648
        data = np.array([1., np.nan, 3], dtype=np.float32)
        arr = SparseArray(data, dtype=np.float32)

        assert arr.dtype == SparseDtype(np.float32)
        tm.assert_numpy_array_equal(arr.sp_values,
                                    np.array([1, 3], dtype=np.float32))
        # Behavior change: np.asarray densifies.
        # tm.assert_numpy_array_equal(arr.sp_values, np.asarray(arr))
        tm.assert_numpy_array_equal(arr.sp_index.indices,
                                    np.array([0, 2], dtype=np.int32))

        for dense in [arr.to_dense(), arr.values]:
            assert dense.dtype == np.float32
            tm.assert_numpy_array_equal(dense, data)

    def test_astype(self):
        # float -> float
        arr = SparseArray([None, None, 0, 2])
        result = arr.astype("Sparse[float32]")
        expected = SparseArray([None, None, 0, 2], dtype=np.dtype('float32'))
        tm.assert_sp_array_equal(result, expected)

        dtype = SparseDtype("float64", fill_value=0)
        result = arr.astype(dtype)
        expected = SparseArray._simple_new(np.array([0., 2.],
                                                    dtype=dtype.subtype),
                                           IntIndex(4, [2, 3]),
                                           dtype)
        tm.assert_sp_array_equal(result, expected)

        dtype = SparseDtype("int64", 0)
        result = arr.astype(dtype)
        expected = SparseArray._simple_new(np.array([0, 2], dtype=np.int64),
                                           IntIndex(4, [2, 3]),
                                           dtype)
        tm.assert_sp_array_equal(result, expected)

        arr = SparseArray([0, np.nan, 0, 1], fill_value=0)
        with tm.assert_raises_regex(ValueError, 'NA'):
            arr.astype('Sparse[i8]')

    def test_astype_bool(self):
        a = pd.SparseArray([1, 0, 0, 1], dtype=SparseDtype(int, 0))
        result = a.astype(bool)
        expected = SparseArray([True, 0, 0, True],
                               dtype=SparseDtype(bool, 0))
        tm.assert_sp_array_equal(result, expected)

        # update fill value
        result = a.astype(SparseDtype(bool, False))
        expected = SparseArray([True, False, False, True],
                               dtype=SparseDtype(bool, False))
        tm.assert_sp_array_equal(result, expected)

    def test_astype_all(self, any_real_dtype):
        vals = np.array([1, 2, 3])
        arr = SparseArray(vals, fill_value=1)
        typ = np.dtype(any_real_dtype)
        res = arr.astype(typ)
        assert res.dtype == SparseDtype(typ, 1)
        assert res.sp_values.dtype == typ

        tm.assert_numpy_array_equal(np.asarray(res.values),
                                    vals.astype(typ))

    def test_set_fill_value(self):
        arr = SparseArray([1., np.nan, 2.], fill_value=np.nan)
        arr.fill_value = 2
        assert arr.fill_value == 2

        arr = SparseArray([1, 0, 2], fill_value=0, dtype=np.int64)
        arr.fill_value = 2
        assert arr.fill_value == 2

        # XXX: this seems fine? You can construct an integer
        # sparsearray with NaN fill value, why not update one?
        # coerces to int
        # msg = "unable to set fill_value 3\\.1 to int64 dtype"
        # with tm.assert_raises_regex(ValueError, msg):
        arr.fill_value = 3.1
        assert arr.fill_value == 3.1

        # msg = "unable to set fill_value nan to int64 dtype"
        # with tm.assert_raises_regex(ValueError, msg):
        arr.fill_value = np.nan
        assert np.isnan(arr.fill_value)

        arr = SparseArray([True, False, True], fill_value=False, dtype=np.bool)
        arr.fill_value = True
        assert arr.fill_value

        # coerces to bool
        # msg = "unable to set fill_value 0 to bool dtype"
        # with tm.assert_raises_regex(ValueError, msg):
        arr.fill_value = 0
        assert arr.fill_value == 0

        # msg = "unable to set fill_value nan to bool dtype"
        # with tm.assert_raises_regex(ValueError, msg):
        arr.fill_value = np.nan
        assert np.isnan(arr.fill_value)

    @pytest.mark.parametrize("val", [[1, 2, 3], np.array([1, 2]), (1, 2, 3)])
    def test_set_fill_invalid_non_scalar(self, val):
        arr = SparseArray([True, False, True], fill_value=False, dtype=np.bool)
        msg = "fill_value must be a scalar"

        with tm.assert_raises_regex(ValueError, msg):
            arr.fill_value = val

    def test_copy_shallow(self):
        arr2 = self.arr.copy(deep=False)
        assert arr2.sp_values is self.arr.sp_values
        assert arr2.sp_index is self.arr.sp_index

    def test_values_asarray(self):
        assert_almost_equal(self.arr.values, self.arr_data)
        assert_almost_equal(self.arr.to_dense(), self.arr_data)

    @pytest.mark.parametrize('data,shape,dtype', [
        ([0, 0, 0, 0, 0], (5,), None),
        ([], (0,), None),
        ([0], (1,), None),
        (['A', 'A', np.nan, 'B'], (4,), np.object)
    ])
    def test_shape(self, data, shape, dtype):
        # GH 21126
        out = SparseArray(data, dtype=dtype)
        assert out.shape == shape

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

        # see gh-14647
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            SparseArray(vals).to_dense(fill=2)

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

        with pytest.raises(IndexError):
            sparse[4:, :]

        with pytest.raises(IndexError):
            # check numpy compat
            dense[4:, :]

    def test_boolean_slice_empty(self):
        arr = pd.SparseArray([0, 1, 2])
        res = arr[[False, False, False]]
        assert res.dtype == arr.dtype

    @pytest.mark.parametrize("op", ["add", "sub", "mul",
                                    "truediv", "floordiv", "pow"])
    def test_binary_operators(self, op):
        op = getattr(operator, op)
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
            assert isinstance(res, SparseArray)
            assert_almost_equal(res.values, exp.values)

            res2 = op(first, second.values)
            assert isinstance(res2, SparseArray)
            tm.assert_sp_array_equal(res, res2)

            res3 = op(first.values, second)
            assert isinstance(res3, SparseArray)
            tm.assert_sp_array_equal(res, res3)

            res4 = op(first, 4)
            assert isinstance(res4, SparseArray)

            # Ignore this if the actual op raises (e.g. pow).
            try:
                exp = op(first.values, 4)
                exp_fv = op(first.fill_value, 4)
            except ValueError:
                pass
            else:
                assert_almost_equal(res4.fill_value, exp_fv)
                assert_almost_equal(res4.values, exp)

        with np.errstate(all="ignore"):
            for first_arr, second_arr in [(arr1, arr2), (farr1, farr2)]:
                _check_op(op, first_arr, second_arr)

    def test_pickle(self):
        def _check_roundtrip(obj):
            unpickled = tm.round_trip_pickle(obj)
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
        exp = SparseArray([1, -1, -1, 3, -1], fill_value=-1, dtype=np.float64)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([1, np.nan, np.nan, 3, np.nan], fill_value=0)
        res = s.fillna(-1)
        exp = SparseArray([1, -1, -1, 3, -1], fill_value=0, dtype=np.float64)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([1, np.nan, 0, 3, 0])
        res = s.fillna(-1)
        exp = SparseArray([1, -1, 0, 3, 0], fill_value=-1, dtype=np.float64)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([1, np.nan, 0, 3, 0], fill_value=0)
        res = s.fillna(-1)
        exp = SparseArray([1, -1, 0, 3, 0], fill_value=0, dtype=np.float64)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([np.nan, np.nan, np.nan, np.nan])
        res = s.fillna(-1)
        exp = SparseArray([-1, -1, -1, -1], fill_value=-1, dtype=np.float64)
        tm.assert_sp_array_equal(res, exp)

        s = SparseArray([np.nan, np.nan, np.nan, np.nan], fill_value=0)
        res = s.fillna(-1)
        exp = SparseArray([-1, -1, -1, -1], fill_value=0, dtype=np.float64)
        tm.assert_sp_array_equal(res, exp)

        # float dtype's fill_value is np.nan, replaced by -1
        s = SparseArray([0., 0., 0., 0.])
        res = s.fillna(-1)
        exp = SparseArray([0., 0., 0., 0.], fill_value=-1)
        tm.assert_sp_array_equal(res, exp)

        # int dtype shouldn't have missing. No changes.
        s = SparseArray([0, 0, 0, 0])
        assert s.dtype == SparseDtype(np.int64)
        assert s.fill_value == 0
        res = s.fillna(-1)
        tm.assert_sp_array_equal(res, s)

        s = SparseArray([0, 0, 0, 0], fill_value=0)
        assert s.dtype == SparseDtype(np.int64)
        assert s.fill_value == 0
        res = s.fillna(-1)
        exp = SparseArray([0, 0, 0, 0], fill_value=0)
        tm.assert_sp_array_equal(res, exp)

        # fill_value can be nan if there is no missing hole.
        # only fill_value will be changed
        s = SparseArray([0, 0, 0, 0], fill_value=np.nan)
        assert s.dtype == SparseDtype(np.int64, fill_value=np.nan)
        assert np.isnan(s.fill_value)
        res = s.fillna(-1)
        exp = SparseArray([0, 0, 0, 0], fill_value=-1)
        tm.assert_sp_array_equal(res, exp)

    def test_fillna_overlap(self):
        s = SparseArray([1, np.nan, np.nan, 3, np.nan])
        # filling with existing value doesn't replace existing value with
        # fill_value, i.e. existing 3 remains in sp_values
        res = s.fillna(3)
        exp = np.array([1, 3, 3, 3, 3], dtype=np.float64)
        tm.assert_numpy_array_equal(res.to_dense(), exp)

        s = SparseArray([1, np.nan, np.nan, 3, np.nan], fill_value=0)
        res = s.fillna(3)
        exp = SparseArray([1, 3, 3, 3, 3], fill_value=0, dtype=np.float64)
        tm.assert_sp_array_equal(res, exp)


class TestSparseArrayAnalytics(object):

    @pytest.mark.parametrize('data,pos,neg', [
        ([True, True, True], True, False),
        ([1, 2, 1], 1, 0),
        ([1.0, 2.0, 1.0], 1.0, 0.0)
    ])
    def test_all(self, data, pos, neg):
        # GH 17570
        out = SparseArray(data).all()
        assert out

        out = SparseArray(data, fill_value=pos).all()
        assert out

        data[1] = neg
        out = SparseArray(data).all()
        assert not out

        out = SparseArray(data, fill_value=pos).all()
        assert not out

    @pytest.mark.parametrize('data,pos,neg', [
        ([True, True, True], True, False),
        ([1, 2, 1], 1, 0),
        ([1.0, 2.0, 1.0], 1.0, 0.0)
    ])
    @td.skip_if_np_lt_115  # prior didn't dispatch
    def test_numpy_all(self, data, pos, neg):
        # GH 17570
        out = np.all(SparseArray(data))
        assert out

        out = np.all(SparseArray(data, fill_value=pos))
        assert out

        data[1] = neg
        out = np.all(SparseArray(data))
        assert not out

        out = np.all(SparseArray(data, fill_value=pos))
        assert not out

        # raises with a different message on py2.
        msg = "the \'out\' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.all,
                               SparseArray(data), out=np.array([]))

    @pytest.mark.parametrize('data,pos,neg', [
        ([False, True, False], True, False),
        ([0, 2, 0], 2, 0),
        ([0.0, 2.0, 0.0], 2.0, 0.0)
    ])
    def test_any(self, data, pos, neg):
        # GH 17570
        out = SparseArray(data).any()
        assert out

        out = SparseArray(data, fill_value=pos).any()
        assert out

        data[1] = neg
        out = SparseArray(data).any()
        assert not out

        out = SparseArray(data, fill_value=pos).any()
        assert not out

    @pytest.mark.parametrize('data,pos,neg', [
        ([False, True, False], True, False),
        ([0, 2, 0], 2, 0),
        ([0.0, 2.0, 0.0], 2.0, 0.0)
    ])
    @td.skip_if_np_lt_115  # prior didn't dispatch
    def test_numpy_any(self, data, pos, neg):
        # GH 17570
        out = np.any(SparseArray(data))
        assert out

        out = np.any(SparseArray(data, fill_value=pos))
        assert out

        data[1] = neg
        out = np.any(SparseArray(data))
        assert not out

        out = np.any(SparseArray(data, fill_value=pos))
        assert not out

        msg = "the \'out\' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.any,
                               SparseArray(data), out=out)

    def test_sum(self):
        data = np.arange(10).astype(float)
        out = SparseArray(data).sum()
        assert out == 45.0

        data[5] = np.nan
        out = SparseArray(data, fill_value=2).sum()
        assert out == 40.0

        out = SparseArray(data, fill_value=np.nan).sum()
        assert out == 40.0

    def test_numpy_sum(self):
        data = np.arange(10).astype(float)
        out = np.sum(SparseArray(data))
        assert out == 45.0

        data[5] = np.nan
        out = np.sum(SparseArray(data, fill_value=2))
        assert out == 40.0

        out = np.sum(SparseArray(data, fill_value=np.nan))
        assert out == 40.0

        msg = "the 'dtype' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.sum,
                               SparseArray(data), dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.sum,
                               SparseArray(data), out=out)

    @pytest.mark.parametrize("data,expected", [
        (np.array([1, 2, 3, 4, 5], dtype=float),  # non-null data
         SparseArray(np.array([1.0, 3.0, 6.0, 10.0, 15.0]))),
        (np.array([1, 2, np.nan, 4, 5], dtype=float),  # null data
         SparseArray(np.array([1.0, 3.0, np.nan, 7.0, 12.0])))
    ])
    @pytest.mark.parametrize("numpy", [True, False])
    def test_cumsum(self, data, expected, numpy):
        cumsum = np.cumsum if numpy else lambda s: s.cumsum()

        out = cumsum(SparseArray(data))
        tm.assert_sp_array_equal(out, expected)

        out = cumsum(SparseArray(data, fill_value=np.nan))
        tm.assert_sp_array_equal(out, expected)

        out = cumsum(SparseArray(data, fill_value=2))
        tm.assert_sp_array_equal(out, expected)

        if numpy:  # numpy compatibility checks.
            msg = "the 'dtype' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, np.cumsum,
                                   SparseArray(data), dtype=np.int64)

            msg = "the 'out' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, np.cumsum,
                                   SparseArray(data), out=out)
        else:
            axis = 1  # SparseArray currently 1-D, so only axis = 0 is valid.
            msg = "axis\\(={axis}\\) out of bounds".format(axis=axis)
            with tm.assert_raises_regex(ValueError, msg):
                SparseArray(data).cumsum(axis=axis)

    def test_mean(self):
        data = np.arange(10).astype(float)
        out = SparseArray(data).mean()
        assert out == 4.5

        data[5] = np.nan
        out = SparseArray(data).mean()
        assert out == 40.0 / 9

    def test_numpy_mean(self):
        data = np.arange(10).astype(float)
        out = np.mean(SparseArray(data))
        assert out == 4.5

        data[5] = np.nan
        out = np.mean(SparseArray(data))
        assert out == 40.0 / 9

        msg = "the 'dtype' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.mean,
                               SparseArray(data), dtype=np.int64)

        msg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.mean,
                               SparseArray(data), out=out)

    def test_ufunc(self):
        # GH 13853 make sure ufunc is applied to fill_value
        sparse = SparseArray([1, np.nan, 2, np.nan, -2])
        result = SparseArray([1, np.nan, 2, np.nan, 2])
        tm.assert_sp_array_equal(abs(sparse), result)
        tm.assert_sp_array_equal(np.abs(sparse), result)

        sparse = SparseArray([1, -1, 2, -2], fill_value=1)
        result = SparseArray([1, 2, 2], sparse_index=sparse.sp_index,
                             fill_value=1)
        tm.assert_sp_array_equal(abs(sparse), result)
        tm.assert_sp_array_equal(np.abs(sparse), result)

        sparse = SparseArray([1, -1, 2, -2], fill_value=-1)
        result = SparseArray([1, 2, 2], sparse_index=sparse.sp_index,
                             fill_value=1)
        tm.assert_sp_array_equal(abs(sparse), result)
        tm.assert_sp_array_equal(np.abs(sparse), result)

        sparse = SparseArray([1, np.nan, 2, np.nan, -2])
        result = SparseArray(np.sin([1, np.nan, 2, np.nan, -2]))
        tm.assert_sp_array_equal(np.sin(sparse), result)

        sparse = SparseArray([1, -1, 2, -2], fill_value=1)
        result = SparseArray(np.sin([1, -1, 2, -2]), fill_value=np.sin(1))
        tm.assert_sp_array_equal(np.sin(sparse), result)

        sparse = SparseArray([1, -1, 0, -2], fill_value=0)
        result = SparseArray(np.sin([1, -1, 0, -2]), fill_value=np.sin(0))
        tm.assert_sp_array_equal(np.sin(sparse), result)

    def test_ufunc_args(self):
        # GH 13853 make sure ufunc is applied to fill_value, including its arg
        sparse = SparseArray([1, np.nan, 2, np.nan, -2])
        result = SparseArray([2, np.nan, 3, np.nan, -1])
        tm.assert_sp_array_equal(np.add(sparse, 1), result)

        sparse = SparseArray([1, -1, 2, -2], fill_value=1)
        result = SparseArray([2, 0, 3, -1], fill_value=2)
        tm.assert_sp_array_equal(np.add(sparse, 1), result)

        sparse = SparseArray([1, -1, 0, -2], fill_value=0)
        result = SparseArray([2, 0, 1, -1], fill_value=1)
        tm.assert_sp_array_equal(np.add(sparse, 1), result)

    def test_nbytes_integer(self):
        arr = SparseArray([1, 0, 0, 0, 2], kind='integer')
        result = arr.nbytes
        # (2 * 8) + 2 * 4
        assert result == 24

    def test_nbytes_block(self):
        arr = SparseArray([1, 2, 0, 0, 0], kind='block')
        result = arr.nbytes
        # (2 * 8) + 4 + 4
        # sp_values, blocs, blenghts
        assert result == 24

    def test_asarray_datetime64(self):
        s = pd.SparseArray(
            pd.to_datetime(['2012', None, None, '2013'])
        )
        np.asarray(s)


def test_setting_fill_value_fillna_still_works():
    # This is why letting users update fill_value / dtype is bad
    # astype has the same problem.
    arr = SparseArray([1., np.nan, 1.0], fill_value=0.0)
    arr.fill_value = np.nan
    result = arr.isna()
    expected = np.array([False, True, False])
    tm.assert_numpy_array_equal(result, expected)


def test_setting_fill_value_updates():
    arr = SparseArray([0.0, np.nan], fill_value=0)
    arr.fill_value = np.nan
    # use private constructor to get the index right
    # otherwise both nans would be un-stored.
    expected = SparseArray._simple_new(
        sparse_array=np.array([np.nan]),
        sparse_index=IntIndex(2, [1]),
        dtype=SparseDtype(float, np.nan),
    )
    tm.assert_sp_array_equal(arr, expected)


@pytest.mark.parametrize("arr, loc", [
    ([None, 1, 2], 0),
    ([0, None, 2], 1),
    ([0, 1, None], 2),
    ([0, 1, 1, None, None], 3),
    ([1, 1, 1, 2], -1),
    ([], -1),
])
def test_first_fill_value_loc(arr, loc):
    result = SparseArray(arr)._first_fill_value_loc()
    assert result == loc


@pytest.mark.parametrize('arr', [
    [1, 2, np.nan, np.nan],
    [1, np.nan, 2, np.nan],
    [1, 2, np.nan],
])
@pytest.mark.parametrize("fill_value", [
    np.nan, 0, 1
])
def test_unique_na_fill(arr, fill_value):
    a = pd.SparseArray(arr, fill_value=fill_value).unique()
    b = pd.Series(arr).unique()
    assert isinstance(a, SparseArray)
    a = np.asarray(a)
    tm.assert_numpy_array_equal(a, b)
