import operator

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core import roperator
from pandas.core.arrays.sparse import SparseArray


@pytest.fixture(params=["integer", "block"])
def kind(request):
    """kind kwarg to pass to SparseArray"""
    return request.param


@pytest.fixture(params=[True, False])
def mix(request):
    """
    Fixture returning True or False, determining whether to operate
    op(sparse, dense) instead of op(sparse, sparse)
    """
    return request.param


class TestSparseArrayArithmetics:
    def _assert(self, a, b):
        # We have to use tm.assert_sp_array_equal. See GH #45126
        tm.assert_numpy_array_equal(a, b)

    def _check_numeric_ops(self, a, b, a_dense, b_dense, mix: bool, op):
        # Check that arithmetic behavior matches non-Sparse Series arithmetic

        if isinstance(a_dense, np.ndarray):
            expected = op(pd.Series(a_dense), b_dense).values
        elif isinstance(b_dense, np.ndarray):
            expected = op(a_dense, pd.Series(b_dense)).values
        else:
            raise NotImplementedError

        with np.errstate(invalid="ignore", divide="ignore"):
            if mix:
                result = op(a, b_dense).to_dense()
            else:
                result = op(a, b).to_dense()

        self._assert(result, expected)

    def _check_bool_result(self, res):
        assert isinstance(res, SparseArray)
        assert isinstance(res.dtype, pd.SparseDtype)
        assert res.dtype.subtype == np.bool_
        assert isinstance(res.fill_value, bool)

    def _check_comparison_ops(self, a, b, a_dense, b_dense):
        with np.errstate(invalid="ignore"):
            # Unfortunately, trying to wrap the computation of each expected
            # value is with np.errstate() is too tedious.
            #
            # sparse & sparse
            self._check_bool_result(a == b)
            self._assert((a == b).to_dense(), a_dense == b_dense)

            self._check_bool_result(a != b)
            self._assert((a != b).to_dense(), a_dense != b_dense)

            self._check_bool_result(a >= b)
            self._assert((a >= b).to_dense(), a_dense >= b_dense)

            self._check_bool_result(a <= b)
            self._assert((a <= b).to_dense(), a_dense <= b_dense)

            self._check_bool_result(a > b)
            self._assert((a > b).to_dense(), a_dense > b_dense)

            self._check_bool_result(a < b)
            self._assert((a < b).to_dense(), a_dense < b_dense)

            # sparse & dense
            self._check_bool_result(a == b_dense)
            self._assert((a == b_dense).to_dense(), a_dense == b_dense)

            self._check_bool_result(a != b_dense)
            self._assert((a != b_dense).to_dense(), a_dense != b_dense)

            self._check_bool_result(a >= b_dense)
            self._assert((a >= b_dense).to_dense(), a_dense >= b_dense)

            self._check_bool_result(a <= b_dense)
            self._assert((a <= b_dense).to_dense(), a_dense <= b_dense)

            self._check_bool_result(a > b_dense)
            self._assert((a > b_dense).to_dense(), a_dense > b_dense)

            self._check_bool_result(a < b_dense)
            self._assert((a < b_dense).to_dense(), a_dense < b_dense)

    def _check_logical_ops(self, a, b, a_dense, b_dense):
        # sparse & sparse
        self._check_bool_result(a & b)
        self._assert((a & b).to_dense(), a_dense & b_dense)

        self._check_bool_result(a | b)
        self._assert((a | b).to_dense(), a_dense | b_dense)
        # sparse & dense
        self._check_bool_result(a & b_dense)
        self._assert((a & b_dense).to_dense(), a_dense & b_dense)

        self._check_bool_result(a | b_dense)
        self._assert((a | b_dense).to_dense(), a_dense | b_dense)

    @pytest.mark.parametrize("scalar", [0, 1, 3])
    @pytest.mark.parametrize("fill_value", [None, 0, 2])
    def test_float_scalar(
        self, kind, mix, all_arithmetic_functions, fill_value, scalar, request
    ):
        op = all_arithmetic_functions
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        a = SparseArray(values, kind=kind, fill_value=fill_value)
        self._check_numeric_ops(a, scalar, values, scalar, mix, op)

    def test_float_scalar_comparison(self, kind):
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])

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

    def test_float_same_index_without_nans(self, kind, mix, all_arithmetic_functions):
        # when sp_index are the same
        op = all_arithmetic_functions

        values = np.array([0.0, 1.0, 2.0, 6.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0])
        rvalues = np.array([0.0, 2.0, 3.0, 4.0, 0.0, 0.0, 1.0, 3.0, 2.0, 0.0])

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind, fill_value=0)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

    def test_float_same_index_with_nans(
        self, kind, mix, all_arithmetic_functions, request
    ):
        # when sp_index are the same
        op = all_arithmetic_functions
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([np.nan, 2, 3, 4, np.nan, 0, 1, 3, 2, np.nan])

        a = SparseArray(values, kind=kind)
        b = SparseArray(rvalues, kind=kind)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

    def test_float_same_index_comparison(self, kind):
        # when sp_index are the same
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([np.nan, 2, 3, 4, np.nan, 0, 1, 3, 2, np.nan])

        a = SparseArray(values, kind=kind)
        b = SparseArray(rvalues, kind=kind)
        self._check_comparison_ops(a, b, values, rvalues)

        values = np.array([0.0, 1.0, 2.0, 6.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0])
        rvalues = np.array([0.0, 2.0, 3.0, 4.0, 0.0, 0.0, 1.0, 3.0, 2.0, 0.0])

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind, fill_value=0)
        self._check_comparison_ops(a, b, values, rvalues)

    def test_float_array(self, kind, mix, all_arithmetic_functions):
        op = all_arithmetic_functions

        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, np.nan, 2, 3, np.nan, 0, 1, 5, 2, np.nan])

        a = SparseArray(values, kind=kind)
        b = SparseArray(rvalues, kind=kind)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)
        self._check_numeric_ops(a, b * 0, values, rvalues * 0, mix, op)

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind, fill_value=0)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, kind=kind, fill_value=1)
        b = SparseArray(rvalues, kind=kind, fill_value=2)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

    def test_float_array_different_kind(self, mix, all_arithmetic_functions):
        op = all_arithmetic_functions

        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, np.nan, 2, 3, np.nan, 0, 1, 5, 2, np.nan])

        a = SparseArray(values, kind="integer")
        b = SparseArray(rvalues, kind="block")
        self._check_numeric_ops(a, b, values, rvalues, mix, op)
        self._check_numeric_ops(a, b * 0, values, rvalues * 0, mix, op)

        a = SparseArray(values, kind="integer", fill_value=0)
        b = SparseArray(rvalues, kind="block")
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, kind="integer", fill_value=0)
        b = SparseArray(rvalues, kind="block", fill_value=0)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, kind="integer", fill_value=1)
        b = SparseArray(rvalues, kind="block", fill_value=2)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

    def test_float_array_comparison(self, kind):
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, np.nan, 2, 3, np.nan, 0, 1, 5, 2, np.nan])

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

    def test_int_array(self, kind, mix, all_arithmetic_functions):
        op = all_arithmetic_functions

        # have to specify dtype explicitly until fixing GH 667
        dtype = np.int64

        values = np.array([0, 1, 2, 0, 0, 0, 1, 2, 1, 0], dtype=dtype)
        rvalues = np.array([2, 0, 2, 3, 0, 0, 1, 5, 2, 0], dtype=dtype)

        a = SparseArray(values, dtype=dtype, kind=kind)
        assert a.dtype == pd.SparseDtype(dtype)
        b = SparseArray(rvalues, dtype=dtype, kind=kind)
        assert b.dtype == pd.SparseDtype(dtype)

        self._check_numeric_ops(a, b, values, rvalues, mix, op)
        self._check_numeric_ops(a, b * 0, values, rvalues * 0, mix, op)

        a = SparseArray(values, fill_value=0, dtype=dtype, kind=kind)
        assert a.dtype == pd.SparseDtype(dtype)
        b = SparseArray(rvalues, dtype=dtype, kind=kind)
        assert b.dtype == pd.SparseDtype(dtype)

        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, fill_value=0, dtype=dtype, kind=kind)
        assert a.dtype == pd.SparseDtype(dtype)
        b = SparseArray(rvalues, fill_value=0, dtype=dtype, kind=kind)
        assert b.dtype == pd.SparseDtype(dtype)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, fill_value=1, dtype=dtype, kind=kind)
        assert a.dtype == pd.SparseDtype(dtype, fill_value=1)
        b = SparseArray(rvalues, fill_value=2, dtype=dtype, kind=kind)
        assert b.dtype == pd.SparseDtype(dtype, fill_value=2)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

    def test_int_array_comparison(self, kind):
        dtype = "int64"
        # int32 NI ATM

        values = np.array([0, 1, 2, 0, 0, 0, 1, 2, 1, 0], dtype=dtype)
        rvalues = np.array([2, 0, 2, 3, 0, 0, 1, 5, 2, 0], dtype=dtype)

        a = SparseArray(values, dtype=dtype, kind=kind)
        b = SparseArray(rvalues, dtype=dtype, kind=kind)
        self._check_comparison_ops(a, b, values, rvalues)
        self._check_comparison_ops(a, b * 0, values, rvalues * 0)

        a = SparseArray(values, dtype=dtype, kind=kind, fill_value=0)
        b = SparseArray(rvalues, dtype=dtype, kind=kind)
        self._check_comparison_ops(a, b, values, rvalues)

        a = SparseArray(values, dtype=dtype, kind=kind, fill_value=0)
        b = SparseArray(rvalues, dtype=dtype, kind=kind, fill_value=0)
        self._check_comparison_ops(a, b, values, rvalues)

        a = SparseArray(values, dtype=dtype, kind=kind, fill_value=1)
        b = SparseArray(rvalues, dtype=dtype, kind=kind, fill_value=2)
        self._check_comparison_ops(a, b, values, rvalues)

    @pytest.mark.parametrize("fill_value", [True, False, np.nan])
    def test_bool_same_index(self, kind, fill_value):
        # GH 14000
        # when sp_index are the same
        values = np.array([True, False, True, True], dtype=np.bool_)
        rvalues = np.array([True, False, True, True], dtype=np.bool_)

        a = SparseArray(values, kind=kind, dtype=np.bool_, fill_value=fill_value)
        b = SparseArray(rvalues, kind=kind, dtype=np.bool_, fill_value=fill_value)
        self._check_logical_ops(a, b, values, rvalues)

    @pytest.mark.parametrize("fill_value", [True, False, np.nan])
    def test_bool_array_logical(self, kind, fill_value):
        # GH 14000
        # when sp_index are the same
        values = np.array([True, False, True, False, True, True], dtype=np.bool_)
        rvalues = np.array([True, False, False, True, False, True], dtype=np.bool_)

        a = SparseArray(values, kind=kind, dtype=np.bool_, fill_value=fill_value)
        b = SparseArray(rvalues, kind=kind, dtype=np.bool_, fill_value=fill_value)
        self._check_logical_ops(a, b, values, rvalues)

    def test_mixed_array_float_int(self, kind, mix, all_arithmetic_functions, request):
        op = all_arithmetic_functions
        rdtype = "int64"
        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, 0, 2, 3, 0, 0, 1, 5, 2, 0], dtype=rdtype)

        a = SparseArray(values, kind=kind)
        b = SparseArray(rvalues, kind=kind)
        assert b.dtype == pd.SparseDtype(rdtype)

        self._check_numeric_ops(a, b, values, rvalues, mix, op)
        self._check_numeric_ops(a, b * 0, values, rvalues * 0, mix, op)

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind)
        assert b.dtype == pd.SparseDtype(rdtype)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind, fill_value=0)
        assert b.dtype == pd.SparseDtype(rdtype)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

        a = SparseArray(values, kind=kind, fill_value=1)
        b = SparseArray(rvalues, kind=kind, fill_value=2)
        assert b.dtype == pd.SparseDtype(rdtype, fill_value=2)
        self._check_numeric_ops(a, b, values, rvalues, mix, op)

    def test_mixed_array_comparison(self, kind):
        rdtype = "int64"
        # int32 NI ATM

        values = np.array([np.nan, 1, 2, 0, np.nan, 0, 1, 2, 1, np.nan])
        rvalues = np.array([2, 0, 2, 3, 0, 0, 1, 5, 2, 0], dtype=rdtype)

        a = SparseArray(values, kind=kind)
        b = SparseArray(rvalues, kind=kind)
        assert b.dtype == pd.SparseDtype(rdtype)

        self._check_comparison_ops(a, b, values, rvalues)
        self._check_comparison_ops(a, b * 0, values, rvalues * 0)

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind)
        assert b.dtype == pd.SparseDtype(rdtype)
        self._check_comparison_ops(a, b, values, rvalues)

        a = SparseArray(values, kind=kind, fill_value=0)
        b = SparseArray(rvalues, kind=kind, fill_value=0)
        assert b.dtype == pd.SparseDtype(rdtype)
        self._check_comparison_ops(a, b, values, rvalues)

        a = SparseArray(values, kind=kind, fill_value=1)
        b = SparseArray(rvalues, kind=kind, fill_value=2)
        assert b.dtype == pd.SparseDtype(rdtype, fill_value=2)
        self._check_comparison_ops(a, b, values, rvalues)

    def test_xor(self):
        s = SparseArray([True, True, False, False])
        t = SparseArray([True, False, True, False])
        result = s ^ t
        sp_index = pd.core.arrays.sparse.IntIndex(4, np.array([0, 1, 2], dtype="int32"))
        expected = SparseArray([False, True, True], sparse_index=sp_index)
        tm.assert_sp_array_equal(result, expected)


@pytest.mark.parametrize("op", [operator.eq, operator.add])
def test_with_list(op):
    arr = SparseArray([0, 1], fill_value=0)
    result = op(arr, [0, 1])
    expected = op(arr, SparseArray([0, 1]))
    tm.assert_sp_array_equal(result, expected)


def test_with_dataframe():
    # GH#27910
    arr = SparseArray([0, 1], fill_value=0)
    df = pd.DataFrame([[1, 2], [3, 4]])
    result = arr.__add__(df)
    assert result is NotImplemented


def test_with_zerodim_ndarray():
    # GH#27910
    arr = SparseArray([0, 1], fill_value=0)

    result = arr * np.array(2)
    expected = arr * 2
    tm.assert_sp_array_equal(result, expected)


@pytest.mark.parametrize("ufunc", [np.abs, np.exp])
@pytest.mark.parametrize(
    "arr", [SparseArray([0, 0, -1, 1]), SparseArray([None, None, -1, 1])]
)
def test_ufuncs(ufunc, arr):
    result = ufunc(arr)
    fill_value = ufunc(arr.fill_value)
    expected = SparseArray(ufunc(np.asarray(arr)), fill_value=fill_value)
    tm.assert_sp_array_equal(result, expected)


@pytest.mark.parametrize(
    "a, b",
    [
        (SparseArray([0, 0, 0]), np.array([0, 1, 2])),
        (SparseArray([0, 0, 0], fill_value=1), np.array([0, 1, 2])),
    ],
)
@pytest.mark.parametrize("ufunc", [np.add, np.greater])
def test_binary_ufuncs(ufunc, a, b):
    # can't say anything about fill value here.
    result = ufunc(a, b)
    expected = ufunc(np.asarray(a), np.asarray(b))
    assert isinstance(result, SparseArray)
    tm.assert_numpy_array_equal(np.asarray(result), expected)


def test_ndarray_inplace():
    sparray = SparseArray([0, 2, 0, 0])
    ndarray = np.array([0, 1, 2, 3])
    ndarray += sparray
    expected = np.array([0, 3, 2, 3])
    tm.assert_numpy_array_equal(ndarray, expected)


def test_sparray_inplace():
    sparray = SparseArray([0, 2, 0, 0])
    ndarray = np.array([0, 1, 2, 3])
    sparray += ndarray
    expected = SparseArray([0, 3, 2, 3], fill_value=0)
    tm.assert_sp_array_equal(sparray, expected)


@pytest.mark.parametrize("cons", [list, np.array, SparseArray])
def test_mismatched_length_cmp_op(cons):
    left = SparseArray([True, True])
    right = cons([True, True, True])
    with pytest.raises(ValueError, match="operands have mismatched length"):
        left & right


@pytest.mark.parametrize("op", [operator.and_, operator.or_, operator.xor])
def test_logical_op_uneven_length_series(op):
    # GH#32119 a logical op between a sparse Series and a differently-indexed
    #  Series aligns (introducing NA) and must match the non-sparse result
    #  instead of raising an uninformative error.
    sparse = pd.Series(SparseArray(np.arange(10)))
    dense = pd.Series(np.arange(11))

    result = op(sparse == 5, dense == 5)
    expected = op(pd.Series(np.arange(10)) == 5, dense == 5)

    assert isinstance(result.dtype, pd.SparseDtype)
    tm.assert_series_equal(result.astype(bool), expected)


@pytest.mark.parametrize("op", [operator.and_, operator.or_, operator.xor])
@pytest.mark.parametrize("subtype", [bool, object])
def test_logical_op_masked_other(op, subtype):
    # GH#68483 a masked operand keeps its own Kleene semantics; densifying it
    #  lost them -- raising on the NA for a bool subtype, silently resolving it
    #  to False for an object one
    values = np.array([True, True, False, False], dtype=subtype)
    other = pd.array([True, pd.NA, True, False], dtype="boolean")

    result = op(SparseArray(values), other)
    expected = op(values, other)
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("op", [operator.and_, operator.or_, operator.xor])
def test_logical_op_masked_other_without_na(op):
    # GH#68483 the operand's dtype decides, not whether it holds NA, so the
    #  result is the masked dtype either way -- as it is for the dense operand
    values = np.array([True, True, False, False])
    other = pd.array([True, False, True, False], dtype="boolean")

    result = op(SparseArray(values), other)
    expected = op(values, other)
    assert result.dtype == pd.BooleanDtype()
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("op", [operator.and_, operator.or_, operator.xor])
def test_logical_op_non_boolean_masked_other(op):
    # GH#68483 only a masked *boolean* operand is deferred to; the other masked
    #  dtypes reach _arith_method, which cannot consume a SparseArray
    # int64 explicitly: on 32-bit the default int would not match the int64 result
    values = np.array([1, 0, 3, 0], dtype="int64")
    other = pd.array([1, 2, 3, 4], dtype="Int64")

    result = op(SparseArray(values), other)
    expected = op(values, other.to_numpy(dtype="int64"))
    assert isinstance(result.dtype, pd.SparseDtype)
    tm.assert_numpy_array_equal(result.to_dense(), expected)


@pytest.mark.parametrize("op", [operator.and_, operator.or_, operator.xor])
@pytest.mark.parametrize("n_sparse, n_other", [(11, 10), (10, 11)])
def test_logical_op_masked_other_uneven_length_series(op, n_sparse, n_other):
    # GH#68483 the alignment path onto the above; either operand can be the one
    #  that gains the NA, and only the shorter-sparse case upcasts it to object
    sparse = pd.Series(SparseArray(np.arange(n_sparse)))
    other = pd.Series(np.arange(n_other) == 5, dtype="boolean")

    result = op(sparse == 5, other)
    expected = op(pd.Series(np.arange(n_sparse)) == 5, other)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "a, b",
    [
        ([0, 1, 2], [0, 1, 2, 3]),
        ([0, 1, 2, 3], [0, 1, 2]),
    ],
)
def test_mismatched_length_arith_op(a, b, all_arithmetic_functions):
    op = all_arithmetic_functions
    with pytest.raises(AssertionError, match=f"length mismatch: {len(a)} vs. {len(b)}"):
        op(SparseArray(a, fill_value=0), np.array(b))


@pytest.mark.parametrize("op", ["add", "sub", "mul", "truediv", "floordiv", "pow"])
@pytest.mark.parametrize("fill_value", [np.nan, 3])
def test_binary_operators(op, fill_value):
    op = getattr(operator, op)
    data1 = np.random.default_rng(2).standard_normal(20)
    data2 = np.random.default_rng(2).standard_normal(20)

    data1[::2] = fill_value
    data2[::3] = fill_value

    first = SparseArray(data1, fill_value=fill_value)
    second = SparseArray(data2, fill_value=fill_value)

    with np.errstate(all="ignore"):
        res = op(first, second)
        exp = SparseArray(
            op(first.to_dense(), second.to_dense()), fill_value=first.fill_value
        )
        assert isinstance(res, SparseArray)
        tm.assert_almost_equal(res.to_dense(), exp.to_dense())

        res2 = op(first, second.to_dense())
        assert isinstance(res2, SparseArray)
        tm.assert_sp_array_equal(res, res2)

        res3 = op(first.to_dense(), second)
        assert isinstance(res3, SparseArray)
        tm.assert_sp_array_equal(res, res3)

        res4 = op(first, 4)
        assert isinstance(res4, SparseArray)

        # Ignore this if the actual op raises (e.g. pow).
        try:
            exp = op(first.to_dense(), 4)
            exp_fv = op(first.fill_value, 4)
        except ValueError:
            pass
        else:
            tm.assert_almost_equal(res4.fill_value, exp_fv)
            tm.assert_almost_equal(res4.to_dense(), exp)


@pytest.mark.parametrize("values", [[1, 2, 3], [True, False, True]])
@pytest.mark.parametrize("op", [operator.add, operator.gt])
@pytest.mark.parametrize("other_kind", ["datetime64", "timedelta64"])
def test_datetimelike_operand_raises_typeerror(values, op, other_kind):
    # GH#68466 the operand was coerced to self.fill_value first, so SparseDtype
    #  rejected the fill_value before the op could reject the operand
    if other_kind == "datetime64":
        other = np.asarray(pd.date_range("2016", periods=3))
    else:
        other = np.asarray(pd.timedelta_range("1 Day", periods=3))
    with pytest.raises(TypeError):
        op(SparseArray(values), other)


def test_mul_timedelta64_operand():
    # GH#68466 valid densely, but the fill_value coercion raised
    tda = np.asarray(pd.timedelta_range("1 Day", periods=3))
    result = SparseArray([1, 2, 3]) * tda
    tm.assert_numpy_array_equal(result.to_dense(), np.arange(1, 4) * tda)


_NO_COMMON_SUBTYPE_VALUES = {
    "i8": np.array([1, 2, 3]),
    "f8": np.array([1.0, np.nan, 3.0]),
    "bool": np.array([True, False, True]),
    "tda": np.asarray(pd.timedelta_range("1 Day", periods=3)),
    "dti": np.asarray(pd.date_range("2016", periods=3)),
}


@pytest.mark.parametrize(
    "lkind, rkind, op",
    [
        ("i8", "tda", operator.mul),
        ("f8", "tda", operator.mul),
        ("i8", "tda", roperator.rmul),
        ("i8", "bool", operator.mul),
        ("i8", "bool", operator.or_),
        ("i8", "bool", operator.xor),
        ("tda", "i8", operator.truediv),
        ("tda", "f8", operator.floordiv),
        ("tda", "dti", operator.add),
        ("dti", "tda", operator.sub),
    ],
)
def test_no_common_subtype_matches_dense(lkind, rkind, op):
    # GH#68562 int64 and m8[us] have no common subtype, so both operands were
    #  cast to object and the result kept that dtype
    left = _NO_COMMON_SUBTYPE_VALUES[lkind]
    right = _NO_COMMON_SUBTYPE_VALUES[rkind]

    result = op(SparseArray(left), right)
    expected = op(pd.Series(left), pd.Series(right))
    assert result.dtype == pd.SparseDtype(expected.dtype)
    tm.assert_numpy_array_equal(result.to_dense(), expected.to_numpy())


@pytest.mark.parametrize("op", [operator.eq, operator.ne])
def test_no_common_subtype_comparison_matches_dense(op):
    # GH#68562 a comparison keeps its own fill value, so only the subtype is
    #  pinned to the dense result
    tda = _NO_COMMON_SUBTYPE_VALUES["tda"]
    i8 = _NO_COMMON_SUBTYPE_VALUES["i8"]

    result = op(SparseArray(tda), i8)
    expected = op(pd.Series(tda), pd.Series(i8))
    assert result.dtype.subtype == expected.dtype
    tm.assert_numpy_array_equal(result.to_dense(), expected.to_numpy())


@pytest.mark.parametrize("op", [divmod, roperator.rdivmod])
def test_no_common_subtype_divmod(op):
    # GH#68562 divmod returns a 2-tuple, which the reflected name used to miss
    left = _NO_COMMON_SUBTYPE_VALUES["i8" if op is roperator.rdivmod else "tda"]
    right = _NO_COMMON_SUBTYPE_VALUES["tda" if op is roperator.rdivmod else "i8"]

    result = op(SparseArray(left), right)
    expected = op(pd.Series(left), pd.Series(right))
    for res, exp in zip(result, expected, strict=True):
        assert res.dtype == pd.SparseDtype(exp.dtype)
        tm.assert_numpy_array_equal(res.to_dense(), exp.to_numpy())


def test_no_common_subtype_fill_value_matches_subtype():
    # GH#68562 the float operand contributes an np.nan fill value, which is not
    #  the flavor of NA a timedelta64 result holds
    tda = _NO_COMMON_SUBTYPE_VALUES["tda"]
    result = SparseArray([1.0, np.nan, 3.0]) * tda

    assert isinstance(result.fill_value, np.timedelta64)
    assert isinstance(result[1], np.timedelta64)


def test_no_common_subtype_preserves_index_kind():
    # GH#68562 the dense path re-sparsifies from scratch, so it has to be told
    #  which kind of index to rebuild
    tda = _NO_COMMON_SUBTYPE_VALUES["tda"]
    result = SparseArray([1, 2, 3], kind="block") * tda
    assert result.kind == "block"


def test_no_common_subtype_na_fill_value_keeps_object():
    # GH#68562 a pd.NA fill escapes _get_fill's ValueError fallback, so such an
    #  operand keeps the object-cast path rather than raising
    arr = SparseArray([1, 2, 3], fill_value=pd.NA)

    result = arr * np.array([True, False, True])
    assert result.dtype.subtype == np.dtype(object)
    assert result.fill_value is pd.NA
    tm.assert_numpy_array_equal(result.to_dense(), np.array([1, 0, 3], dtype=object))


def test_no_common_subtype_both_with_gaps():
    # GH#68562 neither operand is dense, so this used to reach the splib kernels
    #  and fail on the missing sparse_mul_object
    left = SparseArray([1, 0, 3, 4])
    right = SparseArray(np.array([1, 2, "NaT", 4], dtype="m8[us]"))

    result = left * right
    expected = pd.Series(left.to_dense()) * pd.Series(right.to_dense())
    assert result.dtype == pd.SparseDtype(expected.dtype)
    tm.assert_numpy_array_equal(result.to_dense(), expected.to_numpy())
