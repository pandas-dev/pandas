"""
Tests for ExtensionArray.__array_function__ beyond the extension base tests.
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.arrays import (
    NumpyExtensionArray,
    SparseArray,
)
from pandas.tests.extension.array_with_attr import FloatAttrArray


def test_delete_with_extension_array_obj():
    # GH#26380 the ExtensionArray is the indexer, not the array
    result = np.delete(np.arange(3), pd.array([0]))
    tm.assert_numpy_array_equal(result, np.array([1, 2]))

    result = pd.Index([1, 2, 3]).delete(pd.array([0]))
    tm.assert_index_equal(result, pd.Index([2, 3]))


@pytest.mark.parametrize("axis", [0, -1, None])
def test_axis_spellings(axis):
    arr = pd.array([1, None, 3], dtype="Int64")
    result = np.delete(arr, 0, axis=axis)
    tm.assert_extension_array_equal(result, arr[1:])

    result = np.concatenate([arr, arr[:1]], axis=axis)
    tm.assert_extension_array_equal(result, arr.take([0, 1, 2, 0]))


@pytest.mark.parametrize(
    "kwargs",
    [{"dtype": object}, {"axis": None, "dtype": object}, {"casting": "unsafe"}],
)
def test_concatenate_with_kwargs_falls_back(kwargs):
    arr = pd.array([1, None], dtype="Int64")
    result = np.concatenate([arr, arr], **kwargs)
    expected = np.concatenate([np.asarray(arr)] * 2, **kwargs)
    tm.assert_numpy_array_equal(result, expected)


def test_concatenate_out_falls_back():
    arr = pd.array([1, None], dtype="Int64")
    out = np.empty(4, dtype=object)
    result = np.concatenate([arr, arr], out=out)
    assert result is out
    expected = np.concatenate([np.asarray(arr)] * 2, out=np.empty(4, dtype=object))
    tm.assert_numpy_array_equal(out, expected)


@pytest.mark.parametrize(
    "other",
    [
        pd.array([1.5], dtype="Float64"),
        pd.Categorical([1]),
        np.array([1]),
    ],
)
def test_concatenate_mismatched_falls_back(other):
    arr = pd.array([1, None], dtype="Int64")
    result = np.concatenate([arr, other])
    expected = np.concatenate([np.asarray(arr), np.asarray(other)])
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize(
    "left, right",
    [
        (pd.Categorical(["a"]), pd.Categorical(["b"])),
        (
            pd.array(pd.date_range("2000", periods=2)),
            pd.array(pd.date_range("2000", periods=2, tz="UTC")),
        ),
    ],
)
def test_concatenate_same_type_different_dtype_falls_back(left, right):
    result = np.concatenate([left, right])
    expected = np.concatenate([np.asarray(left), np.asarray(right)])
    tm.assert_numpy_array_equal(result, expected)


def test_2d_falls_back():
    arr = pd.date_range("2000", periods=4)._data.reshape(2, 2)
    expected = np.asarray(arr)

    result = np.delete(arr, 0, axis=0)
    tm.assert_numpy_array_equal(result, np.delete(expected, 0, axis=0))

    result = np.concatenate([arr, arr])
    tm.assert_numpy_array_equal(result, np.concatenate([expected, expected]))

    cond = np.array([[True, False], [False, True]])
    result = np.where(cond, arr, pd.NaT)
    tm.assert_numpy_array_equal(result, np.where(cond, expected, pd.NaT))


class ArrayFunctionOverride:
    def __array_function__(self, func, types, args, kwargs):
        return "override"


def test_defers_to_other_array_function():
    arr = pd.array([1, 2], dtype="Int64")
    other = ArrayFunctionOverride()
    assert np.concatenate([arr, other]) == "override"
    assert np.concatenate([other, arr]) == "override"
    assert np.where(np.array([True, False]), arr, other) == "override"
    assert np.where(np.array([True, False]), other, arr) == "override"


def test_defers_to_extension_array_subclass_override():
    class SubclassOverride(NumpyExtensionArray):
        def __array_function__(self, func, types, args, kwargs):
            return "override"

    arr = pd.array([1, 2], dtype="Int64")
    other = SubclassOverride(np.array([3, 4]))
    assert np.concatenate([arr, other]) == "override"


def test_subclass_override_delegating_to_super():
    class DelegatingOverride(NumpyExtensionArray):
        def __array_function__(self, func, types, args, kwargs):
            return super().__array_function__(func, types, args, kwargs)

    arr = DelegatingOverride(np.array([1, 2, 3]))
    assert np.mean(arr) == 2.0
    result = np.delete(arr, 0)
    tm.assert_extension_array_equal(result, DelegatingOverride(np.array([2, 3])))


@pytest.mark.parametrize(
    "left, right",
    [
        (pd.array([3, None, 2], dtype="Int64"), pd.array([2, None], dtype="Int64")),
        (SparseArray([3, 1, 2]), SparseArray([2, 5])),
    ],
)
def test_nested_numpy_calls_unchanged(left, right):
    # GH#26380 np.setxor1d calls np.concatenate internally, which must keep
    #  returning an ndarray there
    result = np.setxor1d(left, right, assume_unique=True)
    expected = np.setxor1d(np.asarray(left), np.asarray(right), assume_unique=True)
    tm.assert_numpy_array_equal(result, expected)


def test_concat_same_type_calling_np_concatenate():
    # GH#26380 a subclass implementing _concat_same_type with np.concatenate
    #  must not recurse endlessly, whether called via numpy or directly
    class ConcatWithNumpy(NumpyExtensionArray):
        @classmethod
        def _concat_same_type(cls, to_concat):
            return cls(np.concatenate(to_concat))

    arr = ConcatWithNumpy(np.arange(3))
    expected = ConcatWithNumpy(np.array([0, 1, 2, 0, 1, 2]))

    result = np.concatenate([arr, arr])
    tm.assert_extension_array_equal(result, expected)

    result = ConcatWithNumpy._concat_same_type([arr, arr])
    tm.assert_extension_array_equal(result, expected)


def test_where_series_condition():
    # GH#46594
    left = pd.Series(pd.array([1, 2, 3, 4], dtype="Int64"))
    right = pd.Series(pd.array([2, 1, None, 8], dtype="Int64"))
    cond = (left > right).fillna(False)
    expected = pd.array([2, 2, None, 8], dtype="Int64")

    result = np.where(cond, left.array, right.array)
    tm.assert_extension_array_equal(result, expected)

    # a Series other is aligned by position, not label
    result = np.where(cond, left.array, right.set_axis([3, 2, 1, 0]))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("other", [None, np.nan, pd.NA])
def test_where_na_scalar(other):
    arr = pd.array([1, 2, 3], dtype="Int64")
    cond = np.array([True, False, True])
    expected = pd.array([1, None, 3], dtype="Int64")

    result = np.where(cond, arr, other)
    tm.assert_extension_array_equal(result, expected)

    result = np.where(~cond, other, arr)
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "x, y",
    [
        # value not holdable by the dtype
        (pd.array([1, 2, 3], dtype="Int64"), 1.5),
        (pd.Categorical(["a", "b", "a"]), "c"),
        (pd.array([1, 2, 3], dtype="Int8"), 500.0),
        # numpy-typed values that numpy would not cast safely
        (pd.array([1, 2, 3], dtype="Int8"), np.int64(500)),
        (pd.array([1, 2, 3], dtype="Float32"), np.array([0.1, 0.2, 0.3])),
        (pd.array([1, 2, 3], dtype="Int64"), np.array([4.0, 5.0, 6.0])),
        # extension dtypes differ
        (pd.array([1, 2, 3], dtype="Int64"), pd.array([1.5, 2, 3], dtype="Float64")),
        (
            pd.array([1, 2, 3], dtype="Int64"),
            pd.Series(pd.array([1.5, 2, 3], dtype="Float64")),
        ),
        # broadcasting
        (pd.array([1, 2, 3], dtype="Int64"), pd.array([4], dtype="Int64")),
        # list-likes other than arrays
        (pd.array([1, 2, 3], dtype="Int64"), [4, 5, 6]),
        # _where upcasts instead of raising
        (SparseArray([1, 0, 2]), 1.5),
        # no __setitem__
        (FloatAttrArray(np.array([1.0, 2.0, 3.0])), 0.0),
    ],
)
def test_where_falls_back(x, y):
    cond = np.array([True, False, True])
    result = np.where(cond, x, y)
    expected = np.where(cond, np.asarray(x), np.asarray(y))
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type, values, other",
    [
        ("bool_", [True, False, True], np.array([0.5, 0.5, 0.5])),
        ("string", ["a", "b", "c"], np.array([7, 8, 9])),
    ],
)
def test_where_arrow_with_ndarray_falls_back(pa_type, values, other):
    # setitem would cast other lossily
    pa = pytest.importorskip("pyarrow")
    arr = pd.array(values, dtype=pd.ArrowDtype(getattr(pa, pa_type)()))
    cond = np.array([True, False, True])
    result = np.where(cond, arr, other)
    expected = np.where(cond, np.asarray(arr), other)
    tm.assert_numpy_array_equal(result, expected)


def test_where_boolean_extension_array_condition():
    arr = pd.array([1, 2, 3], dtype="Int64")
    cond = pd.array([True, False, True], dtype="boolean")
    result = np.where(cond, arr, pd.array([4, 5, 6], dtype="Int64"))
    tm.assert_extension_array_equal(result, pd.array([1, 5, 3], dtype="Int64"))


@pytest.mark.parametrize("cond", [np.array([1, 0, 1]), np.bool_(True)])
def test_where_non_bool_or_broadcast_condition_falls_back(cond):
    arr = pd.array([1, 2, 3], dtype="Int64")
    other = pd.array([4, 5, 6], dtype="Int64")
    result = np.where(cond, arr, other)
    expected = np.where(cond, np.asarray(arr), np.asarray(other))
    tm.assert_numpy_array_equal(result, expected)


def test_where_condition_with_na_raises():
    arr = pd.array([1, 2, 3], dtype="Int64")
    cond = pd.array([True, None, False], dtype="boolean")
    with pytest.raises(TypeError, match="boolean value of NA is ambiguous"):
        np.where(cond, arr, arr)


def test_where_condition_only():
    # np.where(cond) is np.nonzero, unchanged
    result = np.where(pd.array([True, False, True], dtype="boolean"))
    tm.assert_numpy_array_equal(result[0], np.array([0, 2]))


@pytest.mark.parametrize(
    "other",
    [SparseArray([4, 5, 6], fill_value=5), pd.array([4, 5, 6], dtype="Int64")],
)
def test_sparse_where_extension_array_other(other):
    # GH#26380 SparseArray._where calls np.where on other internally
    ser = pd.Series(SparseArray([1, 0, 2]))
    result = ser.where(np.array([True, False, True]), pd.Series(other))
    expected = pd.Series(SparseArray([1, 5, 2]))
    tm.assert_series_equal(result, expected)


def test_sparse_fillna_extension_array_value():
    # GH#26380 SparseArray.fillna calls np.where on value internally
    arr = SparseArray([np.nan, 1.0, np.nan], fill_value=0)
    result = arr.fillna(pd.array([5.0, 6.0, 7.0], dtype="Float64"))
    tm.assert_numpy_array_equal(result.to_dense(), np.array([5.0, 1.0, 7.0]))


def test_numpy_extension_array_where_extension_array_value():
    # GH#26380 NDArrayBackedExtensionArray._where calls np.where on value
    arr = NumpyExtensionArray(np.arange(4, dtype=np.int64))
    mask = np.array([True, False, True, False])
    result = arr._where(mask, pd.array([4, 5, 6, 7], dtype="Int64"))
    expected = NumpyExtensionArray(np.array([0, 5, 2, 7], dtype=np.int64))
    tm.assert_extension_array_equal(result, expected)
