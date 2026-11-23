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


class ArrayFunctionOverride:
    def __array_function__(self, func, types, args, kwargs):
        return "override"


def test_defers_to_other_array_function():
    arr = pd.array([1, 2], dtype="Int64")
    other = ArrayFunctionOverride()
    assert np.concatenate([arr, other]) == "override"
    assert np.concatenate([other, arr]) == "override"


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
