"""
Additional tests for PandasArray that aren't covered by
the interface tests.
"""
import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas.core.arrays import PandasArray
import pandas.util.testing as tm


def test_to_numpy():
    arr = PandasArray(np.array([1, 2, 3]))
    result = arr.to_numpy()
    assert result is arr._ndarray

    result = arr.to_numpy(copy=True)
    assert result is not arr._ndarray

    result = arr.to_numpy(dtype='f8')
    expected = np.array([1, 2, 3], dtype='f8')
    tm.assert_numpy_array_equal(result, expected)


def test_setitem():
    ser = pd.Series([1, 2, 3])
    ser.array[0] = 10
    expected = pd.Series([10, 2, 3])
    tm.assert_series_equal(ser, expected)


def test_bad_reduce_raises():
    arr = np.array([1, 2, 3], dtype='int64')
    arr = PandasArray(arr)
    msg = "cannot perform not_a_method with type int"
    with pytest.raises(TypeError, match=msg):
        arr._reduce(msg)


def test_from_sequence_dtype():
    arr = np.array([1, 2, 3], dtype='int64')
    result = PandasArray._from_sequence(arr, dtype='uint64')
    expected = PandasArray(np.array([1, 2, 3], dtype='uint64'))
    tm.assert_extension_array_equal(result, expected)


def test_validate_reduction_keyword_args():
    arr = PandasArray(np.array([1, 2, 3]))
    msg = "the 'keepdims' parameter is not supported .*all"
    with pytest.raises(ValueError, match=msg):
        arr.all(keepdims=True)


@td.skip_if_no("numpy", min_version="1.13.0")
def test_ufunc():
    arr = PandasArray([-1.0, 0.0, 1.0])
    result = np.abs(arr)
    expected = PandasArray(np.abs(arr._ndarray))
    tm.assert_extension_array_equal(result, expected)

    r1, r2 = np.divmod(arr, np.add(arr, 2))
    e1, e2 = np.divmod(arr._ndarray, np.add(arr._ndarray, 2))
    e1 = PandasArray(e1)
    e2 = PandasArray(e2)
    tm.assert_extension_array_equal(r1, e1)
    tm.assert_extension_array_equal(r2, e2)


@td.skip_if_no("numpy", min_version="1.13.0")
def test_basic_binop():
    # Just a basic smoke test. The EA interface tests exercise this
    # more thoroughly.
    x = PandasArray([1, 2, 3])
    result = x + x
    expected = PandasArray([2, 4, 6])
    tm.assert_extension_array_equal(result, expected)


def test_series_constructor_with_copy():
    ndarray = np.array([1, 2, 3])
    ser = pd.Series(PandasArray(ndarray), copy=True)

    assert ser.values is not ndarray


def test_series_constructor_with_astype():
    ndarray = np.array([1, 2, 3])
    result = pd.Series(PandasArray(ndarray), dtype="float64")
    expected = pd.Series([1.0, 2.0, 3.0], dtype="float64")
    tm.assert_series_equal(result, expected)
