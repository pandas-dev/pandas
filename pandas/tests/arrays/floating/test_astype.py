import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_astype():
    # with missing values
    arr = pd.array([0.1, 0.2, None], dtype="Float64")

    with pytest.raises(ValueError, match="cannot convert to 'int64'-dtype NumPy"):
        arr.astype("int64")

    with pytest.raises(ValueError, match="cannot convert to 'bool'-dtype NumPy"):
        arr.astype("bool")

    result = arr.astype("float64")
    expected = np.array([0.1, 0.2, np.nan], dtype="float64")
    tm.assert_numpy_array_equal(result, expected)

    # no missing values
    arr = pd.array([0.0, 1.0, 0.5], dtype="Float64")
    result = arr.astype("int64")
    expected = np.array([0, 1, 0], dtype="int64")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype("bool")
    expected = np.array([False, True, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)


def test_astype_to_floating_array():
    # astype to FloatingArray
    arr = pd.array([0.0, 1.0, None], dtype="Float64")

    result = arr.astype("Float64")
    tm.assert_extension_array_equal(result, arr)
    result = arr.astype(pd.Float64Dtype())
    tm.assert_extension_array_equal(result, arr)
    result = arr.astype("Float32")
    expected = pd.array([0.0, 1.0, None], dtype="Float32")
    tm.assert_extension_array_equal(result, expected)


def test_astype_to_boolean_array():
    # astype to BooleanArray
    arr = pd.array([0.0, 1.0, None], dtype="Float64")

    result = arr.astype("boolean")
    expected = pd.array([False, True, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)
    result = arr.astype(pd.BooleanDtype())
    tm.assert_extension_array_equal(result, expected)


def test_astype_to_integer_array():
    # astype to IntegerArray
    arr = pd.array([0.0, 1.5, None], dtype="Float64")

    result = arr.astype("Int64")
    expected = pd.array([0, 1, None], dtype="Int64")
    tm.assert_extension_array_equal(result, expected)


def test_astype_str():
    a = pd.array([0.1, 0.2, None], dtype="Float64")
    expected = np.array(["0.1", "0.2", "<NA>"], dtype=object)

    tm.assert_numpy_array_equal(a.astype(str), expected)
    tm.assert_numpy_array_equal(a.astype("str"), expected)
