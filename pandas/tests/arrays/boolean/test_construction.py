from decimal import Decimal

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.arrays import BooleanArray
from pandas.core.arrays.boolean import coerce_to_array


def test_boolean_array_constructor():
    values = np.array([True, False, True, False], dtype="bool")
    mask = np.array([False, False, False, True], dtype="bool")

    result = BooleanArray(values, mask)
    expected = pd.array([True, False, True, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    with pytest.raises(TypeError, match="values should be boolean numpy array"):
        BooleanArray(values.tolist(), mask)

    with pytest.raises(TypeError, match="mask should be boolean numpy array"):
        BooleanArray(values, mask.tolist())

    with pytest.raises(TypeError, match="values should be boolean numpy array"):
        BooleanArray(values.astype(int), mask)

    with pytest.raises(TypeError, match="mask should be boolean numpy array"):
        BooleanArray(values, None)

    with pytest.raises(ValueError, match="values.shape must match mask.shape"):
        BooleanArray(values.reshape(1, -1), mask)

    with pytest.raises(ValueError, match="values.shape must match mask.shape"):
        BooleanArray(values, mask.reshape(1, -1))


def test_boolean_array_constructor_copy():
    values = np.array([True, False, True, False], dtype="bool")
    mask = np.array([False, False, False, True], dtype="bool")

    result = BooleanArray(values, mask)
    assert result._data is values
    assert result._mask is mask

    result = BooleanArray(values, mask, copy=True)
    assert result._data is not values
    assert result._mask is not mask


def test_to_boolean_array():
    expected = BooleanArray(
        np.array([True, False, True]), np.array([False, False, False])
    )

    result = pd.array([True, False, True], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)
    result = pd.array(np.array([True, False, True]), dtype="boolean")
    tm.assert_extension_array_equal(result, expected)
    result = pd.array(np.array([True, False, True], dtype=object), dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    # with missing values
    expected = BooleanArray(
        np.array([True, False, True]), np.array([False, False, True])
    )

    result = pd.array([True, False, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)
    result = pd.array(np.array([True, False, None], dtype=object), dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_to_boolean_array_all_none():
    expected = BooleanArray(np.array([True, True, True]), np.array([True, True, True]))

    result = pd.array([None, None, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)
    result = pd.array(np.array([None, None, None], dtype=object), dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "a, b",
    [
        ([True, False, None, np.nan, pd.NA], [True, False, None, None, None]),
        ([True, np.nan], [True, None]),
        ([True, pd.NA], [True, None]),
        ([np.nan, np.nan], [None, None]),
        (np.array([np.nan, np.nan], dtype=float), [None, None]),
    ],
)
def test_to_boolean_array_missing_indicators(a, b):
    result = pd.array(a, dtype="boolean")
    expected = pd.array(b, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "values, expected",
    [
        (["foo", "", "False", None], [True, False, True, None]),
        ([0, 1, 2, -1, None], [False, True, True, True, None]),
        ([0.0, 0.5, np.nan], [False, True, None]),
        ([0, "a", 1.5, pd.NA, pd.NaT], [False, True, True, None, None]),
        ([np.nan, {"a": 1}, {}], [None, True, False]),
        # a list is treated as object dtype, where every Timestamp is truthy
        ([pd.Timestamp(0), pd.NaT], [True, None]),
        (np.array(["foo", ""]), [True, False]),
        (np.array([0, 1, 2]), [False, True, True]),
        (np.array([0.0, 0.5, np.nan]), [False, True, None]),
        (np.array([0j, 1j, complex(np.nan, 0)]), [False, True, None]),
        (np.array([0, 5, "NaT"], dtype="m8[ns]"), [False, True, None]),
        (
            np.array(["1970-01-01", "2000-01-01", "NaT"], dtype="M8[ns]"),
            [False, True, None],
        ),
        (pd.date_range("1970-01-01", periods=2, tz="UTC"), [False, True]),
        (pd.period_range("2000", periods=2), [True, True]),
    ],
)
def test_to_boolean_array_matches_astype_bool(values, expected):
    # GH#37614 values are converted as by astype(bool), with missing values masked
    result = pd.array(values, dtype="boolean")
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "values, expected",
    [
        (
            [0.0, np.nan, None, pd.NA, pd.NaT, Decimal("nan")],
            [False, True, None, None, None, True],
        ),
        (np.array([0.0, np.nan]), [False, True]),
        (np.array([0j, complex(np.nan, 0)]), [False, True]),
    ],
)
def test_to_boolean_array_distinguish_nan_and_na(values, expected):
    # with NaN distinct from NA, NaN is a value and converts as bool(np.nan)
    with pd.option_context("future.distinguish_nan_and_na", True):
        result = pd.array(values, dtype="boolean")
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_float_ea_astype_boolean_distinguish_nan_and_na():
    with pd.option_context("future.distinguish_nan_and_na", True):
        arr = pd.array([0.0, np.nan, None], dtype="Float64")
        result = arr.astype("boolean")
    expected = pd.array([False, True, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_to_boolean_array_2d_raises():
    with pytest.raises(TypeError, match="values must be a 1D list-like"):
        pd.array([[True, False]], dtype="boolean")


def test_to_boolean_array_from_integer_array():
    result = pd.array(np.array([1, 0, 1, 0]), dtype="boolean")
    expected = pd.array([True, False, True, False], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    # with missing values
    result = pd.array(np.array([1, 0, 1, None]), dtype="boolean")
    expected = pd.array([True, False, True, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_to_boolean_array_from_float_array():
    result = pd.array(np.array([1.0, 0.0, 1.0, 0.0]), dtype="boolean")
    expected = pd.array([True, False, True, False], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    # with missing values
    result = pd.array(np.array([1.0, 0.0, 1.0, np.nan]), dtype="boolean")
    expected = pd.array([True, False, True, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_to_boolean_array_integer_like():
    # integers of 0's and 1's
    result = pd.array([1, 0, 1, 0], dtype="boolean")
    expected = pd.array([True, False, True, False], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

    # with missing values
    result = pd.array([1, 0, 1, None], dtype="boolean")
    expected = pd.array([True, False, True, None], dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_coerce_to_array():
    # TODO this is currently not public API
    values = np.array([True, False, True, False], dtype="bool")
    mask = np.array([False, False, False, True], dtype="bool")
    result = BooleanArray(*coerce_to_array(values, mask=mask))
    expected = BooleanArray(values, mask)
    tm.assert_extension_array_equal(result, expected)
    assert result._data is values
    assert result._mask is mask
    result = BooleanArray(*coerce_to_array(values, mask=mask, copy=True))
    expected = BooleanArray(values, mask)
    tm.assert_extension_array_equal(result, expected)
    assert result._data is not values
    assert result._mask is not mask

    # mixed missing from values and mask
    values = [True, False, None, False]
    mask = np.array([False, False, False, True], dtype="bool")
    result = BooleanArray(*coerce_to_array(values, mask=mask))
    expected = BooleanArray(
        np.array([True, False, True, True]), np.array([False, False, True, True])
    )
    tm.assert_extension_array_equal(result, expected)
    result = BooleanArray(*coerce_to_array(np.array(values, dtype=object), mask=mask))
    tm.assert_extension_array_equal(result, expected)
    result = BooleanArray(*coerce_to_array(values, mask=mask.tolist()))
    tm.assert_extension_array_equal(result, expected)

    # raise errors for wrong dimension
    values = np.array([True, False, True, False], dtype="bool")
    mask = np.array([False, False, False, True], dtype="bool")

    # passing 2D values is OK as long as no mask
    coerce_to_array(values.reshape(1, -1))

    with pytest.raises(ValueError, match="values.shape and mask.shape must match"):
        coerce_to_array(values.reshape(1, -1), mask=mask)

    with pytest.raises(ValueError, match="values.shape and mask.shape must match"):
        coerce_to_array(values, mask=mask.reshape(1, -1))


def test_coerce_to_array_from_boolean_array():
    # passing BooleanArray to coerce_to_array
    values = np.array([True, False, True, False], dtype="bool")
    mask = np.array([False, False, False, True], dtype="bool")
    arr = BooleanArray(values, mask)
    result = BooleanArray(*coerce_to_array(arr))
    tm.assert_extension_array_equal(result, arr)
    # no copy
    assert result._data is arr._data
    assert result._mask is arr._mask

    result = BooleanArray(*coerce_to_array(arr), copy=True)
    tm.assert_extension_array_equal(result, arr)
    assert result._data is not arr._data
    assert result._mask is not arr._mask

    with pytest.raises(ValueError, match="cannot pass mask for BooleanArray input"):
        coerce_to_array(arr, mask=mask)


def test_coerce_to_numpy_array():
    # with missing values -> object dtype
    arr = pd.array([True, False, None], dtype="boolean")
    result = np.array(arr)
    expected = np.array([True, False, pd.NA], dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    # also with no missing values -> object dtype
    arr = pd.array([True, False, True], dtype="boolean")
    result = np.array(arr)
    expected = np.array([True, False, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)

    # force bool dtype
    result = np.array(arr, dtype="bool")
    expected = np.array([True, False, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)
    # with missing values will raise error
    arr = pd.array([True, False, None], dtype="boolean")
    msg = (
        "cannot convert to 'bool'-dtype NumPy array with missing values. "
        "Specify an appropriate 'na_value' for this dtype."
    )
    with pytest.raises(ValueError, match=msg):
        np.array(arr, dtype="bool")


def test_to_boolean_array_from_strings():
    result = BooleanArray._from_sequence_of_strings(
        np.array(["True", "False", "1", "1.0", "0", "0.0", np.nan], dtype=object),
        dtype=pd.BooleanDtype(),
    )
    expected = BooleanArray(
        np.array([True, False, True, True, False, False, False]),
        np.array([False, False, False, False, False, False, True]),
    )

    tm.assert_extension_array_equal(result, expected)


def test_to_boolean_array_from_strings_invalid_string():
    with pytest.raises(ValueError, match="cannot be cast"):
        BooleanArray._from_sequence_of_strings(["donkey"], dtype=pd.BooleanDtype())


@pytest.mark.parametrize("box", [True, False], ids=["series", "array"])
def test_to_numpy(box):
    con = pd.Series if box else pd.array
    # default (with or without missing values) -> object dtype
    arr = con([True, False, True], dtype="boolean")
    result = arr.to_numpy()
    expected = np.array([True, False, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)

    arr = con([True, False, None], dtype="boolean")
    result = arr.to_numpy()
    expected = np.array([True, False, pd.NA], dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    arr = con([True, False, None], dtype="boolean")
    result = arr.to_numpy(dtype="str")
    expected = np.array([True, False, pd.NA], dtype=f"{tm.ENDIAN}U5")
    tm.assert_numpy_array_equal(result, expected)

    # no missing values -> can convert to bool, otherwise raises
    arr = con([True, False, True], dtype="boolean")
    result = arr.to_numpy(dtype="bool")
    expected = np.array([True, False, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)

    arr = con([True, False, None], dtype="boolean")
    with pytest.raises(ValueError, match="cannot convert to 'bool'-dtype"):
        result = arr.to_numpy(dtype="bool")

    # specify dtype and na_value
    arr = con([True, False, None], dtype="boolean")
    result = arr.to_numpy(dtype=object, na_value=None)
    expected = np.array([True, False, None], dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.to_numpy(dtype=bool, na_value=False)
    expected = np.array([True, False, False], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.to_numpy(dtype="int64", na_value=-99)
    expected = np.array([1, 0, -99], dtype="int64")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.to_numpy(dtype="float64", na_value=np.nan)
    expected = np.array([1, 0, np.nan], dtype="float64")
    tm.assert_numpy_array_equal(result, expected)

    # converting to int or float without specifying na_value raises
    with pytest.raises(ValueError, match="cannot convert to 'int64'-dtype"):
        arr.to_numpy(dtype="int64")


def test_to_numpy_copy():
    # to_numpy can be zero-copy if no missing values
    arr = pd.array([True, False, True], dtype="boolean")
    result = arr.to_numpy(dtype=bool)
    result[0] = False
    tm.assert_extension_array_equal(
        arr, pd.array([False, False, True], dtype="boolean")
    )

    arr = pd.array([True, False, True], dtype="boolean")
    result = arr.to_numpy(dtype=bool, copy=True)
    result[0] = False
    tm.assert_extension_array_equal(arr, pd.array([True, False, True], dtype="boolean"))


def test_to_numpy_readonly():
    arr = pd.array([True, False, True], dtype="boolean")
    arr._readonly = True
    result = arr.to_numpy(dtype=bool)
    assert not result.flags.writeable

    result = arr.to_numpy(dtype="int64")
    assert result.flags.writeable
