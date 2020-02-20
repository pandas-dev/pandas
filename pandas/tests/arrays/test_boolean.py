import operator

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm
from pandas.arrays import BooleanArray
from pandas.core.arrays.boolean import coerce_to_array
from pandas.tests.extension.base import BaseOpsUtil


def make_data():
    return [True, False] * 4 + [np.nan] + [True, False] * 44 + [np.nan] + [True, False]


@pytest.fixture
def dtype():
    return pd.BooleanDtype()


@pytest.fixture
def data(dtype):
    return pd.array(make_data(), dtype=dtype)


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

    with pytest.raises(ValueError, match="values must be a 1D array"):
        BooleanArray(values.reshape(1, -1), mask)

    with pytest.raises(ValueError, match="mask must be a 1D array"):
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
    "values",
    [
        ["foo", "bar"],
        ["1", "2"],
        # "foo",
        [1, 2],
        [1.0, 2.0],
        pd.date_range("20130101", periods=2),
        np.array(["foo"]),
        np.array([1, 2]),
        np.array([1.0, 2.0]),
        [np.nan, {"a": 1}],
    ],
)
def test_to_boolean_array_error(values):
    # error in converting existing arrays to BooleanArray
    msg = "Need to pass bool-like value"
    with pytest.raises(TypeError, match=msg):
        pd.array(values, dtype="boolean")


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

    with pytest.raises(ValueError, match="values must be a 1D list-like"):
        coerce_to_array(values.reshape(1, -1))

    with pytest.raises(ValueError, match="mask must be a 1D list-like"):
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
    expected = np.array([True, False, True], dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    # force bool dtype
    result = np.array(arr, dtype="bool")
    expected = np.array([True, False, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)
    # with missing values will raise error
    arr = pd.array([True, False, None], dtype="boolean")
    with pytest.raises(ValueError):
        np.array(arr, dtype="bool")


def test_to_boolean_array_from_strings():
    result = BooleanArray._from_sequence_of_strings(
        np.array(["True", "False", np.nan], dtype=object)
    )
    expected = BooleanArray(
        np.array([True, False, False]), np.array([False, False, True])
    )

    tm.assert_extension_array_equal(result, expected)


def test_to_boolean_array_from_strings_invalid_string():
    with pytest.raises(ValueError, match="cannot be cast"):
        BooleanArray._from_sequence_of_strings(["donkey"])


def test_repr():
    df = pd.DataFrame({"A": pd.array([True, False, None], dtype="boolean")})
    expected = "       A\n0   True\n1  False\n2   <NA>"
    assert repr(df) == expected

    expected = "0     True\n1    False\n2     <NA>\nName: A, dtype: boolean"
    assert repr(df.A) == expected

    expected = "<BooleanArray>\n[True, False, <NA>]\nLength: 3, dtype: boolean"
    assert repr(df.A.array) == expected


@pytest.mark.parametrize("box", [True, False], ids=["series", "array"])
def test_to_numpy(box):
    con = pd.Series if box else pd.array
    # default (with or without missing values) -> object dtype
    arr = con([True, False, True], dtype="boolean")
    result = arr.to_numpy()
    expected = np.array([True, False, True], dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    arr = con([True, False, None], dtype="boolean")
    result = arr.to_numpy()
    expected = np.array([True, False, pd.NA], dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    arr = con([True, False, None], dtype="boolean")
    result = arr.to_numpy(dtype="str")
    expected = np.array([True, False, pd.NA], dtype="<U5")
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
    with pytest.raises(ValueError, match="cannot convert to 'float64'-dtype"):
        arr.to_numpy(dtype="float64")


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


def test_astype():
    # with missing values
    arr = pd.array([True, False, None], dtype="boolean")

    with pytest.raises(ValueError, match="cannot convert NA to integer"):
        arr.astype("int64")

    with pytest.raises(ValueError, match="cannot convert float NaN to"):
        arr.astype("bool")

    result = arr.astype("float64")
    expected = np.array([1, 0, np.nan], dtype="float64")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype("str")
    expected = np.array(["True", "False", "<NA>"], dtype="object")
    tm.assert_numpy_array_equal(result, expected)

    # no missing values
    arr = pd.array([True, False, True], dtype="boolean")
    result = arr.astype("int64")
    expected = np.array([1, 0, 1], dtype="int64")
    tm.assert_numpy_array_equal(result, expected)

    result = arr.astype("bool")
    expected = np.array([True, False, True], dtype="bool")
    tm.assert_numpy_array_equal(result, expected)


def test_astype_to_boolean_array():
    # astype to BooleanArray
    arr = pd.array([True, False, None], dtype="boolean")

    result = arr.astype("boolean")
    tm.assert_extension_array_equal(result, arr)
    result = arr.astype(pd.BooleanDtype())
    tm.assert_extension_array_equal(result, arr)


def test_astype_to_integer_array():
    # astype to IntegerArray
    arr = pd.array([True, False, None], dtype="boolean")

    result = arr.astype("Int64")
    expected = pd.array([1, 0, None], dtype="Int64")
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("na", [None, np.nan, pd.NA])
def test_setitem_missing_values(na):
    arr = pd.array([True, False, None], dtype="boolean")
    expected = pd.array([True, None, None], dtype="boolean")
    arr[1] = na
    tm.assert_extension_array_equal(arr, expected)


@pytest.mark.parametrize(
    "ufunc", [np.add, np.logical_or, np.logical_and, np.logical_xor]
)
def test_ufuncs_binary(ufunc):
    # two BooleanArrays
    a = pd.array([True, False, None], dtype="boolean")
    result = ufunc(a, a)
    expected = pd.array(ufunc(a._data, a._data), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_extension_array_equal(result, expected)

    s = pd.Series(a)
    result = ufunc(s, a)
    expected = pd.Series(ufunc(a._data, a._data), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_series_equal(result, expected)

    # Boolean with numpy array
    arr = np.array([True, True, False])
    result = ufunc(a, arr)
    expected = pd.array(ufunc(a._data, arr), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_extension_array_equal(result, expected)

    result = ufunc(arr, a)
    expected = pd.array(ufunc(arr, a._data), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_extension_array_equal(result, expected)

    # BooleanArray with scalar
    result = ufunc(a, True)
    expected = pd.array(ufunc(a._data, True), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_extension_array_equal(result, expected)

    result = ufunc(True, a)
    expected = pd.array(ufunc(True, a._data), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_extension_array_equal(result, expected)

    # not handled types
    with pytest.raises(TypeError):
        ufunc(a, "test")


@pytest.mark.parametrize("ufunc", [np.logical_not])
def test_ufuncs_unary(ufunc):
    a = pd.array([True, False, None], dtype="boolean")
    result = ufunc(a)
    expected = pd.array(ufunc(a._data), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_extension_array_equal(result, expected)

    s = pd.Series(a)
    result = ufunc(s)
    expected = pd.Series(ufunc(a._data), dtype="boolean")
    expected[a._mask] = np.nan
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("values", [[True, False], [True, None]])
def test_ufunc_reduce_raises(values):
    a = pd.array(values, dtype="boolean")
    with pytest.raises(NotImplementedError):
        np.add.reduce(a)


class TestUnaryOps:
    def test_invert(self):
        a = pd.array([True, False, None], dtype="boolean")
        expected = pd.array([False, True, None], dtype="boolean")
        tm.assert_extension_array_equal(~a, expected)

        expected = pd.Series(expected, index=["a", "b", "c"], name="name")
        result = ~pd.Series(a, index=["a", "b", "c"], name="name")
        tm.assert_series_equal(result, expected)

        df = pd.DataFrame({"A": a, "B": [True, False, False]}, index=["a", "b", "c"])
        result = ~df
        expected = pd.DataFrame(
            {"A": expected, "B": [False, True, True]}, index=["a", "b", "c"]
        )
        tm.assert_frame_equal(result, expected)


class TestLogicalOps(BaseOpsUtil):
    def test_numpy_scalars_ok(self, all_logical_operators):
        a = pd.array([True, False, None], dtype="boolean")
        op = getattr(a, all_logical_operators)

        tm.assert_extension_array_equal(op(True), op(np.bool(True)))
        tm.assert_extension_array_equal(op(False), op(np.bool(False)))

    def get_op_from_name(self, op_name):
        short_opname = op_name.strip("_")
        short_opname = short_opname if "xor" in short_opname else short_opname + "_"
        try:
            op = getattr(operator, short_opname)
        except AttributeError:
            # Assume it is the reverse operator
            rop = getattr(operator, short_opname[1:])
            op = lambda x, y: rop(y, x)

        return op

    def test_empty_ok(self, all_logical_operators):
        a = pd.array([], dtype="boolean")
        op_name = all_logical_operators
        result = getattr(a, op_name)(True)
        tm.assert_extension_array_equal(a, result)

        result = getattr(a, op_name)(False)
        tm.assert_extension_array_equal(a, result)

        # TODO: pd.NA
        # result = getattr(a, op_name)(pd.NA)
        # tm.assert_extension_array_equal(a, result)

    def test_logical_length_mismatch_raises(self, all_logical_operators):
        op_name = all_logical_operators
        a = pd.array([True, False, None], dtype="boolean")
        msg = "Lengths must match to compare"

        with pytest.raises(ValueError, match=msg):
            getattr(a, op_name)([True, False])

        with pytest.raises(ValueError, match=msg):
            getattr(a, op_name)(np.array([True, False]))

        with pytest.raises(ValueError, match=msg):
            getattr(a, op_name)(pd.array([True, False], dtype="boolean"))

    def test_logical_nan_raises(self, all_logical_operators):
        op_name = all_logical_operators
        a = pd.array([True, False, None], dtype="boolean")
        msg = "Got float instead"

        with pytest.raises(TypeError, match=msg):
            getattr(a, op_name)(np.nan)

    @pytest.mark.parametrize("other", ["a", 1])
    def test_non_bool_or_na_other_raises(self, other, all_logical_operators):
        a = pd.array([True, False], dtype="boolean")
        with pytest.raises(TypeError, match=str(type(other).__name__)):
            getattr(a, all_logical_operators)(other)

    def test_kleene_or(self):
        # A clear test of behavior.
        a = pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        b = pd.array([True, False, None] * 3, dtype="boolean")
        result = a | b
        expected = pd.array(
            [True, True, True, True, False, None, True, None, None], dtype="boolean"
        )
        tm.assert_extension_array_equal(result, expected)

        result = b | a
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_extension_array_equal(
            a, pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        )
        tm.assert_extension_array_equal(
            b, pd.array([True, False, None] * 3, dtype="boolean")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (pd.NA, [True, None, None]),
            (True, [True, True, True]),
            (np.bool_(True), [True, True, True]),
            (False, [True, False, None]),
            (np.bool_(False), [True, False, None]),
        ],
    )
    def test_kleene_or_scalar(self, other, expected):
        # TODO: test True & False
        a = pd.array([True, False, None], dtype="boolean")
        result = a | other
        expected = pd.array(expected, dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

        result = other | a
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_extension_array_equal(
            a, pd.array([True, False, None], dtype="boolean")
        )

    def test_kleene_and(self):
        # A clear test of behavior.
        a = pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        b = pd.array([True, False, None] * 3, dtype="boolean")
        result = a & b
        expected = pd.array(
            [True, False, None, False, False, False, None, False, None], dtype="boolean"
        )
        tm.assert_extension_array_equal(result, expected)

        result = b & a
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_extension_array_equal(
            a, pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        )
        tm.assert_extension_array_equal(
            b, pd.array([True, False, None] * 3, dtype="boolean")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (pd.NA, [None, False, None]),
            (True, [True, False, None]),
            (False, [False, False, False]),
            (np.bool_(True), [True, False, None]),
            (np.bool_(False), [False, False, False]),
        ],
    )
    def test_kleene_and_scalar(self, other, expected):
        a = pd.array([True, False, None], dtype="boolean")
        result = a & other
        expected = pd.array(expected, dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

        result = other & a
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_extension_array_equal(
            a, pd.array([True, False, None], dtype="boolean")
        )

    def test_kleene_xor(self):
        a = pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        b = pd.array([True, False, None] * 3, dtype="boolean")
        result = a ^ b
        expected = pd.array(
            [False, True, None, True, False, None, None, None, None], dtype="boolean"
        )
        tm.assert_extension_array_equal(result, expected)

        result = b ^ a
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_extension_array_equal(
            a, pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        )
        tm.assert_extension_array_equal(
            b, pd.array([True, False, None] * 3, dtype="boolean")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (pd.NA, [None, None, None]),
            (True, [False, True, None]),
            (np.bool_(True), [False, True, None]),
            (np.bool_(False), [True, False, None]),
        ],
    )
    def test_kleene_xor_scalar(self, other, expected):
        a = pd.array([True, False, None], dtype="boolean")
        result = a ^ other
        expected = pd.array(expected, dtype="boolean")
        tm.assert_extension_array_equal(result, expected)

        result = other ^ a
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_extension_array_equal(
            a, pd.array([True, False, None], dtype="boolean")
        )

    @pytest.mark.parametrize(
        "other", [True, False, pd.NA, [True, False, None] * 3],
    )
    def test_no_masked_assumptions(self, other, all_logical_operators):
        # The logical operations should not assume that masked values are False!
        a = pd.arrays.BooleanArray(
            np.array([True, True, True, False, False, False, True, False, True]),
            np.array([False] * 6 + [True, True, True]),
        )
        b = pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        if isinstance(other, list):
            other = pd.array(other, dtype="boolean")

        result = getattr(a, all_logical_operators)(other)
        expected = getattr(b, all_logical_operators)(other)
        tm.assert_extension_array_equal(result, expected)

        if isinstance(other, BooleanArray):
            other._data[other._mask] = True
            a._data[a._mask] = False

            result = getattr(a, all_logical_operators)(other)
            expected = getattr(b, all_logical_operators)(other)
            tm.assert_extension_array_equal(result, expected)


class TestComparisonOps(BaseOpsUtil):
    def _compare_other(self, data, op_name, other):
        op = self.get_op_from_name(op_name)

        # array
        result = pd.Series(op(data, other))
        expected = pd.Series(op(data._data, other), dtype="boolean")
        # propagate NAs
        expected[data._mask] = pd.NA

        tm.assert_series_equal(result, expected)

        # series
        s = pd.Series(data)
        result = op(s, other)

        expected = pd.Series(data._data)
        expected = op(expected, other)
        expected = expected.astype("boolean")
        # propagate NAs
        expected[data._mask] = pd.NA

        tm.assert_series_equal(result, expected)

    def test_compare_scalar(self, data, all_compare_operators):
        op_name = all_compare_operators
        self._compare_other(data, op_name, True)

    def test_compare_array(self, data, all_compare_operators):
        op_name = all_compare_operators
        other = pd.array([True] * len(data), dtype="boolean")
        self._compare_other(data, op_name, other)
        other = np.array([True] * len(data))
        self._compare_other(data, op_name, other)
        other = pd.Series([True] * len(data))
        self._compare_other(data, op_name, other)

    @pytest.mark.parametrize("other", [True, False, pd.NA])
    def test_scalar(self, other, all_compare_operators):
        op = self.get_op_from_name(all_compare_operators)
        a = pd.array([True, False, None], dtype="boolean")

        result = op(a, other)

        if other is pd.NA:
            expected = pd.array([None, None, None], dtype="boolean")
        else:
            values = op(a._data, other)
            expected = BooleanArray(values, a._mask, copy=True)
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        result[0] = None
        tm.assert_extension_array_equal(
            a, pd.array([True, False, None], dtype="boolean")
        )

    def test_array(self, all_compare_operators):
        op = self.get_op_from_name(all_compare_operators)
        a = pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        b = pd.array([True, False, None] * 3, dtype="boolean")

        result = op(a, b)

        values = op(a._data, b._data)
        mask = a._mask | b._mask
        expected = BooleanArray(values, mask)
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        result[0] = None
        tm.assert_extension_array_equal(
            a, pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")
        )
        tm.assert_extension_array_equal(
            b, pd.array([True, False, None] * 3, dtype="boolean")
        )


class TestArithmeticOps(BaseOpsUtil):
    def test_error(self, data, all_arithmetic_operators):
        # invalid ops

        op = all_arithmetic_operators
        s = pd.Series(data)
        ops = getattr(s, op)
        opa = getattr(data, op)

        # invalid scalars
        with pytest.raises(TypeError):
            ops("foo")
        with pytest.raises(TypeError):
            ops(pd.Timestamp("20180101"))

        # invalid array-likes
        if op not in ("__mul__", "__rmul__"):
            # TODO(extension) numpy's mul with object array sees booleans as numbers
            with pytest.raises(TypeError):
                ops(pd.Series("foo", index=s.index))

        # 2d
        result = opa(pd.DataFrame({"A": s}))
        assert result is NotImplemented

        with pytest.raises(NotImplementedError):
            opa(np.arange(len(s)).reshape(-1, len(s)))


@pytest.mark.parametrize("dropna", [True, False])
def test_reductions_return_types(dropna, data, all_numeric_reductions):
    op = all_numeric_reductions
    s = pd.Series(data)
    if dropna:
        s = s.dropna()

    if op in ("sum", "prod"):
        assert isinstance(getattr(s, op)(), np.int64)
    elif op in ("min", "max"):
        assert isinstance(getattr(s, op)(), np.bool_)
    else:
        # "mean", "std", "var", "median", "kurt", "skew"
        assert isinstance(getattr(s, op)(), np.float64)


@pytest.mark.parametrize(
    "values, exp_any, exp_all, exp_any_noskip, exp_all_noskip",
    [
        ([True, pd.NA], True, True, True, pd.NA),
        ([False, pd.NA], False, False, pd.NA, False),
        ([pd.NA], False, True, pd.NA, pd.NA),
        ([], False, True, False, True),
    ],
)
def test_any_all(values, exp_any, exp_all, exp_any_noskip, exp_all_noskip):
    # the methods return numpy scalars
    exp_any = pd.NA if exp_any is pd.NA else np.bool_(exp_any)
    exp_all = pd.NA if exp_all is pd.NA else np.bool_(exp_all)
    exp_any_noskip = pd.NA if exp_any_noskip is pd.NA else np.bool_(exp_any_noskip)
    exp_all_noskip = pd.NA if exp_all_noskip is pd.NA else np.bool_(exp_all_noskip)

    for con in [pd.array, pd.Series]:
        a = con(values, dtype="boolean")
        assert a.any() is exp_any
        assert a.all() is exp_all
        assert a.any(skipna=False) is exp_any_noskip
        assert a.all(skipna=False) is exp_all_noskip

        assert np.any(a.any()) is exp_any
        assert np.all(a.all()) is exp_all


# TODO when BooleanArray coerces to object dtype numpy array, need to do conversion
# manually in the indexing code
# def test_indexing_boolean_mask():
#     arr = pd.array([1, 2, 3, 4], dtype="Int64")
#     mask = pd.array([True, False, True, False], dtype="boolean")
#     result = arr[mask]
#     expected = pd.array([1, 3], dtype="Int64")
#     tm.assert_extension_array_equal(result, expected)

#     # missing values -> error
#     mask = pd.array([True, False, True, None], dtype="boolean")
#     with pytest.raises(IndexError):
#         result = arr[mask]


@td.skip_if_no("pyarrow", min_version="0.15.0")
def test_arrow_array(data):
    # protocol added in 0.15.0
    import pyarrow as pa

    arr = pa.array(data)

    # TODO use to_numpy(na_value=None) here
    data_object = np.array(data, dtype=object)
    data_object[data.isna()] = None
    expected = pa.array(data_object, type=pa.bool_(), from_pandas=True)
    assert arr.equals(expected)


@td.skip_if_no("pyarrow", min_version="0.15.1.dev")
def test_arrow_roundtrip():
    # roundtrip possible from arrow 1.0.0
    import pyarrow as pa

    data = pd.array([True, False, None], dtype="boolean")
    df = pd.DataFrame({"a": data})
    table = pa.table(df)
    assert table.field("a").type == "bool"
    result = table.to_pandas()
    assert isinstance(result["a"].dtype, pd.BooleanDtype)
    tm.assert_frame_equal(result, df)


def test_value_counts_na():
    arr = pd.array([True, False, pd.NA], dtype="boolean")
    result = arr.value_counts(dropna=False)
    expected = pd.Series([1, 1, 1], index=[True, False, pd.NA], dtype="Int64")
    tm.assert_series_equal(result, expected)

    result = arr.value_counts(dropna=True)
    expected = pd.Series([1, 1], index=[True, False], dtype="Int64")
    tm.assert_series_equal(result, expected)


def test_diff():
    a = pd.array(
        [True, True, False, False, True, None, True, None, False], dtype="boolean"
    )
    result = pd.core.algorithms.diff(a, 1)
    expected = pd.array(
        [None, False, True, False, True, None, None, None, None], dtype="boolean"
    )
    tm.assert_extension_array_equal(result, expected)

    s = pd.Series(a)
    result = s.diff()
    expected = pd.Series(expected)
    tm.assert_series_equal(result, expected)
