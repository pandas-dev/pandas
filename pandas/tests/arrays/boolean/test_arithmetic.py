import operator

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def data():
    """Fixture returning boolean array with valid and missing values."""
    return pd.array(
        [True, False] * 2 + [np.nan] + [True, False] + [np.nan] + [True, False],
        dtype="boolean",
    )


@pytest.fixture
def left_array():
    """Fixture returning boolean array with valid and missing values."""
    return pd.array([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean")


@pytest.fixture
def right_array():
    """Fixture returning boolean array with valid and missing values."""
    return pd.array([True, False, None] * 3, dtype="boolean")


# Basic test for the arithmetic array ops
# -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "opname, exp",
    [
        ("add", [True, True, None, True, False, None, None, None, None]),
        ("mul", [True, False, None, False, False, None, None, None, None]),
    ],
    ids=["add", "mul"],
)
def test_add_mul(left_array, right_array, opname, exp):
    op = getattr(operator, opname)
    result = op(left_array, right_array)
    expected = pd.array(exp, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)


def test_sub(left_array, right_array):
    msg = (
        r"numpy boolean subtract, the `-` operator, is (?:deprecated|not supported), "
        r"use the bitwise_xor, the `\^` operator, or the logical_xor function instead\."
    )
    with pytest.raises(TypeError, match=msg):
        left_array - right_array


def test_div(left_array, right_array):
    msg = "operator '.*' not implemented for bool dtypes"
    with pytest.raises(NotImplementedError, match=msg):
        # check that we are matching the non-masked Series behavior
        pd.Series(left_array._data) / pd.Series(right_array._data)

    with pytest.raises(NotImplementedError, match=msg):
        left_array / right_array


@pytest.mark.parametrize(
    "opname",
    [
        "floordiv",
        "mod",
        "pow",
    ],
)
def test_op_int8(left_array, right_array, opname):
    op = getattr(operator, opname)
    if opname != "mod":
        msg = "operator '.*' not implemented for bool dtypes"
        with pytest.raises(NotImplementedError, match=msg):
            result = op(left_array, right_array)
        return
    result = op(left_array, right_array)
    expected = op(left_array.astype("Int8"), right_array.astype("Int8"))
    tm.assert_extension_array_equal(result, expected)


# Test generic characteristics / errors
# -----------------------------------------------------------------------------


def test_error_invalid_values(data, all_arithmetic_operators):
    # invalid ops
    op = all_arithmetic_operators
    s = pd.Series(data)
    ops = getattr(s, op)

    # invalid scalars
    msg = (
        "did not contain a loop with signature matching types|"
        "BooleanArray cannot perform the operation|"
        "not supported for the input types, and the inputs could not be safely coerced "
        "to any supported types according to the casting rule ''safe''|"
        "not supported for dtype"
    )
    with pytest.raises(TypeError, match=msg):
        ops("foo")
    msg = "|".join(
        [
            r"unsupported operand type\(s\) for",
            "Concatenation operation is not implemented for NumPy arrays",
            "has no kernel",
            "not supported for dtype",
        ]
    )
    with pytest.raises(TypeError, match=msg):
        ops(pd.Timestamp("20180101"))

    # invalid array-likes
    if op not in ("__mul__", "__rmul__"):
        # TODO(extension) numpy's mul with object array sees booleans as numbers
        msg = "|".join(
            [
                r"unsupported operand type\(s\) for",
                "can only concatenate str",
                "not all arguments converted during string formatting",
                "has no kernel",
                "not implemented",
                "not supported for dtype",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            ops(pd.Series("foo", index=s.index))

# Tests for BooleanArray logical and comparison operations
# -------------------------------------------------------------------------

@pytest.mark.parametrize(
    "opname, expected",[
        ("and_", [True, False, None, False, False, None, None, None, None]),
        ("or_", [True, True, None, True, False, None, None, None, None]),
        ("xor", [False, True, None, True, False, None, None, None, None]),
    ],ids=["and", "or", "xor"],
)

def test_logical_ops(left_array, right_array, opname, expected):
    """Test bitwise logical operations (&, |, ^) for BooleanArray."""
    op = getattr(operator, opname)
    result = op(left_array, right_array)
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

def test_invert(left_array):
    """Test inversion (~) for BooleanArray."""
    result = operator.invert(left_array)
    expected = pd.array(
        [False, False, False, True, True, True, None, None, None],
        dtype="boolean",
    )
    tm.assert_extension_array_equal(result, expected)

@pytest.mark.parametrize(
    "opname, expected",[
        ("eq", [True, False, None, False, True, None, None, None, None]),
        ("ne", [False, True, None, True, False, None, None, None, None]),
        ("gt", [False, True, None, False, False, None, None, None, None]),
        ("lt", [False, False, None, False, False, None, None, None, None]),
        ("ge", [True, False, None, True, True, None, None, None, None]),
        ("le", [True, True, None, False, True, None, None, None, None]),
    ],ids=["==", "!=", ">", "<", ">=", "<="],
)
def test_comparison_ops(left_array, right_array, opname, expected):
    """Test comparison operators (==, !=, >, <, >=, <=) for BooleanArray."""
    op = getattr(operator, opname)
    result = op(left_array, right_array)
    expected = pd.array(expected, dtype="boolean")
    tm.assert_extension_array_equal(result, expected)

def test_logical_with_scalar(left_array):
    """Test logical ops with scalar values (True/False/NA)."""
    result_true = left_array & True
    result_false = left_array | False
    result_na = left_array ^ pd.NA

    expected_true = left_array  # '&' True→Same
    expected_false = left_array  # '|' False→same
    expected_na = pd.array([None] * len(left_array), dtype="boolean")  # XOR NA→All NA

    tm.assert_extension_array_equal(result_true, expected_true)
    tm.assert_extension_array_equal(result_false, expected_false)
    tm.assert_extension_array_equal(result_na, expected_na)
