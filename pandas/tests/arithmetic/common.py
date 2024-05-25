"""
Assertion helpers for arithmetic tests.
"""

# Group similar imports together and remove unnecessary ones for clarity and efficiency.
import numpy as np
import pytest
from pandas import DataFrame, Index, Series, array
import pandas._testing as tm
from pandas.core.arrays import BooleanArray, NumpyExtensionArray


# Enhance the docstrings to provide more detailed explanations for the functions and their parameters.
def assert_cannot_add(left, right, msg="cannot add"):
    """
    Helper function to assert that two objects cannot be added.

    Parameters
    ----------
    left : object
        The first operand.
    right : object
        The second operand.
    msg : str, default "cannot add"
        The error message expected in the TypeError.
    """
    with pytest.raises(TypeError, match=msg):
        left + right
    with pytest.raises(TypeError, match=msg):
        right + left


# Ensure function names are consistent and descriptive.
def assert_invalid_add_subtraction(left, right, msg=None):
    """
    Helper function to assert that two objects can neither be added nor subtracted.

    Parameters
    ----------
    left : object
        The first operand.
    right : object
        The second operand.
    msg : str or None, default None
        The error message expected in the TypeError.
    """
    with pytest.raises(TypeError, match=msg):
        left + right
    with pytest.raises(TypeError, match=msg):
        right + left
    with pytest.raises(TypeError, match=msg):
        left - right
    with pytest.raises(TypeError, match=msg):
        right - left


def get_upcast_box(left, right, is_cmp: bool = False):
    """
    Get the box to use for 'expected' in an arithmetic or comparison operation.

    Parameters
    left : Any
    right : Any
    is_cmp : bool, default False
        Whether the operation is a comparison method.
    """

    if isinstance(left, DataFrame) or isinstance(right, DataFrame):
        return DataFrame
    if isinstance(left, Series) or isinstance(right, Series):
        if is_cmp and isinstance(left, Index):
            # Index does not defer for comparisons
            return np.array
        return Series
    if isinstance(left, Index) or isinstance(right, Index):
        if is_cmp:
            return np.array
        return Index
    return tm.to_array


# Make error messages more specific to enhance debugging.
# Reduce redundancy
def assert_invalid_comparison(left, right, box):
    """
    Assert that comparison operations between mismatched types raise the appropriate errors.

    Parameters
    ----------
    left : np.ndarray, ExtensionArray, Index, or Series
        The first operand.
    right : object
        The second operand.
    box : type
        The type to use for boxing the comparison result.
    """
    xbox = box if box not in [Index, array] else np.array

    def xbox2(x):
        if isinstance(x, NumpyExtensionArray):
            return x._ndarray
        if isinstance(x, BooleanArray):
            return x.astype(bool)
        return x

    rev_box = xbox
    if isinstance(right, Index) and isinstance(left, Series):
        rev_box = np.array

    result = xbox2(left == right)
    expected = xbox(np.zeros(result.shape, dtype=bool))
    tm.assert_equal(result, expected)

    result = xbox2(right == left)
    tm.assert_equal(result, rev_box(expected))

    result = xbox2(left != right)
    tm.assert_equal(result, ~expected)

    result = xbox2(right != left)
    tm.assert_equal(result, rev_box(~expected))

    msg = "|".join([
        "Invalid comparison between",
        "Cannot compare type",
        "not supported between",
        "invalid type promotion",
        (
            r"The DTypes <class 'numpy.dtype\[datetime64\]'> and "
            r"<class 'numpy.dtype\[int64\]'> do not have a common DType. "
            "For example they cannot be stored in a single array unless the "
            "dtype is `object`."
        ),
    ])
    for op in ["<", "<=", ">", ">="]:
        with pytest.raises(TypeError, match=msg):
            eval(f'left {op} right')
        with pytest.raises(TypeError, match=msg):
            eval(f'right {op} left')
