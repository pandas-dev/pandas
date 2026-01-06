"""
Assertion helpers for arithmetic tests.
"""

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Index,
    Series,
    array,
)
import pandas._testing as tm
from pandas.core.arrays import (
    BooleanArray,
    NumpyExtensionArray,
)


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


def assert_invalid_addsub_type(left, right, msg=None):
    """
    Helper function to assert that two objects can
    neither be added nor subtracted.

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


def assert_invalid_comparison(left, right, box):
    """
    Assert that comparison operations with mismatched types behave correctly.

    Parameters
    ----------
    left : np.ndarray, ExtensionArray, Index, or Series
    right : object
    box : {pd.DataFrame, pd.Series, pd.Index, pd.array, tm.to_array}
    """
    # Not for tznaive-tzaware comparison

    # Note: not quite the same as how we do this for tm.box_expected
    xbox = box if box not in [Index, array] else np.array

    def xbox2(x):
        # Eventually we'd like this to be tighter, but for now we'll
        #  just exclude NumpyExtensionArray[bool]
        if isinstance(x, NumpyExtensionArray):
            return x._ndarray
        if isinstance(x, BooleanArray):
            # NB: we are assuming no pd.NAs for now
            return x.astype(bool)
        return x

    result = xbox2(left == right)
    expected = xbox(np.zeros(result.shape, dtype=np.bool_))

    tm.assert_equal(result, expected)

    result = xbox2(right == left)
    tm.assert_equal(result, xbox(expected))

    result = xbox2(left != right)
    tm.assert_equal(result, ~expected)

    result = xbox2(right != left)
    tm.assert_equal(result, xbox(~expected))

    msg = "|".join(
        [
            "Invalid comparison between",
            "Cannot compare type",
            "not supported between",
            "invalid type promotion",
            (
                # GH#36706 npdev 1.20.0 2020-09-28
                r"The DTypes <class 'numpy.dtype\[datetime64\]'> and "
                r"<class 'numpy.dtype\[int64\]'> do not have a common DType. "
                "For example they cannot be stored in a single array unless the "
                "dtype is `object`."
            ),
        ]
    )
    with pytest.raises(TypeError, match=msg):
        left < right
    with pytest.raises(TypeError, match=msg):
        left <= right
    with pytest.raises(TypeError, match=msg):
        left > right
    with pytest.raises(TypeError, match=msg):
        left >= right
    with pytest.raises(TypeError, match=msg):
        right < left
    with pytest.raises(TypeError, match=msg):
        right <= left
    with pytest.raises(TypeError, match=msg):
        right > left
    with pytest.raises(TypeError, match=msg):
        right >= left
