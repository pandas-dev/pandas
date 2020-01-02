"""
Assertion helpers for arithmetic tests.
"""
import numpy as np
import pytest

from pandas import DataFrame, Index, Series
import pandas._testing as tm


def assert_invalid_addsub_type(left, right, msg=None):
    """
    Helper to assert that left and right can be neither added nor subtracted.

    Parameters
    ---------
    left : object
    right : object
    msg : str or None, default None
    """
    with pytest.raises(TypeError, match=msg):
        left + right
    with pytest.raises(TypeError, match=msg):
        right + left
    with pytest.raises(TypeError, match=msg):
        left - right
    with pytest.raises(TypeError, match=msg):
        right - left


def get_upcast_box(box, vector):
    """
    Given two box-types, find the one that takes priority
    """
    if box is DataFrame or isinstance(vector, DataFrame):
        return DataFrame
    if box is Series or isinstance(vector, Series):
        return Series
    if box is Index or isinstance(vector, Index):
        return Index
    return box


def assert_invalid_comparison(left, right, box):
    """
    Assert that comparison operations with mismatched types behave correctly.

    Parameters
    ----------
    left : np.ndarray, ExtensionArray, Index, or Series
    right : object
    box : {pd.DataFrame, pd.Series, pd.Index, tm.to_array}
    """
    # Not for tznaive-tzaware comparison

    # Note: not quite the same as how we do this for tm.box_expected
    xbox = box if box is not Index else np.array

    result = left == right
    expected = xbox(np.zeros(result.shape, dtype=np.bool_))

    tm.assert_equal(result, expected)

    result = right == left
    tm.assert_equal(result, expected)

    result = left != right
    tm.assert_equal(result, ~expected)

    result = right != left
    tm.assert_equal(result, ~expected)

    msg = "Invalid comparison between"
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
