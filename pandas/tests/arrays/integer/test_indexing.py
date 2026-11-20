import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_array_setitem_nullable_boolean_mask():
    # GH 31446
    ser = pd.Series([1, 2], dtype="Int64")
    result = ser.where(ser > 1)
    expected = pd.Series([pd.NA, 2], dtype="Int64")
    tm.assert_series_equal(result, expected)


def test_array_setitem():
    # GH 31446
    arr = pd.array([1, 2], dtype="Int64")
    arr[arr > 1] = 1

    expected = pd.array([1, 1], dtype="Int64")
    tm.assert_extension_array_equal(arr, expected)


def test_setitem_series_non_default_index(any_int_ea_dtype):
    # GH#72193 the lossless-cast check looked values up by label, not position
    arr = pd.array([1, 2], dtype=any_int_ea_dtype)
    arr[:] = pd.Series([3.0, np.nan], index=[5, 6])

    expected = pd.array([3, None], dtype=any_int_ea_dtype)
    tm.assert_extension_array_equal(arr, expected)


def test_setitem_series_keeps_precision():
    # GH#72193 looking up by label compared 2**53 + 1 with another element that
    #  equals it as a float, so the lost precision went undetected
    pytest.importorskip("pyarrow")
    big = 2**53 + 1
    arr = pd.array([0, 0, 0], dtype="Int64")
    value = pd.Series([big, 2**53, None], index=[1, 0, 2], dtype="int64[pyarrow]")
    arr[:] = value

    expected = pd.array([big, 2**53, None], dtype="Int64")
    tm.assert_extension_array_equal(arr, expected)


def test_setitem_masked_series_out_of_bounds_raises():
    # GH#72193 converting through float clipped 2**63 to 2**63 - 1
    arr = pd.array([0, 0], dtype="Int64")
    value = pd.Series([2**63, None], dtype="UInt64")
    with pytest.raises(TypeError, match="cannot safely cast"):
        arr[:] = value
