import decimal

import numpy as np
import pytest

from pandas.core.dtypes.cast import maybe_downcast_to_dtype

from pandas import (
    Series,
)
import pandas._testing as tm


@pytest.mark.parametrize(
    "arr,dtype,expected",
    [
        (
            # This is a judgement call, but we do _not_ downcast Decimal
            #  objects
            np.array([decimal.Decimal("0.0")]),
            np.dtype("int64"),
            np.array([decimal.Decimal("0.0")]),
        ),
    ],
)
def test_downcast(arr, expected, dtype):
    result = maybe_downcast_to_dtype(arr, dtype)
    tm.assert_numpy_array_equal(result, expected)


def test_downcast_booleans():
    # see gh-16875: coercing of booleans.
    ser = Series([True, True, False])
    result = maybe_downcast_to_dtype(ser, np.dtype(np.float64))

    expected = ser.values
    tm.assert_numpy_array_equal(result, expected)


def test_downcast_conversion_empty(any_real_numpy_dtype):
    dtype = any_real_numpy_dtype
    arr = np.array([], dtype=dtype)
    result = maybe_downcast_to_dtype(arr, np.dtype("int64"))
    tm.assert_numpy_array_equal(result, np.array([], dtype=np.int64))


@pytest.mark.parametrize("klass", [np.datetime64, np.timedelta64])
def test_datetime_likes_nan(klass):
    dtype = np.dtype(klass.__name__ + "[ns]")
    arr = np.array([1, 2, np.nan])

    exp = np.array([1, 2, klass("NaT")], dtype)
    res = maybe_downcast_to_dtype(arr, dtype)
    tm.assert_numpy_array_equal(res, exp)
