import decimal
import warnings

import numpy as np
import pytest

from pandas.core.dtypes.cast import (
    maybe_downcast_numeric,
    maybe_downcast_to_dtype,
)

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


@pytest.mark.parametrize(
    "value, dtype",
    [
        # float64(np.iinfo(np.int64).max) == 2**63, so the top edge needs a
        #  strict inequality to be rejected.
        (2.0**63, "int64"),
        (-(2.0**63) - 2048, "int64"),
        (2.0**64, "uint64"),
        (-1.0, "uint64"),
        (2.0**31, "int32"),
        (np.inf, "int64"),
        (-np.inf, "int64"),
    ],
)
def test_downcast_float_out_of_integer_range(value, dtype):
    # GH#66394 casting an out-of-range float to an integer dtype is undefined
    #  behavior, so we must not downcast; previously on aarch64 the saturated
    #  cast compared equal in float64 space and was silently kept.
    arr = np.array([value])

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = maybe_downcast_numeric(arr, np.dtype(dtype))

    tm.assert_numpy_array_equal(result, arr)


@pytest.mark.parametrize(
    "value, dtype",
    [
        # Not exactly representable as float64, but still in range and
        #  round-tripping exactly.
        (2.0**62, "int64"),
        (2.0**63, "uint64"),
        (2.0**31 - 1, "int32"),
        (1.0, "int64"),
    ],
)
def test_downcast_float_within_integer_range(value, dtype):
    # GH#66394 the bounds check must not block downcasts that are exact
    arr = np.array([value])

    result = maybe_downcast_numeric(arr, np.dtype(dtype))

    tm.assert_numpy_array_equal(result, np.array([value], dtype=dtype))


@pytest.mark.parametrize("value", [2**63, 2**64, -(2**64), 2.0**64, np.inf])
def test_downcast_object_out_of_integer_range(value):
    # GH#66394 object entries that don't fit the target dtype previously raised
    #  OverflowError out of the astype
    arr = np.array([value], dtype=object)

    result = maybe_downcast_numeric(arr, np.dtype("int64"))

    tm.assert_numpy_array_equal(result, arr)


def test_downcast_object_within_integer_range():
    # GH#66394 in-range object entries still downcast
    arr = np.array([1, 2], dtype=object)

    result = maybe_downcast_numeric(arr, np.dtype("int64"))

    tm.assert_numpy_array_equal(result, np.array([1, 2], dtype="int64"))


def test_downcast_do_round_out_of_range():
    # GH#66394 rounding can push a value out of range, so the bounds check
    #  has to happen after rounding
    arr = np.array([np.iinfo(np.int32).max + 0.6])

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = maybe_downcast_numeric(arr, np.dtype("int32"), do_round=True)

    tm.assert_numpy_array_equal(result, arr)


@pytest.mark.parametrize(
    "klass, nat",
    [
        (np.datetime64, np.datetime64("NaT", "ns")),
        (np.timedelta64, np.timedelta64("NaT", "ns")),
    ],
)
def test_datetime_likes_nan(klass, nat):
    dtype = np.dtype(klass.__name__ + "[ns]")
    arr = np.array([1, 2, np.nan])

    exp = np.array([1, 2, nat], dtype)
    res = maybe_downcast_to_dtype(arr, dtype)
    tm.assert_numpy_array_equal(res, exp)
