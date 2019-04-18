# -*- coding: utf-8 -*-

import numpy as np
import pytest

from pandas.core.dtypes.cast import maybe_downcast_to_dtype

from pandas import DatetimeIndex, Series, Timestamp
from pandas.util import testing as tm


@pytest.mark.parametrize("arr,dtype,expected", [
    (np.array([8.5, 8.6, 8.7, 8.8, 8.9999999999995]), "infer",
     np.array([8.5, 8.6, 8.7, 8.8, 8.9999999999995])),

    (np.array([8., 8., 8., 8., 8.9999999999995]), "infer",
     np.array([8, 8, 8, 8, 9], dtype=np.int64)),

    (np.array([8., 8., 8., 8., 9.0000000000005]), "infer",
     np.array([8, 8, 8, 8, 9], dtype=np.int64)),
])
def test_downcast(arr, expected, dtype):
    result = maybe_downcast_to_dtype(arr, dtype)
    tm.assert_numpy_array_equal(result, expected)


def test_downcast_booleans():
    # see gh-16875: coercing of booleans.
    ser = Series([True, True, False])
    result = maybe_downcast_to_dtype(ser, np.dtype(np.float64))

    expected = ser
    tm.assert_series_equal(result, expected)


def test_downcast_conversion_no_nan(any_real_dtype):
    dtype = any_real_dtype
    expected = np.array([1, 2])
    arr = np.array([1.0, 2.0], dtype=dtype)

    result = maybe_downcast_to_dtype(arr, "infer")
    tm.assert_almost_equal(result, expected, check_dtype=False)


def test_downcast_conversion_nan(float_dtype):
    dtype = float_dtype
    data = [1.0, 2.0, np.nan]

    expected = np.array(data, dtype=dtype)
    arr = np.array(data, dtype=dtype)

    result = maybe_downcast_to_dtype(arr, "infer")
    tm.assert_almost_equal(result, expected)


def test_downcast_conversion_empty(any_real_dtype):
    dtype = any_real_dtype
    arr = np.array([], dtype=dtype)
    result = maybe_downcast_to_dtype(arr, "int64")
    tm.assert_numpy_array_equal(result, np.array([], dtype=np.int64))


@pytest.mark.parametrize("klass", [np.datetime64, np.timedelta64])
def test_datetime_likes_nan(klass):
    dtype = klass.__name__ + "[ns]"
    arr = np.array([1, 2, np.nan])

    exp = np.array([1, 2, klass("NaT")], dtype)
    res = maybe_downcast_to_dtype(arr, dtype)
    tm.assert_numpy_array_equal(res, exp)


@pytest.mark.parametrize("as_asi", [True, False])
def test_datetime_with_timezone(as_asi):
    # see gh-15426
    ts = Timestamp("2016-01-01 12:00:00", tz="US/Pacific")
    exp = DatetimeIndex([ts, ts])

    obj = exp.asi8 if as_asi else exp
    res = maybe_downcast_to_dtype(obj, exp.dtype)

    tm.assert_index_equal(res, exp)
