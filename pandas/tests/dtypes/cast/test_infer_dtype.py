from datetime import date, datetime, timedelta

import numpy as np
import pytest

from pandas.core.dtypes.cast import (
    cast_scalar_to_array,
    infer_dtype_from_array,
    infer_dtype_from_scalar,
)
from pandas.core.dtypes.common import is_dtype_equal

from pandas import Categorical, Period, Series, Timedelta, Timestamp, date_range
from pandas.util import testing as tm


@pytest.fixture(params=[True, False])
def pandas_dtype(request):
    return request.param


def test_infer_dtype_from_int_scalar(any_int_dtype):
    # Test that infer_dtype_from_scalar is
    # returning correct dtype for int and float.
    data = np.dtype(any_int_dtype).type(12)
    dtype, val = infer_dtype_from_scalar(data)
    assert dtype == type(data)


def test_infer_dtype_from_float_scalar(float_dtype):
    float_dtype = np.dtype(float_dtype).type
    data = float_dtype(12)

    dtype, val = infer_dtype_from_scalar(data)
    assert dtype == float_dtype


@pytest.mark.parametrize("data,exp_dtype", [(12, np.int64), (np.float(12), np.float64)])
def test_infer_dtype_from_python_scalar(data, exp_dtype):
    dtype, val = infer_dtype_from_scalar(data)
    assert dtype == exp_dtype


@pytest.mark.parametrize("bool_val", [True, False])
def test_infer_dtype_from_boolean(bool_val):
    dtype, val = infer_dtype_from_scalar(bool_val)
    assert dtype == np.bool_


def test_infer_dtype_from_complex(complex_dtype):
    data = np.dtype(complex_dtype).type(1)
    dtype, val = infer_dtype_from_scalar(data)
    assert dtype == np.complex_


@pytest.mark.parametrize(
    "data", [np.datetime64(1, "ns"), Timestamp(1), datetime(2000, 1, 1, 0, 0)]
)
def test_infer_dtype_from_datetime(data):
    dtype, val = infer_dtype_from_scalar(data)
    assert dtype == "M8[ns]"


@pytest.mark.parametrize("data", [np.timedelta64(1, "ns"), Timedelta(1), timedelta(1)])
def test_infer_dtype_from_timedelta(data):
    dtype, val = infer_dtype_from_scalar(data)
    assert dtype == "m8[ns]"


@pytest.mark.parametrize("freq", ["M", "D"])
def test_infer_dtype_from_period(freq, pandas_dtype):
    p = Period("2011-01-01", freq=freq)
    dtype, val = infer_dtype_from_scalar(p, pandas_dtype=pandas_dtype)

    if pandas_dtype:
        exp_dtype = "period[{0}]".format(freq)
        exp_val = p.ordinal
    else:
        exp_dtype = np.object_
        exp_val = p

    assert dtype == exp_dtype
    assert val == exp_val


@pytest.mark.parametrize(
    "data", [date(2000, 1, 1), "foo", Timestamp(1, tz="US/Eastern")]
)
def test_infer_dtype_misc(data):
    dtype, val = infer_dtype_from_scalar(data)
    assert dtype == np.object_


@pytest.mark.parametrize("tz", ["UTC", "US/Eastern", "Asia/Tokyo"])
def test_infer_from_scalar_tz(tz, pandas_dtype):
    dt = Timestamp(1, tz=tz)
    dtype, val = infer_dtype_from_scalar(dt, pandas_dtype=pandas_dtype)

    if pandas_dtype:
        exp_dtype = "datetime64[ns, {0}]".format(tz)
        exp_val = dt.value
    else:
        exp_dtype = np.object_
        exp_val = dt

    assert dtype == exp_dtype
    assert val == exp_val


def test_infer_dtype_from_scalar_errors():
    msg = "invalid ndarray passed to infer_dtype_from_scalar"

    with pytest.raises(ValueError, match=msg):
        infer_dtype_from_scalar(np.array([1]))


@pytest.mark.parametrize(
    "arr, expected, pandas_dtype",
    [
        ("foo", np.object_, False),
        (b"foo", np.object_, False),
        (1, np.int_, False),
        (1.5, np.float_, False),
        ([1], np.int_, False),
        (np.array([1], dtype=np.int64), np.int64, False),
        ([np.nan, 1, ""], np.object_, False),
        (np.array([[1.0, 2.0]]), np.float_, False),
        (Categorical(list("aabc")), np.object_, False),
        (Categorical([1, 2, 3]), np.int64, False),
        (Categorical(list("aabc")), "category", True),
        (Categorical([1, 2, 3]), "category", True),
        (Timestamp("20160101"), np.object_, False),
        (np.datetime64("2016-01-01"), np.dtype("=M8[D]"), False),
        (date_range("20160101", periods=3), np.dtype("=M8[ns]"), False),
        (
            date_range("20160101", periods=3, tz="US/Eastern"),
            "datetime64[ns, US/Eastern]",
            True,
        ),
        (Series([1.0, 2, 3]), np.float64, False),
        (Series(list("abc")), np.object_, False),
        (
            Series(date_range("20160101", periods=3, tz="US/Eastern")),
            "datetime64[ns, US/Eastern]",
            True,
        ),
    ],
)
def test_infer_dtype_from_array(arr, expected, pandas_dtype):
    dtype, _ = infer_dtype_from_array(arr, pandas_dtype=pandas_dtype)
    assert is_dtype_equal(dtype, expected)


@pytest.mark.parametrize(
    "obj,dtype",
    [
        (1, np.int64),
        (1.1, np.float64),
        (Timestamp("2011-01-01"), "datetime64[ns]"),
        (Timestamp("2011-01-01", tz="US/Eastern"), np.object),
        (Period("2011-01-01", freq="D"), np.object),
    ],
)
def test_cast_scalar_to_array(obj, dtype):
    shape = (3, 2)

    exp = np.empty(shape, dtype=dtype)
    exp.fill(obj)

    arr = cast_scalar_to_array(shape, obj, dtype=dtype)
    tm.assert_numpy_array_equal(arr, exp)
