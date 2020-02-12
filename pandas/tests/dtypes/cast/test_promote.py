"""
These test the method maybe_promote from core/dtypes/cast.py
"""

import datetime

import numpy as np
import pytest

from pandas._libs.tslibs import NaT

from pandas.core.dtypes.cast import maybe_promote
from pandas.core.dtypes.common import (
    is_complex_dtype,
    is_datetime64_dtype,
    is_datetime_or_timedelta_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
    is_scalar,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.dtypes import DatetimeTZDtype
from pandas.core.dtypes.missing import isna

import pandas as pd


@pytest.fixture(
    params=[
        bool,
        "uint8",
        "int32",
        "uint64",
        "float32",
        "float64",
        "complex64",
        "complex128",
        "M8[ns]",
        "m8[ns]",
        str,
        bytes,
        object,
    ]
)
def any_numpy_dtype_reduced(request):
    """
    Parameterized fixture for numpy dtypes, reduced from any_numpy_dtype.

    * bool
    * 'int32'
    * 'uint64'
    * 'float32'
    * 'float64'
    * 'complex64'
    * 'complex128'
    * 'M8[ns]'
    * 'M8[ns]'
    * str
    * bytes
    * object
    """
    return request.param


def _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar=None):
    """
    Auxiliary function to unify testing of scalar/array promotion.

    Parameters
    ----------
    dtype : dtype
        The value to pass on as the first argument to maybe_promote.
    fill_value : scalar
        The value to pass on as the second argument to maybe_promote as
        a scalar.
    expected_dtype : dtype
        The expected dtype returned by maybe_promote (by design this is the
        same regardless of whether fill_value was passed as a scalar or in an
        array!).
    exp_val_for_scalar : scalar
        The expected value for the (potentially upcast) fill_value returned by
        maybe_promote.
    """
    assert is_scalar(fill_value)

    # here, we pass on fill_value as a scalar directly; the expected value
    # returned from maybe_promote is fill_value, potentially upcast to the
    # returned dtype.
    result_dtype, result_fill_value = maybe_promote(dtype, fill_value)
    expected_fill_value = exp_val_for_scalar

    assert result_dtype == expected_dtype
    _assert_match(result_fill_value, expected_fill_value)


def _assert_match(result_fill_value, expected_fill_value):
    # GH#23982/25425 require the same type in addition to equality/NA-ness
    res_type = type(result_fill_value)
    ex_type = type(expected_fill_value)
    if res_type.__name__ == "uint64":
        # No idea why, but these (sometimes) do not compare as equal
        assert ex_type.__name__ == "uint64"
    elif res_type.__name__ == "ulonglong":
        # On some builds we get this instead of np.uint64
        # Note: cant check res_type.dtype.itemsize directly on numpy 1.18
        assert res_type(0).itemsize == 8
        assert ex_type == res_type or ex_type == np.uint64
    else:
        # On some builds, type comparison fails, e.g. np.int32 != np.int32
        assert res_type == ex_type or res_type.__name__ == ex_type.__name__

    match_value = result_fill_value == expected_fill_value

    # Note: type check above ensures that we have the _same_ NA value
    # for missing values, None == None (which is checked
    # through match_value above), but np.nan != np.nan and pd.NaT != pd.NaT
    match_missing = isna(result_fill_value) and isna(expected_fill_value)

    assert match_value or match_missing


@pytest.mark.parametrize(
    "dtype, fill_value, expected_dtype",
    [
        # size 8
        ("int8", 1, "int8"),
        ("int8", np.iinfo("int8").max + 1, "int16"),
        ("int8", np.iinfo("int16").max + 1, "int32"),
        ("int8", np.iinfo("int32").max + 1, "int64"),
        ("int8", np.iinfo("int64").max + 1, "object"),
        ("int8", -1, "int8"),
        ("int8", np.iinfo("int8").min - 1, "int16"),
        ("int8", np.iinfo("int16").min - 1, "int32"),
        ("int8", np.iinfo("int32").min - 1, "int64"),
        ("int8", np.iinfo("int64").min - 1, "object"),
        # keep signed-ness as long as possible
        ("uint8", 1, "uint8"),
        ("uint8", np.iinfo("int8").max + 1, "uint8"),
        ("uint8", np.iinfo("uint8").max + 1, "uint16"),
        ("uint8", np.iinfo("int16").max + 1, "uint16"),
        ("uint8", np.iinfo("uint16").max + 1, "uint32"),
        ("uint8", np.iinfo("int32").max + 1, "uint32"),
        ("uint8", np.iinfo("uint32").max + 1, "uint64"),
        ("uint8", np.iinfo("int64").max + 1, "uint64"),
        ("uint8", np.iinfo("uint64").max + 1, "object"),
        # max of uint8 cannot be contained in int8
        ("uint8", -1, "int16"),
        ("uint8", np.iinfo("int8").min - 1, "int16"),
        ("uint8", np.iinfo("int16").min - 1, "int32"),
        ("uint8", np.iinfo("int32").min - 1, "int64"),
        ("uint8", np.iinfo("int64").min - 1, "object"),
        # size 16
        ("int16", 1, "int16"),
        ("int16", np.iinfo("int8").max + 1, "int16"),
        ("int16", np.iinfo("int16").max + 1, "int32"),
        ("int16", np.iinfo("int32").max + 1, "int64"),
        ("int16", np.iinfo("int64").max + 1, "object"),
        ("int16", -1, "int16"),
        ("int16", np.iinfo("int8").min - 1, "int16"),
        ("int16", np.iinfo("int16").min - 1, "int32"),
        ("int16", np.iinfo("int32").min - 1, "int64"),
        ("int16", np.iinfo("int64").min - 1, "object"),
        ("uint16", 1, "uint16"),
        ("uint16", np.iinfo("int8").max + 1, "uint16"),
        ("uint16", np.iinfo("uint8").max + 1, "uint16"),
        ("uint16", np.iinfo("int16").max + 1, "uint16"),
        ("uint16", np.iinfo("uint16").max + 1, "uint32"),
        ("uint16", np.iinfo("int32").max + 1, "uint32"),
        ("uint16", np.iinfo("uint32").max + 1, "uint64"),
        ("uint16", np.iinfo("int64").max + 1, "uint64"),
        ("uint16", np.iinfo("uint64").max + 1, "object"),
        ("uint16", -1, "int32"),
        ("uint16", np.iinfo("int8").min - 1, "int32"),
        ("uint16", np.iinfo("int16").min - 1, "int32"),
        ("uint16", np.iinfo("int32").min - 1, "int64"),
        ("uint16", np.iinfo("int64").min - 1, "object"),
        # size 32
        ("int32", 1, "int32"),
        ("int32", np.iinfo("int8").max + 1, "int32"),
        ("int32", np.iinfo("int16").max + 1, "int32"),
        ("int32", np.iinfo("int32").max + 1, "int64"),
        ("int32", np.iinfo("int64").max + 1, "object"),
        ("int32", -1, "int32"),
        ("int32", np.iinfo("int8").min - 1, "int32"),
        ("int32", np.iinfo("int16").min - 1, "int32"),
        ("int32", np.iinfo("int32").min - 1, "int64"),
        ("int32", np.iinfo("int64").min - 1, "object"),
        ("uint32", 1, "uint32"),
        ("uint32", np.iinfo("int8").max + 1, "uint32"),
        ("uint32", np.iinfo("uint8").max + 1, "uint32"),
        ("uint32", np.iinfo("int16").max + 1, "uint32"),
        ("uint32", np.iinfo("uint16").max + 1, "uint32"),
        ("uint32", np.iinfo("int32").max + 1, "uint32"),
        ("uint32", np.iinfo("uint32").max + 1, "uint64"),
        ("uint32", np.iinfo("int64").max + 1, "uint64"),
        ("uint32", np.iinfo("uint64").max + 1, "object"),
        ("uint32", -1, "int64"),
        ("uint32", np.iinfo("int8").min - 1, "int64"),
        ("uint32", np.iinfo("int16").min - 1, "int64"),
        ("uint32", np.iinfo("int32").min - 1, "int64"),
        ("uint32", np.iinfo("int64").min - 1, "object"),
        # size 64
        ("int64", 1, "int64"),
        ("int64", np.iinfo("int8").max + 1, "int64"),
        ("int64", np.iinfo("int16").max + 1, "int64"),
        ("int64", np.iinfo("int32").max + 1, "int64"),
        ("int64", np.iinfo("int64").max + 1, "object"),
        ("int64", -1, "int64"),
        ("int64", np.iinfo("int8").min - 1, "int64"),
        ("int64", np.iinfo("int16").min - 1, "int64"),
        ("int64", np.iinfo("int32").min - 1, "int64"),
        ("int64", np.iinfo("int64").min - 1, "object"),
        ("uint64", 1, "uint64"),
        ("uint64", np.iinfo("int8").max + 1, "uint64"),
        ("uint64", np.iinfo("uint8").max + 1, "uint64"),
        ("uint64", np.iinfo("int16").max + 1, "uint64"),
        ("uint64", np.iinfo("uint16").max + 1, "uint64"),
        ("uint64", np.iinfo("int32").max + 1, "uint64"),
        ("uint64", np.iinfo("uint32").max + 1, "uint64"),
        ("uint64", np.iinfo("int64").max + 1, "uint64"),
        ("uint64", np.iinfo("uint64").max + 1, "object"),
        ("uint64", -1, "object"),
        ("uint64", np.iinfo("int8").min - 1, "object"),
        ("uint64", np.iinfo("int16").min - 1, "object"),
        ("uint64", np.iinfo("int32").min - 1, "object"),
        ("uint64", np.iinfo("int64").min - 1, "object"),
    ],
)
def test_maybe_promote_int_with_int(dtype, fill_value, expected_dtype):
    dtype = np.dtype(dtype)
    expected_dtype = np.dtype(expected_dtype)

    # output is not a generic int, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_int_with_float(any_int_dtype, float_dtype):
    dtype = np.dtype(any_int_dtype)
    fill_dtype = np.dtype(float_dtype)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling int with float always upcasts to float64
    expected_dtype = np.float64
    # fill_value can be different float type
    exp_val_for_scalar = np.float64(fill_value)

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_float_with_int(float_dtype, any_int_dtype):

    dtype = np.dtype(float_dtype)
    fill_dtype = np.dtype(any_int_dtype)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling float with int always keeps float dtype
    # because: np.finfo('float32').max > np.iinfo('uint64').max
    expected_dtype = dtype
    # output is not a generic float, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


@pytest.mark.parametrize(
    "dtype, fill_value, expected_dtype",
    [
        # float filled with float
        ("float32", 1, "float32"),
        ("float32", np.finfo("float32").max * 1.1, "float64"),
        ("float64", 1, "float64"),
        ("float64", np.finfo("float32").max * 1.1, "float64"),
        # complex filled with float
        ("complex64", 1, "complex64"),
        ("complex64", np.finfo("float32").max * 1.1, "complex128"),
        ("complex128", 1, "complex128"),
        ("complex128", np.finfo("float32").max * 1.1, "complex128"),
        # float filled with complex
        ("float32", 1 + 1j, "complex64"),
        ("float32", np.finfo("float32").max * (1.1 + 1j), "complex128"),
        ("float64", 1 + 1j, "complex128"),
        ("float64", np.finfo("float32").max * (1.1 + 1j), "complex128"),
        # complex filled with complex
        ("complex64", 1 + 1j, "complex64"),
        ("complex64", np.finfo("float32").max * (1.1 + 1j), "complex128"),
        ("complex128", 1 + 1j, "complex128"),
        ("complex128", np.finfo("float32").max * (1.1 + 1j), "complex128"),
    ],
)
def test_maybe_promote_float_with_float(dtype, fill_value, expected_dtype):

    dtype = np.dtype(dtype)
    expected_dtype = np.dtype(expected_dtype)

    # output is not a generic float, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_bool_with_any(any_numpy_dtype_reduced):
    dtype = np.dtype(bool)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling bool with anything but bool casts to object
    expected_dtype = np.dtype(object) if fill_dtype != bool else fill_dtype
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_any_with_bool(any_numpy_dtype_reduced):
    dtype = np.dtype(any_numpy_dtype_reduced)
    fill_value = True

    # filling anything but bool with bool casts to object
    expected_dtype = np.dtype(object) if dtype != bool else dtype
    # output is not a generic bool, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_bytes_with_any(bytes_dtype, any_numpy_dtype_reduced):
    dtype = np.dtype(bytes_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # we never use bytes dtype internally, always promote to object
    expected_dtype = np.dtype(np.object_)
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_any_with_bytes(any_numpy_dtype_reduced, bytes_dtype):
    dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype
    fill_value = b"abc"

    # we never use bytes dtype internally, always promote to object
    expected_dtype = np.dtype(np.object_)
    # output is not a generic bytes, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_datetime64_with_any(datetime64_dtype, any_numpy_dtype_reduced):
    dtype = np.dtype(datetime64_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling datetime with anything but datetime casts to object
    if is_datetime64_dtype(fill_dtype):
        expected_dtype = dtype
        # for datetime dtypes, scalar values get cast to to_datetime64
        exp_val_for_scalar = pd.Timestamp(fill_value).to_datetime64()
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


@pytest.mark.parametrize(
    "fill_value",
    [
        pd.Timestamp("now"),
        np.datetime64("now"),
        datetime.datetime.now(),
        datetime.date.today(),
    ],
    ids=["pd.Timestamp", "np.datetime64", "datetime.datetime", "datetime.date"],
)
def test_maybe_promote_any_with_datetime64(
    any_numpy_dtype_reduced, datetime64_dtype, fill_value
):
    dtype = np.dtype(any_numpy_dtype_reduced)

    # filling datetime with anything but datetime casts to object
    if is_datetime64_dtype(dtype):
        expected_dtype = dtype
        # for datetime dtypes, scalar values get cast to pd.Timestamp.value
        exp_val_for_scalar = pd.Timestamp(fill_value).to_datetime64()
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_datetimetz_with_any_numpy_dtype(
    tz_aware_fixture, any_numpy_dtype_reduced
):
    dtype = DatetimeTZDtype(tz=tz_aware_fixture)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling datetimetz with any numpy dtype casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_datetimetz_with_datetimetz(tz_aware_fixture, tz_aware_fixture2):
    dtype = DatetimeTZDtype(tz=tz_aware_fixture)
    fill_dtype = DatetimeTZDtype(tz=tz_aware_fixture2)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = pd.Series([10 ** 9], dtype=fill_dtype)[0]

    # filling datetimetz with datetimetz casts to object, unless tz matches
    exp_val_for_scalar = fill_value
    if dtype.tz == fill_dtype.tz:
        expected_dtype = dtype
    else:
        expected_dtype = np.dtype(object)

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


@pytest.mark.parametrize("fill_value", [None, np.nan, NaT])
def test_maybe_promote_datetimetz_with_na(tz_aware_fixture, fill_value):

    dtype = DatetimeTZDtype(tz=tz_aware_fixture)

    expected_dtype = dtype
    exp_val_for_scalar = NaT

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


@pytest.mark.parametrize(
    "fill_value",
    [
        pd.Timestamp("now"),
        np.datetime64("now"),
        datetime.datetime.now(),
        datetime.date.today(),
    ],
    ids=["pd.Timestamp", "np.datetime64", "datetime.datetime", "datetime.date"],
)
def test_maybe_promote_any_numpy_dtype_with_datetimetz(
    any_numpy_dtype_reduced, tz_aware_fixture, fill_value
):
    dtype = np.dtype(any_numpy_dtype_reduced)
    fill_dtype = DatetimeTZDtype(tz=tz_aware_fixture)

    fill_value = pd.Series([fill_value], dtype=fill_dtype)[0]

    # filling any numpy dtype with datetimetz casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_timedelta64_with_any(timedelta64_dtype, any_numpy_dtype_reduced):
    dtype = np.dtype(timedelta64_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling timedelta with anything but timedelta casts to object
    if is_timedelta64_dtype(fill_dtype):
        expected_dtype = dtype
        # for timedelta dtypes, scalar values get cast to pd.Timedelta.value
        exp_val_for_scalar = pd.Timedelta(fill_value).to_timedelta64()
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


@pytest.mark.parametrize(
    "fill_value",
    [pd.Timedelta(days=1), np.timedelta64(24, "h"), datetime.timedelta(1)],
    ids=["pd.Timedelta", "np.timedelta64", "datetime.timedelta"],
)
def test_maybe_promote_any_with_timedelta64(
    any_numpy_dtype_reduced, timedelta64_dtype, fill_value
):
    dtype = np.dtype(any_numpy_dtype_reduced)

    # filling anything but timedelta with timedelta casts to object
    if is_timedelta64_dtype(dtype):
        expected_dtype = dtype
        # for timedelta dtypes, scalar values get cast to pd.Timedelta.value
        exp_val_for_scalar = pd.Timedelta(fill_value).to_timedelta64()
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_string_with_any(string_dtype, any_numpy_dtype_reduced):
    dtype = np.dtype(string_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling string with anything casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_any_with_string(any_numpy_dtype_reduced, string_dtype):
    dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype
    fill_value = "abc"

    # filling anything with a string casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_object_with_any(object_dtype, any_numpy_dtype_reduced):
    dtype = np.dtype(object_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling object with anything stays object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


def test_maybe_promote_any_with_object(any_numpy_dtype_reduced, object_dtype):
    dtype = np.dtype(any_numpy_dtype_reduced)

    # create array of object dtype from a scalar value (i.e. passing
    # dtypes.common.is_scalar), which can however not be cast to int/float etc.
    fill_value = pd.DateOffset(1)

    # filling object with anything stays object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


@pytest.mark.parametrize("fill_value", [None, np.nan, NaT])
def test_maybe_promote_any_numpy_dtype_with_na(any_numpy_dtype_reduced, fill_value):
    dtype = np.dtype(any_numpy_dtype_reduced)

    if is_integer_dtype(dtype) and fill_value is not NaT:
        # integer + other missing value (np.nan / None) casts to float
        expected_dtype = np.float64
        exp_val_for_scalar = np.nan
    elif is_object_dtype(dtype) and fill_value is NaT:
        # inserting into object does not cast the value
        # but *does* cast None to np.nan
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value
    elif is_datetime_or_timedelta_dtype(dtype):
        # datetime / timedelta cast all missing values to dtyped-NaT
        expected_dtype = dtype
        exp_val_for_scalar = dtype.type("NaT", "ns")
    elif fill_value is NaT:
        # NaT upcasts everything that's not datetime/timedelta to object
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = NaT
    elif is_float_dtype(dtype) or is_complex_dtype(dtype):
        # float / complex + missing value (!= NaT) stays the same
        expected_dtype = dtype
        exp_val_for_scalar = np.nan
    else:
        # all other cases cast to object, and use np.nan as missing value
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = np.nan

    _check_promote(dtype, fill_value, expected_dtype, exp_val_for_scalar)


@pytest.mark.parametrize("dim", [0, 2, 3])
def test_maybe_promote_dimensions(any_numpy_dtype_reduced, dim):
    dtype = np.dtype(any_numpy_dtype_reduced)

    # create 0-dim array of given dtype; casts "1" to correct dtype
    fill_array = np.array(1, dtype=dtype)

    # expand to desired dimension:
    for _ in range(dim):
        fill_array = np.expand_dims(fill_array, 0)

    if dtype != object:
        # test against 1-dimensional case
        with pytest.raises(ValueError, match="fill_value must be a scalar"):
            maybe_promote(dtype, np.array([1], dtype=dtype))

        with pytest.raises(ValueError, match="fill_value must be a scalar"):
            maybe_promote(dtype, fill_array)

    else:
        expected_dtype, expected_missing_value = maybe_promote(
            dtype, np.array([1], dtype=dtype)
        )
        result_dtype, result_missing_value = maybe_promote(dtype, fill_array)
        assert result_dtype == expected_dtype
        _assert_match(result_missing_value, expected_missing_value)
