"""
These test the method maybe_promote from core/dtypes/cast.py
"""

import datetime

import numpy as np
import pytest

from pandas._libs.tslibs import NaT, iNaT
from pandas.compat import is_platform_windows

from pandas.core.dtypes.cast import (
    _maybe_promote_with_array,
    _maybe_promote_with_scalar,
    maybe_promote,
)
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
from pandas.core.dtypes.dtypes import DatetimeTZDtype, PandasExtensionDtype

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


@pytest.fixture(
    params=[(True, None), (True, object), (False, None)],
    ids=["True-None", "True-object", "False-None"],
)
def box(request):
    """
    Parametrized fixture determining whether/how to transform fill_value.

    Since fill_value is defined on a per-test basis, the actual transformation
    (based on this fixture) is executed in _check_promote.

    Returns
    -------
    boxed : Boolean
        Whether fill_value should be wrapped in an np.array.
    box_dtype : dtype
        The dtype to pass to np.array([fill_value], dtype=box_dtype). If None,
        then this is passed on unmodified, and corresponds to the numpy default
        dtype for the given fill_value.

    * (True, None)    # fill_value wrapped in array with default dtype
    * (True, object)  # fill_value wrapped in array with object dtype
    * (False, None)   # fill_value passed on as scalar
    """
    return request.param


def _safe_dtype_assert(left_dtype, right_dtype):
    """
    Compare two dtypes without raising TypeError.
    """
    if isinstance(right_dtype, PandasExtensionDtype):
        # switch order of equality check because numpy dtypes (e.g. if
        # left_dtype is np.object_) do not know some expected dtypes (e.g.
        # DatetimeTZDtype) and would raise a TypeError in their __eq__-method.
        assert right_dtype == left_dtype
    else:
        assert left_dtype == right_dtype


def _check_promote(
    dtype,
    fill_value,
    boxed,
    box_dtype,
    expected_dtype,
    exp_val_for_scalar=None,
    exp_val_for_array=None,
):
    """
    Auxiliary function to unify testing of scalar/array promotion.

    Parameters
    ----------
    dtype : dtype
        The value to pass on as the first argument to maybe_promote.
    fill_value : scalar
        The value to pass on as the second argument to maybe_promote, either as
        a scalar, or boxed into an array (depending on the parameter `boxed`).
    boxed : Boolean
        Parameter whether fill_value should be passed to maybe_promote
        directly, or wrapped in an array (of dtype box_dtype).
    box_dtype : dtype
        The dtype to enforce when wrapping fill_value into an np.array.
    expected_dtype : dtype
        The expected dtype returned by maybe_promote (by design this is the
        same regardless of whether fill_value was passed as a scalar or in an
        array!).
    exp_val_for_scalar : scalar
        The expected value for the (potentially upcast) fill_value returned by
        maybe_promote.
    exp_val_for_array : scalar
        The expected missing value marker for the expected_dtype (which is
        returned by maybe_promote when it receives an array).
    """
    assert is_scalar(fill_value)

    if boxed:
        # in this case, we pass on fill_value wrapped in an array of specified
        # box_dtype; the expected value returned from maybe_promote is the
        # missing value marker for the returned dtype.
        fill_array = np.array([fill_value], dtype=box_dtype)
        result_dtype, result_fill_value = maybe_promote(dtype, fill_array)
        expected_fill_value = exp_val_for_array
    else:
        # here, we pass on fill_value as a scalar directly; the expected value
        # returned from maybe_promote is fill_value, potentially upcast to the
        # returned dtype.
        result_dtype, result_fill_value = maybe_promote(dtype, fill_value)
        expected_fill_value = exp_val_for_scalar

    _safe_dtype_assert(result_dtype, expected_dtype)

    # for equal values, also check type (relevant e.g. for int vs float, resp.
    # for different datetimes and timedeltas)
    match_value = result_fill_value == expected_fill_value and type(
        result_fill_value
    ) == type(expected_fill_value)

    # for missing values, None == None and iNaT == iNaT (which is checked
    # through match_value above), but np.nan != np.nan and pd.NaT != pd.NaT
    match_missing = (result_fill_value is np.nan and expected_fill_value is np.nan) or (
        result_fill_value is NaT and expected_fill_value is NaT
    )

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
def test_maybe_promote_int_with_int(dtype, fill_value, expected_dtype, box):
    dtype = np.dtype(dtype)
    expected_dtype = np.dtype(expected_dtype)
    boxed, box_dtype = box  # read from parametrized fixture

    # output is not a python int, but a numpy int of expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]
    # no missing value marker for integers
    exp_val_for_array = None if expected_dtype != "object" else np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_int_with_float(any_int_dtype, float_dtype, box):
    dtype = np.dtype(any_int_dtype)
    fill_dtype = np.dtype(float_dtype)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling int with float always upcasts to float64
    expected_dtype = np.float64
    # fill_value can be different float type
    exp_val_for_scalar = np.float64(fill_value)
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_float_with_int(float_dtype, any_int_dtype, box):

    dtype = np.dtype(float_dtype)
    fill_dtype = np.dtype(any_int_dtype)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling float with int always keeps float dtype
    # because: np.finfo('float32').max > np.iinfo('uint64').max
    expected_dtype = dtype
    # output is not a generic float, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


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
def test_maybe_promote_float_with_float(dtype, fill_value, expected_dtype, box):

    dtype = np.dtype(dtype)
    expected_dtype = np.dtype(expected_dtype)
    boxed, box_dtype = box  # read from parametrized fixture

    # output is not a generic float, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_bool_with_any(any_numpy_dtype_reduced, box):
    dtype = np.dtype(bool)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling bool with anything but bool casts to object
    expected_dtype = np.dtype(object) if fill_dtype != bool else fill_dtype
    exp_val_for_scalar = fill_value
    exp_val_for_array = np.nan if fill_dtype != bool else None

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_any_with_bool(any_numpy_dtype_reduced, box):
    dtype = np.dtype(any_numpy_dtype_reduced)
    fill_value = True
    boxed, box_dtype = box  # read from parametrized fixture

    # filling anything but bool with bool casts to object
    expected_dtype = np.dtype(object) if dtype != bool else dtype
    # output is not a generic bool, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]
    exp_val_for_array = np.nan if dtype != bool else None

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_bytes_with_any(bytes_dtype, any_numpy_dtype_reduced, box):
    dtype = np.dtype(bytes_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling bytes with anything but bytes casts to object
    expected_dtype = (
        dtype if issubclass(fill_dtype.type, np.bytes_) else np.dtype(object)
    )
    exp_val_for_scalar = fill_value
    exp_val_for_array = None if issubclass(fill_dtype.type, np.bytes_) else np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


# override parametrization of box to add special case for bytes
@pytest.mark.parametrize(
    "box",
    [
        (True, None),  # fill_value wrapped in array with default dtype
        (True, "bytes"),  # fill_value in array with generic bytes dtype
        (True, object),  # fill_value wrapped in array with object dtype
        (False, None),  # fill_value passed on as scalar
    ],
    ids=["True-None", "True-bytes", "True-object", "False-None"],
)
def test_maybe_promote_any_with_bytes(any_numpy_dtype_reduced, bytes_dtype, box):
    dtype = np.dtype(any_numpy_dtype_reduced)
    fill_dtype = np.dtype(bytes_dtype)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype
    fill_value = b"abc"

    # special case for box_dtype (cannot use fixture in parametrization)
    box_dtype = fill_dtype if box_dtype == "bytes" else box_dtype

    # filling bytes with anything but bytes casts to object
    expected_dtype = dtype if issubclass(dtype.type, np.bytes_) else np.dtype(object)
    # output is not a generic bytes, but corresponds to expected_dtype
    exp_val_for_scalar = np.array([fill_value], dtype=expected_dtype)[0]
    exp_val_for_array = None if issubclass(dtype.type, np.bytes_) else np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_datetime64_with_any(
    datetime64_dtype, any_numpy_dtype_reduced, box
):
    dtype = np.dtype(datetime64_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling datetime with anything but datetime casts to object
    if is_datetime64_dtype(fill_dtype):
        expected_dtype = dtype
        # for datetime dtypes, scalar values get cast to pd.Timestamp.value
        exp_val_for_scalar = pd.Timestamp(fill_value).value
        exp_val_for_array = iNaT
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value
        exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


# override parametrization of box to add special case for dt_dtype
@pytest.mark.parametrize(
    "box",
    [
        (True, None),  # fill_value wrapped in array with default dtype
        (True, "dt_dtype"),  # fill_value in array with explicit datetime dtype
        (True, object),  # fill_value wrapped in array with object dtype
        (False, None),  # fill_value passed on as scalar
    ],
    ids=["True-None", "True-dt_dtype", "True-object", "False-None"],
)
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
    any_numpy_dtype_reduced, datetime64_dtype, fill_value, box
):
    dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # special case for box_dtype
    box_dtype = np.dtype(datetime64_dtype) if box_dtype == "dt_dtype" else box_dtype

    # filling datetime with anything but datetime casts to object
    if is_datetime64_dtype(dtype):
        expected_dtype = dtype
        # for datetime dtypes, scalar values get cast to pd.Timestamp.value
        exp_val_for_scalar = pd.Timestamp(fill_value).value
        exp_val_for_array = iNaT
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value
        exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_datetimetz_with_any_numpy_dtype(
    tz_aware_fixture, any_numpy_dtype_reduced, box
):
    dtype = DatetimeTZDtype(tz=tz_aware_fixture)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling datetimetz with any numpy dtype casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_datetimetz_with_datetimetz(
    tz_aware_fixture, tz_aware_fixture2, box
):
    dtype = DatetimeTZDtype(tz=tz_aware_fixture)
    fill_dtype = DatetimeTZDtype(tz=tz_aware_fixture2)
    boxed, box_dtype = box  # read from parametrized fixture

    from dateutil.tz import tzlocal

    if is_platform_windows() and tz_aware_fixture2 == tzlocal():
        pytest.xfail("Cannot process fill_value with this dtype, see GH 24310")
    if dtype.tz == fill_dtype.tz:
        # here we should keep the datetime64tz dtype, but since that cannot be
        # inferred correctly for fill_value, the calling dtype ends up being
        # compared to a tz-naive datetime64-dtype, and must therefore upcast
        pytest.xfail("cannot infer datetime64tz dtype, see GH 23554")

    # create array of given dtype; casts "1" to correct dtype
    fill_value = pd.Series([10 ** 9], dtype=fill_dtype)[0]

    # filling datetimetz with datetimetz casts to object, unless tz matches
    exp_val_for_scalar = fill_value
    if dtype.tz == fill_dtype.tz:
        expected_dtype = dtype
        exp_val_for_array = NaT
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


@pytest.mark.parametrize(
    "fill_value", [None, np.nan, NaT, iNaT], ids=["None", "np.nan", "pd.NaT", "iNaT"]
)
def test_maybe_promote_datetimetz_with_na(tz_aware_fixture, fill_value, box):

    dtype = DatetimeTZDtype(tz=tz_aware_fixture)
    boxed, box_dtype = box  # read from parametrized fixture

    expected_dtype = dtype
    # DatetimeTZDtype does not use iNaT as missing value marker
    exp_val_for_scalar = NaT
    exp_val_for_array = NaT

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


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
    any_numpy_dtype_reduced, tz_aware_fixture, fill_value, box
):
    dtype = np.dtype(any_numpy_dtype_reduced)
    fill_dtype = DatetimeTZDtype(tz=tz_aware_fixture)
    boxed, box_dtype = box  # read from parametrized fixture

    if is_datetime64_dtype(dtype):
        # fill_dtype does not get inferred correctly to datetime64tz but to
        # datetime64, which then falsely matches with datetime64 dtypes.
        pytest.xfail("cannot infer datetime64tz dtype, see GH 23554")

    fill_value = pd.Series([fill_value], dtype=fill_dtype)[0]

    # filling any numpy dtype with datetimetz casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_timedelta64_with_any(
    timedelta64_dtype, any_numpy_dtype_reduced, box
):
    dtype = np.dtype(timedelta64_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling timedelta with anything but timedelta casts to object
    if is_timedelta64_dtype(fill_dtype):
        expected_dtype = dtype
        # for timedelta dtypes, scalar values get cast to pd.Timedelta.value
        exp_val_for_scalar = pd.Timedelta(fill_value).value
        exp_val_for_array = iNaT
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value
        exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


@pytest.mark.parametrize(
    "fill_value",
    [pd.Timedelta(days=1), np.timedelta64(24, "h"), datetime.timedelta(1)],
    ids=["pd.Timedelta", "np.timedelta64", "datetime.timedelta"],
)
# override parametrization of box to add special case for td_dtype
@pytest.mark.parametrize(
    "box",
    [
        (True, None),  # fill_value wrapped in array with default dtype
        (True, "td_dtype"),  # fill_value in array with explicit timedelta dtype
        (True, object),  # fill_value wrapped in array with object dtype
        (False, None),  # fill_value passed on as scalar
    ],
    ids=["True-None", "True-td_dtype", "True-object", "False-None"],
)
def test_maybe_promote_any_with_timedelta64(
    any_numpy_dtype_reduced, timedelta64_dtype, fill_value, box
):
    dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # special case for box_dtype
    box_dtype = np.dtype(timedelta64_dtype) if box_dtype == "td_dtype" else box_dtype

    # filling anything but timedelta with timedelta casts to object
    if is_timedelta64_dtype(dtype):
        expected_dtype = dtype
        # for timedelta dtypes, scalar values get cast to pd.Timedelta.value
        exp_val_for_scalar = pd.Timedelta(fill_value).value
        exp_val_for_array = iNaT
    else:
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value
        exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_string_with_any(string_dtype, any_numpy_dtype_reduced, box):
    dtype = np.dtype(string_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling string with anything casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


# override parametrization of box to add special case for str
@pytest.mark.parametrize(
    "box",
    [
        # disabled due to too many xfails; see GH 23982 / 25425
        (True, None),  # fill_value wrapped in array with default dtype
        (True, "str"),  # fill_value wrapped in array with generic string-dtype
        (True, object),  # fill_value wrapped in array with object dtype
        (False, None),  # fill_value passed on as scalar
    ],
    ids=["True-None", "True-str", "True-object", "False-None"],
)
def test_maybe_promote_any_with_string(any_numpy_dtype_reduced, string_dtype, box):
    dtype = np.dtype(any_numpy_dtype_reduced)
    fill_dtype = np.dtype(string_dtype)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype
    fill_value = "abc"

    # special case for box_dtype (cannot use fixture in parametrization)
    box_dtype = fill_dtype if box_dtype == "str" else box_dtype

    # filling anything with a string casts to object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_object_with_any(object_dtype, any_numpy_dtype_reduced, box):
    dtype = np.dtype(object_dtype)
    fill_dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of given dtype; casts "1" to correct dtype
    fill_value = np.array([1], dtype=fill_dtype)[0]

    # filling object with anything stays object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


def test_maybe_promote_any_with_object(any_numpy_dtype_reduced, object_dtype, box):
    dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    # create array of object dtype from a scalar value (i.e. passing
    # dtypes.common.is_scalar), which can however not be cast to int/float etc.
    fill_value = pd.DateOffset(1)

    # filling object with anything stays object
    expected_dtype = np.dtype(object)
    exp_val_for_scalar = fill_value
    exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


@pytest.mark.parametrize(
    "fill_value", [None, np.nan, NaT, iNaT], ids=["None", "np.nan", "pd.NaT", "iNaT"]
)
# override parametrization of box, because default dtype for na is always float
@pytest.mark.parametrize(
    "box",
    [
        (True, object),  # fill_value wrapped in array with object dtype
        (False, None),  # fill_value passed on as scalar
    ],
    ids=["True-object", "False-None"],
)
def test_maybe_promote_any_numpy_dtype_with_na(
    any_numpy_dtype_reduced, fill_value, box
):
    dtype = np.dtype(any_numpy_dtype_reduced)
    boxed, box_dtype = box  # read from parametrized fixture

    if is_integer_dtype(dtype) and dtype == "uint64" and fill_value == iNaT:
        # uint64 + negative int casts to object; iNaT is considered as missing
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = np.nan
    elif is_integer_dtype(dtype) and fill_value == iNaT:
        # other integer + iNaT casts to int64
        expected_dtype = np.int64
        exp_val_for_scalar = iNaT
    elif is_integer_dtype(dtype) and fill_value is not NaT:
        # integer + other missing value (np.nan / None) casts to float
        expected_dtype = np.float64
        exp_val_for_scalar = np.nan
    elif is_object_dtype(dtype) and (fill_value == iNaT or fill_value is NaT):
        # inserting into object does not cast the value
        # but *does* cast None to np.nan
        expected_dtype = np.dtype(object)
        exp_val_for_scalar = fill_value
    elif is_datetime_or_timedelta_dtype(dtype):
        # datetime / timedelta cast all missing values to iNaT
        expected_dtype = dtype
        exp_val_for_scalar = iNaT
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

    # array case has same expected_dtype; but returns corresponding na-marker
    if is_integer_dtype(expected_dtype):
        # integers cannot hold NaNs; _maybe_promote_with_array returns None
        exp_val_for_array = None
    elif is_datetime_or_timedelta_dtype(expected_dtype):
        exp_val_for_array = iNaT
    else:  # expected_dtype = float / complex / object
        exp_val_for_array = np.nan

    _check_promote(
        dtype,
        fill_value,
        boxed,
        box_dtype,
        expected_dtype,
        exp_val_for_scalar,
        exp_val_for_array,
    )


@pytest.mark.parametrize("dim", [0, 2, 3])
def test_maybe_promote_dimensions(any_numpy_dtype_reduced, dim):
    dtype = np.dtype(any_numpy_dtype_reduced)

    # create 0-dim array of given dtype; casts "1" to correct dtype
    fill_array = np.array(1, dtype=dtype)

    # expand to desired dimension:
    for _ in range(dim):
        fill_array = np.expand_dims(fill_array, 0)

    # test against 1-dimensional case
    expected_dtype, expected_missing_value = maybe_promote(
        dtype, np.array([1], dtype=dtype)
    )

    result_dtype, result_missing_value = maybe_promote(dtype, fill_array)

    assert result_dtype == expected_dtype
    # None == None, iNaT == iNaT, but np.nan != np.nan
    assert (result_missing_value == expected_missing_value) or (
        result_missing_value is np.nan and expected_missing_value is np.nan
    )

    # same again for _maybe_promote_with_array (for coverage)
    result_dtype, result_missing_value = _maybe_promote_with_array(dtype, fill_array)

    assert result_dtype == expected_dtype
    # None == None, iNaT == iNaT, but np.nan != np.nan
    assert (result_missing_value == expected_missing_value) or (
        result_missing_value is np.nan and expected_missing_value is np.nan
    )


def test_maybe_promote_raises(any_numpy_dtype):
    msg = "fill_value must either be scalar, or a Series / Index / np.ndarra.*"
    with pytest.raises(ValueError, match=msg):
        # something that's neither scalar, nor Series / Index / np.ndarray
        maybe_promote(any_numpy_dtype, [1, 2, 3])

    msg = "fill_value must either be a Series / Index / np.ndarray, received.*"
    with pytest.raises(ValueError, match=msg):
        # something that's not a Series / Index / np.ndarray
        _maybe_promote_with_array(any_numpy_dtype, 1)

    msg = "fill_value must be a scalar, received .*"
    with pytest.raises(ValueError, match=msg):
        # something that's not scalar
        _maybe_promote_with_scalar(any_numpy_dtype, pd.Series([1, 2, 3]))
