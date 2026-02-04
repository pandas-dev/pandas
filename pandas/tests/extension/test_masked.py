"""
This file contains a minimal set of tests for compliance with the extension
array interface test suite, and should contain no other tests.
The test suite for the full functionality of the array is located in
`pandas/tests/arrays/`.

The tests in this file are inherited from the BaseExtensionTests, and only
minimal tweaks should be applied to get the tests passing (by overwriting a
parent method).

Additional tests should either be added to one of the BaseExtensionTests
classes (if they are relevant for the extension interface for all dtypes), or
be added to the array-specific tests in `pandas/tests/arrays/`.

"""

import numpy as np
import pytest

from pandas.compat import (
    IS64,
    is_platform_windows,
)
from pandas.compat.numpy import np_version_gt2

from pandas.core.dtypes.common import (
    is_float_dtype,
    is_signed_integer_dtype,
    is_unsigned_integer_dtype,
)

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays.boolean import BooleanDtype
from pandas.core.arrays.floating import (
    Float32Dtype,
    Float64Dtype,
)
from pandas.core.arrays.integer import (
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    UInt8Dtype,
    UInt16Dtype,
    UInt32Dtype,
    UInt64Dtype,
)
from pandas.tests.extension import base

is_windows_or_32bit = (is_platform_windows() and not np_version_gt2) or not IS64

pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:invalid value encountered in divide:RuntimeWarning"
    ),
    pytest.mark.filterwarnings("ignore:Mean of empty slice:RuntimeWarning"),
    # overflow only relevant for Floating dtype cases cases
    pytest.mark.filterwarnings("ignore:overflow encountered in reduce:RuntimeWarning"),
]


def make_data():
    return [1, 2, 3, 4, pd.NA, 10, 11, pd.NA, 99, 100]


def make_float_data():
    return [0.1, 0.2, 0.3, 0.4, pd.NA, 1.0, 1.1, pd.NA, 9.9, 10.0]


def make_bool_data():
    return [True, False] * 2 + [np.nan] + [True, False] + [np.nan] + [True, False]


@pytest.fixture(
    params=[
        Int8Dtype,
        Int16Dtype,
        Int32Dtype,
        Int64Dtype,
        UInt8Dtype,
        UInt16Dtype,
        UInt32Dtype,
        UInt64Dtype,
        Float32Dtype,
        Float64Dtype,
        BooleanDtype,
    ]
)
def dtype(request):
    return request.param()


@pytest.fixture
def data(dtype):
    if dtype.kind == "f":
        data = make_float_data()
    elif dtype.kind == "b":
        data = make_bool_data()
    else:
        data = make_data()
    return pd.array(data, dtype=dtype)


@pytest.fixture
def data_for_twos(dtype):
    if dtype.kind == "b":
        return pd.array(np.ones(10), dtype=dtype)
    return pd.array(np.ones(10) * 2, dtype=dtype)


@pytest.fixture
def data_missing(dtype):
    if dtype.kind == "f":
        return pd.array([pd.NA, 0.1], dtype=dtype)
    elif dtype.kind == "b":
        return pd.array([np.nan, True], dtype=dtype)
    return pd.array([pd.NA, 1], dtype=dtype)


@pytest.fixture
def data_for_sorting(dtype):
    if dtype.kind == "f":
        return pd.array([0.1, 0.2, 0.0], dtype=dtype)
    elif dtype.kind == "b":
        return pd.array([True, True, False], dtype=dtype)
    return pd.array([1, 2, 0], dtype=dtype)


@pytest.fixture
def data_missing_for_sorting(dtype):
    if dtype.kind == "f":
        return pd.array([0.1, pd.NA, 0.0], dtype=dtype)
    elif dtype.kind == "b":
        return pd.array([True, np.nan, False], dtype=dtype)
    return pd.array([1, pd.NA, 0], dtype=dtype)


@pytest.fixture
def na_cmp():
    # we are pd.NA
    return lambda x, y: x is pd.NA and y is pd.NA


@pytest.fixture
def data_for_grouping(dtype):
    if dtype.kind == "f":
        b = 0.1
        a = 0.0
        c = 0.2
    elif dtype.kind == "b":
        b = True
        a = False
        c = b
    else:
        b = 1
        a = 0
        c = 2

    na = pd.NA
    return pd.array([b, b, na, na, a, a, b, c], dtype=dtype)


class TestMaskedArrays(base.ExtensionTests):
    _combine_le_expected_dtype = "boolean"

    @pytest.fixture(autouse=True)
    def skip_if_doesnt_support_2d(self, dtype, request):
        # Override the fixture so that we run these tests.
        assert not dtype._supports_2d
        # If dtype._supports_2d is ever changed to True, then this fixture
        #  override becomes unnecessary.

    @pytest.mark.parametrize("na_action", [None, "ignore"])
    def test_map(self, data_missing, na_action, using_nan_is_na):
        result = data_missing.map(lambda x: x, na_action=na_action)
        if data_missing.dtype == Float32Dtype() and using_nan_is_na:
            # map roundtrips through objects, which converts to float64
            expected = data_missing.to_numpy(dtype="float64", na_value=np.nan)
        else:
            expected = data_missing.to_numpy()
        tm.assert_numpy_array_equal(result, expected)

    def test_map_na_action_ignore(self, data_missing_for_sorting, using_nan_is_na):
        zero = data_missing_for_sorting[2]
        result = data_missing_for_sorting.map(lambda x: zero, na_action="ignore")
        if data_missing_for_sorting.dtype.kind == "b":
            expected = np.array([False, pd.NA, False], dtype=object)
        elif not using_nan_is_na:
            # TODO: would we prefer to get NaN in this case to get a non-object?
            expected = np.array([zero, pd.NA, zero], dtype=object)
        else:
            expected = np.array([zero, np.nan, zero])
        tm.assert_numpy_array_equal(result, expected)

    def _get_expected_exception(self, op_name, obj, other):
        try:
            dtype = tm.get_dtype(obj)
        except AttributeError:
            # passed arguments reversed
            dtype = tm.get_dtype(other)

        if dtype.kind == "b":
            if op_name.strip("_").lstrip("r") in ["pow", "truediv", "floordiv"]:
                # match behavior with non-masked bool dtype
                return NotImplementedError
            elif op_name in ["__sub__", "__rsub__"]:
                # exception message would include "numpy boolean subtract""
                return TypeError
            return None
        return None

    def _cast_pointwise_result(self, op_name: str, obj, other, pointwise_result):
        sdtype = tm.get_dtype(obj)
        expected = pointwise_result

        if sdtype.kind == "b":
            if op_name in (
                "__mod__",
                "__rmod__",
            ):
                # combine keeps boolean type
                expected = expected.astype("Int8")

        return expected

    def test_divmod_series_array(self, data, data_for_twos, request):
        if data.dtype.kind == "b":
            mark = pytest.mark.xfail(
                reason="Inconsistency between floordiv and divmod; we raise for "
                "floordiv but not for divmod. This matches what we do for "
                "non-masked bool dtype."
            )
            request.applymarker(mark)
        super().test_divmod_series_array(data, data_for_twos)

    def _supports_reduction(self, ser: pd.Series, op_name: str) -> bool:
        if op_name in ["any", "all"] and ser.dtype.kind != "b":
            pytest.skip(reason="Tested in tests/reductions/test_reductions.py")
        return True

    def check_reduce(self, ser: pd.Series, op_name: str, skipna: bool):
        # overwrite to ensure pd.NA is tested instead of np.nan
        # https://github.com/pandas-dev/pandas/issues/30958

        cmp_dtype = "int64"
        if ser.dtype.kind == "f":
            # Item "dtype[Any]" of "Union[dtype[Any], ExtensionDtype]" has
            # no attribute "numpy_dtype"
            cmp_dtype = ser.dtype.numpy_dtype  # type: ignore[union-attr]
        elif ser.dtype.kind == "b":
            if op_name in ["min", "max"]:
                cmp_dtype = "bool"

        # TODO: prod with integer dtypes does *not* match the result we would
        #  get if we used object for cmp_dtype. In that cae the object result
        #  is a large integer while the non-object case overflows and returns 0
        alt = ser.dropna().astype(cmp_dtype)
        if op_name == "count":
            result = getattr(ser, op_name)()
            expected = getattr(alt, op_name)()
        else:
            result = getattr(ser, op_name)(skipna=skipna)
            expected = getattr(alt, op_name)(skipna=skipna)
            if not skipna and ser.isna().any() and op_name not in ["any", "all"]:
                expected = pd.NA
        tm.assert_almost_equal(result, expected)

    def _get_expected_reduction_dtype(self, arr, op_name: str, skipna: bool):
        if is_float_dtype(arr.dtype):
            cmp_dtype = arr.dtype.name
        elif op_name in ["mean", "median", "var", "std", "skew", "kurt", "sem"]:
            cmp_dtype = "Float64"
        elif op_name in ["max", "min"]:
            cmp_dtype = arr.dtype.name
        elif arr.dtype in ["Int64", "UInt64"]:
            cmp_dtype = arr.dtype.name
        elif is_signed_integer_dtype(arr.dtype):
            # TODO: Why does Window Numpy 2.0 dtype depend on skipna?
            cmp_dtype = (
                "Int32"
                if (is_platform_windows() and (not np_version_gt2 or not skipna))
                or not IS64
                else "Int64"
            )
        elif is_unsigned_integer_dtype(arr.dtype):
            cmp_dtype = (
                "UInt32"
                if (is_platform_windows() and (not np_version_gt2 or not skipna))
                or not IS64
                else "UInt64"
            )
        elif arr.dtype.kind == "b":
            if op_name in ["min", "max"]:
                cmp_dtype = "boolean"
            elif op_name in ["sum", "prod"]:
                cmp_dtype = (
                    "Int32"
                    if (is_platform_windows() and (not np_version_gt2 or not skipna))
                    or not IS64
                    else "Int64"
                )
            else:
                raise TypeError("not supposed to reach this")
        else:
            raise TypeError("not supposed to reach this")
        return cmp_dtype

    def _supports_accumulation(self, ser: pd.Series, op_name: str) -> bool:
        return True

    def check_accumulate(self, ser: pd.Series, op_name: str, skipna: bool):
        # overwrite to ensure pd.NA is tested instead of np.nan
        # https://github.com/pandas-dev/pandas/issues/30958
        length = 64
        if is_windows_or_32bit:
            # Item "ExtensionDtype" of "Union[dtype[Any], ExtensionDtype]" has
            # no attribute "itemsize"
            if not ser.dtype.itemsize == 8:  # type: ignore[union-attr]
                length = 32

        if ser.dtype.name.startswith("U"):
            expected_dtype = f"UInt{length}"
        elif ser.dtype.name.startswith("I"):
            expected_dtype = f"Int{length}"
        elif ser.dtype.name.startswith("F"):
            # Incompatible types in assignment (expression has type
            # "Union[dtype[Any], ExtensionDtype]", variable has type "str")
            expected_dtype = ser.dtype  # type: ignore[assignment]
        elif ser.dtype.kind == "b":
            if op_name in ("cummin", "cummax"):
                expected_dtype = "boolean"
            else:
                expected_dtype = f"Int{length}"

        if expected_dtype == "Float32" and op_name == "cumprod" and skipna:
            # TODO: xfail?
            pytest.skip(
                f"Float32 precision lead to large differences with op {op_name} "
                f"and skipna={skipna}"
            )

        if op_name == "cumsum":
            pass
        elif op_name in ["cummax", "cummin"]:
            expected_dtype = ser.dtype  # type: ignore[assignment]
        elif op_name == "cumprod":
            ser = ser[:12]
        else:
            raise NotImplementedError(f"{op_name} not supported")

        result = getattr(ser, op_name)(skipna=skipna)
        expected = pd.Series(
            pd.array(
                getattr(ser.astype("float64"), op_name)(skipna=skipna),
                dtype="Float64",
            )
        )
        expected[np.isnan(expected)] = pd.NA
        expected = expected.astype(expected_dtype)
        tm.assert_series_equal(result, expected)

    def test_loc_setitem_with_expansion_preserves_ea_index_dtype(self, data, request):
        super().test_loc_setitem_with_expansion_preserves_ea_index_dtype(data)


@pytest.mark.parametrize(
    "arr", [pd.array([True, False]), pd.array([1, 2]), pd.array([1.0, 2.0])]
)
def test_cast_pointwise_result_all_na_respects_original_dtype(arr):
    # GH#62344
    values = [pd.NA, pd.NA]
    result = arr._cast_pointwise_result(values)
    assert result.dtype == arr.dtype
    assert all(x is pd.NA for x in result)
