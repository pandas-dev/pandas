import re

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


class TestSetitemValidation:
    def _check_setitem_invalid(self, arr, invalid):
        msg = f"Invalid value '{invalid!s}' for dtype '{arr.dtype}'"
        msg = re.escape(msg)
        with pytest.raises(TypeError, match=msg):
            arr[0] = invalid

        with pytest.raises(TypeError, match=msg):
            arr[:] = invalid

        with pytest.raises(TypeError, match=msg):
            arr[[0]] = invalid

        with pytest.raises(TypeError):
            arr[[0]] = [invalid]

        with pytest.raises(TypeError):
            arr[[0]] = np.array([invalid], dtype=object)

        # Series non-coercion, behavior subject to change
        ser = pd.Series(arr)
        with pytest.raises(TypeError, match=msg):
            ser[0] = invalid
            # TODO: so, so many other variants of this...

    _invalid_scalars = [
        1 + 2j,
        "True",
        "1",
        "1.0",
        pd.NaT,
        np.datetime64("NaT", "ns"),
        np.timedelta64("NaT", "ns"),
    ]

    @pytest.mark.parametrize(
        "invalid", [*_invalid_scalars, 1, 1.0, np.int64(1), np.float64(1)]
    )
    def test_setitem_validation_scalar_bool(self, invalid):
        arr = pd.array([True, False, None], dtype="boolean")
        self._check_setitem_invalid(arr, invalid)

    @pytest.mark.parametrize("invalid", [*_invalid_scalars, True, 1.5, np.float64(1.5)])
    def test_setitem_validation_scalar_int(self, invalid, any_int_ea_dtype):
        arr = pd.array([1, 2, None], dtype=any_int_ea_dtype)
        self._check_setitem_invalid(arr, invalid)

    @pytest.mark.parametrize("invalid", [*_invalid_scalars, True])
    def test_setitem_validation_scalar_float(self, invalid, float_ea_dtype):
        arr = pd.array([1, 2, None], dtype=float_ea_dtype)
        self._check_setitem_invalid(arr, invalid)


@pytest.mark.parametrize(
    "dtype, value",
    [
        ("Int8", 1000),
        ("Int16", 40000),
        ("Int32", 2**31),
        ("UInt8", -1),
        ("UInt16", -1),
        ("UInt64", -1),
    ],
)
def test_setitem_listlike_out_of_bounds_raises(dtype, value):
    # GH#65510 an out-of-bounds list-like value must raise instead of
    # silently wrapping through a raw astype
    arr = pd.array([1, 2, 3], dtype=dtype)
    with pytest.raises(TypeError, match="cannot safely cast"):
        arr[[0]] = [value]
    # an NA alongside the out-of-bounds value still raises
    with pytest.raises(TypeError, match="cannot safely cast"):
        arr[[0, 1]] = [value, None]
    # the array is left unchanged
    assert arr.tolist() == [1, 2, 3]


def test_setitem_listlike_tuple():
    # GH#65510 a tuple value is assigned element-wise rather than raising an
    # AttributeError from isna returning a scalar for a tuple
    arr = pd.array([1, 2, 3], dtype="Int64")
    arr[[0, 1]] = (7, 8)
    assert arr.tolist() == [7, 8, 3]

    barr = pd.array([True, False, None], dtype="boolean")
    barr[[0]] = (False,)
    tm.assert_extension_array_equal(
        barr, pd.array([False, False, None], dtype="boolean")
    )


def test_setitem_listlike_in_bounds():
    # GH#65510 an in-bounds list-like value is still assigned
    arr = pd.array([1, 2, 3], dtype="Int8")
    arr[[0]] = [100]
    assert arr.tolist() == [100, 2, 3]


@pytest.mark.parametrize(
    "dtype, value_dtype, value",
    [
        ("Int8", "Int64", 1000),
        ("Int8", "Int16", 1000),
        ("UInt8", "Int64", -1),
        ("UInt64", "Int64", -1),
    ],
)
def test_setitem_masked_array_out_of_bounds_raises(dtype, value_dtype, value):
    # GH#65510 assigning a masked array whose values do not fit the target
    # dtype must raise instead of silently wrapping through a raw astype
    arr = pd.array([1, 2, 3], dtype=dtype)
    value = pd.array([value], dtype=value_dtype)
    with pytest.raises(TypeError, match="cannot safely cast"):
        arr[[0]] = value
    assert arr.tolist() == [1, 2, 3]


def test_setitem_masked_array_in_bounds():
    # GH#65510 an in-bounds masked array of a wider dtype is assigned, and an
    # NA in the value is preserved
    arr = pd.array([1, 2, 3], dtype="Int8")
    arr[[0, 1]] = pd.array([100, None], dtype="Int64")
    tm.assert_extension_array_equal(arr, pd.array([100, None, 3], dtype="Int8"))


@pytest.mark.parametrize(
    "dtype",
    [
        "Float64",
        pytest.param("float64[pyarrow]", marks=td.skip_if_no("pyarrow")),
    ],
)
@pytest.mark.parametrize("indexer", [1, [1], [False, True, False]])
def test_setitem_nan_in_float64_array(dtype, indexer, using_nan_is_na):
    arr = pd.array([0, pd.NA, 1], dtype=dtype)

    arr[indexer] = np.nan
    if not using_nan_is_na:
        assert np.isnan(arr[1])
    else:
        assert arr[1] is pd.NA


@pytest.mark.parametrize(
    "dtype",
    [
        "Int64",
        pytest.param("int64[pyarrow]", marks=td.skip_if_no("pyarrow")),
    ],
)
@pytest.mark.parametrize("indexer", [1, [1], [False, True, False]])
def test_setitem_nan_in_int64_array(dtype, indexer, using_nan_is_na):
    arr = pd.array([0, 1, 2], dtype=dtype)
    if not using_nan_is_na:
        err = TypeError
        msg = "Invalid value 'nan' for dtype 'Int64'"
        if dtype == "int64[pyarrow]":
            import pyarrow as pa

            err = pa.lib.ArrowInvalid
            msg = "Could not convert nan with type float"
        with pytest.raises(err, match=msg):
            arr[indexer] = np.nan
        assert arr[1] == 1
    else:
        arr[indexer] = np.nan
        assert arr[1] is pd.NA
