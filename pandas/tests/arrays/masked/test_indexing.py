import re

import numpy as np
import pytest

from pandas.errors import Pandas4Warning
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


def test_setitem_masked_array_out_of_bounds_behind_mask():
    # GH#55232 a value behind the mask is meaningless, so an NA whose
    # underlying payload does not fit the target dtype is still assigned
    value = pd.array([300, 5], dtype="Int64")
    value[0] = pd.NA
    arr = pd.array([1, 2, 3], dtype="UInt8")
    arr[[0, 1]] = value
    tm.assert_extension_array_equal(arr, pd.array([None, 5, 3], dtype="UInt8"))


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


@pytest.mark.parametrize(
    "dtype, value",
    [
        ("Int8", 1000),
        ("Int16", 40000),
        ("Int32", 2**31),
        ("Int64", 2**64),
        ("UInt8", -1),
        ("UInt8", 256),
        ("UInt16", -1),
        ("UInt64", -1),
    ],
)
@pytest.mark.parametrize("box", [int, np.int64, np.float64])
def test_setitem_scalar_out_of_bounds_raises(dtype, value, box):
    # GH#48867 an out-of-bounds scalar must raise whichever box it arrives in;
    # numpy range-checks a python int but C-casts a numpy scalar for an
    # unsigned target and for the fancy keys below.
    if box is np.int64 and abs(value) >= 2**63:
        pytest.skip(f"{value} does not fit in int64")
    arr = pd.array([1, 2, 3], dtype=dtype)
    msg = re.escape(f"Invalid value '{box(value)!s}' for dtype '{dtype}'")

    for key in [0, slice(None), [0], np.array([True, False, False])]:
        with pytest.raises(TypeError, match=msg):
            arr[key] = box(value)

    ser = pd.Series(arr)
    with pytest.raises(TypeError, match=msg):
        ser[0] = box(value)
    with pytest.raises(TypeError, match=msg):
        ser.iloc[0] = box(value)

    assert arr.tolist() == [1, 2, 3]


@pytest.mark.parametrize("value", [np.int64(-1), -1, np.float64(-1.0)])
def test_fill_scalar_out_of_bounds_raises(value):
    # GH#48867 the other entry points that fill a scalar into a masked array
    msg = re.escape(f"Invalid value '{value!s}' for dtype 'UInt8'")
    ser = pd.Series([1, None, 3], dtype="UInt8")

    with pytest.raises(TypeError, match=msg):
        ser.fillna(value)
    with pytest.raises(TypeError, match=msg):
        ser.where(np.array([True, True, False]), value)
    with pytest.raises(TypeError, match=msg):
        ser.mask(np.array([False, False, True]), value)
    with pytest.raises(TypeError, match=msg):
        ser.array.insert(0, value)


@pytest.mark.parametrize(
    "dtype, value",
    [
        ("Int8", np.int64(-128)),
        ("Int8", np.int64(127)),
        ("UInt8", np.int64(0)),
        ("UInt8", np.int64(255)),
        ("UInt8", np.float64(255.0)),
        ("Int64", np.int64(-(2**63))),
        ("UInt64", np.uint64(2**64 - 1)),
        ("UInt64", 2**64 - 1),
        ("UInt64", np.float64(2**64 - 2048)),
    ],
)
def test_setitem_scalar_in_bounds(dtype, value):
    # GH#48867 the bounds check must not reject a value that does fit,
    # including the endpoints and the largest float below an unsigned maximum
    arr = pd.array([1, 2, 3], dtype=dtype)
    arr[0] = value
    assert arr[0] == value


@pytest.mark.parametrize(
    "dtype, value",
    [
        ("UInt64", np.float64(2**64)),
        ("UInt64", np.float32(2**64)),
        ("UInt32", np.float32(2**32)),
        ("Int64", np.float64(2**63)),
        ("Int32", np.float32(2**31)),
        ("Int16", np.float16(2**15)),
    ],
)
def test_setitem_float_just_out_of_bounds_raises(dtype, value):
    # GH#48867 the bounds are python ints, so comparing a float against them
    # demotes them and lets the first value above the maximum through
    arr = pd.array([1, 2, 3], dtype=dtype)
    msg = re.escape(f"Invalid value '{value!s}' for dtype '{dtype}'")

    with pytest.raises(TypeError, match=msg):
        arr[0] = value

    assert arr.tolist() == [1, 2, 3]


def test_index_fill_out_of_bounds_upcasts():
    # GH#48867 Index finds a common dtype where the array raises, so the
    # bounds check turns a silently wrapped 255 into the value asked for
    idx = pd.Index([1, 2, 3], dtype="UInt8")

    result = idx.putmask(np.array([True, False, False]), np.int64(-1))
    tm.assert_index_equal(result, pd.Index([-1, 2, 3], dtype="Int64"))

    result = idx.where(np.array([False, True, True]), np.int64(256))
    tm.assert_index_equal(result, pd.Index([256, 2, 3], dtype="Int64"))

    msg = "'int64' is not supported as a fill value for UInt8 dtype"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = pd.Index([1, None, 3], dtype="UInt8").fillna(np.int64(-1))
    tm.assert_index_equal(result, pd.Index([1, -1, 3], dtype="Int64"))
