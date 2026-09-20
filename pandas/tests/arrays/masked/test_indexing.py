import re
import warnings

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


def _out_of_bounds_msg(value, dtype):
    return re.escape(
        f"Setting the out-of-bounds value {value!s} into an array of dtype "
        f"{dtype} is deprecated"
    )


@pytest.mark.parametrize(
    "dtype, value, wrapped",
    [
        ("UInt8", np.int64(-1), 255),
        ("UInt8", np.int8(-1), 255),
        ("UInt8", np.int64(256), 0),
        ("UInt16", np.int64(-1), 65535),
        ("UInt32", np.int64(-1), 4294967295),
        ("UInt64", np.int64(-1), 2**64 - 1),
    ],
)
@pytest.mark.parametrize("key", [0, slice(None), [0], np.array([True, False, False])])
def test_setitem_unsigned_out_of_bounds_deprecated(dtype, value, wrapped, key):
    # GH#48867 numpy casts a numpy scalar into an unsigned dtype unchecked,
    #  whatever the key, so the wrong value is stored silently
    arr = pd.array([1, 2, 3], dtype=dtype)
    with tm.assert_produces_warning(
        Pandas4Warning, match=_out_of_bounds_msg(value, dtype)
    ):
        arr[key] = value
    assert arr[0] == wrapped


@pytest.mark.parametrize("box", [int, float])
@pytest.mark.parametrize("value", [-1, 256])
def test_setitem_python_scalar_out_of_bounds_still_raises(box, value):
    # GH#48867 numpy range-checks a python scalar, so nothing is deprecated here
    arr = pd.array([1, 2, 3], dtype="UInt8")
    with tm.assert_produces_warning(None):
        with pytest.raises(OverflowError, match="out of bounds"):
            arr[0] = box(value)
    assert arr.tolist() == [1, 2, 3]


@pytest.mark.parametrize(
    "dtype, value",
    [
        ("Int8", np.int64(1000)),
        ("Int16", np.int64(40000)),
        ("Int32", np.int64(2**31)),
        ("Int64", np.float64(2**63)),
    ],
)
@pytest.mark.parametrize("key", [0, slice(0, 1)])
def test_setitem_signed_scalar_key_still_raises(dtype, value, key):
    # GH#48867 for a signed dtype numpy range-checks an integer or slice key, so
    #  nothing is stored and the deprecation must not fire ahead of the raise
    arr = pd.array([1, 2, 3], dtype=dtype)
    with tm.assert_produces_warning(None):
        # "too large"/"too big": CPython's float-to-int message, spelled
        #  differently on Windows and 32-bit
        msg = "|".join(["out of bounds", "too large", "too big"])
        with pytest.raises(OverflowError, match=msg):
            arr[key] = value
    assert arr.tolist() == [1, 2, 3]


@pytest.mark.parametrize(
    "dtype, value, wrapped",
    [
        ("Int8", np.int64(128), -128),
        ("Int8", np.int64(1000), -24),
        ("Int16", np.int64(40000), -25536),
        ("Int64", np.uint64(2**63), -(2**63)),
    ],
)
@pytest.mark.parametrize("key", [[0], np.array([0]), np.array([True, False, False])])
def test_setitem_signed_fancy_key_out_of_bounds_deprecated(dtype, value, wrapped, key):
    # GH#48867 for a signed dtype numpy range-checks only an integer or slice key;
    #  a fancy or mask key is unchecked, so signed wraps too -- and
    #  where/mask/replace/putmask all reach __setitem__ with a mask key
    arr = pd.array([1, 2, 3], dtype=dtype)
    with tm.assert_produces_warning(
        Pandas4Warning, match=_out_of_bounds_msg(value, dtype)
    ):
        arr[key] = value
    assert arr[0] == wrapped


@pytest.mark.parametrize(
    "dtype, value",
    [
        ("UInt8", np.int64(0)),
        ("UInt8", np.int64(255)),
        ("UInt8", np.float64(255.0)),
        ("UInt8", 255),
        ("UInt64", np.uint64(2**64 - 1)),
        ("UInt64", np.float64(2**64 - 2048)),
        ("Int8", np.int64(-128)),
        ("Int8", np.int64(127)),
    ],
)
def test_setitem_in_bounds_not_deprecated(dtype, value):
    # GH#48867 the endpoints and the largest float below an unsigned maximum fit
    arr = pd.array([1, 2, 3], dtype=dtype)
    with tm.assert_produces_warning(None):
        arr[0] = value
    assert arr[0] == value


@pytest.mark.parametrize(
    "dtype, value",
    [
        ("UInt64", np.float64(2**64)),
        ("UInt64", np.float32(2**64)),
        ("UInt32", np.float32(2**32)),
    ],
)
def test_setitem_float_just_out_of_bounds_deprecated(dtype, value):
    # GH#48867 the bounds are python ints, so comparing a float against them
    #  demotes them and would let the first value above the maximum through.
    #  What numpy stores is platform-dependent, so only the warning is checked.
    arr = pd.array([1, 2, 3], dtype=dtype)
    with tm.assert_produces_warning(
        Pandas4Warning, match=_out_of_bounds_msg(value, dtype)
    ):
        with np.errstate(invalid="ignore"):
            arr[0] = value


@pytest.mark.parametrize(
    "dtype, value, wrapped",
    [("UInt8", np.int64(-1), 255), ("Int8", np.int64(1000), -24)],
)
def test_series_fill_paths_out_of_bounds_deprecated(dtype, value, wrapped):
    # GH#48867 the Series entry points that fill a scalar through __setitem__
    msg = _out_of_bounds_msg(value, dtype)
    ser = pd.Series([1, 2, 3], dtype=dtype)

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.where(np.array([False, True, True]), value)
    tm.assert_series_equal(result, pd.Series([wrapped, 2, 3], dtype=dtype))

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.mask(np.array([True, False, False]), value)
    tm.assert_series_equal(result, pd.Series([wrapped, 2, 3], dtype=dtype))

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.replace(1, value)
    tm.assert_series_equal(result, pd.Series([wrapped, 2, 3], dtype=dtype))


@pytest.mark.parametrize(
    "dtype, value, wrapped",
    [("UInt8", np.int64(-1), 255), ("Int8", np.int64(1000), -24)],
)
def test_index_fill_paths_out_of_bounds_deprecated(dtype, value, wrapped):
    # GH#48867 Index.fillna pre-validates the fill value and then routes through
    #  putmask; the count below is what pins that to one warning
    msg = _out_of_bounds_msg(value, dtype)
    idx = pd.Index([1, 2, 3], dtype=dtype)

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = idx.putmask(np.array([True, False, False]), value)
    tm.assert_index_equal(result, pd.Index([wrapped, 2, 3], dtype=dtype))

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = idx.where(np.array([False, True, True]), value)
    tm.assert_index_equal(result, pd.Index([wrapped, 2, 3], dtype=dtype))

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        result = pd.Index([1, None, 3], dtype=dtype).fillna(value)
    tm.assert_index_equal(result, pd.Index([1, wrapped, 3], dtype=dtype))
    assert sum(issubclass(w.category, Pandas4Warning) for w in record) == 1


def test_paths_that_already_raise_are_not_deprecated():
    # GH#48867 Series.fillna and ExtensionArray.insert reject the value before it
    #  reaches __setitem__, so they must not gain a warning before the raise
    value = np.int64(-1)
    with tm.assert_produces_warning(None):
        with pytest.raises(OverflowError, match="can't convert negative value"):
            pd.Series([1, None, 3], dtype="UInt8").fillna(value)

    with tm.assert_produces_warning(None):
        with pytest.raises(TypeError, match="cannot safely cast non-equivalent"):
            pd.array([1, 2, 3], dtype="UInt8").insert(0, value)


@pytest.mark.parametrize("how", ["iloc", "loc", "at", "where", "mask"])
def test_dataframe_out_of_bounds_deprecated(how):
    # GH#48867 DataFrame setitem/where/mask reach the same __setitem__
    value = np.int64(-1)
    df = pd.DataFrame({"a": pd.array([1, 2, 3], dtype="UInt8")})
    with tm.assert_produces_warning(
        Pandas4Warning, match=_out_of_bounds_msg(value, "UInt8")
    ):
        if how == "iloc":
            df.iloc[0, 0] = value
        elif how == "loc":
            df.loc[0, "a"] = value
        elif how == "at":
            df.at[0, "a"] = value
        elif how == "where":
            df = df.where(pd.DataFrame({"a": [False, True, True]}), value)
        else:
            df = df.mask(pd.DataFrame({"a": [True, False, False]}), value)
    tm.assert_series_equal(df["a"], pd.Series([255, 2, 3], dtype="UInt8", name="a"))


@pytest.mark.parametrize(
    "key",
    [[], slice(0, 0), np.array([], dtype=np.intp), np.array([False, False, False])],
)
def test_setitem_empty_key_still_deprecated(key):
    # GH#48867 validation ignores the key, so this fires even though nothing is
    #  stored; an invalid value already raises for these keys
    arr = pd.array([1, 2, 3], dtype="UInt8")
    with tm.assert_produces_warning(
        Pandas4Warning, match=_out_of_bounds_msg(np.int64(-1), "UInt8")
    ):
        arr[key] = np.int64(-1)
    assert arr.tolist() == [1, 2, 3]


@pytest.mark.parametrize(
    "dtype, data, value",
    [
        ("Float64", [1.0, 2.0], np.float64(1e300)),
        ("Float32", [1.0, 2.0], np.float32(1.5)),
        ("boolean", [True, False], np.bool_(True)),
    ],
)
def test_setitem_non_integer_dtype_not_deprecated(dtype, data, value):
    # GH#48867 np.iinfo rejects a non-integer dtype, so the kind guard in
    #  _warn_if_out_of_bounds is load bearing
    arr = pd.array(data, dtype=dtype)
    with tm.assert_produces_warning(None):
        arr[0] = value
    assert arr[0] == value
