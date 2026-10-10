from __future__ import annotations

from datetime import (
    date,
    time,
)
from decimal import Decimal
import re

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray


def test_setitem_null_slice(data):
    # GH50248
    orig = data.copy()

    result = orig.copy()
    result[:] = data[0]
    expected = ArrowExtensionArray._from_sequence(
        [data[0]] * len(data),
        dtype=data.dtype,
    )
    tm.assert_extension_array_equal(result, expected)

    result = orig.copy()
    result[:] = data[::-1]
    expected = data[::-1]
    tm.assert_extension_array_equal(result, expected)

    result = orig.copy()
    result[:] = data.tolist()
    expected = data
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "dtype, np_dtype",
    [
        ("int64[pyarrow]", "int64"),
        ("uint8[pyarrow]", "uint8"),
        ("double[pyarrow]", "float64"),
        ("timestamp[ns][pyarrow]", "M8[ns]"),
        ("duration[ns][pyarrow]", "m8[ns]"),
    ],
)
def test_setitem_null_slice_no_alias(dtype, np_dtype):
    # GH#67990 the null-slice fast path must not adopt a buffer the caller owns
    expected = pd.array(np.array([10, 20, 30]).astype(np_dtype), dtype=dtype)

    np_values = np.array([10, 20, 30]).astype(np_dtype)
    arr = pd.array([None] * 3, dtype=dtype)
    arr[:] = np_values
    np_values[0] = np_values[1]
    tm.assert_extension_array_equal(arr, expected)

    ser = pd.Series(np.array([10, 20, 30]).astype(np_dtype))
    arr = pd.array([None] * 3, dtype=dtype)
    arr[:] = ser
    ser.iloc[0] = ser.iloc[1]
    tm.assert_extension_array_equal(arr, expected)


def test_setitem_null_slice_no_alias_masked():
    # GH#67990 masked arrays are zero-copy through __arrow_array__
    arr = pd.array([None] * 3, dtype="int64[pyarrow]")
    values = pd.array([10, 20, 30], dtype="Int64")
    arr[:] = values
    values[0] = -1
    expected = pd.array([10, 20, 30], dtype="int64[pyarrow]")
    tm.assert_extension_array_equal(arr, expected)


@pytest.mark.parametrize(
    "wrap",
    [
        lambda np_values: pa.array(np_values),
        lambda np_values: pa.chunked_array([pa.array(np_values)]),
        lambda np_values: pa.chunked_array(
            [pa.array(np_values[:1]), pa.array(np_values[1:])]
        ),
        lambda np_values: ArrowExtensionArray(pa.array(np_values)),
        lambda np_values: pd.Series(
            ArrowExtensionArray(pa.array(np_values)), copy=False
        ),
    ],
    ids=["array", "chunked", "chunked_multi", "extension_array", "series"],
)
def test_setitem_null_slice_no_alias_pyarrow(wrap):
    # GH#67990 a pyarrow array is immutable, but its buffers can still be
    #  zero-copy over a numpy array the caller owns
    np_values = np.array([10, 20, 30], dtype="int64")
    arr = pd.array([None] * 3, dtype="int64[pyarrow]")
    value = wrap(np_values)
    arr[:] = value
    np_values[0] = -1
    expected = pd.array([10, 20, 30], dtype="int64[pyarrow]")
    tm.assert_extension_array_equal(arr, expected)
    # assert_extension_array_equal ignores chunking, so pin it separately
    assert arr._pa_array.num_chunks == getattr(value, "num_chunks", 1)


@pytest.mark.parametrize(
    "dtype",
    [
        pd.StringDtype("pyarrow", na_value=np.nan),
        "string[pyarrow]",
        "binary[pyarrow]",
        ArrowDtype(pa.large_string()),
        ArrowDtype(pa.large_binary()),
    ],
)
def test_setitem_null_slice_string_stays_zero_copy(dtype):
    # GH#67990 the aliasing fix must not undo GH#64529/GH#64530: pa.array never
    #  packs character data into caller-owned memory, so it is not copied here
    is_binary = "binary" in str(dtype)
    values = pd.array(
        [b"a", b"bb", b"ccc"] if is_binary else ["a", "bb", "ccc"], dtype=dtype
    )
    arr = pd.array([None] * 3, dtype=dtype)
    arr[:] = values
    # the character data is the last buffer for every one of these layouts
    assert (
        arr._pa_array.chunks[0].buffers()[-1].address
        == values._pa_array.chunks[0].buffers()[-1].address
    )


def test_setitem_null_slice_no_alias_dictionary():
    # GH#67990 pa.concat_arrays reuses the dictionary child, so the values half
    #  of a dictionary type needs copying too
    dtype = ArrowDtype(pa.dictionary(pa.int32(), pa.int64()))
    np_values = np.array([100, 200], dtype="int64")
    indices = pa.array(np.array([0, 1, 0], dtype="int32"))
    value = pa.DictionaryArray.from_arrays(indices, pa.array(np_values))

    arr = pd.array(value, dtype=dtype)
    arr[:] = value
    np_values[0] = -1
    expected = pd.array(
        pa.DictionaryArray.from_arrays(indices, pa.array([100, 200], type=pa.int64())),
        dtype=dtype,
    )
    tm.assert_extension_array_equal(arr, expected)


@pytest.mark.parametrize("kind", ["list", "struct"])
def test_setitem_null_slice_no_alias_nested(kind):
    # GH#67990 pa.concat_arrays does copy a nested type's children, which is why
    #  only dictionary needs the special case above
    np_values = np.array([1, 2, 3, 4], dtype="int64")
    child = pa.array(np_values)
    if kind == "list":
        value = pa.ListArray.from_arrays(pa.array([0, 2, 4], type=pa.int32()), child)
        expected_data = [[1, 2], [3, 4]]
    else:
        value = pa.StructArray.from_arrays([child], names=["x"])
        expected_data = [{"x": 1}, {"x": 2}, {"x": 3}, {"x": 4}]
    dtype = ArrowDtype(value.type)

    arr = pd.array(value, dtype=dtype)
    arr[:] = value
    np_values[0] = -1
    assert arr.tolist() == expected_data


def test_setitem_null_slice_cow():
    # GH#67990 full-slice assignment must not tie the two frames together
    df = pd.DataFrame({"a": pd.array([1, 2, 3], dtype="int64[pyarrow]")})
    other = pd.DataFrame({"b": [10, 20, 30]})
    df.loc[:, "a"] = other["b"]
    other.iloc[0, 0] = -777
    expected = pd.DataFrame({"a": pd.array([10, 20, 30], dtype="int64[pyarrow]")})
    tm.assert_frame_equal(df, expected)


def test_setitem_invalid_dtype(data):
    # GH50248
    pa_type = data._pa_array.type
    if pa.types.is_string(pa_type) or pa.types.is_binary(pa_type):
        fill_value = 123
        err = TypeError
        msg = "Invalid value '123' for dtype"
    elif (
        pa.types.is_integer(pa_type)
        or pa.types.is_floating(pa_type)
        or pa.types.is_boolean(pa_type)
    ):
        fill_value = "foo"
        err = pa.ArrowInvalid
        msg = "Could not convert"
    else:
        fill_value = "foo"
        err = TypeError
        msg = "Invalid value 'foo' for dtype"
    with pytest.raises(err, match=msg):
        data[:] = fill_value


@pytest.mark.parametrize(
    "target_tz, value_tz", [(None, "US/Eastern"), ("US/Eastern", None)]
)
@pytest.mark.parametrize("as_pydatetime", [False, True])
def test_setitem_timestamp_tz_mismatch_raises(target_tz, value_tz, as_pydatetime):
    # GH#69029 the value was stored as its UTC offset, which under a target of the
    #  other tz-awareness names a different instant. A plain datetime.datetime is
    #  reinterpreted the same way, so the guard is on datetime, not on Timestamp
    arr = pd.array([1, None], dtype=ArrowDtype(pa.timestamp("ns", tz=target_tz)))
    value = pd.Timestamp("2016-01-01", tz=value_tz)
    if as_pydatetime:
        value = value.to_pydatetime()
    msg = re.escape(f"Invalid value '{value}' for dtype '{arr.dtype}'")

    with pytest.raises(TypeError, match=msg):
        arr[0] = value
    with pytest.raises(TypeError, match=msg):
        arr[:1] = value
    with pytest.raises(TypeError, match=msg):
        pd.Series(arr).fillna(value)
    with pytest.raises(TypeError, match=msg):
        pd.Series(arr).where([False, True], value)
    with pytest.raises(TypeError, match=msg):
        pd.Series(arr).mask([True, False], value)
    with pytest.raises(TypeError, match=msg):
        # to_replace has to match, or nothing reaches the boxing step
        pd.Series(arr).replace(arr[0], value)
    with pytest.raises(TypeError, match=msg):
        # take() fills through _validate_setitem_value, so reindex is guarded too
        pd.Series(arr).reindex([0, 1, 2], fill_value=value)


@pytest.mark.parametrize(
    "target_tz, value_tz", [(None, "US/Eastern"), ("US/Eastern", None)]
)
def test_setitem_pa_scalar_tz_mismatch_raises(target_tz, value_tz):
    # GH#69029 a pa.Scalar skips the datetime guard and reaches the boundary at the
    #  trailing cast instead; a null still crosses, since it names no instant
    arr = pd.array([1, None], dtype=ArrowDtype(pa.timestamp("ns", tz=target_tz)))
    value = pa.scalar(pd.Timestamp("2016-01-01", tz=value_tz))

    msg = re.escape(f"Invalid value '{value}' for dtype '{arr.dtype}'")
    with pytest.raises(TypeError, match=msg):
        arr[0] = value

    arr[0] = pa.scalar(None, type=pa.timestamp("ns", tz=value_tz))
    assert arr[0] is pd.NA


def test_setitem_timestamp_other_tz_converts():
    # GH#69029 two tz-aware operands name the same instant, so this is not a mismatch
    arr = pd.array([1, None], dtype=ArrowDtype(pa.timestamp("ns", tz="US/Eastern")))
    arr[0] = pd.Timestamp("2016-01-01 00:00", tz="UTC")
    assert arr[0] == pd.Timestamp("2015-12-31 19:00", tz="US/Eastern")


@pytest.mark.parametrize(
    "value",
    [
        pd.date_range("2016-01-01", periods=3, tz="UTC")._data,
        pd.date_range("2016-01-01", periods=3)._data,
        pd.timedelta_range("1D", periods=3)._data,
        pd.period_range("2016-01-01", periods=3, freq="D")._data,
        np.array(["2016-01-01"] * 3, dtype="M8[ns]"),
        np.array([1, 2, 3], dtype="m8[s]"),
        pd.array([date(2016, 1, 1)] * 3, dtype=ArrowDtype(pa.date32())),
        pd.array([time(1, 2)] * 3, dtype=ArrowDtype(pa.time64("us"))),
        pa.array([1, 2, 3], type=pa.timestamp("us")),
        pd.Timestamp("2016-01-01"),
        pd.Timedelta("1D"),
        date(2016, 1, 1),
        time(1, 2),
    ],
)
def test_setitem_temporal_into_numeric_raises(value):
    # GH#68419 pyarrow would convert these into the integer storage instead of
    #  raising the way every other dtype does
    arr = pd.array([1, 2, 3], dtype="int64[pyarrow]")
    with pytest.raises(TypeError, match="Invalid value"):
        arr[:] = value

    ser = pd.Series([1, 2, 3], dtype="int64[pyarrow]")
    with pytest.raises(TypeError, match="Invalid value"):
        ser.iloc[:] = value


@pytest.mark.parametrize("dtype", ["timestamp[ns][pyarrow]", "duration[ns][pyarrow]"])
@pytest.mark.parametrize("value", [np.array([1, 2, 3]), np.array([1.0, 2.0, 3.0]), 1])
def test_setitem_numeric_into_temporal_raises(dtype, value):
    # GH#68419 mirror of test_setitem_temporal_into_numeric_raises
    arr = pd.array([1, 2, 3], dtype=dtype)
    with pytest.raises(TypeError, match="Invalid value"):
        arr[:] = value


@pytest.mark.parametrize("dtype", ["timestamp[ns][pyarrow]", "duration[ns][pyarrow]"])
def test_setitem_oversized_int_scalar_into_temporal_raises(dtype):
    # GH#68419 an int too wide for any integer dtype infers to object, which would
    #  leave it unsettled; numpy M8/m8 raise this same TypeError for it
    arr = pd.array([1, 2, 3], dtype=dtype)
    with pytest.raises(TypeError, match="Invalid value"):
        arr[0] = 2**70


def test_setitem_decimal_scalar_into_temporal_raises():
    # GH#68419 infer_dtype_from_scalar maps a Decimal to object, which would leave
    #  the scalar reinterpreting while the decimal128 array form raises
    arr = pd.array([1, 2, 3], dtype="timestamp[ns][pyarrow]")
    with pytest.raises(TypeError, match="Invalid value"):
        arr[0] = Decimal(1)


def test_setitem_decimal_scalar_into_decimal_self_accepted():
    # GH#68419 counterpart of test_setitem_decimal_scalar_into_temporal_raises
    dtype = ArrowDtype(pa.decimal128(10, 2))
    arr = pd.array([Decimal("1.00")] * 2, dtype=dtype)
    arr[0] = Decimal("2.00")
    tm.assert_extension_array_equal(
        arr, pd.array([Decimal("2.00"), Decimal("1.00")], dtype=dtype)
    )


@pytest.mark.parametrize(
    "value",
    [
        pd.DatetimeIndex(["2017-01-01"] * 3)._data,
        np.array(["2017-01-01"] * 3, dtype="M8[ns]"),
        np.array(["2017-01-01"] * 3, dtype="U10"),
        np.array(["2017-01-01"] * 3, dtype=object),
        pd.arrays.SparseArray(np.array(["2017-01-01"] * 3, dtype="M8[ns]")),
        pd.Categorical(pd.DatetimeIndex(["2017-01-01"] * 3)),
        pd.array(["2017-01-01"] * 3, dtype="string"),
        pd.array(["2017-01-01"] * 3, dtype=ArrowDtype(pa.string())),
        pd.array(["2017-01-01"] * 3, dtype=ArrowDtype(pa.large_string())),
        pd.array(
            [pd.Timestamp("2017-01-01")] * 3, dtype="timestamp[ns][pyarrow]"
        ).astype(ArrowDtype(pa.dictionary(pa.int32(), pa.timestamp("ns")))),
        pd.Timestamp("2017-01-01"),
    ],
)
def test_setitem_temporal_still_accepted(value):
    # GH#68419 the temporal check must not reject what already worked; a dtype
    #  it cannot classify has to fall through rather than count as non-temporal
    arr = pd.array([pd.Timestamp("2016-01-01")] * 3, dtype="timestamp[ns][pyarrow]")
    arr[:] = value
    expected = pd.array(
        [pd.Timestamp("2017-01-01")] * 3, dtype="timestamp[ns][pyarrow]"
    )
    tm.assert_extension_array_equal(arr, expected)


@pytest.mark.parametrize(
    "pa_type", [pa.dictionary(pa.int32(), pa.timestamp("ns")), pa.string()]
)
def test_setitem_temporal_into_unsettled_self_accepted(pa_type):
    # GH#68419 the self side is three-state too: a dictionary self is temporal
    #  via its value_type, and a string self settles nothing
    dtype = ArrowDtype(pa_type)
    arr = pd.array([pd.Timestamp("2016-01-01")] * 3, dtype="timestamp[ns][pyarrow]")
    arr = arr.astype(dtype)
    value = pd.DatetimeIndex(["2017-01-01"] * 3)._data
    if pa.types.is_dictionary(pa_type):
        # pyarrow cannot build a dictionary array from a numpy-backed value
        value = pd.array(value, dtype="timestamp[ns][pyarrow]").astype(dtype)
    arr[:] = value
    assert arr.dtype == dtype
    assert not arr.isna().any()
    expected = pd.Timestamp("2017-01-01")
    assert arr.astype("timestamp[ns][pyarrow]")[0] == expected


@pytest.mark.parametrize(
    "pa_type, value",
    [
        (pa.string(), pd.Timestamp("2016-01-01")),
        (pa.large_string(), pd.Timestamp("2016-01-01")),
        (pa.binary(), pd.Timestamp("2016-01-01")),
        (pa.duration("ns"), pd.Timestamp("2016-01-01")),
        (pa.timestamp("ns"), pd.Timedelta("1s")),
        (pa.time64("us"), pd.Timedelta("1s")),
        (pa.time32("s"), pd.Timedelta("1s")),
    ],
)
def test_setitem_temporal_scalar_into_mismatched_self_raises(pa_type, value):
    # GH#68419 _box_pa_scalar read pa_type.unit whatever the target was, so a
    #  string self raised AttributeError and a mismatched temporal self silently
    #  stored the integer
    arr = pd.array([None, None], dtype=ArrowDtype(pa_type))
    with pytest.raises(TypeError, match="Invalid value"):
        arr[0] = value


@pytest.mark.parametrize("pa_type", [pa.date32(), pa.date64()])
@pytest.mark.parametrize(
    "value", [pd.Timestamp("2016-01-05"), pd.Timestamp("2016-01-05 12:30:45")]
)
def test_setitem_timestamp_into_date_self(pa_type, value):
    # GH#68419 a Timestamp is a valid date value; reaching for date32's
    #  nonexistent .unit used to make this an AttributeError. A time component is
    #  dropped, matching pd.array([value], dtype=ArrowDtype(pa_type))
    arr = pd.array([date(2016, 1, 1)] * 2, dtype=ArrowDtype(pa_type))
    arr[0] = value
    expected = pd.array([date(2016, 1, 5), date(2016, 1, 1)], dtype=ArrowDtype(pa_type))
    tm.assert_extension_array_equal(arr, expected)


def test_setitem_numeric_still_accepted():
    # GH#68419 counterpart of test_setitem_temporal_still_accepted
    arr = pd.array([1, 2, 3], dtype="int64[pyarrow]")
    arr[:] = np.array([4, 5, 6])
    tm.assert_extension_array_equal(arr, pd.array([4, 5, 6], dtype="int64[pyarrow]"))


@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "timestamp[ns][pyarrow]"])
@pytest.mark.parametrize(
    "value",
    [
        None,
        pd.NA,
        pd.NaT,
        np.nan,
        np.datetime64("NaT", "ns"),
        np.timedelta64("NaT", "ns"),
        [None, None, None],
        pa.array([None] * 3),
        pd.array([None] * 3, dtype=ArrowDtype(pa.null())),
    ],
)
def test_setitem_na_not_rejected(dtype, value):
    # GH#68419 the temporal check must not reject NA of any flavor; a null-typed or
    #  NaT-scalar value carries a dtype but still settles nothing. Int64 rejects NaT
    #  here; that divergence is pre-existing
    arr = pd.array([1, 2, 3], dtype=dtype)
    arr[:] = value
    assert arr.isna().all()


@pytest.mark.parametrize(
    "pa_type", [pa.duration("ns"), pa.timestamp("ns"), pa.time64("us"), pa.int64()]
)
def test_setitem_typed_null_scalar_not_rejected(pa_type):
    # GH#68419 a typed null pa.Scalar is just "assign NA"; isna() does not
    #  recognize it, so the check has to look at is_valid
    arr = pd.array([1, 2, 3], dtype="int64[pyarrow]")
    arr[:] = pa.scalar(None, type=pa_type)
    assert arr.isna().all()


@pytest.mark.parametrize(
    "value",
    [
        pd.array([None] * 3, dtype="timestamp[ns][pyarrow]"),
        pa.array([None] * 3, type=pa.timestamp("ns")),
        pd.DatetimeIndex([pd.NaT] * 3)._data,
    ],
)
def test_setitem_all_na_temporal_array_still_raises(value):
    # GH#68419 an all-NA array still carries a temporal dtype, and numpy int64
    #  and Int64 both reject it; only a *scalar* NA means "just assign NA"
    arr = pd.array([1, 2, 3], dtype="int64[pyarrow]")
    with pytest.raises(TypeError, match="Invalid value"):
        arr[:] = value


@pytest.mark.parametrize("dtype", ["timestamp[ns][pyarrow]", "duration[ns][pyarrow]"])
@pytest.mark.parametrize(
    "value",
    [
        np.array([np.nan] * 3),
        np.array([np.nan] * 3, dtype="float32"),
        pd.array([None] * 3, dtype="Float64"),
        pd.array([None] * 3, dtype="double[pyarrow]"),
    ],
)
def test_setitem_all_na_float_array_not_rejected(dtype, value):
    # GH#68419 an all-NaN float array carries no values to reinterpret. Rejecting it
    #  would break where/fillna against an alignment-produced NaN column, which numpy
    #  M8 upcasts to object; DatetimeArray setitem is stricter and is not followed here.
    #  The masked and pyarrow spellings work only because the escape returns typed nulls
    arr = pd.array([1, 2, 3], dtype=dtype)
    arr[:] = value
    assert arr.isna().all()


@pytest.mark.parametrize("dtype", ["timestamp[ns][pyarrow]", "duration[ns][pyarrow]"])
def test_setitem_empty_float_array_not_rejected(dtype):
    # GH#68419 an empty value is vacuously all-NA
    arr = pd.array([1, 2, 3], dtype=dtype)
    arr[[]] = np.array([], dtype=float)
    assert not arr.isna().any()


def test_setitem_2d_all_na_float_array_still_raises():
    # GH#68419 the all-NA escape sizes its nulls with len(), which reads only axis 0,
    #  so a 2-D value has to stay with _box_pa rather than be flattened
    arr = pd.array([1, 2, 3], dtype="timestamp[ns][pyarrow]")
    with pytest.raises(ValueError, match="Mask must be 1D"):
        arr[:] = np.full((3, 1), np.nan)


def test_setitem_partial_na_float_array_still_raises():
    # GH#68419 the all-NA escape must not widen to "contains NA": the non-NA
    #  entry of [nan, 1.0, nan] is reinterpreted as 1ns past the epoch
    arr = pd.array([1, 2, 3], dtype="timestamp[ns][pyarrow]")
    with pytest.raises(TypeError, match="Invalid value"):
        arr[:] = np.array([np.nan, 1.0, np.nan])


def test_setitem_boolean_replace_with_mask_segfault():
    # GH#52059
    N = 145_000
    arr = ArrowExtensionArray(pa.chunked_array([np.ones((N,), dtype=np.bool_)]))
    expected = arr.copy()
    arr[np.zeros((N,), dtype=np.bool_)] = False
    assert arr._pa_array == expected._pa_array


def test_setitem_null_dtype_replace_with_mask_abort():
    # GH#66703 pc.replace_with_mask aborts for null dtype (apache/arrow#47447).
    # Operations routed through _replace_with_mask must not crash the process.
    ser = pd.Series(
        pa.array([None, None, None], type=pa.null()), dtype=ArrowDtype(pa.null())
    )
    mask = np.array([True, False, False])

    # boolean-mask setitem
    result = ser.copy()
    result[mask] = None
    tm.assert_series_equal(result, ser)

    # where / mask
    tm.assert_series_equal(ser.where(~mask), ser)
    tm.assert_series_equal(ser.mask(mask), ser)

    # combine_first
    tm.assert_series_equal(ser.combine_first(ser), ser)

    # DataFrame.loc setitem
    df = ser.to_frame("a")
    df.loc[mask, "a"] = None
    tm.assert_frame_equal(df, ser.to_frame("a"))


def test_setitem_na_chunked_string_if_else():
    # GH#64320
    df = pd.concat(
        [
            pd.DataFrame({"a": ["x"] * 5, "b": ["x"] * 5}),
            pd.DataFrame({"a": ["x"] * 5, "b": ["x"] * 5}),
        ],
        ignore_index=True,
    )
    for _ in range(5):
        df.loc[[0], "a"] = pd.NA
    assert pd.isna(df["a"].iloc[0])
    assert (df["a"].iloc[1:] == "x").all()
    assert (df["b"] == "x").all()


@pytest.mark.parametrize(
    "pa_type", [pa.binary(), pa.large_binary(), pa.string(), pa.large_string()]
)
@pytest.mark.parametrize("extra_chunk", [True, False])
def test_setitem_na_sliced_chunk_if_else(pa_type, extra_chunk):
    # GH#64320
    values = ["a", "bb", "ccc", "dddd", "eeeee"]
    if pa.types.is_binary(pa_type) or pa.types.is_large_binary(pa_type):
        values = [val.encode() for val in values]
    chunk = pa.array(values, type=pa_type)
    # the first chunk carries a non-zero offset, which pc.if_else mishandles
    chunks = [chunk.slice(3), chunk] if extra_chunk else [chunk.slice(3)]
    arr = ArrowExtensionArray(pa.chunked_array(chunks))
    expected_values = [None, values[4]] + (values if extra_chunk else [])
    expected = ArrowExtensionArray(pa.array(expected_values, type=pa_type))

    arr[[0]] = None

    arr._pa_array.validate(full=True)
    tm.assert_extension_array_equal(arr, expected)


def test_setitem_float_nan_is_na(using_nan_is_na):
    # GH#61732
    ser = pd.Series([-1, 0, 1], dtype="int64[pyarrow]")

    if using_nan_is_na:
        ser[1] = np.nan
        assert ser.isna()[1]
    else:
        msg = "Could not convert nan with type float: tried to convert to int64"
        with pytest.raises(pa.lib.ArrowInvalid, match=msg):
            ser[1] = np.nan

    ser = pd.Series([-1, np.nan, 1], dtype="float64[pyarrow]")
    if using_nan_is_na:
        assert ser.isna()[1]
        assert ser[1] is pd.NA

        ser[1] = np.nan
        assert ser[1] is pd.NA

    else:
        assert not ser.isna()[1]
        assert isinstance(ser[1], float)
        assert np.isnan(ser[1])

        ser[2] = np.nan
        assert isinstance(ser[2], float)
        assert np.isnan(ser[2])
