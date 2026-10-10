from __future__ import annotations

from decimal import Decimal

import numpy as np
import pytest

from pandas._libs import lib
from pandas.compat import pa_version_under18p0

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm
from pandas.api.types import (
    is_numeric_dtype,
)

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray
from pandas.tests.arrays.arrow.common import _require_timezone_database


def test_from_sequence_pa_array(data):
    # https://github.com/pandas-dev/pandas/pull/47034#discussion_r955500784
    # data._pa_array = pa.ChunkedArray
    result = type(data)._from_sequence(data._pa_array, dtype=data.dtype)
    tm.assert_extension_array_equal(result, data)
    assert isinstance(result._pa_array, pa.ChunkedArray)

    result = type(data)._from_sequence(
        data._pa_array.combine_chunks(), dtype=data.dtype
    )
    tm.assert_extension_array_equal(result, data)
    assert isinstance(result._pa_array, pa.ChunkedArray)


def test_from_sequence_pa_array_notimplemented():
    dtype = ArrowDtype(pa.month_day_nano_interval())
    with pytest.raises(NotImplementedError, match="Converting strings to"):
        ArrowExtensionArray._from_sequence_of_strings(["12-1"], dtype=dtype)


def test_from_sequence_of_strings_pa_array(data, request):
    pa_dtype = data.dtype.pyarrow_dtype
    if pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
        _require_timezone_database(request)

    pa_array = data._pa_array.cast(pa.string())
    result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
    tm.assert_extension_array_equal(result, data)

    pa_array = pa_array.combine_chunks()
    result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
    tm.assert_extension_array_equal(result, data)


def test_arrowdtype_construct_from_string_type_with_unsupported_parameters():
    with pytest.raises(NotImplementedError, match="Passing pyarrow type"):
        ArrowDtype.construct_from_string("not_a_real_dype[s, tz=UTC][pyarrow]")

    with pytest.raises(NotImplementedError, match="Passing pyarrow type"):
        ArrowDtype.construct_from_string("decimal(7, 2)[pyarrow]")


def test_arrowdtype_construct_from_string_supports_dt64tz():
    # as of GH#50689, timestamptz is supported
    dtype = ArrowDtype.construct_from_string("timestamp[s, tz=UTC][pyarrow]")
    expected = ArrowDtype(pa.timestamp("s", "UTC"))
    assert dtype == expected


def test_arrowdtype_construct_from_string_type_only_one_pyarrow():
    # GH#51225
    invalid = "int64[pyarrow]foobar[pyarrow]"
    msg = (
        r"Passing pyarrow type specific parameters \(\[pyarrow\]\) in the "
        r"string is not supported\."
    )
    with pytest.raises(NotImplementedError, match=msg):
        pd.Series(range(3), dtype=invalid)


def test_from_sequence_of_strings_boolean():
    true_strings = ["true", "TRUE", "True", "1", "1.0"]
    false_strings = ["false", "FALSE", "False", "0", "0.0"]
    nulls = [None]
    strings = true_strings + false_strings + nulls
    bools = (
        [True] * len(true_strings) + [False] * len(false_strings) + [None] * len(nulls)
    )

    dtype = ArrowDtype(pa.bool_())
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)
    expected = pd.array(bools, dtype="boolean[pyarrow]")
    tm.assert_extension_array_equal(result, expected)

    strings = ["True", "foo"]
    with pytest.raises(pa.ArrowInvalid, match="Failed to parse"):
        ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)


def test_from_sequence_of_strings_empty_string_float():
    # GH#66834 match numpy float64: empty string is not a valid float
    strings = ["1.5", "", "2.0"]
    dtype = ArrowDtype(pa.float64())
    with pytest.raises(ValueError, match=r"could not convert string to double"):
        ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)


def test_from_sequence_of_strings_empty_string_with_na():
    # GH#66834 pd.NA must not make the empty-string check raise TypeError
    strings = np.array(["1.5", pd.NA, ""], dtype=object)
    dtype = ArrowDtype(pa.float64())
    with pytest.raises(ValueError, match=r"could not convert string to double"):
        ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)


def test_from_sequence_of_strings_empty_string_int():
    # GH#66834
    strings = ["1", "", "2"]
    dtype = ArrowDtype(pa.int64())
    with pytest.raises(ValueError, match=r"could not convert string to int64"):
        ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)


@pytest.mark.parametrize("pa_type", [pa.int64(), pa.uint64()])
@pytest.mark.parametrize("box", [np.array, list, pd.Series, pd.Index])
def test_from_sequence_of_strings_int_precision_with_na(pa_type, box):
    # GH#56135 an NA must not route the integers through float64
    strings = box(np.array(["1582218195625938945", None], dtype=object))
    dtype = ArrowDtype(pa_type)
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)
    expected = ArrowExtensionArray(pa.array([1582218195625938945, None], type=pa_type))
    tm.assert_extension_array_equal(result, expected)


def test_from_sequence_of_strings_pa_array_rejects_hex():
    # GH#56135 pyarrow's own string cast reads "0x1F" as 31; to_numeric is what keeps
    #  pa.Array input as strict as list input
    strings = pa.array(["0x1F"], type=pa.string())
    dtype = ArrowDtype(pa.int64())
    with pytest.raises(ValueError, match="Unable to parse string"):
        ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)


def test_from_sequence_of_strings_int_above_uint64():
    # GH#56135 pins the unchanged fallback: too large for the nullable backend,
    #  so the default one still reports the value pyarrow cannot hold
    strings = np.array(["184467440737095516150", None], dtype=object)
    dtype = ArrowDtype(pa.int64())
    with pytest.raises(pa.ArrowInvalid, match="truncated converting to int64"):
        ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)


def test_from_sequence_of_strings_none_float():
    # GH#66834
    strings = ["1.5", None, "2.0"]
    dtype = ArrowDtype(pa.float64())
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)
    expected = ArrowExtensionArray(pa.array([1.5, None, 2.0], type=pa.float64()))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("box", [list, np.array, pd.Series, pa.array])
def test_from_sequence_of_strings_decimal(box):
    # GH#69838 float64 cannot hold every decimal, so parse without it
    strings = box(
        np.array(["1.23456789012345678901", None, " -2.5", "3 "], dtype=object)
    )
    dtype = ArrowDtype(pa.decimal128(30, 20))
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)
    expected = pd.array(
        [Decimal("1.23456789012345678901"), None, Decimal("-2.5"), Decimal("3")],
        dtype=dtype,
    )
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type",
    [
        pa.string(),
        pa.large_string(),
        pytest.param(
            pa.string_view(),
            marks=pytest.mark.xfail(
                pa_version_under18p0,
                reason="string_view cast to string added in pyarrow 18",
                raises=pa.ArrowNotImplementedError,
                strict=True,
            ),
        ),
        pa.binary(),
    ],
)
def test_from_sequence_of_strings_decimal_pa_string_types(pa_type):
    # GH#69838
    strings = pa.array([" 1.5", None], type=pa_type)
    dtype = ArrowDtype(pa.decimal128(10, 2))
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)
    expected = pd.array([Decimal("1.5"), None], dtype=dtype)
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("value", ["", "0x1F", "1.005"])
def test_from_sequence_of_strings_decimal_invalid(value):
    # GH#69838 the last one needs more digits than the scale allows
    dtype = ArrowDtype(pa.decimal128(10, 2))
    with pytest.raises(pa.ArrowInvalid, match="[Dd]ecimal"):
        ArrowExtensionArray._from_sequence_of_strings(["1.5", value], dtype=dtype)


@pytest.mark.parametrize("pa_type", [pa.string(), pa.large_string()])
@pytest.mark.parametrize("chunked", [True, False])
def test_from_sequence_of_strings_duration_sliced(chunked, pa_type):
    # GH#64320: the non-ns duration path used to route strings through pc.if_else,
    # which silently truncated values at a non-zero array offset
    values = ["11", "22", "33", "444444444", None]
    # seconds, not the nanoseconds to_timedelta would infer from a bare integer
    seconds = [11, 22, 33, 444444444, None]
    strings = pa.array(values, type=pa_type)
    if chunked:
        strings = pa.chunked_array([strings.slice(3), strings])
        expected_seconds = seconds[3:] + seconds
    else:
        strings = strings.slice(3)
        expected_seconds = seconds[3:]

    dtype = ArrowDtype(pa.duration("s"))
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=dtype)
    expected = ArrowExtensionArray(pa.array(expected_seconds, type=pa.duration("s")))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type", [pa.null(), pa.dictionary(pa.int32(), pa.string()), pa.int64()]
)
def test_from_sequence_of_strings_duration_non_varbinary(pa_type):
    # GH#64320: only the four affected layouts may route through
    # replace_with_mask. Diverting the rest is at best pointless and at worst
    # wrong: it aborts on the null type, and its numpy fallback drops nulls for
    # dictionary. int64 is unharmed, and is here to pin that it stays that way.
    if pa.types.is_null(pa_type):
        values, expected_seconds = [None, None, None], [None, None, None]
    elif pa.types.is_integer(pa_type):
        values, expected_seconds = [1, 2, None], [1, 2, None]
    else:
        values, expected_seconds = ["1", "2", None], [1, 2, None]
    strings = pa.chunked_array([pa.array(values, type=pa_type)])

    result = ArrowExtensionArray._from_sequence_of_strings(
        strings, dtype=ArrowDtype(pa.duration("s"))
    )

    expected = ArrowExtensionArray(pa.array(expected_seconds, type=pa.duration("s")))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("dtype", ["string", "string[pyarrow]"])
def test_series_from_string_array(dtype):
    arr = pa.array("the quick brown fox".split())
    ser = pd.Series(arr, dtype=dtype)
    expected = pd.Series(ArrowExtensionArray(arr), dtype=dtype)
    tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize(
    "data, arrow_dtype",
    [
        ([b"a", b"b"], pa.large_binary()),
        (["a", "b"], pa.large_string()),
    ],
)
def test_conversion_large_dtypes_from_numpy_array(data, arrow_dtype):
    dtype = ArrowDtype(arrow_dtype)
    result = pd.array(np.array(data), dtype=dtype)
    expected = pd.array(data, dtype=dtype)
    tm.assert_extension_array_equal(result, expected)


def test_astype_from_non_pyarrow(data):
    # GH49795
    np_arr = data.to_numpy()
    pd_array = pd.array(np_arr, dtype=np_arr.dtype)
    result = pd_array.astype(data.dtype)
    assert not isinstance(pd_array.dtype, ArrowDtype)
    assert isinstance(result.dtype, ArrowDtype)
    tm.assert_extension_array_equal(result, data)


def test_astype_float_from_non_pyarrow_str():
    # GH50430
    ser = pd.Series(["1.0"])
    result = ser.astype("float64[pyarrow]")
    expected = pd.Series([1.0], dtype="float64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_astype_errors_ignore():
    # GH 55399
    expected = pd.DataFrame({"col": [17000000]}, dtype="int32[pyarrow]")
    result = expected.astype("float[pyarrow]", errors="ignore")
    tm.assert_frame_equal(result, expected)


def test_to_numpy_with_defaults(data, using_nan_is_na):
    # GH49973
    result = data.to_numpy()

    pa_type = data._pa_array.type
    if pa.types.is_duration(pa_type) or pa.types.is_timestamp(pa_type):
        pytest.skip("Tested in test_to_numpy_temporal")
    elif pa.types.is_date(pa_type):
        expected = np.array(list(data))
    else:
        expected = np.array(data._pa_array)

    if data._hasna and (not is_numeric_dtype(data.dtype) or not using_nan_is_na):
        expected = expected.astype(object)
        expected[pd.isna(data)] = pd.NA

    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_int_with_na(using_nan_is_na):
    # GH51227: ensure to_numpy does not convert int to float
    data = [1, None]
    arr = pd.array(data, dtype="int64[pyarrow]")
    result = arr.to_numpy()
    if not using_nan_is_na:
        expected = np.array([1, pd.NA], dtype=object)
    else:
        expected = np.array([1, np.nan])
        assert isinstance(result[0], float)
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("na_val, exp", [(lib.no_default, np.nan), (1, 1)])
def test_to_numpy_null_array(na_val, exp):
    # GH#52443
    arr = pd.array([pd.NA, pd.NA], dtype="null[pyarrow]")
    result = arr.to_numpy(dtype="float64", na_value=na_val)
    expected = np.array([exp] * 2, dtype="float64")
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_null_array_no_dtype():
    # GH#52443
    arr = pd.array([pd.NA, pd.NA], dtype="null[pyarrow]")
    result = arr.to_numpy(dtype=None)
    expected = np.array([pd.NA] * 2, dtype="object")
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_without_dtype():
    # GH 54808
    arr = pd.array([True, pd.NA], dtype="boolean[pyarrow]")
    result = arr.to_numpy(na_value=False)
    expected = np.array([True, False], dtype=np.bool_)
    tm.assert_numpy_array_equal(result, expected)

    arr = pd.array([1.0, pd.NA], dtype="float32[pyarrow]")
    result = arr.to_numpy(na_value=0.0)
    expected = np.array([1.0, 0.0], dtype=np.float32)
    tm.assert_numpy_array_equal(result, expected)


def test_from_arrow_respecting_given_dtype():
    date_array = pa.array(
        [pd.Timestamp("2019-12-31"), pd.Timestamp("2019-12-31")], type=pa.date32()
    )
    result = date_array.to_pandas(
        types_mapper={pa.date32(): ArrowDtype(pa.date64())}.get
    )
    expected = pd.Series(
        [pd.Timestamp("2019-12-31"), pd.Timestamp("2019-12-31")],
        dtype=ArrowDtype(pa.date64()),
    )
    tm.assert_series_equal(result, expected)


def test_from_arrow_respecting_given_dtype_unsafe():
    array = pa.array([1.5, 2.5], type=pa.float64())
    with tm.external_error_raised(pa.ArrowInvalid):
        array.to_pandas(types_mapper={pa.float64(): ArrowDtype(pa.int64())}.get)


def test_from_arrow_list_of_extension_struct():
    # GH#69869 element access used to segfault after a same-type pyarrow cast
    intervals = pd.arrays.IntervalArray.from_tuples([(0, 1), (2, 3)])
    storage = pa.ListArray.from_arrays(
        pa.array([0, 2], pa.int32()), intervals.__arrow_array__()
    )
    table = pa.table({"x": storage})
    result = table.to_pandas(types_mapper=ArrowDtype)
    assert result["x"].dtype == ArrowDtype(storage.type)
    assert result["x"].iloc[0] == [{"left": 0, "right": 1}, {"left": 2, "right": 3}]


def test_from_arrow_renames_list_field():
    # GH#69869 an equal type with a different list field name is still cast
    arr = pa.array([[1]], type=pa.list_(pa.field("element", pa.int64())))
    dtype = ArrowDtype(pa.list_(pa.int64()))
    result = pa.table({"x": arr}).to_pandas(types_mapper=lambda _: dtype)["x"]
    assert str(result.dtype) == str(dtype)
    assert hash(result.dtype) == hash(dtype)


def test_astype_duration_from_sliced_arrow_strings():
    # GH#64320
    ser = pd.Series(["11", "22", "33", "444444444"], dtype="string[pyarrow]")[2:]

    result = ser.astype("duration[s][pyarrow]")

    expected = pd.Series(
        [pd.Timedelta(seconds=33), pd.Timedelta(seconds=444444444)],
        dtype="duration[s][pyarrow]",
        index=[2, 3],
    )
    tm.assert_series_equal(result, expected)


def test_null_astype_categorical():
    # GH#54908
    ser = pd.Series([None, None], dtype="null[pyarrow]")
    result = pd.Categorical(ser)
    dtype = pd.CategoricalDtype(categories=pd.Index([], dtype="null[pyarrow]"))
    expected = pd.Categorical([None, None], dtype=dtype)
    tm.assert_categorical_equal(result, expected)

    result_astype = ser.astype("category")
    expected_ser = pd.Series(expected)
    tm.assert_series_equal(result_astype, expected_ser)


def test_dictionary_astype_categorical():
    # GH#56672
    arrs = [
        pa.array(np.array(["a", "x", "c", "a"])).dictionary_encode(),
        pa.array(np.array(["a", "d", "c"])).dictionary_encode(),
    ]
    ser = pd.Series(ArrowExtensionArray(pa.chunked_array(arrs)))
    result = ser.astype("category")
    categories = pd.Index(["a", "x", "c", "d"], dtype=ArrowDtype(pa.string()))
    expected = pd.Series(
        ["a", "x", "c", "a", "a", "d", "c"],
        dtype=pd.CategoricalDtype(categories=categories),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("dtype", ["float64", "datetime64[ns]", "timedelta64[ns]"])
def test_astype_int_with_null_to_numpy_dtype(dtype):
    # GH 57093
    ser = pd.Series([1, None], dtype="int64[pyarrow]")
    result = ser.astype(dtype)
    expected = pd.Series([1, None], dtype=dtype)
    tm.assert_series_equal(result, expected)


def test_to_numpy_float():
    # GH#56267
    ser = pd.Series([32, 40, None], dtype="float[pyarrow]")
    result = ser.astype("float64")
    expected = pd.Series([32, 40, np.nan], dtype="float64")
    tm.assert_series_equal(result, expected)


def test_to_numpy_timestamp_to_int():
    # GH 55997
    ser = pd.Series(["2020-01-01 04:30:00"], dtype="timestamp[ns][pyarrow]")
    result = ser.to_numpy(dtype=np.int64)
    expected = np.array([1577853000000000000])
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("arrow_type", [pa.large_string(), pa.string()])
def test_cast_dictionary_different_value_dtype(arrow_type):
    df = pd.DataFrame({"a": ["x", "y"]}, dtype="string[pyarrow]")
    data_type = ArrowDtype(pa.dictionary(pa.int32(), arrow_type))
    result = df.astype({"a": data_type})
    assert result.dtypes.iloc[0] == data_type


def test_astype_dictionary_ordered():
    # GH#58152
    dtype = ArrowDtype(pa.dictionary(pa.int8(), pa.string(), ordered=True))
    ser = pd.Series(["foo", "bar", "foo"]).astype(dtype)
    assert ser.dtype == dtype
    assert ser.tolist() == ["foo", "bar", "foo"]


def test_constructor_dictionary_ordered():
    # GH#58152 pyarrow<25 drops ordered=True in pa.array
    # TODO(pyarrow>=25): remove, pyarrow then keeps the flag itself
    dtype = ArrowDtype(pa.dictionary(pa.int8(), pa.string(), ordered=True))
    ser = pd.Series(["foo", "bar", "foo"], dtype=dtype)
    assert ser.dtype == dtype
    assert ser.tolist() == ["foo", "bar", "foo"]


def test_categorical_from_arrow_dictionary():
    # GH 60563
    df = pd.DataFrame(
        {"A": ["a1", "a2"]}, dtype=ArrowDtype(pa.dictionary(pa.int32(), pa.utf8()))
    )
    result = df.value_counts(dropna=False)
    expected = pd.Series(
        [1, 1],
        index=pd.MultiIndex.from_arrays(
            [pd.Index(["a1", "a2"], dtype=ArrowDtype(pa.string()), name="A")]
        ),
        name="count",
        dtype="int64",
    )
    tm.assert_series_equal(result, expected)


def test_decimal_parse_raises():
    # GH 56984
    ser = pd.Series(["1.2345"], dtype=ArrowDtype(pa.string()))
    with pytest.raises(
        pa.lib.ArrowInvalid, match="Rescaling Decimal(128)? value would cause data loss"
    ):
        ser.astype(ArrowDtype(pa.decimal128(1, 0)))


def test_decimal_parse_succeeds():
    # GH 56984
    ser = pd.Series(["1.2345"], dtype=ArrowDtype(pa.string()))
    dtype = ArrowDtype(pa.decimal128(5, 4))
    result = ser.astype(dtype)
    expected = pd.Series([Decimal("1.2345")], dtype=dtype)
    tm.assert_series_equal(result, expected)
