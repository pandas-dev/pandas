from __future__ import annotations

from datetime import (
    datetime,
    time,
    timedelta,
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


def test_fillna_temporal_into_string_self_accepted():
    # GH#68419 filling a string column with datetimes is a string conversion,
    #  not an integer reinterpretation, and must keep working
    ser = pd.Series(["a", None], dtype=ArrowDtype(pa.string()))
    result = ser.fillna(pd.Series(pd.date_range("2016-01-01", periods=2)))
    assert result.dtype == ArrowDtype(pa.string())
    assert result[0] == "a"
    assert result[1] == "2016-01-02 00:00:00.000000"


def test_fillna_temporal_scalar_into_string_self_raises():
    # GH#68419 the scalar spelling used to raise AttributeError. It stays stricter
    #  than the array spelling above, which pyarrow converts
    ser = pd.Series(["a", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(TypeError, match="Invalid value"):
        ser.fillna(pd.Timestamp("2016-01-01"))


@pytest.mark.parametrize("other_dtype", ["float64", "Float64", "double[pyarrow]"])
def test_where_fillna_all_na_float_other_not_rejected(other_dtype):
    # GH#68419 alignment routinely produces an all-NaN float column, and both
    #  reach _validate_setitem_value
    ser = pd.Series(
        pd.date_range("2016-01-01", periods=3), dtype="timestamp[ns][pyarrow]"
    )
    # None, not np.nan: a masked or Arrow float built from np.nan holds NaN rather
    #  than NA once future.distinguish_nan_and_na is on, which _is_all_na rejects
    result = ser.where(
        np.array([True, False, True]), pd.Series([None] * 3, dtype=other_dtype)
    )
    assert result.dtype == "timestamp[ns][pyarrow]"
    assert result.isna().tolist() == [False, True, False]

    ser = pd.Series([pd.Timestamp("2016-01-01"), None], dtype="timestamp[ns][pyarrow]")
    result = ser.fillna(pd.Series([None] * 2, dtype=other_dtype))
    tm.assert_series_equal(result, ser)


def test_fillna_temporal_reinterpretation_raises():
    # GH#68419 fillna shares _validate_setitem_value, so the same
    #  reinterpretation is rejected there
    ser = pd.Series([1, None], dtype="int64[pyarrow]")
    with pytest.raises(TypeError, match="Invalid value"):
        ser.fillna(pd.Series(pd.date_range("2016-01-01", periods=2)))

    ser = pd.Series([pd.Timestamp("2016-01-01"), None], dtype="timestamp[ns][pyarrow]")
    with pytest.raises(TypeError, match="Invalid value"):
        ser.fillna(1)

    result = ser.fillna(pd.Timestamp("2017-01-01"))
    expected = pd.Series(
        [pd.Timestamp("2016-01-01"), pd.Timestamp("2017-01-01")],
        dtype="timestamp[ns][pyarrow]",
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("meth", ["where", "mask"])
def test_where_mask_temporal_reinterpretation_raises(meth):
    # GH#68419 where/mask reach the guard through setitem
    ser = pd.Series([1, 2, 3], dtype="int64[pyarrow]")
    other = pd.Series(pd.date_range("2016-01-01", periods=3))
    with pytest.raises(TypeError, match="Invalid value"):
        getattr(ser, meth)(np.array([True, False, True]), other)


@pytest.mark.parametrize("limit", [None, 1])
def test_fillna_string_self_agrees_with_limit_path(limit):
    # GH#68419 fillna shares _validate_setitem_value with __setitem__, so the
    #  limit=None and limit=1 paths agree
    arr = pd.array(["a", None], dtype=pd.StringDtype("pyarrow", na_value=np.nan))
    with pytest.raises(TypeError, match="Invalid value for dtype"):
        arr.fillna(np.array([1, 2]), limit=limit)


def test_round():
    dtype = "float64[pyarrow]"

    ser = pd.Series([0.0, 1.23, 2.56, pd.NA], dtype=dtype)
    result = ser.round(1)
    expected = pd.Series([0.0, 1.2, 2.6, pd.NA], dtype=dtype)
    tm.assert_series_equal(result, expected)

    ser = pd.Series([123.4, pd.NA, 56.78], dtype=dtype)
    result = ser.round(-1)
    expected = pd.Series([120.0, pd.NA, 60.0], dtype=dtype)
    tm.assert_series_equal(result, expected)


def test_searchsorted_with_na_raises(data_for_sorting, as_series):
    # GH50447
    b, c, a = data_for_sorting
    arr = data_for_sorting.take([2, 0, 1])  # to get [a, b, c]
    arr[-1] = pd.NA

    if as_series:
        arr = pd.Series(arr)

    msg = (
        "searchsorted requires array to be sorted, "
        "which is impossible with NAs present."
    )
    with pytest.raises(ValueError, match=msg):
        arr.searchsorted(b)


def test_sort_values_dictionary():
    df = pd.DataFrame(
        {
            "a": pd.Series(
                ["x", "y"], dtype=ArrowDtype(pa.dictionary(pa.int32(), pa.string()))
            ),
            "b": [1, 2],
        },
    )
    expected = df.copy()
    result = df.sort_values(by=["a", "b"])
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "ascending, expected_b, expected_index",
    [
        (True, [1, 2], [1, 0]),
        ([False, True], [1, 2], [1, 0]),
        ([True, False], [2, 1], [0, 1]),
    ],
)
def test_sort_values_null(ascending, expected_b, expected_index):
    # GH#54908
    df = pd.DataFrame(
        {
            "a": pd.Series([None, None], dtype="null[pyarrow]"),
            "b": [2, 1],
        }
    )
    result = df.sort_values(["a", "b"], ascending=ascending)
    expected = pd.DataFrame(
        {
            "a": pd.Series([None, None], dtype="null[pyarrow]"),
            "b": expected_b,
        },
        index=expected_index,
    )
    tm.assert_frame_equal(result, expected)


def test_sort_values_null_empty():
    # GH#54908
    df = pd.DataFrame({"a": pd.Series([], dtype="null[pyarrow]")})
    result = df.sort_values(by="a")
    tm.assert_frame_equal(result, df)


def test_sort_values_null_series():
    # GH#54908
    ser = pd.Series([None, None], dtype="null[pyarrow]")
    result = ser.sort_values()
    tm.assert_series_equal(result, ser)


def test_fillna_zero():
    # https://github.com/pandas-dev/pandas/issues/62878 - specific pyarrow bug
    ser = pd.Series([1, 2, 3, 4, pd.NA, 6], dtype="int64[pyarrow]")
    result = ser.fillna(0)
    expected = pd.Series([1, 2, 3, 4, 0, 6], dtype="int64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_sort_readonly():
    arr = pd.array([3, 1, 2], dtype="int64[pyarrow]")
    arr._readonly = True
    with pytest.raises(ValueError, match="Cannot modify read-only array"):
        arr.sort()
    # the array must be left unchanged
    tm.assert_extension_array_equal(arr, pd.array([3, 1, 2], dtype="int64[pyarrow]"))


@pytest.mark.parametrize("pa_type", tm.ALL_INT_PYARROW_DTYPES + tm.FLOAT_PYARROW_DTYPES)
def test_describe_numeric_data(pa_type):
    # GH 52470
    data = pd.Series([1, 2, 3], dtype=ArrowDtype(pa_type))
    result = data.describe()
    expected = pd.Series(
        [3, 2, 1, 1, 1.5, 2.0, 2.5, 3],
        dtype=ArrowDtype(pa.float64()),
        index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("pa_type", tm.TIMEDELTA_PYARROW_DTYPES)
def test_describe_timedelta_data(pa_type):
    # GH53001
    data = pd.Series(range(1, 10), dtype=ArrowDtype(pa_type))
    result = data.describe()
    expected = pd.Series(
        [9, *pd.to_timedelta([5, 2, 1, 3, 5, 7, 9], unit=pa_type.unit).tolist()],
        dtype=object,
        index=["count", "mean", "std", "min", "25%", "50%", "75%", "max"],
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("pa_type", tm.DATETIME_PYARROW_DTYPES)
def test_describe_datetime_data(pa_type):
    # GH53001
    data = pd.Series(range(1, 10), dtype=ArrowDtype(pa_type))
    result = data.describe()
    expected = pd.Series(
        [9]
        + [
            pd.Timestamp(v, tz=pa_type.tz, unit=pa_type.unit)
            for v in [5, 1, 3, 5, 7, 9]
        ],
        dtype=object,
        index=["count", "mean", "min", "25%", "50%", "75%", "max"],
    )
    tm.assert_series_equal(result, expected)


def test_ufunc_retains_missing():
    # GH#62800
    ser = pd.Series([0.1, pd.NA], dtype="float64[pyarrow]")

    result = np.sin(ser)

    expected = pd.Series([np.sin(0.1), pd.NA], dtype="float64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_np_ufunc_pyarrow_distinguish_nan_na():
    # GH#62506 - ufuncs on pyarrow arrays with distinguish_nan_and_na=True
    # should work instead of raising TypeError from object dtype conversion.
    with pd.option_context("future.distinguish_nan_and_na", True):
        ser = pd.Series([1.0, float("nan"), None], dtype="double[pyarrow]")

        result = np.isnan(ser)
        expected = pd.Series([False, True, pd.NA], dtype="bool[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = np.isfinite(ser)
        expected = pd.Series([True, False, pd.NA], dtype="bool[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = np.sqrt(pd.Series([1.0, 4.0, None], dtype="double[pyarrow]"))
        expected = pd.Series([1.0, 2.0, pd.NA], dtype="double[pyarrow]")
        tm.assert_series_equal(result, expected)

        # multi-return ufunc (tuple path)
        ser = pd.Series([1.5, 2.7, None], dtype="double[pyarrow]")
        frac, integ = np.modf(ser)
        expected_frac = pd.Series([0.5, 0.7, pd.NA], dtype="double[pyarrow]")
        expected_integ = pd.Series([1.0, 2.0, pd.NA], dtype="double[pyarrow]")
        tm.assert_series_equal(frac, expected_frac)
        tm.assert_series_equal(integ, expected_integ)


def test_map_numeric_na_action():
    # GH#62164 - _cast_pointwise_result retains Arrow dtype
    ser = pd.Series([32, 40, None], dtype="int64[pyarrow]")
    result = ser.map(lambda x: 42, na_action="ignore")
    expected = pd.Series([42, 42, None], dtype="int64[pyarrow]")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "dtype, func, expected_dtype",
    [
        (
            ArrowDtype(pa.timestamp("s")),
            lambda x: x + timedelta(microseconds=123456),
            ArrowDtype(pa.timestamp("us")),
        ),
        (
            ArrowDtype(pa.time32("s")),
            lambda x: x.replace(microsecond=123456),
            ArrowDtype(pa.time64("us")),
        ),
        (
            ArrowDtype(pa.decimal128(3, 1)),
            lambda x: x + Decimal("0.01"),
            ArrowDtype(pa.decimal128(3, 2)),
        ),
    ],
)
def test_map_finer_resolution_no_arrowinvalid(dtype, func, expected_dtype):
    # GH#62523 mapping to a finer resolution/precision that cannot be cast
    #  down to the original dtype should keep the finer dtype rather than
    #  raising ArrowInvalid
    if pa.types.is_timestamp(dtype.pyarrow_dtype):
        data = [datetime(2020, 1, 1), datetime(2020, 1, 2)]
    elif pa.types.is_time(dtype.pyarrow_dtype):
        data = [time(1, 2, 3), time(4, 5, 6)]
    else:
        data = [Decimal("1.5"), Decimal("2.5")]

    ser = pd.Series(data, dtype=dtype)
    result = ser.map(func)
    assert result.dtype == expected_dtype
    expected = pd.Series([func(x) for x in data], dtype=expected_dtype)
    tm.assert_series_equal(result, expected)


def test_cast_pontwise_result_decimal_nan():
    # GH#62522 we don't want to get back null[pyarrow] here
    ser = pd.Series([], dtype="float64[pyarrow]")
    arr = ser.array
    item = Decimal("NaN")

    result = arr._cast_pointwise_result([item])

    pa_type = result.dtype.pyarrow_dtype
    assert pa.types.is_decimal(pa_type)


@pytest.mark.parametrize("pa_type", tm.TIMEDELTA_PYARROW_DTYPES)
def test_duration_fillna_numpy(pa_type):
    # GH 54707
    ser1 = pd.Series([None, 2], dtype=ArrowDtype(pa_type))
    ser2 = pd.Series(np.array([1, 3], dtype=f"m8[{pa_type.unit}]"))
    result = ser1.fillna(ser2)
    expected = pd.Series([1, 2], dtype=ArrowDtype(pa_type))
    tm.assert_series_equal(result, expected)


def test_factorize_chunked_dictionary():
    # GH 54844
    pa_array = pa.chunked_array(
        [pa.array(["a"]).dictionary_encode(), pa.array(["b"]).dictionary_encode()]
    )
    ser = pd.Series(ArrowExtensionArray(pa_array))
    res_indices, res_uniques = ser.factorize()
    exp_indices = np.array([0, 1], dtype=np.intp)
    exp_uniques = pd.Index(ArrowExtensionArray(pa_array.combine_chunks()))
    tm.assert_numpy_array_equal(res_indices, exp_indices)
    tm.assert_index_equal(res_uniques, exp_uniques)


def test_factorize_dictionary_with_na():
    # GH#60567
    arr = pd.array(
        ["a1", pd.NA], dtype=ArrowDtype(pa.dictionary(pa.int32(), pa.utf8()))
    )
    indices, uniques = arr.factorize(use_na_sentinel=False)
    expected_indices = np.array([0, 1], dtype=np.intp)
    expected_uniques = pd.array(["a1", None], dtype=ArrowDtype(pa.string()))
    tm.assert_numpy_array_equal(indices, expected_indices)
    tm.assert_extension_array_equal(uniques, expected_uniques)


def _reencodable(dict_arr) -> bool:
    # GH#69024 factorize re-encodes a stored dictionary in index space or by
    #  decoding it; which value types support either varies across pyarrow versions
    values = dict_arr.dictionary
    try:
        values.dictionary_encode()
        pa.compute.take(values, pa.array([0], type=pa.int32()))
        return True
    except pa.ArrowNotImplementedError:
        pass
    try:
        dict_arr.cast(values.type).dictionary_encode()
        return True
    except pa.ArrowNotImplementedError:
        return False


@pytest.mark.parametrize("use_na_sentinel", [True, False])
@pytest.mark.parametrize(
    "values",
    [
        pa.array(["a", "b"], type=pa.string_view()),
        pa.array([b"a", b"b"], type=pa.binary_view()),
        pa.array([[1], [2]], type=pa.list_(pa.int64())),
        pa.array(np.array([1, 2], dtype=np.float16)),
    ],
)
def test_factorize_dictionary_unsupported_value_type(values, use_na_sentinel):
    # GH#69024 factorize keeps the stored dictionary for a value type it cannot
    #  re-encode, rather than raising
    dict_arr = pa.DictionaryArray.from_arrays(
        pa.array([0, 1, 0], type=pa.int32()), values
    )
    if _reencodable(dict_arr):
        pytest.skip(f"pyarrow can re-encode a {values.type} dictionary")
    arr = pd.array(dict_arr, dtype=ArrowDtype(dict_arr.type))

    indices, uniques = arr.factorize(use_na_sentinel=use_na_sentinel)
    tm.assert_numpy_array_equal(indices, np.array([0, 1, 0], dtype=np.intp))
    # compared as pyarrow, since ArrowDtype.type raises for the view types
    assert uniques._pa_array.combine_chunks().equals(values)


@pytest.mark.parametrize(
    "values",
    [
        pa.array(["a", "b"], type=pa.string_view()),
        pa.array([b"a", b"b"], type=pa.binary_view()),
        pa.array([[1], [2]], type=pa.list_(pa.int64())),
        pa.array(np.array([1, 2], dtype=np.float16)),
    ],
)
def test_factorize_dictionary_unsupported_value_type_with_na(values):
    # GH#69024 keeping the stored dictionary leaves the null in index space, which
    #  the sentinel can express and use_na_sentinel=False cannot, so that leg raises
    #  rather than handing back a -1 it promised not to use
    dict_arr = pa.DictionaryArray.from_arrays(
        pa.array([0, None, 1], type=pa.int32()), values
    )
    if _reencodable(dict_arr):
        pytest.skip(f"pyarrow can re-encode a {values.type} dictionary")
    arr = pd.array(dict_arr, dtype=ArrowDtype(dict_arr.type))

    indices, uniques = arr.factorize()
    tm.assert_numpy_array_equal(indices, np.array([0, -1, 1], dtype=np.intp))
    assert uniques._pa_array.combine_chunks().equals(values)

    with pytest.raises(pa.ArrowNotImplementedError):
        arr.factorize(use_na_sentinel=False)


@pytest.mark.parametrize(
    "use_na_sentinel",
    [
        pytest.param(
            True,
            marks=pytest.mark.xfail(
                reason="the null stays in the kept dictionary and in uniques",
                strict=True,
            ),
        ),
        False,
    ],
)
@pytest.mark.parametrize(
    "values",
    [
        pa.array(["a", None], type=pa.string_view()),
        pa.array([b"a", None], type=pa.binary_view()),
        pa.array([[1], None], type=pa.list_(pa.int64())),
        pa.array(np.array([1, 0], dtype=np.float16), mask=np.array([False, True])),
    ],
)
def test_factorize_dictionary_unsupported_value_type_null_in_dictionary(
    values, use_na_sentinel
):
    # GH#69024 a null stored only in a kept dictionary is already the separate
    #  code use_na_sentinel=False asks for
    dict_arr = pa.DictionaryArray.from_arrays(
        pa.array([0, 1, 0], type=pa.int32()), values
    )
    if _reencodable(dict_arr):
        pytest.skip(f"pyarrow can re-encode a {values.type} dictionary")
    arr = pd.array(dict_arr, dtype=ArrowDtype(dict_arr.type))

    indices, uniques = arr.factorize(use_na_sentinel=use_na_sentinel)
    if use_na_sentinel:
        tm.assert_numpy_array_equal(indices, np.array([0, -1, 0], dtype=np.intp))
        assert uniques._pa_array.combine_chunks().equals(values[:1])
    else:
        tm.assert_numpy_array_equal(indices, np.array([0, 1, 0], dtype=np.intp))
        assert uniques._pa_array.combine_chunks().equals(values)


@pytest.mark.parametrize("use_na_sentinel", [True, False])
@pytest.mark.parametrize("empty", [False, True])
def test_factorize_dictionary_chunked_null_in_dictionary(empty, use_na_sentinel):
    # GH#69024 pyarrow refuses to unify chunk dictionaries when one holds a null
    enc = pa.array(["a", None, "a"]).dictionary_encode(null_encoding="encode")
    dtype = ArrowDtype(enc.type)
    pieces = [
        pd.Series(pd.array(enc, dtype=dtype)),
        pd.Series(pd.array(pa.array(["b"]).dictionary_encode(), dtype=dtype)),
    ]
    if empty:
        pieces = [piece.iloc[:0] for piece in pieces]
    arr = pd.concat(pieces, ignore_index=True).array
    assert arr._pa_array.num_chunks == 2

    indices, uniques = arr.factorize(use_na_sentinel=use_na_sentinel)
    if empty:
        expected_indices = np.array([], dtype=np.intp)
        expected_uniques = []
    elif use_na_sentinel:
        expected_indices = np.array([0, -1, 0, 1], dtype=np.intp)
        expected_uniques = ["a", "b"]
    else:
        expected_indices = np.array([0, 1, 0, 2], dtype=np.intp)
        expected_uniques = ["a", None, "b"]
    tm.assert_numpy_array_equal(indices, expected_indices)
    tm.assert_extension_array_equal(
        uniques, pd.array(expected_uniques, dtype=ArrowDtype(pa.string()))
    )


def test_factorize_dictionary_null_in_dictionary():
    # GH#69024 a null stored in the dictionary rather than the indices must still
    #  reach the sentinel, and must not show up in uniques
    dict_arr = pa.DictionaryArray.from_arrays(
        pa.array([0, 1, 0], type=pa.int32()), pa.array(["a", None])
    )
    arr = pd.array(dict_arr, dtype=ArrowDtype(dict_arr.type))

    indices, uniques = arr.factorize()
    tm.assert_numpy_array_equal(indices, np.array([0, -1, 0], dtype=np.intp))
    tm.assert_extension_array_equal(
        uniques, pd.array(["a"], dtype=ArrowDtype(pa.string()))
    )

    indices, uniques = arr.factorize(use_na_sentinel=False)
    tm.assert_numpy_array_equal(indices, np.array([0, 1, 0], dtype=np.intp))
    tm.assert_extension_array_equal(
        uniques, pd.array(["a", None], dtype=ArrowDtype(pa.string()))
    )


def test_factorize_dictionary_null_in_dictionary_and_indices():
    # GH#69024 a null entry and a null index are the same missing value, so they
    #  share one code rather than splitting NA across two groups
    dict_arr = pa.DictionaryArray.from_arrays(
        pa.array([0, None, 1], type=pa.int32()), pa.array(["a", None])
    )
    arr = pd.array(dict_arr, dtype=ArrowDtype(dict_arr.type))

    indices, uniques = arr.factorize(use_na_sentinel=False)
    tm.assert_numpy_array_equal(indices, np.array([0, 1, 1], dtype=np.intp))
    tm.assert_extension_array_equal(
        uniques, pd.array(["a", None], dtype=ArrowDtype(pa.string()))
    )


@pytest.mark.parametrize("use_na_sentinel", [True, False])
def test_factorize_dictionary_of_dictionary(use_na_sentinel):
    # GH#69024 the inner dictionary is re-encoded too: "b", which no index points
    #  at, is not a unique, and a null index gets its own code when asked
    inner = pa.array(["a", "b"]).dictionary_encode()
    dict_arr = pa.DictionaryArray.from_arrays(
        pa.array([0, None, 0], type=pa.int32()), inner
    )
    arr = pd.array(dict_arr, dtype=ArrowDtype(dict_arr.type))

    indices, uniques = arr.factorize(use_na_sentinel=use_na_sentinel)
    if use_na_sentinel:
        expected_indices = np.array([0, -1, 0], dtype=np.intp)
        expected_uniques = ["a"]
    else:
        expected_indices = np.array([0, 1, 0], dtype=np.intp)
        expected_uniques = ["a", None]
    tm.assert_numpy_array_equal(indices, expected_indices)
    tm.assert_extension_array_equal(
        uniques, pd.array(expected_uniques, dtype=ArrowDtype(pa.string()))
    )


def test_factorize_dictionary_unreferenced_entries():
    # GH#69024 a dictionary may keep entries no index points at, e.g. after a
    #  filter, and those are not uniques of the values that remain
    enc = pa.array(["a", "b", "c", "a"]).dictionary_encode()
    arr = pd.array(enc, dtype=ArrowDtype(enc.type))[2:]

    indices, uniques = arr.factorize()
    tm.assert_numpy_array_equal(indices, np.array([0, 1], dtype=np.intp))
    tm.assert_extension_array_equal(
        uniques, pd.array(["c", "a"], dtype=ArrowDtype(pa.string()))
    )

    df = pd.DataFrame({"k": arr, "v": [3, 4]})
    result = df.groupby("k", sort=False)["v"].sum()
    expected = pd.Series(
        [3, 4],
        name="v",
        index=pd.Index(pd.array(["c", "a"], dtype=arr.dtype), name="k"),
    )
    tm.assert_series_equal(result, expected)


def test_factorize_dictionary_duplicate_entries():
    # GH#69024 a dictionary may hold the same value twice, which would otherwise
    #  split one group in two
    dict_arr = pa.DictionaryArray.from_arrays(
        pa.array([0, 1, 0], type=pa.int32()), pa.array(["a", "a", "b"])
    )
    arr = pd.array(dict_arr, dtype=ArrowDtype(dict_arr.type))

    indices, uniques = arr.factorize()
    tm.assert_numpy_array_equal(indices, np.array([0, 0, 0], dtype=np.intp))
    tm.assert_extension_array_equal(
        uniques, pd.array(["a"], dtype=ArrowDtype(pa.string()))
    )


def test_factorize_dictionary_chunks_overflow_index_type():
    # GH#69024 the chunks' dictionaries unify to more entries than int8 indices hold
    dtype = pa.dictionary(pa.int8(), pa.string())
    chunks = [
        pa.array([str(i) for i in range(start, start + 100)])
        .dictionary_encode()
        .cast(dtype)
        for start in [0, 100]
    ]
    arr = ArrowExtensionArray(pa.chunked_array(chunks))

    indices, uniques = arr.factorize()
    tm.assert_numpy_array_equal(indices, np.arange(200, dtype=np.intp))
    tm.assert_extension_array_equal(
        uniques,
        pd.array([str(i) for i in range(200)], dtype=ArrowDtype(pa.string())),
    )


def test_factorize_null():
    # GH#54908
    arr = ArrowExtensionArray(pa.array([None, None], type=pa.null()))
    indices, uniques = arr.factorize(use_na_sentinel=True)
    expected_indices = np.array([-1, -1], dtype=np.intp)
    expected_uniques = ArrowExtensionArray(pa.chunked_array([], type=pa.null()))
    tm.assert_numpy_array_equal(indices, expected_indices)
    tm.assert_extension_array_equal(uniques, expected_uniques)

    indices, uniques = arr.factorize(use_na_sentinel=False)
    expected_indices = np.array([0, 0], dtype=np.intp)
    expected_uniques = ArrowExtensionArray(pa.array([None], type=pa.null()))
    tm.assert_numpy_array_equal(indices, expected_indices)
    tm.assert_extension_array_equal(uniques, expected_uniques)


def test_factorize_null_empty():
    # GH#54908
    arr = ArrowExtensionArray(pa.array([], type=pa.null()))
    indices, uniques = arr.factorize(use_na_sentinel=True)
    expected_indices = np.array([], dtype=np.intp)
    expected_uniques = ArrowExtensionArray(pa.chunked_array([], type=pa.null()))
    tm.assert_numpy_array_equal(indices, expected_indices)
    tm.assert_extension_array_equal(uniques, expected_uniques)

    indices, uniques = arr.factorize(use_na_sentinel=False)
    expected_indices = np.array([], dtype=np.intp)
    expected_uniques = ArrowExtensionArray(pa.chunked_array([], type=pa.null()))
    tm.assert_numpy_array_equal(indices, expected_indices)
    tm.assert_extension_array_equal(uniques, expected_uniques)


def test_value_counts_returns_pyarrow_int64(data):
    # GH 51462
    data = data[:10]
    result = data.value_counts()
    assert result.dtype == ArrowDtype(pa.int64())


def test_interpolate_not_numeric(data):
    if not data.dtype._is_numeric:
        ser = pd.Series(data)
        msg = re.escape(f"Cannot interpolate with {ser.dtype} dtype")
        with pytest.raises(TypeError, match=msg):
            pd.Series(data).interpolate()


@pytest.mark.parametrize("dtype", ["int64[pyarrow]", "float64[pyarrow]"])
def test_interpolate_linear(dtype):
    # GH#65345 results should match the masked (e.g. Int64) dtypes:
    # upcast to float, and fill the trailing NA going forward
    ser = pd.Series([None, 1, 2, None, 4, None], dtype=dtype)
    result = ser.interpolate()
    expected = pd.Series([None, 1.0, 2.0, 3.0, 4.0, 4.0], dtype="float64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_interpolate_linear_consecutive_na():
    # GH#65345 consecutive interior NAs were left unfilled
    ser = pd.Series([1, 2, 3, None, None, 6, 7], dtype="int64[pyarrow]")
    result = ser.interpolate(method="linear", limit_direction="forward")
    expected = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], dtype="float64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_interpolate_linear_int_fractional():
    # GH#65345 result should not truncate the interpolated value (1 instead of 1.5)
    ser = pd.Series([1, None, 2], dtype="int64[pyarrow]")
    result = ser.interpolate(method="linear")
    expected = pd.Series([1.0, 1.5, 2.0], dtype="float64[pyarrow]")
    tm.assert_series_equal(result, expected)
