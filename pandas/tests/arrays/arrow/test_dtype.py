from __future__ import annotations

import numpy as np
import pytest

from pandas._libs import lib
from pandas.compat import (
    pa_version_under19p0,
)

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
    CategoricalDtypeType,
)

import pandas as pd
import pandas._testing as tm
from pandas.api.types import (
    is_bool_dtype,
    is_datetime64_any_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_signed_integer_dtype,
    is_unsigned_integer_dtype,
)

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray
from pandas.core.arrays.arrow.extension_types import ArrowPeriodType


@pytest.mark.parametrize(
    "arrow_dtype, expected_type",
    [
        [pa.binary(), bytes],
        [pa.binary(16), bytes],
        [pa.large_binary(), bytes],
        [pa.large_string(), str],
        [pa.list_(pa.int64()), list],
        [pa.large_list(pa.int64()), list],
        [pa.map_(pa.string(), pa.int64()), list],
        [pa.struct([("f1", pa.int8()), ("f2", pa.string())]), dict],
        [pa.dictionary(pa.int64(), pa.int64()), CategoricalDtypeType],
    ],
)
def test_arrow_dtype_type(arrow_dtype, expected_type):
    # GH 51845
    # TODO: Redundant with test_getitem_scalar once arrow_dtype exists in data fixture
    assert ArrowDtype(arrow_dtype).type == expected_type


def test_is_bool_dtype():
    # GH 22667
    data = ArrowExtensionArray(pa.array([True, False, True]))
    assert is_bool_dtype(data)
    assert pd.core.common.is_bool_indexer(data)
    s = pd.Series(range(len(data)))
    result = s[data]
    expected = s[np.asarray(data)]
    tm.assert_series_equal(result, expected)


def test_is_numeric_dtype(data):
    # GH 50563
    pa_type = data.dtype.pyarrow_dtype
    if (
        pa.types.is_floating(pa_type)
        or pa.types.is_integer(pa_type)
        or pa.types.is_decimal(pa_type)
    ):
        assert is_numeric_dtype(data)
    else:
        assert not is_numeric_dtype(data)


def test_is_integer_dtype(data):
    # GH 50667
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_integer(pa_type):
        assert is_integer_dtype(data)
    else:
        assert not is_integer_dtype(data)


def test_is_signed_integer_dtype(data):
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_signed_integer(pa_type):
        assert is_signed_integer_dtype(data)
    else:
        assert not is_signed_integer_dtype(data)


def test_is_unsigned_integer_dtype(data):
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_unsigned_integer(pa_type):
        assert is_unsigned_integer_dtype(data)
    else:
        assert not is_unsigned_integer_dtype(data)


def test_is_datetime64_any_dtype(data):
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_timestamp(pa_type) or pa.types.is_date(pa_type):
        assert is_datetime64_any_dtype(data)
    else:
        assert not is_datetime64_any_dtype(data)


def test_is_float_dtype(data):
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_floating(pa_type):
        assert is_float_dtype(data)
    else:
        assert not is_float_dtype(data)


def test_infer_dtype_pyarrow_dtype(data, request):
    res = lib.infer_dtype(data)
    assert res != "unknown-array"

    if res in ["datetime64", "timedelta64"]:
        # infer_dtype on the pyarrow-backed array returns datetime64/timedelta64
        # via _TYPE_MAP, but infer_dtype on list(data) returns datetime/timedelta
        # because the elements are pd.Timestamp/pd.Timedelta (PyDateTime/PyDelta).
        mark = pytest.mark.xfail(
            reason="infer_dtype(arrow_array) vs infer_dtype(list) naming mismatch"
        )
        request.applymarker(mark)

    assert res == lib.infer_dtype(list(data), skipna=True)


def test_fixed_size_list():
    # GH#55000
    ser = pd.Series(
        [[1, 2], [3, 4]], dtype=ArrowDtype(pa.list_(pa.int64(), list_size=2))
    )
    result = ser.dtype.type
    assert result == list


@pytest.mark.skipif(
    pa_version_under19p0, reason="pa.json_ was introduced in pyarrow v19.0"
)
def test_arrow_json_type():
    # GH 60958
    dtype = ArrowDtype(pa.json_(pa.string()))
    result = dtype.type
    assert result == str


@pytest.mark.parametrize(
    "type_name, expected_size",
    [
        # Integer types
        ("int8", 1),
        ("int16", 2),
        ("int32", 4),
        ("int64", 8),
        ("uint8", 1),
        ("uint16", 2),
        ("uint32", 4),
        ("uint64", 8),
        # Floating point types
        ("float16", 2),
        ("float32", 4),
        ("float64", 8),
        # Boolean
        ("bool_", 1),
        # Date and timestamp types
        ("date32", 4),
        ("date64", 8),
        ("timestamp", 8),
        # Time types
        ("time32", 4),
        ("time64", 8),
        # Decimal types
        ("decimal128", 16),
        ("decimal256", 32),
    ],
)
def test_arrow_dtype_itemsize_fixed_width(type_name, expected_size):
    # GH 57948

    parametric_type_map = {
        "timestamp": pa.timestamp("ns"),
        "time32": pa.time32("s"),
        "time64": pa.time64("ns"),
        "decimal128": pa.decimal128(38, 10),
        "decimal256": pa.decimal256(76, 10),
    }

    if type_name in parametric_type_map:
        arrow_type = parametric_type_map.get(type_name)
    else:
        arrow_type = getattr(pa, type_name)()
    dtype = ArrowDtype(arrow_type)

    if type_name == "bool_":
        expected_size = dtype.numpy_dtype.itemsize

    assert dtype.itemsize == expected_size, (
        f"{type_name} expected {expected_size}, got {dtype.itemsize} "
        f"(bit_width={getattr(dtype.pyarrow_dtype, 'bit_width', 'N/A')})"
    )


@pytest.mark.parametrize("type_name", ["string", "binary", "large_string"])
def test_arrow_dtype_itemsize_variable_width(type_name):
    # GH 57948

    arrow_type = getattr(pa, type_name)()
    dtype = ArrowDtype(arrow_type)

    assert dtype.itemsize == dtype.numpy_dtype.itemsize


def test_arrowextensiondtype_dataframe_repr():
    # GH 54062
    df = pd.DataFrame(
        pd.period_range("2012", periods=3),
        columns=["col"],
        dtype=ArrowDtype(ArrowPeriodType("D")),
    )
    result = repr(df)
    # TODO: repr value may not be expected; address how
    # pyarrow.ExtensionType values are displayed
    expected = "     col\n0  15340\n1  15341\n2  15342"
    assert result == expected
