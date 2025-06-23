import numpy as np
import pytest

from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    DatetimeTZDtype,
    IntervalDtype,
    PeriodDtype,
    SparseDtype,
)
from pandas.core.dtypes.factory import (
    boolean,
    categorical,
    date,
    datetime,
    decimal,
    duration,
    floating,
    integer,
    interval,
    list,
    map,
    period,
    sparse,
    string,
    struct,
)

from pandas import Series
from pandas.core.api import (
    ArrowDtype,
    BooleanDtype,
    Float32Dtype,
    Float64Dtype,
    Int32Dtype,
    Int64Dtype,
    StringDtype,
)


# String
def test_string_default():
    result = string()
    assert result == StringDtype()
    assert str(result) == "string"


def test_string_with_mode():
    pa = pytest.importorskip("pyarrow")
    result = string(mode="binary", backend="pyarrow")
    assert result == ArrowDtype(pa.binary())
    assert str(result) == "binary[pyarrow]"


def test_string_invalid_mode():
    with pytest.raises(ValueError, match="mode must be one of"):
        string(mode="invalid", backend="pyarrow")


# Datetime
def test_datetime_default():
    result = datetime()
    assert result == np.dtype("datetime64[ns]")
    assert isinstance(result, np.dtype)


def test_datetime_with_tz():
    result = datetime(tz="UTC")
    assert isinstance(result, DatetimeTZDtype)
    assert str(result) == "datetime64[ns, UTC]"


def test_datetime_pyarrow():
    pytest.importorskip("pyarrow")
    result = datetime(backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "timestamp[ns][pyarrow]"


def test_datetime_invalid_unit():
    with pytest.raises(ValueError, match="unit must be one of"):
        datetime(unit="invalid", backend="numpy")


# Integer
def test_integer_default():
    result = integer()
    assert result == Int64Dtype()
    assert str(result) == "Int64"


def test_integer_with_bits():
    result = integer(bits=32)
    assert result == Int32Dtype()
    assert str(result) == "Int32"


def test_integer_numpy():
    result = integer(bits=64, backend="numpy")
    assert result == np.dtype("int64")
    assert str(result) == "int64"


def test_integer_pyarrow():
    pytest.importorskip("pyarrow")
    result = integer(bits=64, backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "int64[pyarrow]"


# Floating
def test_floating_default():
    result = floating()
    assert result == Float64Dtype()
    assert str(result) == "Float64"


def test_floating_with_bits():
    result = floating(bits=32)
    assert result == Float32Dtype()
    assert str(result) == "Float32"


def test_floating_numpy():
    result = floating(bits=64, backend="numpy")
    assert result == np.dtype("float64")
    assert str(result) == "float64"


def test_floating_pyarrow():
    pytest.importorskip("pyarrow")
    result = floating(bits=64, backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "double[pyarrow]"


# Decimal
def test_decimal_default():
    pytest.importorskip("pyarrow")
    result = decimal(precision=38, scale=10)
    assert isinstance(result, ArrowDtype)
    assert str(result) == "decimal128(38, 10)[pyarrow]"


def test_decimal_with_precision_scale():
    pytest.importorskip("pyarrow")
    result = decimal(precision=10, scale=2)
    assert isinstance(result, ArrowDtype)
    assert str(result) == "decimal128(10, 2)[pyarrow]"


# Boolean
def test_boolean_default():
    result = boolean()
    assert result == BooleanDtype()
    assert str(result) == "boolean"


def test_boolean_pyarrow():
    pytest.importorskip("pyarrow")
    result = boolean(backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "bool[pyarrow]"


# List
def test_list_default():
    result = list()
    assert result == np.dtype("object")
    assert isinstance(result, np.dtype)


def test_list_pyarrow():
    pa = pytest.importorskip("pyarrow")
    result = list(backend="pyarrow", value_type=pa.int64())
    assert isinstance(result, ArrowDtype)
    assert str(result) == "list<item: int64>[pyarrow]"


def test_list_large():
    pa = pytest.importorskip("pyarrow")
    result = list(backend="pyarrow", value_type=pa.string(), large=True)
    assert isinstance(result, ArrowDtype)
    assert str(result) == "large_list<item: string>[pyarrow]"


# Categorical
def test_categorical_default():
    result = categorical()
    assert isinstance(result, CategoricalDtype)
    assert str(result) == "category"


def test_categorical_pyarrow():
    pytest.importorskip("pyarrow")
    result = categorical(backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "dictionary<values=string, indices=int32, ordered=0>[pyarrow]"


# Interval
def test_interval_default():
    result = interval()
    assert isinstance(result, IntervalDtype)
    assert str(result) == "interval"


def test_interval_pyarrow():
    pytest.importorskip("pyarrow")
    result = interval(backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "struct<left: double, right: double>[pyarrow]"


# Period
def test_period_default():
    result = period()
    assert isinstance(result, PeriodDtype)
    assert str(result) == "period[D]"


def test_period_pyarrow():
    pytest.importorskip("pyarrow")
    result = period(backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "month_day_nano_interval[pyarrow]"


# Date
def test_date_default():
    pytest.importorskip("pyarrow")
    result = date()
    assert isinstance(result, ArrowDtype)
    assert str(result) == "date32[day][pyarrow]"


def test_date_pyarrow():
    pytest.importorskip("pyarrow")
    result = date(backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "date32[day][pyarrow]"


# Duration
def test_duration_default():
    pytest.importorskip("pyarrow")
    result = duration()
    assert isinstance(result, ArrowDtype)
    assert str(result) == "duration[ns][pyarrow]"


def test_duration_pyarrow():
    pytest.importorskip("pyarrow")
    result = duration(backend="pyarrow")
    assert isinstance(result, ArrowDtype)
    assert str(result) == "duration[ns][pyarrow]"


# Map
def test_map_default():
    pa = pytest.importorskip("pyarrow")
    result = map(index_type=pa.string(), value_type=pa.int64())
    assert isinstance(result, ArrowDtype)
    assert str(result) == "map<string, int64>[pyarrow]"


def test_map_custom_types():
    pa = pytest.importorskip("pyarrow")
    result = map(index_type=pa.string(), value_type=pa.float64())
    assert isinstance(result, ArrowDtype)
    assert str(result) == "map<string, double>[pyarrow]"


# Struct
def test_struct_default():
    pa = pytest.importorskip("pyarrow")
    result = struct(fields=[("a", pa.int64()), ("b", pa.string())])
    assert isinstance(result, ArrowDtype)
    assert str(result) == "struct<a: int64, b: string>[pyarrow]"


def test_struct_custom_fields():
    pa = pytest.importorskip("pyarrow")
    fields = [("x", pa.float32()), ("y", pa.int16())]
    result = struct(fields=fields)
    assert isinstance(result, ArrowDtype)
    assert str(result) == "struct<x: float, y: int16>[pyarrow]"


# Sparse
def test_sparse_default():
    result = sparse()
    assert result == SparseDtype(np.float64, fill_value=np.nan)
    assert isinstance(result, SparseDtype)
    assert str(result) == "Sparse[float64, nan]"


def test_sparse_with_dtype():
    result = sparse(dtype=np.int64)
    assert result == SparseDtype(np.int64, fill_value=0)
    assert str(result) == "Sparse[int64, 0]"


def test_sparse_with_fill_value():
    result = sparse(fill_value=-1)
    assert result == SparseDtype(np.float64, fill_value=-1)
    assert str(result) == "Sparse[float64, -1]"


def test_sparse_backend_invalid():
    with pytest.raises(
        ValueError, match="Sparse types are only supported with NumPy backend"
    ):
        sparse(backend="pyarrow")


def test_sparse_series_creation():
    data = [1.0, None, None, 3.0, None]
    s_sparse = Series(data, dtype=sparse())
    assert s_sparse.dtype == SparseDtype(np.float64, fill_value=np.nan)
    assert s_sparse.memory_usage() < Series(data, dtype=np.float64).memory_usage()
