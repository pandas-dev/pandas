from decimal import Decimal

import numpy as np
import pytest

import pandas as pd

pa = pytest.importorskip("pyarrow")


def test_int32_scalar_arithmetic_preserves_dtype():
    ser = pd.Series([1, 2, 3], dtype="int32[pyarrow]")

    result = ser + 1
    assert result.tolist() == [2, 3, 4]
    assert result.dtype == "int32[pyarrow]"

    result = ser * 2
    assert result.tolist() == [2, 4, 6]
    assert result.dtype == "int32[pyarrow]"

    result = ser + np.int32(1)
    assert result.tolist() == [2, 3, 4]
    assert result.dtype == "int32[pyarrow]"

    result = ser + ser[0]
    assert result.tolist() == [2, 3, 4]
    assert result.dtype == "int32[pyarrow]"

    promoted = ser + 2**31
    assert promoted.dtype == "int64[pyarrow]"

    overflowing = pd.Series([2**31 - 1], dtype="int32[pyarrow]") + 1
    assert overflowing.tolist() == [2**31]
    assert overflowing.dtype == "int64[pyarrow]"


def test_float32_scalar_arithmetic_preserves_dtype():
    ser = pd.Series([1, 2, 3], dtype="float[pyarrow]")

    result = ser + 1.5
    assert result.tolist() == [2.5, 3.5, 4.5]
    assert result.dtype == "float[pyarrow]"

    result = ser * 1.5
    assert result.tolist() == [1.5, 3.0, 4.5]
    assert result.dtype == "float[pyarrow]"

    result = ser + np.float32(1.5)
    assert result.tolist() == [2.5, 3.5, 4.5]
    assert result.dtype == "float[pyarrow]"

    result = ser + ser[0]
    assert result.tolist() == [2.0, 3.0, 4.0]
    assert result.dtype == "float[pyarrow]"

    large = pd.Series([1e38], dtype="float[pyarrow]")
    small_product = large * 1e-50
    assert small_product.iloc[0] == pytest.approx(1e-12, rel=1e-6)
    assert small_product.dtype == "float[pyarrow]"

    large_product = large * 1e100
    assert np.isfinite(large_product.iloc[0])
    assert large_product.dtype == "double[pyarrow]"


def test_decimal_scalar_arithmetic_dtype_promotion():
    dtype = pd.ArrowDtype(pa.decimal128(10, 2))
    ser = pd.Series([Decimal("1.00"), Decimal("2.00")], dtype=dtype)

    added = ser + 1
    multiplied = ser * 2

    assert added.tolist() == [Decimal("2.00"), Decimal("3.00")]
    assert added.dtype == dtype
    assert multiplied.tolist() == [Decimal("2.00"), Decimal("4.00")]
    assert multiplied.dtype == dtype

    narrow_dtype = pd.ArrowDtype(pa.decimal128(2, 0))
    narrow = pd.Series([Decimal("99")], dtype=narrow_dtype)
    overflow = narrow + 1

    assert overflow.tolist() == [Decimal("100")]
    assert overflow.dtype != narrow_dtype
    assert overflow.dtype.pyarrow_dtype.precision >= 3
