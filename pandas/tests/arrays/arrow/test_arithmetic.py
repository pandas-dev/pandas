from decimal import Decimal

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


def test_int_scalar_arithmetic_keeps_dtype():
    # GH#69436 arithmetic with a scalar that fits should keep the array's
    #  dtype, matching numpy, where a Python scalar adopts the array's dtype
    ser = pd.Series([1, 2, 3], dtype="int32[pyarrow]")

    tm.assert_series_equal(ser + 1, pd.Series([2, 3, 4], dtype="int32[pyarrow]"))
    tm.assert_series_equal(ser * 2, pd.Series([2, 4, 6], dtype="int32[pyarrow]"))
    tm.assert_series_equal(ser - 1, pd.Series([0, 1, 2], dtype="int32[pyarrow]"))
    tm.assert_series_equal(1 - ser, pd.Series([0, -1, -2], dtype="int32[pyarrow]"))
    tm.assert_series_equal(ser**2, pd.Series([1, 4, 9], dtype="int32[pyarrow]"))
    tm.assert_series_equal(
        ser + np.int32(1), pd.Series([2, 3, 4], dtype="int32[pyarrow]")
    )
    tm.assert_series_equal(ser + ser[0], pd.Series([2, 3, 4], dtype="int32[pyarrow]"))

    # division by an integer still promotes to float, as in numpy
    result = ser / 2
    assert result.dtype == "double[pyarrow]"
    tm.assert_series_equal(result, pd.Series([0.5, 1.0, 1.5], dtype="double[pyarrow]"))

    # a float scalar on an integer array promotes to float, as in numpy
    assert (ser + 1.5).dtype == "double[pyarrow]"


def test_int_scalar_arithmetic_upcasts_when_result_does_not_fit():
    # GH#69436 a result that does not fit in the array's dtype keeps the
    #  wider dtype instead of raising
    ser = pd.Series([1, 2, 3], dtype="int32[pyarrow]")

    assert (ser + 2**31).dtype == "int64[pyarrow]"
    assert (ser * 2**30).dtype == "int64[pyarrow]"

    overflowing = pd.Series([2**31 - 1], dtype="int32[pyarrow]") + 1
    assert overflowing.tolist() == [2**31]
    assert overflowing.dtype == "int64[pyarrow]"

    # an unsigned array going negative also keeps the wider dtype
    unsigned = pd.Series([1, 2, 3], dtype="uint8[pyarrow]")
    assert (unsigned - 300).dtype == "int64[pyarrow]"
    assert (unsigned + 300).dtype == "int64[pyarrow]"

    # unsigned dtypes are not narrowed at all: subtracting can go
    #  negative, so whether the result fits would depend on the data
    assert (unsigned + 1).dtype == "int64[pyarrow]"


def test_float_scalar_arithmetic_keeps_dtype():
    # GH#69436
    ser = pd.Series([1.0, 2.0, 3.0], dtype="float[pyarrow]")

    tm.assert_series_equal(
        ser + 1.5, pd.Series([2.5, 3.5, 4.5], dtype="float[pyarrow]")
    )
    tm.assert_series_equal(
        ser * 1.5, pd.Series([1.5, 3.0, 4.5], dtype="float[pyarrow]")
    )
    tm.assert_series_equal(
        ser + np.float32(1.5), pd.Series([2.5, 3.5, 4.5], dtype="float[pyarrow]")
    )
    tm.assert_series_equal(
        ser + ser[0], pd.Series([2.0, 3.0, 4.0], dtype="float[pyarrow]")
    )

    # the operation is computed in double precision before narrowing, so a
    #  tiny product is not lost to float32 rounding
    large = pd.Series([1e38], dtype="float[pyarrow]")
    small_product = large * 1e-50
    assert small_product.iloc[0] == pytest.approx(1e-12, rel=1e-6)
    assert small_product.dtype == "float[pyarrow]"


def test_float_scalar_arithmetic_upcasts_on_overflow():
    # GH#69436 a result that overflows the array's dtype keeps the wider
    #  dtype instead of silently becoming infinite
    large = pd.Series([1e38], dtype="float[pyarrow]")

    result = large * 1e100
    assert np.isfinite(result.iloc[0])
    assert result.dtype == "double[pyarrow]"


def test_decimal_scalar_arithmetic_keeps_dtype():
    # GH#69436
    dtype = pd.ArrowDtype(pa.decimal128(10, 2))
    ser = pd.Series([Decimal("1.00"), Decimal("2.00")], dtype=dtype)

    tm.assert_series_equal(
        ser + 1, pd.Series([Decimal("2.00"), Decimal("3.00")], dtype=dtype)
    )
    tm.assert_series_equal(
        ser * 2, pd.Series([Decimal("2.00"), Decimal("4.00")], dtype=dtype)
    )
    tm.assert_series_equal(
        ser + Decimal("1.50"),
        pd.Series([Decimal("2.50"), Decimal("3.50")], dtype=dtype),
    )

    # division changes the scale, so the wider result is kept
    assert (ser / 2).dtype == pd.ArrowDtype(pa.decimal128(30, 22))

    # a result that does not fit in the array's precision keeps the wider dtype
    narrow_dtype = pd.ArrowDtype(pa.decimal128(2, 0))
    narrow = pd.Series([Decimal("99")], dtype=narrow_dtype)
    overflow = narrow + 1
    assert overflow.tolist() == [Decimal("100")]
    assert overflow.dtype != narrow_dtype
    assert isinstance(overflow.dtype, pd.ArrowDtype)
    assert overflow.dtype.pyarrow_dtype.precision >= 3


def test_scalar_arithmetic_with_missing_values_keeps_dtype():
    # GH#69436
    ser = pd.Series([1, None, 3], dtype="int32[pyarrow]")
    result = ser + 1
    assert result.dtype == "int32[pyarrow]"
    tm.assert_series_equal(result, pd.Series([2, None, 4], dtype="int32[pyarrow]"))
