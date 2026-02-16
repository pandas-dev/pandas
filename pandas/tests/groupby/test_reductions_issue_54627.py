"""
Test for Issue #54627: groupby.var() does not return arrow types

https://github.com/pandas-dev/pandas/issues/54627
"""

import decimal

import numpy as np
import pytest

import pandas as pd
import pyarrow as pa
from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


@pytest.mark.filterwarnings("ignore::FutureWarning")
class TestGroupByVarArrowDtype:
    """Test groupby.var() preserves Arrow dtype."""

    def test_groupby_var_arrow_dtype(self):
        """Test that groupby.var() returns Arrow dtype when input has Arrow dtype."""
        # Create DataFrame with Arrow types
        df = DataFrame({
            "A": Series([True, True], dtype="bool[pyarrow]"),
            "B": Series(
                [decimal.Decimal(123), decimal.Decimal(12)],
                dtype=pd.ArrowDtype(pa.decimal128(6, 3))
            ),
        })

        # Groupby and compute variance
        result = df.groupby("A").var()

        # Expected: Column B should have float64[pyarrow] dtype
        # Bug: Column B has float64 dtype
        assert "pyarrow" in str(result["B"].dtype), (
            f"Expected Arrow dtype, got {result['B'].dtype}"
        )

    def test_groupby_var_float_pyarrow(self):
        """Test groupby.var() with float64[pyarrow] dtype."""
        df = DataFrame({
            "key": Series([1, 1, 2], dtype="int64[pyarrow]"),
            "value": Series([1.0, 2.0, 3.0], dtype="float64[pyarrow]"),
        })

        result = df.groupby("key").var()

        # The result should preserve Arrow dtype
        assert "pyarrow" in str(result["value"].dtype), (
            f"Expected Arrow dtype, got {result['value'].dtype}"
        )

    def test_groupby_var_decimal_pyarrow(self):
        """Test groupby.var() with decimal[pyarrow] dtype."""
        df = DataFrame({
            "key": Series([1, 1, 2], dtype="int64[pyarrow]"),
            "value": Series(
                [decimal.Decimal(1.5), decimal.Decimal(2.5), decimal.Decimal(3.5)],
                dtype=pd.ArrowDtype(pa.decimal128(10, 2))
            ),
        })

        result = df.groupby("key").var()

        # The result should preserve Arrow dtype
        assert "pyarrow" in str(result["value"].dtype), (
            f"Expected Arrow dtype, got {result['value'].dtype}"
        )


if __name__ == "__main__":
    # Run the tests
    test = TestGroupByVarArrowDtype()
    try:
        test.test_groupby_var_arrow_dtype()
        print("✅ test_groupby_var_arrow_dtype passed")
    except Exception as e:
        print(f"❌ test_groupby_var_arrow_dtype failed: {e}")
    
    try:
        test.test_groupby_var_float_pyarrow()
        print("✅ test_groupby_var_float_pyarrow passed")
    except Exception as e:
        print(f"❌ test_groupby_var_float_pyarrow failed: {e}")
    
    try:
        test.test_groupby_var_decimal_pyarrow()
        print("✅ test_groupby_var_decimal_pyarrow passed")
    except Exception as e:
        print(f"❌ test_groupby_var_decimal_pyarrow failed: {e}")
