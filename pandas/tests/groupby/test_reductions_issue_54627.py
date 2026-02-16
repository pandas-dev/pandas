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


pytestmark = pytest.mark.filterwarnings("ignore::FutureWarning")


class TestGroupByVarArrowDtype:
    """Test groupby.var() preserves Arrow dtype."""

    def test_groupby_var_arrow_dtype_series(self):
        """Test that groupby.var() returns Arrow dtype for Series with Arrow dtype."""
        # Create Series with Arrow type
        s = Series([1.0, 2.0, 3.0, 4.0, 5.0], dtype="float64[pyarrow]")
        key = Series([1, 1, 2, 2, 2], dtype="int64[pyarrow]"])

        result = s.groupby(key).var()

        # Expected: result should have float64[pyarrow] dtype
        # Bug: result has float64 dtype
        assert "pyarrow" in str(result.dtype), (
            f"Expected Arrow dtype, got {result.dtype}"
        )

    def test_groupby_var_arrow_dtype_dataframe(self):
        """Test that groupby.var() returns Arrow dtype for DataFrame with Arrow dtype."""
        # Create DataFrame with Arrow types
        df = DataFrame({
            "key": Series([1, 1, 2, 2], dtype="int64[pyarrow]"),
            "value": Series([1.0, 2.0, 3.0, 4.0], dtype="float64[pyarrow]"),
        })

        result = df.groupby("key").var()

        # Expected: Column 'value' should have float64[pyarrow] dtype
        assert "pyarrow" in str(result["value"].dtype), (
            f"Expected Arrow dtype for value column, got {result['value'].dtype}"
        )

    def test_groupby_var_arrow_dtype_decimal(self):
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
        # Note: decimal128 converts to float64 in variance calculation
        # but should still be pyarrow float64
        assert "pyarrow" in str(result["value"].dtype), (
            f"Expected Arrow dtype, got {result['value'].dtype}"
        )

    def test_groupby_var_without_arrow_returns_numpy(self):
        """Test that groupby.var() returns numpy dtype when input is numpy."""
        # Create DataFrame with regular numpy types
        df = DataFrame({
            "key": [1, 1, 2, 2],
            "value": [1.0, 2.0, 3.0, 4.0],
        })

        result = df.groupby("key").var()

        # Should be regular numpy float64
        assert result["value"].dtype == np.float64
        assert "pyarrow" not in str(result["value"].dtype)

    def test_groupby_var_mixed_types(self):
        """Test groupby.var() with mixed Arrow and numpy columns."""
        df = DataFrame({
            "key": Series([1, 1, 2, 2], dtype="int64[pyarrow]"),
            "arrow_col": Series([1.0, 2.0, 3.0, 4.0], dtype="float64[pyarrow]"),
            "numpy_col": [1.0, 2.0, 3.0, 4.0],
        })

        result = df.groupby("key").var()

        # Arrow column should preserve Arrow dtype
        assert "pyarrow" in str(result["arrow_col"].dtype), (
            f"Expected Arrow dtype for arrow_col, got {result['arrow_col'].dtype}"
        )
        # Numpy column should remain numpy
        assert result["numpy_col"].dtype == np.float64


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
