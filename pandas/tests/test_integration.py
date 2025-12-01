"""
Integration tests for pandas modules.

These tests verify interactions between multiple modules/components:
- pandas.core.series (Series construction)
- pandas.core.frame (DataFrame construction)
- pandas.core.dtypes (dtype handling)
- pandas.core.internals (internal data management)
- pandas.util._validators (validation utilities)
- pandas.core.missing (missing data handling)
"""
import numpy as np
import pytest

import pandas as pd
from pandas import Series, DataFrame, Index
from pandas.core.missing import clean_fill_method
from pandas._libs import lib
from pandas.util._validators import (
    validate_args_and_kwargs,
    validate_fillna_kwargs,
    check_dtype_backend,
    validate_percentile,
)


class TestSandeepIntegration:
    """Integration tests by Sandeep Ramavath covering Series-DataFrame-dtype interactions."""
    
    def test_series_to_dataframe_dtype_preservation(self):
        """Test Series.to_frame() preserves dtype through internals conversion.
        
        This exercises interaction between:
        - pandas.core.series.Series.to_frame()
        - pandas.core.internals (manager conversion)
        - pandas.core.frame.DataFrame
        - pandas.core.dtypes (dtype preservation)
        """
        # Create Series with specific dtype
        s = Series([1, 2, 3], name="test_col", dtype="int32")
        
        # Convert to DataFrame - should preserve dtype through internal conversion
        df = s.to_frame()
        
        assert isinstance(df, DataFrame)
        assert df.columns[0] == "test_col"
        assert df["test_col"].dtype == np.dtype("int32")
        assert len(df) == 3
        assert (df["test_col"] == s).all()

    def test_dataframe_from_dict_mixed_series_dtypes(self):
        """Test DataFrame construction from dict with mixed Series dtypes.
        
        This exercises interaction between:
        - pandas.core.frame.DataFrame.__init__
        - pandas.core.internals.construction.dict_to_mgr
        - pandas.core.series.Series (multiple instances with different dtypes)
        - pandas.core.dtypes (type coercion and preservation)
        """
        # Create Series with different dtypes
        s1 = Series([1, 2, 3], dtype="int64")
        s2 = Series([1.0, 2.0, 3.0], dtype="float32")
        s3 = Series(["a", "b", "c"], dtype="object")
        
        # Build DataFrame from dict of Series
        df = DataFrame({"col1": s1, "col2": s2, "col3": s3})
        
        # Verify each column maintains its original dtype
        assert df["col1"].dtype == np.dtype("int64")
        assert df["col2"].dtype == np.dtype("float32")
        assert df["col3"].dtype == np.dtype("object")
        assert len(df) == 3


class TestNithikeshIntegration:
    """Integration tests by Nithikesh Bobbili covering validation-missing data interactions."""
    
    def test_validate_fillna_with_clean_method(self):
        """Test validate_fillna_kwargs delegates to clean_fill_method.
        
        This exercises interaction between:
        - pandas.util._validators.validate_fillna_kwargs
        - pandas.core.missing.clean_fill_method
        - method normalization and validation
        """
        # Test method normalization through validate_fillna_kwargs
        value, method = validate_fillna_kwargs(None, "pad")
        assert value is None
        assert method == clean_fill_method("pad")
        
        # Test alternate method names
        value, method = validate_fillna_kwargs(None, "ffill")
        assert method == clean_fill_method("ffill")
        
        # Both None should raise
        with pytest.raises(ValueError, match="Must specify a fill"):
            validate_fillna_kwargs(None, None)
    
    def test_series_fillna_integration(self):
        """Test Series.fillna() and ffill() use validation and missing data modules.
        
        This exercises interaction between:
        - pandas.core.series.Series.fillna() / ffill()
        - pandas.util._validators.validate_fillna_kwargs (internally)
        - pandas.core.missing (fill methods)
        - pandas.core.internals (data modification)
        """
        # Create Series with missing values
        s = Series([1.0, np.nan, 3.0, np.nan, 5.0])
        
        # ffill uses forward fill method - interacts with missing data module
        result = s.ffill()
        expected = Series([1.0, 1.0, 3.0, 3.0, 5.0])
        pd.testing.assert_series_equal(result, expected)
        
        # fillna with value - validation ensures value is acceptable
        result = s.fillna(value=0.0)
        expected = Series([1.0, 0.0, 3.0, 0.0, 5.0])
        pd.testing.assert_series_equal(result, expected)

class TestMallikarjunaIntegration:
    """Integration tests by Mallikarjuna covering dtype_backend-libs interactions."""
    
    def test_check_dtype_backend_with_lib_sentinel(self):
        """Test check_dtype_backend with lib.no_default sentinel.
        
        This exercises interaction between:
        - pandas.util._validators.check_dtype_backend
        - pandas._libs.lib.no_default (sentinel value)
        - validation of backend options
        """
        # Should accept sentinel without exception
        check_dtype_backend(lib.no_default)
        
        # Should accept valid backends
        check_dtype_backend("numpy_nullable")
        check_dtype_backend("pyarrow")
        
        # Should reject unknown backend
        with pytest.raises(ValueError, match="dtype_backend .* is invalid"):
            check_dtype_backend("not_a_backend")
    
    def test_percentile_validation_with_numpy_arrays(self):
        """Test validate_percentile with numpy array interaction.
        
        This exercises interaction between:
        - pandas.util._validators.validate_percentile
        - numpy array conversion and validation
        - pandas statistical methods that use percentiles
        """
        # Single percentile as float
        result = validate_percentile(0.5)
        assert isinstance(result, np.ndarray)
        assert result == 0.5
        
        # Multiple percentiles as list
        result = validate_percentile([0.25, 0.5, 0.75])
        expected = np.array([0.25, 0.5, 0.75])
        np.testing.assert_array_equal(result, expected)
        
        # Invalid percentile should raise
        with pytest.raises(ValueError, match="percentiles should all be"):
            validate_percentile(1.5)
        
        with pytest.raises(ValueError, match="percentiles should all be"):
            validate_percentile([0.25, 1.5, 0.75])


