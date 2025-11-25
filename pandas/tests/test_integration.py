"""
Integration tests for pandas modules.

These tests verify interactions between multiple modules/components:
- pandas.core.series (Series construction)
- pandas.core.frame (DataFrame construction)
- pandas.core.dtypes (dtype handling)

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


