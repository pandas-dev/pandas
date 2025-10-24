"""
Regression tests for issue #61434: Improved error message for incompatible merge types

Tests that:
1. Merging with polars.DataFrame raises TypeError with helpful message
2. Merging with other incompatible types also gets helpful messages
3. Normal pandas merges still work correctly
"""

import pytest
import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm


class TestMergeIncompatibleTypes:
    """Test merge error messages with incompatible DataFrame types."""
    
    def test_merge_with_polars_dataframe(self):
        """
        Test that merging with polars.DataFrame raises helpful TypeError.
        
        Regression test for issue #61434.
        """
        pytest.importorskip("polars")
        import polars as pl
        
        pdf = DataFrame({
            "key": ["a", "b", "c"],
            "value_x": [1, 2, 3]
        })
        
        plf = pl.DataFrame({
            "key": ["a", "b", "c"],
            "value_y": [10, 20, 30]
        })
        
        with pytest.raises(TypeError, match=".*polars.*pandas.*"):
            pd.merge(pdf, plf, on="key")
    
    def test_merge_polars_to_pandas_conversion(self):
        """
        Test that converting polars to pandas works.
        
        Shows the workaround mentioned in error message.
        """
        pytest.importorskip("polars")
        import polars as pl
        
        pdf = DataFrame({
            "key": ["a", "b", "c"],
            "value_x": [1, 2, 3]
        })
        
        plf = pl.DataFrame({
            "key": ["a", "b", "c"],
            "value_y": [10, 20, 30]
        })
        
        # Convert polars to pandas - this should work
        plf_pd = plf.to_pandas()
        result = pd.merge(pdf, plf_pd, on="key")
        
        expected = DataFrame({
            "key": ["a", "b", "c"],
            "value_x": [1, 2, 3],
            "value_y": [10, 20, 30]
        })
        
        tm.assert_frame_equal(result, expected)
    
    def test_merge_with_dict(self):
        """Test that merging with dict raises TypeError with helpful message."""
        df = DataFrame({"key": ["a", "b"], "value": [1, 2]})
        
        dict_obj = {"key": ["a", "b"], "value": [3, 4]}
        
        with pytest.raises(TypeError, match=".*dict.*"):
            pd.merge(df, dict_obj, on="key")
    
    def test_merge_with_list(self):
        """Test that merging with list raises TypeError with helpful message."""
        df = DataFrame({"key": ["a", "b"], "value": [1, 2]})
        
        list_obj = [["a", 1], ["b", 2]]
        
        msg = "Can only merge Series or DataFrame objects"
        
        with pytest.raises(TypeError, match=msg):
            pd.merge(df, list_obj, on="key")
    
    def test_merge_pandas_baseline(self):
        """
        Test that normal pandas merge still works.
        
        Baseline test to ensure fix doesn't break existing functionality.
        """
        df1 = DataFrame({
            "key": ["a", "b", "c"],
            "value_x": [1, 2, 3]
        })
        
        df2 = DataFrame({
            "key": ["a", "b", "c"],
            "value_y": [10, 20, 30]
        })
        
        result = pd.merge(df1, df2, on="key")
        
        expected = DataFrame({
            "key": ["a", "b", "c"],
            "value_x": [1, 2, 3],
            "value_y": [10, 20, 30]
        })
        
        tm.assert_frame_equal(result, expected)
    
    def test_merge_with_series_name(self):
        """Test that merging with named Series works (baseline)."""
        df = DataFrame({"key": ["a", "b", "c"], "value_x": [1, 2, 3]})
        s = Series([10, 20, 30], name="value_y")
        
        result = pd.merge(df, s, left_index=True, right_index=True)
        
        expected = DataFrame({
            "key": ["a", "b", "c"],
            "value_x": [1, 2, 3],
            "value_y": [10, 20, 30]
        })
        
        tm.assert_frame_equal(result, expected)
    
    def test_merge_with_unnamed_series(self):
        """Test that merging with unnamed Series raises helpful error."""
        df = DataFrame({"key": ["a", "b", "c"], "value": [1, 2, 3]})
        s = Series([10, 20, 30])  # No name
        
        msg = "Cannot merge a Series without a name"
        
        with pytest.raises(ValueError, match=msg):
            pd.merge(df, s, left_index=True, right_index=True)
