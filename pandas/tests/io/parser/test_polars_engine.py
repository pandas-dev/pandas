"""
Test polars engine for read_csv
"""

import pytest

from pandas.compat import HAS_POLARS
import pandas as pd
import numpy as np
import pandas._testing as tm


class TestPolarsEngine:
    """Tests for the polars engine."""

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_basic(self):
        """Test basic functionality with polars engine."""
        csv_data = "a,b,c\n1,2,3\n4,5,6"
        
        result = pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars")
        expected = pd.read_csv(pd.io.common.StringIO(csv_data), engine="python")
        
        tm.assert_frame_equal(result, expected)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_with_header(self):
        """Test polars engine with custom header."""
        csv_data = "col1,col2,col3\n1,2,3\n4,5,6"
        
        result = pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", header=0)
        expected = pd.read_csv(pd.io.common.StringIO(csv_data), engine="python", header=0)
        
        tm.assert_frame_equal(result, expected)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_with_names(self):
        """Test polars engine with custom column names."""
        csv_data = "1,2,3\n4,5,6"
        names = ["x", "y", "z"]
        
        result = pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", names=names, header=None)
        expected = pd.read_csv(pd.io.common.StringIO(csv_data), engine="python", names=names, header=None)
        
        tm.assert_frame_equal(result, expected)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_with_usecols_string(self):
        """Test polars engine with usecols as strings."""
        csv_data = "a,b,c\n1,2,3\n4,5,6"
        
        result = pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", usecols=["a", "c"])
        expected = pd.read_csv(pd.io.common.StringIO(csv_data), engine="python", usecols=["a", "c"])
        
        tm.assert_frame_equal(result, expected)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_unsupported_chunksize(self):
        """Test that polars engine raises error for chunksize."""
        csv_data = "a,b,c\n1,2,3\n4,5,6"
        
        with pytest.raises(ValueError, match="not supported with the 'polars' engine"):
            pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", chunksize=1)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_unsupported_iterator(self):
        """Test that polars engine raises error for iterator."""
        csv_data = "a,b,c\n1,2,3\n4,5,6"
        
        with pytest.raises(ValueError, match="not supported with the 'polars' engine"):
            pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", iterator=True)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_with_nrows(self):
        """Test polars engine with nrows parameter."""
        csv_data = "a,b,c\n1,2,3\n4,5,6\n7,8,9"
        
        result = pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", nrows=2)
        expected = pd.read_csv(pd.io.common.StringIO(csv_data), engine="python", nrows=2)
        
        tm.assert_frame_equal(result, expected)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_string_na_values(self):
        """Test polars engine with na_values."""
        csv_data = "a,b,c\n1,NULL,3\n4,5,NULL"
        
        result = pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", na_values=["NULL"])
        expected = pd.read_csv(pd.io.common.StringIO(csv_data), engine="python", na_values=["NULL"])
        
        tm.assert_frame_equal(result, expected)

    @pytest.mark.skipif(not HAS_POLARS, reason="polars not installed")
    def test_polars_engine_dict_na_values_error(self):
        """Test that polars engine raises error for dict na_values."""
        csv_data = "a,b,c\n1,2,3\n4,5,6"
        
        with pytest.raises(ValueError, match="doesn't support passing a dict for na_values"):
            pd.read_csv(pd.io.common.StringIO(csv_data), engine="polars", na_values={"a": ["1"]})
