import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    check,
)
import pandas._testing as tm


class TestCheck:
    def test_basic_functionality(self):
        """Test basic functionality of pd.check()."""
        df = DataFrame({
            'A': [1, 2, None, 4],
            'B': ['x', 'y', 'x', None],
            'C': [1.0, 2.0, 3.0, 4.0]
        })
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [3, 2, 4],
            'non_null': [3, 3, 4],
            'missing': [1, 1, 0],
            'missing_pct': [25.00, 25.00, 0.00]
        }, index=['A', 'B', 'C'])
        
        tm.assert_frame_equal(result, expected)

    def test_empty_dataframe(self):
        """Test check() with empty DataFrame."""
        df = DataFrame()
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [],
            'non_null': [],
            'missing': [],
            'missing_pct': []
        }).astype('int64')
        expected['missing_pct'] = expected['missing_pct'].astype('float64')
        
        tm.assert_frame_equal(result, expected)

    def test_all_null_column(self):
        """Test check() with a column that is all null."""
        df = DataFrame({
            'A': [1, 2, 3],
            'B': [None, None, None],
            'C': ['x', 'y', 'z']
        })
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [3, 0, 3],
            'non_null': [3, 0, 3],
            'missing': [0, 3, 0],
            'missing_pct': [0.00, 100.00, 0.00]
        }, index=['A', 'B', 'C'])
        
        tm.assert_frame_equal(result, expected)

    def test_no_missing_values(self):
        """Test check() with DataFrame that has no missing values."""
        df = DataFrame({
            'A': [1, 2, 3, 4],
            'B': ['w', 'x', 'y', 'z'],
            'C': [1.1, 2.2, 3.3, 4.4]
        })
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [4, 4, 4],
            'non_null': [4, 4, 4],
            'missing': [0, 0, 0],
            'missing_pct': [0.00, 0.00, 0.00]
        }, index=['A', 'B', 'C'])
        
        tm.assert_frame_equal(result, expected)

    def test_round_digits_parameter(self):
        """Test check() with different round_digits parameter."""
        df = DataFrame({
            'A': [1, None, None],  # 2/3 = 66.666... % missing
            'B': [1, 2, 3]
        })
        
        # Test with default round_digits=2
        result_default = check(df)
        expected_default = DataFrame({
            'unique': [1, 3],
            'non_null': [1, 3],
            'missing': [2, 0],
            'missing_pct': [66.67, 0.00]
        }, index=['A', 'B'])
        tm.assert_frame_equal(result_default, expected_default)
        
        # Test with round_digits=0
        result_zero = check(df, round_digits=0)
        expected_zero = DataFrame({
            'unique': [1, 3],
            'non_null': [1, 3],
            'missing': [2, 0],
            'missing_pct': [67.0, 0.0]
        }, index=['A', 'B'])
        tm.assert_frame_equal(result_zero, expected_zero)
        
        # Test with round_digits=4
        result_four = check(df, round_digits=4)
        expected_four = DataFrame({
            'unique': [1, 3],
            'non_null': [1, 3],
            'missing': [2, 0],
            'missing_pct': [66.6667, 0.0000]
        }, index=['A', 'B'])
        tm.assert_frame_equal(result_four, expected_four)

    def test_various_dtypes(self):
        """Test check() with various data types."""
        df = DataFrame({
            'int_col': [1, 2, None],
            'float_col': [1.1, None, 3.3],
            'str_col': ['a', 'b', None],
            'bool_col': [True, False, None],
            'datetime_col': pd.to_datetime(['2020-01-01', '2020-01-02', None])
        })
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [2, 2, 2, 2, 2],
            'non_null': [2, 2, 2, 2, 2],
            'missing': [1, 1, 1, 1, 1],
            'missing_pct': [33.33, 33.33, 33.33, 33.33, 33.33]
        }, index=['int_col', 'float_col', 'str_col', 'bool_col', 'datetime_col'])
        
        tm.assert_frame_equal(result, expected)

    def test_duplicate_values(self):
        """Test check() with columns containing duplicate values."""
        df = DataFrame({
            'A': [1, 1, 2, 2, 2],
            'B': ['x', 'x', 'x', 'y', 'y'],
            'C': [1, 1, 1, 1, 1]  # All same value
        })
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [2, 2, 1],
            'non_null': [5, 5, 5],
            'missing': [0, 0, 0],
            'missing_pct': [0.00, 0.00, 0.00]
        }, index=['A', 'B', 'C'])
        
        tm.assert_frame_equal(result, expected)

    def test_single_row_dataframe(self):
        """Test check() with single row DataFrame."""
        df = DataFrame({
            'A': [1],
            'B': [None],
            'C': ['test']
        })
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [1, 0, 1],
            'non_null': [1, 0, 1],
            'missing': [0, 1, 0],
            'missing_pct': [0.00, 100.00, 0.00]
        }, index=['A', 'B', 'C'])
        
        tm.assert_frame_equal(result, expected)

    def test_single_column_dataframe(self):
        """Test check() with single column DataFrame."""
        df = DataFrame({
            'A': [1, 2, None, 4]
        })
        
        result = check(df)
        
        expected = DataFrame({
            'unique': [3],
            'non_null': [3],
            'missing': [1],
            'missing_pct': [25.00]
        }, index=['A'])
        
        tm.assert_frame_equal(result, expected)

    def test_non_dataframe_raises_error(self):
        """Test that check() raises appropriate error for non-DataFrame input."""
        with pytest.raises(AttributeError):
            check("not a dataframe")
        
        with pytest.raises(AttributeError):
            check([1, 2, 3])

    def test_return_type(self):
        """Test that check() returns a DataFrame."""
        df = DataFrame({'A': [1, 2, 3]})
        result = check(df)
        assert isinstance(result, DataFrame)

    def test_column_order_preserved(self):
        """Test that the order of columns is preserved in the result."""
        df = DataFrame({
            'Z': [1, 2, 3],
            'A': [4, 5, 6],
            'M': [7, 8, 9]
        })
        
        result = check(df)
        
        expected_index = ['Z', 'A', 'M']
        tm.assert_index_equal(result.index, pd.Index(expected_index))