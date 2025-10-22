"""
Tests for DataFrame.agg with numeric_only parameter and list of functions.
This tests the fix for GH#49352.
"""

import numpy as np
import pytest
import pandas as pd
from pandas import DataFrame
import pandas._testing as tm


class TestFrameAggNumericOnly:
    """Tests for DataFrame.agg with numeric_only parameter and list of functions."""

    def test_agg_list_numeric_only_mixed_dtypes(self):
        """GH#49352 - Main test case from the issue."""
        df = DataFrame({
            'A': [1, 2, 3, 4, 5],
            'B': [10.5, 20.5, 30.5, 40.5, 50.5],
            'C': ['a', 'b', 'c', 'd', 'e']
        })
        result = df.agg(['min', 'max', 'mean'], numeric_only=True)
        expected = DataFrame({
            'A': [1.0, 5.0, 3.0],
            'B': [10.5, 50.5, 30.5]
        }, index=['min', 'max', 'mean'])
        tm.assert_frame_equal(result, expected)

    def test_agg_list_numeric_only_all_numeric(self):
        """Should work when all columns are numeric."""
        df = DataFrame({
            'A': [1, 2, 3],
            'B': [10, 20, 30]
        })
        result = df.agg(['sum', 'mean'], numeric_only=True)
        expected = DataFrame({
            'A': [6.0, 2.0],
            'B': [60.0, 20.0]
        }, index=['sum', 'mean'])
        tm.assert_frame_equal(result, expected)

    def test_agg_list_numeric_only_no_numeric(self):
        """Should return empty DataFrame when no numeric columns."""
        df = DataFrame({
            'A': ['a', 'b', 'c'],
            'B': ['x', 'y', 'z']
        })
        result = df.agg(['min', 'max'], numeric_only=True)
        expected = DataFrame(index=['min', 'max'])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("funcs,expected_index", [
        (['sum', 'mean'], ['sum', 'mean']),
        ([np.sum, np.mean], ['sum', 'mean']),
        (['sum', np.mean], ['sum', 'mean']),
        ([np.sum, 'mean'], ['sum', 'mean']),
    ])
    def test_agg_list_numeric_only_various_function_types(self, funcs, expected_index):
        """Test with different combinations of string and numpy functions."""
        df = DataFrame({
            'A': [1, 2, 3],
            'B': [10, 20, 30],
            'C': ['a', 'b', 'c']
        })
        result = df.agg(funcs, numeric_only=True)
        expected = DataFrame({
            'A': [6.0, 2.0],
            'B': [60.0, 20.0]
        }, index=expected_index)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("funcs", [
        ['min', 'max'],
        ['sum', 'mean', 'std'],
        ['min', 'max', 'mean', 'median'],
    ])
    def test_agg_list_numeric_only_different_function_counts(self, funcs):
        """Test with different numbers of functions."""
        df = DataFrame({
            'A': [1, 2, 3, 4, 5],
            'B': [10, 20, 30, 40, 50],
            'C': ['a', 'b', 'c', 'd', 'e']
        })
        result = df.agg(funcs, numeric_only=True)

        # Verify structure
        assert isinstance(result, DataFrame)
        assert list(result.columns) == ['A', 'B']
        assert list(result.index) == funcs
        assert result.shape == (len(funcs), 2)

    @pytest.mark.parametrize("data,expected_cols", [
        # Only integers
        ({'A': [1, 2, 3], 'B': [4, 5, 6], 'C': ['x', 'y', 'z']}, ['A', 'B']),
        # Only floats
        ({'A': [1.1, 2.2], 'B': [3.3, 4.4], 'C': ['x', 'y']}, ['A', 'B']),
        # Mix of int and float
        ({'int': [1, 2], 'float': [1.5, 2.5], 'str': ['a', 'b']}, ['int', 'float']),
        # Single numeric column
        ({'num': [1, 2, 3], 'text': ['a', 'b', 'c']}, ['num']),
    ])
    def test_agg_list_numeric_only_various_dtypes(self, data, expected_cols):
        """Test with various numeric dtype combinations."""
        df = DataFrame(data)
        result = df.agg(['sum', 'mean'], numeric_only=True)

        assert isinstance(result, DataFrame)
        assert list(result.columns) == expected_cols
        assert list(result.index) == ['sum', 'mean']

    @pytest.mark.parametrize("numeric_only", [True, False, None])
    def test_agg_list_numeric_only_parameter_values(self, numeric_only):
        """Test with different numeric_only parameter values."""
        df = DataFrame({
            'A': [1, 2, 3],
            'B': [10, 20, 30]
        })

        if numeric_only is None:
            result = df.agg(['sum', 'mean'])
        else:
            result = df.agg(['sum', 'mean'], numeric_only=numeric_only)

        expected = DataFrame({
            'A': [6, 2.0],
            'B': [60, 20.0]
        }, index=['sum', 'mean'])
        tm.assert_frame_equal(result, expected)

    def test_agg_list_numeric_only_false_with_strings(self):
        """Verify numeric_only=False works with min/max on strings."""
        df = DataFrame({
            'A': [1, 2, 3],
            'B': ['a', 'b', 'c']
        })
        result = df.agg(['min', 'max'], numeric_only=False)
        expected = DataFrame({
            'A': [1, 3],
            'B': ['a', 'c']
        }, index=['min', 'max'])
        tm.assert_frame_equal(result, expected)

    def test_agg_list_numeric_only_preserves_column_order(self):
        """Test that column order is preserved."""
        df = DataFrame({
            'Z': [1, 2, 3],
            'A': [10, 20, 30],
            'M': [100, 200, 300],
            'text': ['a', 'b', 'c']
        })
        result = df.agg(['sum', 'mean'], numeric_only=True)

        assert list(result.columns) == ['Z', 'A', 'M']

    @pytest.mark.parametrize("single_func", ['sum', 'mean', 'min', 'max'])
    def test_agg_single_function_still_works(self, single_func):
        """Verify that single function (not a list) still works."""
        df = DataFrame({
            'A': [1, 2, 3],
            'B': [10, 20, 30],
            'C': ['a', 'b', 'c']
        })
        result = df.agg(single_func, numeric_only=True)

        assert isinstance(result, pd.Series)
        assert 'A' in result.index
        assert 'B' in result.index
        assert 'C' not in result.index

    def test_agg_list_numeric_only_with_int_and_float(self):
        """Test that both int and float columns are included."""
        df = DataFrame({
            'int_col': [1, 2, 3],
            'float_col': [1.5, 2.5, 3.5],
            'str_col': ['a', 'b', 'c']
        })
        result = df.agg(['sum', 'mean'], numeric_only=True)
        expected = DataFrame({
            'int_col': [6.0, 2.0],
            'float_col': [7.5, 2.5]
        }, index=['sum', 'mean'])
        tm.assert_frame_equal(result, expected)

    def test_agg_list_numeric_only_single_row(self):
        """Test with single row DataFrame."""
        df = DataFrame({
            'A': [1],
            'B': [10],
            'C': ['x']
        })
        result = df.agg(['sum', 'mean'], numeric_only=True)
        expected = DataFrame({
            'A': [1.0, 1.0],
            'B': [10.0, 10.0]
        }, index=['sum', 'mean'])
        tm.assert_frame_equal(result, expected)

    # ========== NEW TESTS - Additional Edge Cases ==========

    def test_agg_list_numeric_only_with_nans(self):
        """Test DataFrame with NaN values."""
        df = DataFrame({
            'A': [1, np.nan, 3],
            'B': [10, 20, np.nan],
            'C': ['x', 'y', 'z']
        })
        result = df.agg(['sum', 'mean'], numeric_only=True)
        expected = DataFrame({
            'A': [4.0, 2.0],
            'B': [30.0, 15.0]
        }, index=['sum', 'mean'])
        tm.assert_frame_equal(result, expected)

    def test_agg_list_numeric_only_with_datetime(self):
        """Test that datetime columns are excluded with numeric_only=True."""
        df = DataFrame({
            'num': [1, 2, 3],
            'date': pd.date_range('2020-01-01', periods=3),
            'text': ['a', 'b', 'c']
        })
        result = df.agg(['sum', 'mean'], numeric_only=True)
        expected = DataFrame({
            'num': [6.0, 2.0]
        }, index=['sum', 'mean'])
        tm.assert_frame_equal(result, expected)

    def test_agg_list_numeric_only_large_dataframe(self):
        """Test with a larger DataFrame for performance verification."""
        np.random.seed(42)
        df = DataFrame({
            'A': np.random.randint(1, 100, 1000),
            'B': np.random.randn(1000),
            'C': ['text'] * 1000
        })
        result = df.agg(['sum', 'mean', 'std'], numeric_only=True)

        # Just verify structure, not exact values due to randomness
        assert isinstance(result, DataFrame)
        assert list(result.columns) == ['A', 'B']
        assert list(result.index) == ['sum', 'mean', 'std']
        assert result.shape == (3, 2)
