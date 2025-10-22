"""
Tests for the CoW DataFrame.replace fix for np.nan dictionary replacement bug.

Regression tests for GH#62787: Enabling Copy on Write with DataFrame.replace
Raises Exception with np.nan as replacement value.
"""

import numpy as np

import pandas as pd
from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


class TestReplaceCoWFix:
    """Tests for the CoW replace fix for GH#62787."""

    def test_replace_dict_with_nan_cow_enabled(self):
        """Test that dictionary replacement with np.nan works with CoW enabled."""
        # GH#62787
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame(
                {
                    "A": [1, 2],
                    "B": ["b", "i like pandas"],
                }
            )
            df["Name"] = "I Have a Name"
            df["Name2"] = "i like pandas"

            # This should not raise an error
            replace_mappings = {
                pd.NA: None,
                pd.NaT: None,
                np.nan: None,  # This was causing the bug
            }
            result = df.replace(replace_mappings)

            # Should return a DataFrame without errors
            assert isinstance(result, DataFrame)
            # The original data should remain unchanged since we are
            # replacing values that don't exist
            tm.assert_frame_equal(result, df)

    def test_replace_dict_with_various_na_values_cow(self):
        """Test dictionary replacement with various NA values under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            # Create DataFrame with actual NA values to replace
            df = DataFrame(
                {
                    "A": [1, np.nan, 3],
                    "B": [pd.NA, "test", pd.NaT],
                    "C": ["x", "y", "z"],
                }
            )

            replace_mappings = {
                pd.NA: "replaced_NA",
                pd.NaT: "replaced_NaT",
                np.nan: "replaced_nan",
            }

            result = df.replace(replace_mappings)

            expected = DataFrame(
                {
                    "A": [1, "replaced_nan", 3],
                    "B": ["replaced_NA", "test", "replaced_NaT"],
                    "C": ["x", "y", "z"],
                }
            )

            tm.assert_frame_equal(result, expected)

    def test_replace_dict_nan_series_cow(self):
        """Test Series replace with np.nan in dictionary under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            s = Series([1, np.nan, 3, np.nan])

            replace_mappings = {np.nan: "missing", 1: "one"}

            result = s.replace(replace_mappings)
            expected = Series(["one", "missing", 3, "missing"])

            tm.assert_series_equal(result, expected)

    def test_replace_dict_empty_cow(self):
        """Test empty dictionary replacement under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame({"A": [1, 2], "B": ["a", "b"]})

            # Empty replacement dict should work
            result = df.replace({})
            tm.assert_frame_equal(result, df)

    def test_replace_dict_with_nan_inplace_cow(self):
        """Test inplace dictionary replacement with np.nan under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame({"A": [1, np.nan, 3], "B": ["x", "y", "z"]})

            replace_mappings = {np.nan: -999}
            result = df.replace(replace_mappings, inplace=True)

            # inplace=True should return None
            assert result is None

            expected = DataFrame({"A": [1, -999, 3], "B": ["x", "y", "z"]})

            tm.assert_frame_equal(df, expected)

    def test_replace_mixed_types_with_nan_cow(self):
        """Test mixed type replacement including np.nan under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame(
                {
                    "int_col": [1, 2, 3],
                    "float_col": [1.1, np.nan, 3.3],
                    "str_col": ["a", "b", "c"],
                    "mixed_col": [1, "text", np.nan],
                }
            )

            replace_mappings = {np.nan: "MISSING", 1: "ONE", "a": "LETTER_A"}

            result = df.replace(replace_mappings)

            expected = DataFrame(
                {
                    "int_col": ["ONE", 2, 3],
                    "float_col": [1.1, "MISSING", 3.3],
                    "str_col": ["LETTER_A", "b", "c"],
                    "mixed_col": ["ONE", "text", "MISSING"],
                }
            )

            tm.assert_frame_equal(result, expected)

    def test_replace_cow_vs_no_cow_consistency(self):
        """Test that CoW and non-CoW modes give same results."""
        df_data = {"A": [1, np.nan, 3], "B": ["x", "y", "z"]}
        replace_mappings = {np.nan: "REPLACED"}

        # Test with CoW enabled
        with pd.option_context("mode.copy_on_write", True):
            df_cow = DataFrame(df_data)
            result_cow = df_cow.replace(replace_mappings)

        # Test with CoW disabled
        with pd.option_context("mode.copy_on_write", False):
            df_no_cow = DataFrame(df_data)
            result_no_cow = df_no_cow.replace(replace_mappings)

        # Results should be identical
        tm.assert_frame_equal(result_cow, result_no_cow)

    def test_replace_complex_nested_dict_with_nan_cow(self):
        """Test complex nested dictionary replacements with np.nan under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame(
                {"A": [1, np.nan, 3], "B": [4, 5, np.nan], "C": ["x", "y", "z"]}
            )

            # Column-specific replacements
            replace_mappings = {"A": {np.nan: -1, 1: 100}, "B": {np.nan: -2, 4: 400}}

            result = df.replace(replace_mappings)

            expected = DataFrame(
                {"A": [100, -1, 3], "B": [400, 5, -2], "C": ["x", "y", "z"]}
            )

            tm.assert_frame_equal(result, expected)

    def test_replace_regex_with_nan_cow(self):
        """Test regex replacement combined with np.nan under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame(
                {"text": ["hello world", "foo bar", "test"], "nums": [1, np.nan, 3]}
            )

            # First do dictionary replacement, then regex
            replace_mappings = {np.nan: "MISSING"}
            result = df.replace(replace_mappings)

            # Then regex replacement
            result = result.replace(r"hello.*", "GREETING", regex=True)

            expected = DataFrame(
                {"text": ["GREETING", "foo bar", "test"], "nums": [1, "MISSING", 3]}
            )

            tm.assert_frame_equal(result, expected)

    def test_replace_multiple_nan_types_cow(self):
        """Test replacement of different NaN types in same operation."""
        with pd.option_context("mode.copy_on_write", True):
            # Create DataFrame with different types of missing values
            df = DataFrame(
                {
                    "float_nan": [1.0, np.nan, 3.0],
                    "pd_na": ["a", pd.NA, "c"],
                    "pd_nat": [
                        pd.Timestamp("2020-01-01"),
                        pd.NaT,
                        pd.Timestamp("2020-01-03"),
                    ],
                }
            )

            replace_mappings = {
                np.nan: "float_missing",
                pd.NA: "string_missing",
                pd.NaT: pd.Timestamp("1900-01-01"),
            }

            result = df.replace(replace_mappings)

            expected = DataFrame(
                {
                    "float_nan": [1.0, "float_missing", 3.0],
                    "pd_na": ["a", "string_missing", "c"],
                    "pd_nat": [
                        pd.Timestamp("2020-01-01"),
                        pd.Timestamp("1900-01-01"),
                        pd.Timestamp("2020-01-03"),
                    ],
                }
            )

            tm.assert_frame_equal(result, expected)


class TestReplaceCoWEdgeCases:
    """Edge case tests for the CoW replace fix."""

    def test_replace_nan_with_none_cow(self):
        """Test specific case from bug report: np.nan -> None."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame(
                {
                    "A": [1, 2],
                    "B": ["b", "i like pandas"],
                }
            )
            df["Name"] = "I Have a Name"
            df["Name2"] = "i like pandas"

            # This exact case from the bug report
            replace_mappings = {pd.NA: None, pd.NaT: None, np.nan: None}

            # Should not raise ValueError about weakref
            result = df.replace(replace_mappings)
            assert isinstance(result, DataFrame)

    def test_replace_large_dict_with_nan_cow(self):
        """Test large replacement dictionary including np.nan."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame({"A": range(100), "B": [np.nan] * 100})

            # Large replacement dict to stress test weak reference handling
            replace_dict = {i: f"num_{i}" for i in range(50)}
            replace_dict[np.nan] = "missing"

            result = df.replace(replace_dict)

            # Verify it works without error
            assert len(result) == 100
            assert all(result["B"] == "missing")

    def test_replace_chained_operations_cow(self):
        """Test chained replace operations with np.nan under CoW."""
        with pd.option_context("mode.copy_on_write", True):
            df = DataFrame({"A": [1, np.nan, 3, np.nan], "B": ["a", "b", "c", "d"]})

            # Chain multiple replace operations
            result = (
                df.replace({np.nan: -1}).replace({1: "ONE"}).replace({"a": "LETTER_A"})
            )

            expected = DataFrame(
                {"A": ["ONE", -1, 3, -1], "B": ["LETTER_A", "b", "c", "d"]}
            )

            tm.assert_frame_equal(result, expected)
