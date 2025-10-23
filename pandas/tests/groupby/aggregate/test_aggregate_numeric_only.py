"""
Tests for GroupBy.agg with numeric_only parameter and list of functions.
This tests the GroupBy part of the fix for GH#49352.
"""

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame
import pandas._testing as tm


class TestGroupByAggNumericOnly:
    """Tests for GroupBy.agg with numeric_only parameter and list of functions."""

    def test_groupby_agg_list_numeric_only_basic(self):
        """GH#49352 - Basic GroupBy aggregation with mixed dtypes."""
        df = DataFrame(
            {
                "key": ["A", "B", "A", "B", "A"],
                "num1": [1, 2, 3, 4, 5],
                "num2": [10, 20, 30, 40, 50],
                "text": ["a", "b", "c", "d", "e"],
            }
        )
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)
        expected = DataFrame(
            [[9, 3.0, 90, 30.0], [6, 3.0, 60, 30.0]],
            index=pd.Index(["A", "B"], name="key"),
            columns=pd.MultiIndex.from_product([["num1", "num2"], ["sum", "mean"]]),
        )
        tm.assert_frame_equal(result, expected)

    def test_groupby_agg_list_numeric_only_all_numeric(self):
        """GroupBy with all numeric columns."""
        df = DataFrame(
            {
                "key": ["X", "Y", "X", "Y"],
                "val1": [1, 2, 3, 4],
                "val2": [10, 20, 30, 40],
            }
        )
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)
        expected = DataFrame(
            [[4, 2.0, 40, 20.0], [6, 3.0, 60, 30.0]],
            index=pd.Index(["X", "Y"], name="key"),
            columns=pd.MultiIndex.from_product([["val1", "val2"], ["sum", "mean"]]),
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "funcs,expected_func_names",
        [
            (["sum", "mean"], ["sum", "mean"]),
            ([np.sum, np.mean], ["sum", "mean"]),
            (["sum", np.mean], ["sum", "mean"]),
            (["min", "max", "mean"], ["min", "max", "mean"]),
        ],
    )
    def test_groupby_agg_list_numeric_only_various_functions(
        self, funcs, expected_func_names
    ):
        """Test GroupBy with different function combinations."""
        df = DataFrame(
            {
                "key": ["A", "A", "B", "B"],
                "val": [1, 2, 3, 4],
                "text": ["a", "b", "c", "d"],
            }
        )
        result = df.groupby("key").agg(funcs, numeric_only=True)

        assert isinstance(result, DataFrame)
        assert result.columns.levels[0].tolist() == ["val"]
        assert result.columns.levels[1].tolist() == expected_func_names
        assert result.shape == (2, len(funcs))

    @pytest.mark.parametrize(
        "group_cols",
        [
            ["key1"],
            ["key1", "key2"],
        ],
    )
    def test_groupby_agg_list_numeric_only_multiple_groups(self, group_cols):
        """Test GroupBy with single and multiple grouping columns."""
        df = DataFrame(
            {
                "key1": ["A", "A", "B", "B"],
                "key2": ["X", "Y", "X", "Y"],
                "val": [1, 2, 3, 4],
                "text": ["a", "b", "c", "d"],
            }
        )
        result = df.groupby(group_cols).agg(["sum", "mean"], numeric_only=True)

        assert isinstance(result, DataFrame)
        assert result.columns.levels[0].tolist() == ["val"]
        assert result.columns.levels[1].tolist() == ["sum", "mean"]

    @pytest.mark.parametrize(
        "data,expected_cols",
        [
            # Int and float
            (
                {
                    "key": ["A", "A", "B"],
                    "int": [1, 2, 3],
                    "float": [1.5, 2.5, 3.5],
                    "str": ["x", "y", "z"],
                },
                ["int", "float"],
            ),
            # Only int
            (
                {"key": ["A", "B", "A"], "num": [1, 2, 3], "text": ["a", "b", "c"]},
                ["num"],
            ),
            # Multiple numeric
            (
                {
                    "key": ["A", "A"],
                    "n1": [1, 2],
                    "n2": [3, 4],
                    "n3": [5, 6],
                    "str": ["x", "y"],
                },
                ["n1", "n2", "n3"],
            ),
        ],
    )
    def test_groupby_agg_list_numeric_only_various_dtypes(self, data, expected_cols):
        """Test GroupBy with various numeric column combinations."""
        df = DataFrame(data)
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)

        assert isinstance(result, DataFrame)
        assert result.columns.levels[0].tolist() == expected_cols
        assert result.columns.levels[1].tolist() == ["sum", "mean"]

    def test_groupby_agg_list_numeric_only_mixed_int_float(self):
        """Test that both int and float columns are included in GroupBy."""
        df = DataFrame(
            {
                "key": ["A", "A", "B", "B"],
                "int_col": [1, 2, 3, 4],
                "float_col": [1.5, 2.5, 3.5, 4.5],
                "str_col": ["a", "b", "c", "d"],
            }
        )
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)
        expected = DataFrame(
            [[3, 1.5, 4.0, 2.0], [7, 3.5, 8.0, 4.0]],
            index=pd.Index(["A", "B"], name="key"),
            columns=pd.MultiIndex.from_product(
                [["int_col", "float_col"], ["sum", "mean"]]
            ),
        )
        tm.assert_frame_equal(result, expected)

    def test_groupby_agg_list_numeric_only_preserves_column_order(self):
        """Test that GroupBy preserves column order."""
        df = DataFrame(
            {
                "key": ["A", "A", "B", "B"],
                "z_col": [1, 2, 3, 4],
                "a_col": [10, 20, 30, 40],
                "m_col": [100, 200, 300, 400],
                "text": ["a", "b", "c", "d"],
            }
        )
        result = df.groupby("key").agg(["sum"], numeric_only=True)

        assert result.columns.levels[0].tolist() == ["z_col", "a_col", "m_col"]

    @pytest.mark.parametrize("numeric_only", [True, False])
    def test_groupby_agg_list_numeric_only_parameter_values(self, numeric_only):
        """Test GroupBy with numeric_only=True and False."""
        df = DataFrame({"key": ["A", "A", "B", "B"], "val": [1, 2, 3, 4]})
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=numeric_only)
        expected = DataFrame(
            [[3, 1.5], [7, 3.5]],
            index=pd.Index(["A", "B"], name="key"),
            columns=pd.MultiIndex.from_product([["val"], ["sum", "mean"]]),
        )
        tm.assert_frame_equal(result, expected)

    def test_groupby_agg_list_numeric_only_single_group(self):
        """Test GroupBy with a single group."""
        df = DataFrame(
            {
                "key": ["A", "A", "A"],
                "val1": [1, 2, 3],
                "val2": [10, 20, 30],
                "text": ["x", "y", "z"],
            }
        )
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)
        expected = DataFrame(
            [[6, 2.0, 60, 20.0]],
            index=pd.Index(["A"], name="key"),
            columns=pd.MultiIndex.from_product([["val1", "val2"], ["sum", "mean"]]),
        )
        tm.assert_frame_equal(result, expected)

    def test_groupby_agg_list_numeric_only_many_groups(self):
        """Test GroupBy with many groups."""
        df = DataFrame(
            {
                "key": ["A", "B", "C", "D", "E"],
                "val": [1, 2, 3, 4, 5],
                "text": ["a", "b", "c", "d", "e"],
            }
        )
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)

        assert isinstance(result, DataFrame)
        assert len(result) == 5
        assert result.columns.levels[0].tolist() == ["val"]
        assert result.columns.levels[1].tolist() == ["sum", "mean"]

    # ========== NEW TESTS - Additional Edge Cases ==========

    def test_groupby_agg_list_numeric_only_with_nans(self):
        """Test GroupBy with NaN values."""
        df = DataFrame(
            {
                "key": ["A", "A", "B", "B"],
                "val": [1, np.nan, 3, 4],
                "text": ["a", "b", "c", "d"],
            }
        )
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)

        assert isinstance(result, DataFrame)
        assert result.loc["A", ("val", "sum")] == 1.0
        assert result.loc["B", ("val", "sum")] == 7.0

    @pytest.mark.parametrize("as_index", [True, False])
    def test_groupby_agg_list_numeric_only_as_index(self, as_index):
        """Test GroupBy with as_index parameter."""
        df = DataFrame(
            {
                "key": ["A", "A", "B", "B"],
                "val": [1, 2, 3, 4],
                "text": ["a", "b", "c", "d"],
            }
        )
        result = df.groupby("key", as_index=as_index).agg(
            ["sum", "mean"], numeric_only=True
        )

        if as_index:
            assert result.index.name == "key"
        else:
            assert "key" in result.columns

    def test_groupby_agg_list_numeric_only_datetime_column(self):
        """Test GroupBy with datetime columns excluded."""
        df = DataFrame(
            {
                "key": ["A", "A", "B", "B"],
                "val": [1, 2, 3, 4],
                "date": pd.date_range("2020-01-01", periods=4),
                "text": ["a", "b", "c", "d"],
            }
        )
        result = df.groupby("key").agg(["sum", "mean"], numeric_only=True)
        expected = DataFrame(
            [[3, 1.5], [7, 3.5]],
            index=pd.Index(["A", "B"], name="key"),
            columns=pd.MultiIndex.from_product([["val"], ["sum", "mean"]]),
        )
        tm.assert_frame_equal(result, expected)
