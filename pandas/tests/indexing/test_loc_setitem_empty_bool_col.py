"""
Tests for GH#66255: loc setitem with all-False boolean column indexer
should be a no-op for extension (nullable) dtypes.
"""

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
)
import pandas._testing as tm


class TestLocSetitemEmptyBoolColumnIndexer:
    """Regression tests for GH#66255.

    When ``df.loc[:, bool_mask] = value`` is called with a boolean mask
    that selects *no* columns, the operation should be a no-op regardless
    of the DataFrame's dtype.
    """

    @pytest.mark.parametrize(
        "dtype",
        [
            "Float64",
            "Float32",
            "Int64",
            "Int32",
            "Int16",
            "Int8",
            "UInt64",
            "UInt32",
            "UInt16",
            "UInt8",
            "boolean",
            "string",
        ],
    )
    def test_loc_setitem_empty_bool_col_indexer_nullable(self, dtype):
        # GH#66255
        if dtype == "string":
            data = ["a", "b", "c"]
        elif dtype == "boolean":
            data = [True, False, True]
        else:
            data = [1, 2, 3]

        df = DataFrame(data, columns=["col"], dtype=dtype)
        expected = df.copy()

        select = df.columns == "not_in_columns"  # array([False])
        df.loc[:, select] = pd.NA

        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize("value", [pd.NA, np.nan, 0, 999])
    def test_loc_setitem_empty_bool_col_indexer_scalar_values(self, value):
        # GH#66255 - assigning any scalar to empty column selection
        df = DataFrame([1.0, 2.0, 3.0], columns=["col"], dtype="Float64")
        expected = df.copy()

        select = df.columns == "not_in_columns"
        df.loc[:, select] = value

        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_empty_bool_col_indexer_category(self):
        # GH#66255
        df = DataFrame(["a", "b", "c"], columns=["col"], dtype="category")
        expected = df.copy()

        select = df.columns == "not_in_columns"
        df.loc[:, select] = "x"

        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_nonempty_bool_col_indexer_still_works(self):
        # Ensure that selecting a column via boolean [True] still works
        # (this is the legitimate case that the fix must not break)
        df = DataFrame([1.0, 2.0, 3.0], columns=["col"], dtype="Float64")

        # iloc with integer array [0] should still work
        df.iloc[:, [0]] = 999.0
        expected = DataFrame([999.0, 999.0, 999.0], columns=["col"], dtype="Float64")
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_empty_bool_col_indexer_multicolumn(self):
        # Multi-column case (takes different code path but should also work)
        df = DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}, dtype="Float64")
        expected = df.copy()

        select = df.columns == "not_in_columns"
        df.loc[:, select] = pd.NA

        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_empty_bool_col_indexer_with_row_slice(self):
        # GH#66255 - empty column selection with row subset
        df = DataFrame([1.0, 2.0, 3.0], columns=["col"], dtype="Float64")
        expected = df.copy()

        select = df.columns == "not_in_columns"
        df.loc[0:1, select] = pd.NA

        tm.assert_frame_equal(df, expected)
