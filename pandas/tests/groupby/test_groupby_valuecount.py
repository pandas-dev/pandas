"""
these are systematically testing all of the args to value_counts
with different size combinations. This is to ensure stability of the sorting
and proper parameter handling
"""

from itertools import product

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Grouper,
    MultiIndex,
)
import pandas._testing as tm


@pytest.fixture()
def dataframe1():
    # Create a multi-column, multi-row index
    TEST_DF = pd.DataFrame({
        "c1": ["a", "b", "c"],
        "c2": ["x", "y", "y"]
        }, index=[0, 1, 1])
    return TEST_DF
    """
    c1 c2
    0  a  x
    1  b  y
    1  c  y
    """

def test_groupby_row_subset_on_multiindex(dataframe1):

    # Test value_counts with argumrnt subset
    subset_result = dataframe.groupby(level=0).value_counts(subset=["c2"])

    expected_df = pd.DataFrame({
        "c2": ["x", "y", "y"]
        }, index=[0, 1, 1])
    #value_counts without argument subset
    expected_result = expected_df.groupby(level=0).values_count()
    """
    expected
       c2
    0  x     1
    1  y     2

    actual(before fix)
       c1  c2
    0  a   x     1
    1  b   y     1
       c   y     1

    """
    tm.assert_series_equal(result, expected)

def test_groupby_column_subset_on_multiindex(dataframe1):

    # Test value_counts with argumrnt subset
    subset_result = dataframe.groupby(level=1).value_counts(subset=["c2"])

    expected_df = pd.DataFrame({
        "c2": ["x", "y", "y"]
        }, index=[0, 1, 1])
    #value_counts without argument subset
    expected_result = expected_df.groupby(level=0).values_count()
    """
    expected
       c2
    0  x     1
    1  y     2

    actual(before fix)
       c1  c2
    0  a   x     1
    1  b   y     1
       c   y     1

    """
    tm.assert_series_equal(result, expected)

#run test through parametrize
@pytest.mark.parametrize("subset_col", [
    [], 
    ["c4"], 
    ["c1"], 
    ["c1", "c2", "c3", "c4"], 
    ["c1", "c2", "c3"]])
def test_moderate_groupby_value_counts_subset(subset_col):
    temp_df = pd.DataFrame(
        {
            "c1": ["a", "b", "c", "b"],
            "c2": [3, 2, 4, 2],
            "c3": [3, 1, 1, 2],
            "c4": ["a", "b", "c", "b"],
            "c5": ["x", "y", "y", "x"],
        },
        index=[0, 1, 1, 2],
    )
    subset_result = temp_df.groupby(level=0).value_counts(subset=subset_col)
    subset_expected = temp_df.groupby(level=0)[subset_col].value_counts()

    """
    actual output before fixing(p1)
        c1  c2  c3  c4  c5
    0  a   3   3   a   x     1
    1  b   2   1   b   y     1
       c   4   1   c   y     1
    2  b   2   2   b   x     1
    expected
       c1
    0  a     1
    1  b     1
       c     1
    2  b     1
    """
    tm.assert_series_equal(subset_result, subset_expected)