"""
Test for GH-61123: nrows parameter with adjacent tables in Excel files.
"""
import os
import pytest
import pandas as pd
import pandas._testing as tm


@pytest.mark.skipif(not os.path.exists("pandas/io/excel/_openpyxl.py"), reason="openpyxl not installed")
def test_nrows_with_adjacent_tables(tmp_path):
    """
    Test that nrows parameter correctly handles adjacent tables.

    This test creates two Excel files:
    1. One with a blank row between two tables
    2. One with no blank row between two tables

    Then it verifies that reading with nrows=3 returns only the first table
    in both cases.
    """
    # Create test files
    file1 = tmp_path / "with_blank.xlsx"
    file2 = tmp_path / "no_blank.xlsx"

    # Create test data
    df_upper = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
    df_lower = pd.DataFrame({"A": [7, 8, 9], "B": [10, 11, 12]})

    # Create file with blank row between tables
    with pd.ExcelWriter(file1) as writer:
        df_upper.to_excel(writer, sheet_name="Sheet1", index=False)
        # Add blank row by starting lower table at row 5 (0-based index + header)
        df_lower.to_excel(writer, sheet_name="Sheet1", startrow=5, index=False)

    # Create file with no blank row between tables
    with pd.ExcelWriter(file2) as writer:
        df_upper.to_excel(writer, sheet_name="Sheet1", index=False)
        # No blank row, lower table starts right after (row 4 = header of second table)
        df_lower.to_excel(writer, sheet_name="Sheet1", startrow=4, index=False)

    # Read with nrows=3 (should only get the first table)
    df1 = pd.read_excel(file1, nrows=3)
    df2 = pd.read_excel(file2, nrows=3)

    # Expected result - just the first table
    expected = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

    # Verify results
    tm.assert_frame_equal(df1, expected)
    tm.assert_frame_equal(df2, expected)

    # Verify shapes
    assert df1.shape == (3, 2)
    assert df2.shape == (3, 2)

    # Verify last row doesn't contain headers from second table
    assert df2.iloc[-1, 0] == 3
    assert df2.iloc[-1, 1] == 6
