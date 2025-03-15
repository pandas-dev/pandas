"""
Standalone script to test nrows parameter with adjacent tables in Excel files.
This script can be run directly with Python without using pytest.

Usage:
    python pandas/tests/io/excel/run_nrows_test.py
"""
import os
import tempfile
import pandas as pd


def run_test():
    """
    Test that nrows parameter correctly handles adjacent tables.

    This test creates two Excel files:
    1. One with a blank row between two tables
    2. One with no blank row between two tables

    Then it verifies that reading with nrows=3 returns only the first table
    in both cases.
    """
    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create test files
        file1 = os.path.join(tmp_dir, "with_blank.xlsx")
        file2 = os.path.join(tmp_dir, "no_blank.xlsx")

        # Create test data
        df_upper = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        df_lower = pd.DataFrame({"A": [7, 8, 9], "B": [10, 11, 12]})

        print("Creating Excel files...")

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

        print("Reading Excel files with nrows=3...")

        # Read with nrows=3 (should only get the first table)
        df1 = pd.read_excel(file1, nrows=3)
        df2 = pd.read_excel(file2, nrows=3)

        # Expected result - just the first table
        expected = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

        # Verify results
        print("Verifying results...")
        pd.testing.assert_frame_equal(df1, expected)
        pd.testing.assert_frame_equal(df2, expected)

        # Verify shapes
        assert df1.shape == (3, 2), f"Expected (3, 2) but got {df1.shape}"
        assert df2.shape == (3, 2), f"Expected (3, 2) but got {df2.shape}"

        # Verify last row doesn't contain headers from second table
        assert df2.iloc[-1, 0] == 3, f"Expected 3 but got {df2.iloc[-1, 0]}"
        assert df2.iloc[-1, 1] == 6, f"Expected 6 but got {df2.iloc[-1, 1]}"

        print("All tests passed!")


if __name__ == "__main__":
    run_test()
