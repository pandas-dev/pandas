"""
Tests for reading Excel files with adjacent tables.
"""
import pytest
import pandas as pd
import pandas._testing as tm


class TestExcelAdjacentTables:
    """Tests for reading Excel files with adjacent tables."""

    @pytest.mark.parametrize("engine", ["openpyxl"])
    def test_nrows_with_adjacent_tables(self, engine, tmp_path):
        """
        Test that nrows parameter correctly handles adjacent tables.

        GH-61123: When using nrows to limit the number of rows read from an Excel file,
        the function should correctly handle cases where tables are adjacent (no blank
        row between them).
        """
        # Create test files with tables with and without blank rows between them
        # File 1: Two tables with a blank row between
        file1 = tmp_path / "test1.xlsx"
        df_upper = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        df_lower = pd.DataFrame({"A": [7, 8, 9], "B": [10, 11, 12]})

        with pd.ExcelWriter(file1, engine=engine) as writer:
            df_upper.to_excel(writer, sheet_name="Sheet1", index=False)
            # Add blank row by starting lower table at row 5 (0-based index + header)
            df_lower.to_excel(writer, sheet_name="Sheet1", startrow=5, index=False)

        # File 2: Two tables with no blank row
        file2 = tmp_path / "test2.xlsx"
        with pd.ExcelWriter(file2, engine=engine) as writer:
            df_upper.to_excel(writer, sheet_name="Sheet1", index=False)
            # No blank row, lower table starts right after (row 4 = header of second table)
            df_lower.to_excel(writer, sheet_name="Sheet1", startrow=4, index=False)

        # Read first 3 rows (header + 3 data rows)
        # Using nrows=3 to get exactly the upper table without blank rows
        df1 = pd.read_excel(file1, header=0, nrows=3, engine=engine)
        df2 = pd.read_excel(file2, header=0, nrows=3, engine=engine)

        # Expected data - just the upper table
        expected = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

        # Check content
        tm.assert_frame_equal(df1, expected)
        tm.assert_frame_equal(df2, expected)

        # Verify we didn't read the header of the next table in df2
        # If we did, the last row would contain column headers from the second table
        assert df1.shape == (3, 2), f"Expected (3, 2) but got {df1.shape}"
        assert df2.shape == (3, 2), f"Expected (3, 2) but got {df2.shape}"

        # Check specific values in the last row to ensure we didn't read the header
        assert df2.iloc[-1, 0] == 3, f"Expected 3 but got {df2.iloc[-1, 0]}"
        assert df2.iloc[-1, 1] == 6, f"Expected 6 but got {df2.iloc[-1, 1]}"
