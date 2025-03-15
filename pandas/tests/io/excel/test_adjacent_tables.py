from __future__ import annotations

import pytest
import pandas as pd
import pandas._testing as tm

from pandas.io.excel import ExcelWriter


class TestAdjacentTables:
    """Tests for reading Excel files with adjacent tables."""

    @pytest.mark.parametrize(
        "engine,read_ext",
        [
            pytest.param("openpyxl", ".xlsx", marks=[pytest.mark.skip_if_no("openpyxl")]),
            pytest.param("xlsxwriter", ".xlsx", marks=[pytest.mark.skip_if_no("xlsxwriter")]),
        ],
    )
    def test_excel_read_adjacent_tables_nrows(self, engine, read_ext, tmp_path):
        """
        Test that nrows parameter correctly handles adjacent tables with and without blank rows.

        GH-61123
        """
        # Create test files with tables with and without blank rows between them
        # File 1: Two tables with a blank row between
        file1 = tmp_path / f"test1{read_ext}"
        df_upper = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        df_lower = pd.DataFrame({"A": [7, 8, 9], "B": [10, 11, 12]})

        with ExcelWriter(file1, engine=engine) as writer:
            df_upper.to_excel(writer, sheet_name="Sheet1", index=False)
            # Add blank row by starting lower table at row 5 (0-based index + header)
            df_lower.to_excel(writer, sheet_name="Sheet1", startrow=5, index=False)

        # File 2: Two tables with no blank row
        file2 = tmp_path / f"test2{read_ext}"
        with ExcelWriter(file2, engine=engine) as writer:
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

        # Fix the comparison warning by checking string values properly
        last_row_values = [str(x) for x in df2.iloc[-1].values]
        assert "A" not in last_row_values, "Second table header was incorrectly included"
        assert "B" not in last_row_values, "Second table header was incorrectly included"
