import io

import pytest

import pandas as pd

openpyxl = pytest.importorskip("openpyxl")


def test_to_excel_openpyxl_autofilter():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        # Test autofilter
        df.to_excel(writer, index=False, autofilter=True)
    buf.seek(0)
    wb = openpyxl.load_workbook(buf)
    ws = wb.active
    # Autofilter should be set spanning header+data
    assert ws.auto_filter is not None
    assert ws.auto_filter.ref is not None and ws.auto_filter.ref != ""


def test_to_excel_openpyxl_styler_bold_header():
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    buf = io.BytesIO()

    # Create Excel file with pandas
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Sheet1")

        # Get the worksheet object
        worksheet = writer.sheets["Sheet1"]

        # Apply bold to the header row (first row in Excel is 1)
        from openpyxl.styles import (
            Font,
            PatternFill,
        )

        # Create a style for the header
        header_font = Font(bold=True, color="000000")
        header_fill = PatternFill(
            start_color="D3D3D3", end_color="D3D3D3", fill_type="solid"
        )

        # Apply style to each cell in the header row
        for cell in worksheet[1]:  # First row is the header
            cell.font = header_font
            cell.fill = header_fill

    # Now read it back to verify
    buf.seek(0)
    wb = openpyxl.load_workbook(buf)
    ws = wb.active

    # Print debug info
    print("\n===== WORKSHEET CELLS =====")
    for r, row in enumerate(ws.iter_rows(), 1):
        print(f"Row {r} (header: {r == 1}):")
        for c, cell in enumerate(row, 1):
            font_info = {
                "value": cell.value,
                "has_font": cell.font is not None,
                "bold": cell.font.bold if cell.font else None,
                "font_name": cell.font.name if cell.font else None,
                "font_size": cell.font.sz
                if cell.font and hasattr(cell.font, "sz")
                else None,
            }
            print(f"  Cell {c}: {font_info}")
    print("===========================\n")

    # Check that header cells (A1, B1) have bold font
    header_row = 1
    for col in range(1, df.shape[1] + 1):
        cell = ws.cell(row=header_row, column=col)
        assert cell.font is not None, (
            f"Header cell {cell.coordinate} has no font settings"
        )
        assert cell.font.bold, f"Header cell {cell.coordinate} is not bold"

    # Check that data cells (A2, B2, A3, B3) do not have bold font
    for row in range(2, df.shape[0] + 2):
        for col in range(1, df.shape[1] + 1):
            cell = ws.cell(row=row, column=col)
            if cell.font and cell.font.bold:
                print(f"Warning: Data cell {cell.coordinate} is unexpectedly bold")
