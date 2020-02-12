import warnings

import pytest

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.excel import ExcelWriter

xlsxwriter = pytest.importorskip("xlsxwriter")

pytestmark = pytest.mark.parametrize("ext", [".xlsx"])


def test_column_format(ext):
    # Test that column formats are applied to cells. Test for issue #9167.
    # Applicable to xlsxwriter only.
    with warnings.catch_warnings():
        # Ignore the openpyxl lxml warning.
        warnings.simplefilter("ignore")
        openpyxl = pytest.importorskip("openpyxl")

    with tm.ensure_clean(ext) as path:
        frame = DataFrame({"A": [123456, 123456], "B": [123456, 123456]})

        writer = ExcelWriter(path)
        frame.to_excel(writer)

        # Add a number format to col B and ensure it is applied to cells.
        num_format = "#,##0"
        write_workbook = writer.book
        write_worksheet = write_workbook.worksheets()[0]
        col_format = write_workbook.add_format({"num_format": num_format})
        write_worksheet.set_column("B:B", None, col_format)
        writer.save()

        read_workbook = openpyxl.load_workbook(path)
        try:
            read_worksheet = read_workbook["Sheet1"]
        except TypeError:
            # compat
            read_worksheet = read_workbook.get_sheet_by_name(name="Sheet1")

        # Get the number format from the cell.
        try:
            cell = read_worksheet["B2"]
        except TypeError:
            # compat
            cell = read_worksheet.cell("B2")

        try:
            read_num_format = cell.number_format
        except AttributeError:
            read_num_format = cell.style.number_format._format_code

        assert read_num_format == num_format


def test_write_append_mode_raises(ext):
    msg = "Append mode is not supported with xlsxwriter!"

    with tm.ensure_clean(ext) as f:
        with pytest.raises(ValueError, match=msg):
            ExcelWriter(f, engine="xlsxwriter", mode="a")
