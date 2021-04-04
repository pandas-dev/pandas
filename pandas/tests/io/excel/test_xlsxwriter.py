from datetime import datetime
import warnings

import pytest

from pandas import DataFrame
import pandas._testing as tm

from pandas.io.excel import ExcelWriter

xlsxwriter = pytest.importorskip("xlsxwriter")

pytestmark = pytest.mark.parametrize("ext", [".xlsx"])


@pytest.mark.parametrize(
    "data,num_format,fmt,col,expected",
    [
        (
            DataFrame({"A": [123456, 123456], "B": [123456, 123456]}),
            "#,##0",
            None,
            "B",
            "#,##0",
        ),
        (
            DataFrame(
                {
                    "A": [
                        datetime(2013, 1, 13, 1, 2, 3),
                        datetime(2013, 1, 13, 2, 45, 56),
                    ],
                    "B": [
                        datetime(1998, 5, 26, 23, 33, 4),
                        datetime(2014, 2, 28, 13, 5, 13),
                    ],
                }
            ),
            "DD-MMM-YYYY HH:MM:SS",
            None,
            "B",
            "DD-MMM-YYYY HH:MM:SS",
        ),
        (tm.makeTimeDataFrame()[:5], "MMM", None, "A", "YYYY-MM-DD HH:MM:SS"),
        (
            tm.makeTimeDataFrame()[:5],
            "MMM",
            "DD-MMM-YYYY HH:MM:SS",
            "A",
            "DD-MMM-YYYY HH:MM:SS",
        ),
    ],
)
def test_column_format(ext, data, num_format, fmt, col, expected):
    # Test that column formats are applied to cells. Test for issue #9167.
    # Applicable to xlsxwriter only.
    with warnings.catch_warnings():
        # Ignore the openpyxl lxml warning.
        warnings.simplefilter("ignore")
        openpyxl = pytest.importorskip("openpyxl")

    with tm.ensure_clean(ext) as path:
        frame = data

        with ExcelWriter(path, datetime_format=fmt) as writer:
            frame.to_excel(writer)

            # Add a number format to col
            write_workbook = writer.book
            write_worksheet = write_workbook.worksheets()[0]
            col_format = write_workbook.add_format({"num_format": num_format})
            write_worksheet.set_column(f"{col}:{col}", None, col_format)

        read_workbook = openpyxl.load_workbook(path)
        try:
            read_worksheet = read_workbook["Sheet1"]
        except TypeError:
            # compat
            read_worksheet = read_workbook.get_sheet_by_name(name="Sheet1")

        # Get the number format from the cell.
        try:
            cell = read_worksheet[f"{col}2"]
        except TypeError:
            # compat
            cell = read_worksheet.cell(f"{col}2")

        try:
            read_num_format = cell.number_format
        except AttributeError:
            read_num_format = cell.style.number_format._format_code

        assert read_num_format == expected


def test_write_append_mode_raises(ext):
    msg = "Append mode is not supported with xlsxwriter!"

    with tm.ensure_clean(ext) as f:
        with pytest.raises(ValueError, match=msg):
            ExcelWriter(f, engine="xlsxwriter", mode="a")
