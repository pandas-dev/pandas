import warnings

import pytest

from pandas import DataFrame, read_excel, MultiIndex
from pandas.util.testing import ensure_clean, assert_frame_equal

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

    with ensure_clean(ext) as path:
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
        except Exception:
            read_num_format = cell.style.number_format._format_code

        assert read_num_format == num_format


def test_write_append_mode_raises(ext):
    msg = "Append mode is not supported with xlsxwriter!"

    with ensure_clean(ext) as f:
        with pytest.raises(ValueError, match=msg):
            ExcelWriter(f, engine="xlsxwriter", mode="a")


def test_constant_memory_regularindex(ext):
    # Test if cells are written row by row which is the requirement
    # when exporting excel using xlsxwriter with constant_memory
    # set True, for regular index.
    # Test for issue #15392.
    # Applicable to xlsxwriter only.
    with ensure_clean(ext) as path:
        df = DataFrame({"A": [123456, 123456], "B": [123456, 123456]})

        with ExcelWriter(
            path, engine="xlsxwriter", options=dict(constant_memory=True)
        ) as writer:
            df.to_excel(writer)

        read_df = read_excel(path, header=0, index_col=0)

        assert_frame_equal(df, read_df)


def test_constant_memory_multiindex(ext):
    # Test if cells are written row by row which is the requirement
    # when exporting excel using xlsxwriter with constant_memory set
    # True, for MultiIndex.
    # Test for issue #15392.
    # Applicable to xlsxwriter only.
    with ensure_clean(ext) as path:
        df = DataFrame({"A": [123456, 123456], "B": [123456, 123456]})
        df.index = MultiIndex.from_arrays([["a", "a"], [1, 2]])

        with ExcelWriter(
            path, engine="xlsxwriter", options=dict(constant_memory=True)
        ) as writer:
            df.to_excel(writer)

        read_df = read_excel(path, header=0, index_col=[0, 1])

        assert_frame_equal(df, read_df)

def test_constant_memory_multiheader(ext):
    # Test if cells of a header of MultiIndex are written row by row
    # Test for issue #15392.
    # Applicable to xlsxwriter only.
    with ensure_clean(ext) as path:
        df = DataFrame({"A": [123456, 123456], "B": [123456, 123456]})
        df.columns = MultiIndex.from_arrays([["a", "a"], [1, 2]])

        with ExcelWriter(
            path, engine="xlsxwriter", options=dict(constant_memory=True)
        ) as writer:
            df.to_excel(writer)

        read_df = read_excel(path, header=[0, 1], index_col=0)

        assert_frame_equal(df, read_df)