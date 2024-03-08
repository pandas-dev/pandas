import contextlib
import uuid

import pytest

from pandas import DataFrame

from pandas.io.excel import ExcelWriter

xlsxwriter = pytest.importorskip("xlsxwriter")


@pytest.fixture
def ext():
    return ".xlsx"


@pytest.fixture
def tmp_excel(ext, tmp_path):
    tmp = tmp_path / f"{uuid.uuid4()}{ext}"
    tmp.touch()
    return str(tmp)


def test_column_format(tmp_excel):
    # Test that column formats are applied to cells. Test for issue #9167.
    # Applicable to xlsxwriter only.
    openpyxl = pytest.importorskip("openpyxl")

    frame = DataFrame({"A": [123456, 123456], "B": [123456, 123456]})

    with ExcelWriter(tmp_excel) as writer:
        frame.to_excel(writer)

        # Add a number format to col B and ensure it is applied to cells.
        num_format = "#,##0"
        write_workbook = writer.book
        write_worksheet = write_workbook.worksheets()[0]
        col_format = write_workbook.add_format({"num_format": num_format})
        write_worksheet.set_column("B:B", None, col_format)

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as read_workbook:
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


def test_write_append_mode_raises(tmp_excel):
    msg = "Append mode is not supported with xlsxwriter!"

    with pytest.raises(ValueError, match=msg):
        ExcelWriter(tmp_excel, engine="xlsxwriter", mode="a")


@pytest.mark.parametrize("nan_inf_to_errors", [True, False])
def test_engine_kwargs(tmp_excel, nan_inf_to_errors):
    # GH 42286
    engine_kwargs = {"options": {"nan_inf_to_errors": nan_inf_to_errors}}
    with ExcelWriter(
        tmp_excel, engine="xlsxwriter", engine_kwargs=engine_kwargs
    ) as writer:
        assert writer.book.nan_inf_to_errors == nan_inf_to_errors


def test_book_and_sheets_consistent(tmp_excel):
    # GH#45687 - Ensure sheets is updated if user modifies book
    with ExcelWriter(tmp_excel, engine="xlsxwriter") as writer:
        assert writer.sheets == {}
        sheet = writer.book.add_worksheet("test_name")
        assert writer.sheets == {"test_name": sheet}
