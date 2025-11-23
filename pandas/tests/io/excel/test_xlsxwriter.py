import contextlib
import uuid

import pytest

import pandas as pd
from pandas import DataFrame

from pandas.io.excel import ExcelWriter

xlsxwriter = pytest.importorskip("xlsxwriter")

# xfail marker for pending autofilter feature; see #62994
xfail_autofilter = pytest.mark.xfail(
    reason="Excel header autofilter not yet implemented on main; see #62994",
    strict=False,
)


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


def test_autofilter_empty_dataframe(tmp_excel):
    # GH 61194 - Edge case: empty DataFrame with autofilter
    openpyxl = pytest.importorskip("openpyxl")
    df = DataFrame()
    df.to_excel(tmp_excel, engine="xlsxwriter", autofilter=True)

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        ws = wb.active
        # Empty DataFrame should still set autofilter (even if range is just header)
        assert ws.auto_filter.ref is not None


def test_autofilter_single_row(tmp_excel):
    # GH 61194 - Edge case: single row DataFrame
    openpyxl = pytest.importorskip("openpyxl")
    df = DataFrame({"A": [1], "B": [2]})
    df.to_excel(tmp_excel, engine="xlsxwriter", autofilter=True, index=False)

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        ws = wb.active
        assert ws.auto_filter.ref is not None
        # Should cover header + 1 row: A1:B2
        assert ws.auto_filter.ref == "A1:B2"


def test_autofilter_single_column(tmp_excel):
    # GH 61194 - Edge case: single column DataFrame
    openpyxl = pytest.importorskip("openpyxl")
    df = DataFrame({"A": [1, 2, 3]})
    df.to_excel(tmp_excel, engine="xlsxwriter", autofilter=True, index=False)

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        ws = wb.active
        assert ws.auto_filter.ref is not None
        # Should cover header + 3 rows: A1:A4
        assert ws.auto_filter.ref == "A1:A4"


def test_autofilter_no_header(tmp_excel):
    # GH 61194 - Edge case: autofilter with header=False
    openpyxl = pytest.importorskip("openpyxl")
    df = DataFrame([[1, 2], [3, 4]])
    df.to_excel(
        tmp_excel,
        engine="xlsxwriter",
        autofilter=True,
        header=False,
        index=False,
    )

    with contextlib.closing(openpyxl.load_workbook(tmp_excel)) as wb:
        ws = wb.active
        assert ws.auto_filter.ref is not None
        # Without header, filter should start at first row: A1:B2
        assert ws.auto_filter.ref == "A1:B2"


@xfail_autofilter
def test_to_excel_autofilter_xlsxwriter(tmp_excel):
    openpyxl = pytest.importorskip("openpyxl")

    df = DataFrame({"A": [1, 2], "B": [3, 4]})
    # Write with xlsxwriter, verify via openpyxl that an autofilter exists
    df.to_excel(tmp_excel, engine="xlsxwriter", index=False, autofilter=True)

    wb = openpyxl.load_workbook(tmp_excel)
    try:
        ws = wb[wb.sheetnames[0]]
        assert ws.auto_filter is not None
        assert ws.auto_filter.ref is not None
        # Verify filter covers all columns (A and B)
        assert "A" in ws.auto_filter.ref
        assert "B" in ws.auto_filter.ref
    finally:
        wb.close()


@xfail_autofilter
def test_to_excel_autofilter_startrow_startcol_xlsxwriter(tmp_excel):
    openpyxl = pytest.importorskip("openpyxl")

    df = DataFrame({"A": [1, 2], "B": [3, 4]})
    df.to_excel(
        tmp_excel,
        engine="xlsxwriter",
        index=False,
        autofilter=True,
        startrow=2,
        startcol=1,
    )

    wb = openpyxl.load_workbook(tmp_excel)
    try:
        ws = wb[wb.sheetnames[0]]
        assert ws.auto_filter is not None
        assert ws.auto_filter.ref is not None
        # Filter should be offset by startrow=2 and startcol=1 (B3:D5)
        assert ws.auto_filter.ref.startswith("B")
        assert "3" in ws.auto_filter.ref
    finally:
        wb.close()


@xfail_autofilter
def test_to_excel_autofilter_multiindex_merge_cells_xlsxwriter(tmp_excel):
    openpyxl = pytest.importorskip("openpyxl")

    df = DataFrame(
        [[1, 2, 3, 4], [5, 6, 7, 8]],
        columns=pd.MultiIndex.from_tuples(
            [("A", "a"), ("A", "b"), ("B", "a"), ("B", "b")]
        ),
    )
    df.to_excel(
        tmp_excel,
        engine="xlsxwriter",
        index=False,
        autofilter=True,
        merge_cells=True,
    )

    wb = openpyxl.load_workbook(tmp_excel)
    try:
        ws = wb[wb.sheetnames[0]]
        assert ws.auto_filter is not None
        assert ws.auto_filter.ref is not None
    finally:
        wb.close()


@xfail_autofilter
def test_to_excel_autofilter_multiindex_no_merge_xlsxwriter(tmp_excel):
    openpyxl = pytest.importorskip("openpyxl")

    df = DataFrame(
        [[1, 2, 3, 4], [5, 6, 7, 8]],
        columns=pd.MultiIndex.from_tuples(
            [("A", "a"), ("A", "b"), ("B", "a"), ("B", "b")]
        ),
    )
    df.to_excel(
        tmp_excel,
        engine="xlsxwriter",
        index=False,
        autofilter=True,
        merge_cells=False,
    )

    wb = openpyxl.load_workbook(tmp_excel)
    try:
        ws = wb[wb.sheetnames[0]]
        assert ws.auto_filter is not None
        assert ws.auto_filter.ref is not None
    finally:
        wb.close()
