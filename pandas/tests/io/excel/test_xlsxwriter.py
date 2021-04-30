from datetime import (
    date,
    datetime,
)
import warnings

import numpy as np
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

        with ExcelWriter(path) as writer:
            frame.to_excel(writer)

            # Add a number format to col B and ensure it is applied to cells.
            num_format = "#,##0"
            write_workbook = writer.book
            write_worksheet = write_workbook.worksheets()[0]
            col_format = write_workbook.add_format({"num_format": num_format})
            write_worksheet.set_column("B:B", None, col_format)

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


@pytest.mark.parametrize("c_idx_levels", [1, 2])
@pytest.mark.parametrize("r_idx_levels", [1, 2])
@pytest.mark.parametrize("head", [True, False, ["col1", "col2"]])
@pytest.mark.parametrize("ind", [True, False])
@pytest.mark.parametrize(
    "data,fmt,out,default_out",
    [
        # Format specified should not have any impact on the data
        # as the type of the data is int
        (np.random.randint(1, 100), "0%", "General", None),
        # Format specified should format the data
        (
            datetime(2013, 1, 13, 18, 20, 52),
            "DD/MM/YYYY",
            "DD/MM/YYYY",
            "YYYY-MM-DD HH:MM:SS",
        ),
        # Format specified should format the data
        (date(2014, 1, 31), "MMM", "MMM", "YYYY-MM-DD"),
        # Format specified should not have any impact on the data
        # as the type of the data is bool
        (True, "MMM", "General", None),
    ],
)
def test_num_formats(
    ext, c_idx_levels, r_idx_levels, head, ind, data, fmt, out, default_out
):
    # GH 30275
    # Testing out num_formats with various dtypes
    df = tm.makeCustomDataframe(
        6,
        2,
        c_idx_names=False,
        r_idx_names=False,
        c_idx_nlevels=c_idx_levels,
        r_idx_nlevels=r_idx_levels,
        data_gen_f=lambda r, c: data,
    )
    format_dict = dict(zip(df.columns, [fmt] * len(df.columns)))

    if c_idx_levels > 1:
        # Changing the value of out because num_formats has not been implemented
        # for multindex columns and the format specified will not be format the data
        out = default_out if default_out is not None else "General"
        if not ind:
            pytest.skip("Not implemented")

    out_fmt = [out] * len(df.columns)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        openpyxl = pytest.importorskip("openpyxl")

    with tm.ensure_clean(ext) as path:
        with ExcelWriter(path, num_formats=format_dict) as writer:
            df.to_excel(writer, columns=None, header=head, index=ind)

        # Changing cells to be checked as per the row levels
        if ind:
            if r_idx_levels == 1:
                cells = ["B4", "C4"]
            else:
                cells = ["C4", "D4"]
        else:
            cells = ["A4", "B4"]

        num_formats = []

        read_workbook = openpyxl.load_workbook(path)
        try:
            read_worksheet = read_workbook["Sheet1"]
        except TypeError:
            read_worksheet = read_workbook.get_sheet_by_name(name="Sheet1")

        for cl in cells:
            try:
                cell = read_worksheet[cl]
            except TypeError:
                cell = read_worksheet.cell(cl)

            try:
                read_num_format = cell.number_format
            except AttributeError:
                read_num_format = cell.style.number_format._format_code
            num_formats.append(read_num_format)

        assert num_formats == out_fmt


@pytest.mark.parametrize(
    "df,col_fmts,cols,out_fmt",
    [
        # Checking formats with timeseries data
        (
            tm.makeTimeDataFrame()[:6],
            {"B": "0%", "C": "#,##0"},
            None,
            ["YYYY-MM-DD HH:MM:SS", "General", "General", "General", "General"],
        ),
        # Checking formats with duplicated columns
        (
            DataFrame(
                {
                    "A": [1, 1, 1, 1],
                    "B": [2.34, 2.34, 2.34, 2.34],
                    "C": [
                        datetime(2014, 1, 31),
                        datetime(1999, 9, 24),
                        datetime(2014, 1, 31),
                        datetime(1999, 9, 24),
                    ],
                    "D": ["AA", "BB", "AA", "BB"],
                }
            ),
            {
                "B": "0%",
                "C": "YYYY",
            },
            ["B", "C", "C", "A"],
            ["General", "General", "YYYY", "YYYY", "General"],
        ),
    ],
)
def test_num_formats_others(ext, df, col_fmts, cols, out_fmt):
    # GH 30275
    # Testing out num_formats for cases not covered in test_num_formats
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        openpyxl = pytest.importorskip("openpyxl")

    with tm.ensure_clean(ext) as path:
        with ExcelWriter(path, num_formats=col_fmts) as writer:
            df.to_excel(writer, columns=cols, header=True, index=True)

        num_formats = []
        cells = ["A3", "B3", "C3", "D3", "E3"]

        read_workbook = openpyxl.load_workbook(path)
        try:
            read_worksheet = read_workbook["Sheet1"]
        except TypeError:
            read_worksheet = read_workbook.get_sheet_by_name(name="Sheet1")

        for cl in cells:
            try:
                cell = read_worksheet[cl]
            except TypeError:
                cell = read_worksheet.cell(cl)

            try:
                read_num_format = cell.number_format
            except AttributeError:
                read_num_format = cell.style.number_format._format_code
            num_formats.append(read_num_format)
        assert num_formats == out_fmt


def test_check_exceptions(ext):
    # GH 30275
    with tm.ensure_clean(ext) as path:
        with tm.ensure_clean(ext) as path:
            with pytest.raises(
                TypeError,
                match="Invalid type list, num_formats must be dict.",
            ):
                ExcelWriter(path, num_formats=[1, 2, 3])

            with pytest.raises(TypeError, match="Format for 'B' is not a string."):
                ExcelWriter(path, num_formats={"A": "0%", "B": 234})
