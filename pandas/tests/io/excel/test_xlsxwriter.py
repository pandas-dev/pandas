from datetime import (
    date,
    datetime,
)
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


@pytest.mark.parametrize("head", [True, False, ["col1", "col2"]])
@pytest.mark.parametrize(
    "data,fmt_set,out_fmt",
    [
        (
            datetime(2013, 1, 13, 18, 20, 52),
            "DD/MM/YYYY",
            "DD/MM/YYYY",
        ),
        (date(2014, 1, 31), "MMM", "MMM"),
    ],
)
def test_formatters(ext, head, data, fmt_set, out_fmt):
    # GH 30275
    df = tm.makeCustomDataframe(
        6,
        2,
        c_idx_names=False,
        r_idx_names=False,
        data_gen_f=lambda r, c: data,
    )
    out_fmt_list = [out_fmt] * len(df.columns)
    dict_for_formatters = dict(zip(df.columns, [fmt_set] * len(df.columns)))
    cells_to_check = ["B4", "C4"]

    with tm.ensure_clean(ext) as path:
        with ExcelWriter(path, formatters=dict_for_formatters) as writer:
            df.to_excel(writer, columns=None, header=head, index=True)

        openpyxl = pytest.importorskip("openpyxl")
        read_workbook = openpyxl.load_workbook(path)
        read_worksheet = read_workbook["Sheet1"]

        formats = []
        for cl in cells_to_check:
            cell = read_worksheet[cl]
            read_num_format = cell.number_format
            formats.append(read_num_format)

        assert formats == out_fmt_list


@pytest.mark.parametrize("head", [True, False, ["col1", "col2"]])
@pytest.mark.parametrize(
    "data,fmt_set,out_fmt",
    [
        (
            datetime(2013, 1, 13, 18, 20, 52),
            "DD/MM/YYYY",
            "YYYY-MM-DD HH:MM:SS",
        ),
        (date(2014, 1, 31), "MMM", "YYYY-MM-DD"),
    ],
)
def test_formatters_multiindex_cols(ext, head, data, fmt_set, out_fmt):
    # GH 30275
    df = tm.makeCustomDataframe(
        6,
        2,
        c_idx_names=False,
        r_idx_names=False,
        c_idx_nlevels=2,
        data_gen_f=lambda r, c: data,
    )
    out_fmt_list = [out_fmt] * len(df.columns)
    dict_for_formatters = dict(zip(df.columns, [fmt_set] * len(df.columns)))
    cells_to_check = ["B4", "C4"]

    with tm.ensure_clean(ext) as path:
        with ExcelWriter(path, formatters=dict_for_formatters) as writer:
            df.to_excel(writer, columns=None, header=head, index=True)

        openpyxl = pytest.importorskip("openpyxl")
        read_workbook = openpyxl.load_workbook(path)
        read_worksheet = read_workbook["Sheet1"]

        formats = []
        for cl in cells_to_check:
            cell = read_worksheet[cl]
            read_num_format = cell.number_format
            formats.append(read_num_format)

        assert formats == out_fmt_list


@pytest.mark.parametrize(
    "df,dict_for_formatters,cols,out_fmt_list",
    [
        # Timeseries data
        (
            tm.makeTimeDataFrame()[:6],
            {"B": "0%", "C": "#,##0"},
            None,
            ["YYYY-MM-DD HH:MM:SS", "General", "General", "General", "General"],
        ),
        # Duplicated columns
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
                    "D": ["abc", "abc", "def", "def"],
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
def test_formatters_others(ext, df, dict_for_formatters, cols, out_fmt_list):
    # GH 30275
    cells_to_check = ["A3", "B3", "C3", "D3", "E3"]

    with tm.ensure_clean(ext) as path:
        with ExcelWriter(path, formatters=dict_for_formatters) as writer:
            df.to_excel(writer, columns=cols, header=True, index=True)

        openpyxl = pytest.importorskip("openpyxl")
        read_workbook = openpyxl.load_workbook(path)
        read_worksheet = read_workbook["Sheet1"]

        formats = []
        for cl in cells_to_check:
            cell = read_worksheet[cl]
            read_num_format = cell.number_format
            formats.append(read_num_format)
        assert formats == out_fmt_list


def test_check_exceptions(ext):
    # GH 30275
    with tm.ensure_clean(ext) as path:
        with tm.ensure_clean(ext) as path:
            with pytest.raises(
                TypeError,
                match="Invalid type list, formatters must be dict.",
            ):
                ExcelWriter(path, formatters=[1, 2, 3])

            with pytest.raises(TypeError, match="Format for 'B' is not a string."):
                ExcelWriter(path, formatters={"A": "0%", "B": 234})
