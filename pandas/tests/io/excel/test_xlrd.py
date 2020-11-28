import pytest

import pandas as pd
import pandas._testing as tm

from pandas.io.excel import ExcelFile

xlrd = pytest.importorskip("xlrd")
xlwt = pytest.importorskip("xlwt")


@pytest.fixture(autouse=True)
def skip_ods_and_xlsb_files(read_ext):
    if read_ext == ".ods":
        pytest.skip("Not valid for xlrd")
    if read_ext == ".xlsb":
        pytest.skip("Not valid for xlrd")


def test_read_xlrd_book(read_ext, frame):
    df = frame

    engine = "xlrd"
    sheet_name = "SheetA"

    with tm.ensure_clean(read_ext) as pth:
        df.to_excel(pth, sheet_name)
        book = xlrd.open_workbook(pth)

        with ExcelFile(book, engine=engine) as xl:
            result = pd.read_excel(xl, sheet_name=sheet_name, index_col=0)
            tm.assert_frame_equal(df, result)

        result = pd.read_excel(book, sheet_name=sheet_name, engine=engine, index_col=0)
        tm.assert_frame_equal(df, result)


# TODO: test for openpyxl as well
def test_excel_table_sheet_by_index(datapath, read_ext):
    path = datapath("io", "data", "excel", f"test1{read_ext}")
    with ExcelFile(path, engine="xlrd") as excel:
        with pytest.raises(xlrd.XLRDError):
            pd.read_excel(excel, sheet_name="asdf")


def test_excel_file_warning_with_xlsx_file(datapath):
    # GH 29375
    path = datapath("io", "data", "excel", "test1.xlsx")
    # DeprecationWarning: "This method will be removed in future versions.
    # Use 'tree.iter()' or 'list(tree.iter())' instead."
    with tm.assert_produces_warning(
        DeprecationWarning, check_stacklevel=False, raise_on_extra_warnings=False
    ):
        ExcelFile(path, engine=None)


def test_read_excel_warning_with_xlsx_file(tmpdir, datapath):
    # GH 29375
    path = datapath("io", "data", "excel", "test1.xlsx")
    # DeprecationWarning: "This method will be removed in future versions.
    # Use 'tree.iter()' or 'list(tree.iter())' instead."
    with tm.assert_produces_warning(
        DeprecationWarning, check_stacklevel=False, raise_on_extra_warnings=False
    ):
        pd.read_excel(path, "Sheet1", engine=None)
