import pytest

import pandas as pd
import pandas.util.testing as tm

from pandas.io.excel import ExcelFile

xlrd = pytest.importorskip("xlrd")
xlwt = pytest.importorskip("xlwt")


@pytest.fixture(autouse=True)
def skip_ods_files(read_ext):
    if read_ext == ".ods":
        pytest.skip("Not valid for xlrd")


def test_read_xlrd_book(read_ext, frame):
    df = frame

    engine = "xlrd"
    sheet_name = "SheetA"

    with tm.ensure_clean(read_ext) as pth:
        df.to_excel(pth, sheet_name)
        book = xlrd.open_workbook(pth)

        with ExcelFile(book, engine=engine) as xl:
            result = pd.read_excel(xl, sheet_name, index_col=0)
            tm.assert_frame_equal(df, result)

        result = pd.read_excel(book, sheet_name=sheet_name, engine=engine, index_col=0)
        tm.assert_frame_equal(df, result)


# TODO: test for openpyxl as well
def test_excel_table_sheet_by_index(datapath, read_ext):
    path = datapath("io", "data", "excel", "test1{}".format(read_ext))
    with pd.ExcelFile(path) as excel:
        with pytest.raises(xlrd.XLRDError):
            pd.read_excel(excel, "asdf")


# See issue #29375
def test_excel_file_warning_with_default_engine(datapath):
    path = datapath("io", "data", "test1.xls")
    with tm.assert_produces_warning(FutureWarning):
        pd.ExcelFile(path)


# See issue #29375
def test_read_excel_warning_with_default_engine(tmpdir, datapath):
    path = datapath("io", "data", "test1.xls")
    with tm.assert_produces_warning(FutureWarning):
        pd.read_excel(path, "Sheet1")
