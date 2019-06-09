import pandas.util._test_decorators as td

import pandas as pd
import pandas.util.testing as tm
from pandas.util.testing import ensure_clean

from pandas.io.excel import ExcelFile


@td.skip_if_no('xlrd')
class TestXlrdReader:
    """
    This is the base class for the xlrd tests, and 3 different file formats
    are supported: xls, xlsx, xlsm
    """

    @td.skip_if_no("xlwt")
    def test_read_xlrd_book(self, read_ext, frame):
        import xlrd
        df = frame

        engine = "xlrd"
        sheet_name = "SheetA"

        with ensure_clean(read_ext) as pth:
            df.to_excel(pth, sheet_name)
            book = xlrd.open_workbook(pth)

            with ExcelFile(book, engine=engine) as xl:
                result = pd.read_excel(xl, sheet_name, index_col=0)
                tm.assert_frame_equal(df, result)

            result = pd.read_excel(book, sheet_name=sheet_name,
                                   engine=engine, index_col=0)
            tm.assert_frame_equal(df, result)
