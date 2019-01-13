import contextlib
from functools import partial
import os
import warnings

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame
from pandas.core.config import get_option, set_option
import pandas.util.testing as tm
from pandas.util.testing import ensure_clean

from pandas.io.excel import ExcelFile, read_excel

_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()
_frame = DataFrame(_seriesd)[:10]
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])[:10]
_tsframe = tm.makeTimeDataFrame()[:5]
_mixed_frame = _frame.copy()
_mixed_frame['foo'] = 'bar'


@contextlib.contextmanager
def ignore_xlrd_time_clock_warning():
    """
    Context manager to ignore warnings raised by the xlrd library,
    regarding the deprecation of `time.clock` in Python 3.7.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action='ignore',
            message='time.clock has been deprecated',
            category=DeprecationWarning)
        yield


@td.skip_if_no('xlrd', '1.0.0')
class XlrdSharedItems(object):

    @pytest.fixture(autouse=True)
    def setup_method(self, datapath):
        self.dirpath = datapath("io", "data")
        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()

    def get_excelfile(self, basename, ext):
        """
        Return test data ExcelFile instance.

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        excel : io.excel.ExcelFile
        """
        return ExcelFile(os.path.join(self.dirpath, basename + ext))

    def get_exceldf(self, basename, ext, *args, **kwds):
        """
        Return test data DataFrame.

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        df : DataFrame
        """
        pth = os.path.join(self.dirpath, basename + ext)
        return read_excel(pth, *args, **kwds)


class XlrdReadingTestsBase(XlrdSharedItems):
    # This is based on ExcelWriterBase

    @pytest.fixture(autouse=True, params=['xlrd'])
    def set_engine(self, request):
        func_name = "get_exceldf"
        old_func = getattr(self, func_name)
        new_func = partial(old_func, engine=request.param)
        setattr(self, func_name, new_func)
        yield
        setattr(self, func_name, old_func)

    @td.skip_if_no("xlrd", "1.0.1")  # see gh-22682
    def test_usecols_int(self, ext):

        df_ref = self.get_csv_refdf("test1")
        df_ref = df_ref.reindex(columns=["A", "B", "C"])

        # usecols as int
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            with ignore_xlrd_time_clock_warning():
                df1 = self.get_exceldf("test1", ext, "Sheet1",
                                       index_col=0, usecols=3)

        # usecols as int
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            with ignore_xlrd_time_clock_warning():
                df2 = self.get_exceldf("test1", ext, "Sheet2", skiprows=[1],
                                       index_col=0, usecols=3)

        # parse_cols instead of usecols, usecols as int
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            with ignore_xlrd_time_clock_warning():
                df3 = self.get_exceldf("test1", ext, "Sheet2", skiprows=[1],
                                       index_col=0, parse_cols=3)

        # TODO add index to xls file)
        tm.assert_frame_equal(df1, df_ref, check_names=False)
        tm.assert_frame_equal(df2, df_ref, check_names=False)
        tm.assert_frame_equal(df3, df_ref, check_names=False)

    @td.skip_if_no('xlrd', '1.0.1')  # GH-22682
    def test_usecols_list(self, ext):

        dfref = self.get_csv_refdf('test1')
        dfref = dfref.reindex(columns=['B', 'C'])
        df1 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols=[0, 2, 3])
        df2 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols=[0, 2, 3])

        with tm.assert_produces_warning(FutureWarning):
            with ignore_xlrd_time_clock_warning():
                df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                                       index_col=0, parse_cols=[0, 2, 3])

        # TODO add index to xls file)
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)
        tm.assert_frame_equal(df3, dfref, check_names=False)

    @td.skip_if_no('xlrd', '1.0.1')  # GH-22682
    def test_usecols_str(self, ext):

        dfref = self.get_csv_refdf('test1')

        df1 = dfref.reindex(columns=['A', 'B', 'C'])
        df2 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols='A:D')
        df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols='A:D')

        with tm.assert_produces_warning(FutureWarning):
            with ignore_xlrd_time_clock_warning():
                df4 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                                       index_col=0, parse_cols='A:D')

        # TODO add index to xls, read xls ignores index name ?
        tm.assert_frame_equal(df2, df1, check_names=False)
        tm.assert_frame_equal(df3, df1, check_names=False)
        tm.assert_frame_equal(df4, df1, check_names=False)

        df1 = dfref.reindex(columns=['B', 'C'])
        df2 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols='A,C,D')
        df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols='A,C,D')
        # TODO add index to xls file
        tm.assert_frame_equal(df2, df1, check_names=False)
        tm.assert_frame_equal(df3, df1, check_names=False)

        df1 = dfref.reindex(columns=['B', 'C'])
        df2 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols='A,C:D')
        df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols='A,C:D')
        tm.assert_frame_equal(df2, df1, check_names=False)
        tm.assert_frame_equal(df3, df1, check_names=False)

    @td.skip_if_no('xlrd', '1.0.1')  # GH-22682
    def test_deprecated_sheetname(self, ext):
        # gh-17964
        excel = self.get_excelfile('test1', ext)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            read_excel(excel, sheetname='Sheet1')

        with pytest.raises(TypeError):
            read_excel(excel, sheet='Sheet1')

    @td.skip_if_no('xlrd', '1.0.1')  # GH-22682
    def test_excel_table_sheet_by_index(self, ext):

        excel = self.get_excelfile('test1', ext)
        dfref = self.get_csv_refdf('test1')

        df1 = read_excel(excel, 0, index_col=0)
        df2 = read_excel(excel, 1, skiprows=[1], index_col=0)
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)

        df1 = excel.parse(0, index_col=0)
        df2 = excel.parse(1, skiprows=[1], index_col=0)
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)

        df3 = read_excel(excel, 0, index_col=0, skipfooter=1)
        tm.assert_frame_equal(df3, df1.iloc[:-1])

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            df4 = read_excel(excel, 0, index_col=0, skip_footer=1)
            tm.assert_frame_equal(df3, df4)

        df3 = excel.parse(0, index_col=0, skipfooter=1)
        tm.assert_frame_equal(df3, df1.iloc[:-1])

        import xlrd
        with pytest.raises(xlrd.XLRDError):
            read_excel(excel, 'asdf')

    @td.skip_if_no("xlrd", "1.0.1")  # see gh-22682
    def test_sheet_name_and_sheetname(self, ext):
        # gh-10559: Minor improvement: Change "sheet_name" to "sheetname"
        # gh-10969: DOC: Consistent var names (sheetname vs sheet_name)
        # gh-12604: CLN GH10559 Rename sheetname variable to sheet_name
        # gh-20920: ExcelFile.parse() and pd.read_xlsx() have different
        #           behavior for "sheetname" argument
        filename = "test1"
        sheet_name = "Sheet1"

        df_ref = self.get_csv_refdf(filename)
        df1 = self.get_exceldf(filename, ext,
                               sheet_name=sheet_name, index_col=0)  # doc
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            with ignore_xlrd_time_clock_warning():
                df2 = self.get_exceldf(filename, ext, index_col=0,
                                       sheetname=sheet_name)  # backward compat

        excel = self.get_excelfile(filename, ext)
        df1_parse = excel.parse(sheet_name=sheet_name, index_col=0)  # doc
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            df2_parse = excel.parse(index_col=0,
                                    sheetname=sheet_name)  # backward compat

        tm.assert_frame_equal(df1, df_ref, check_names=False)
        tm.assert_frame_equal(df2, df_ref, check_names=False)
        tm.assert_frame_equal(df1_parse, df_ref, check_names=False)
        tm.assert_frame_equal(df2_parse, df_ref, check_names=False)


@pytest.mark.parametrize("ext", ['.xls', '.xlsx', '.xlsm'])
class TestXlrdReader(XlrdReadingTestsBase):
    """
    This is the base class for the xlrd tests, and 3 different file formats
    are supported: xls, xlsx, xlsm
    """

    @td.skip_if_no("xlwt")
    def test_read_xlrd_book(self, ext):
        import xlrd
        df = self.frame

        engine = "xlrd"
        sheet_name = "SheetA"

        with ensure_clean(ext) as pth:
            df.to_excel(pth, sheet_name)
            book = xlrd.open_workbook(pth)

            with ExcelFile(book, engine=engine) as xl:
                result = read_excel(xl, sheet_name, index_col=0)
                tm.assert_frame_equal(df, result)

            result = read_excel(book, sheet_name=sheet_name,
                                engine=engine, index_col=0)
            tm.assert_frame_equal(df, result)


class _XlrdWriterBase(XlrdSharedItems):

    @pytest.fixture(autouse=True)
    def set_engine_and_path(self, request, merge_cells, engine, ext):
        """Fixture to set engine and open file for use in each test case

        Rather than requiring `engine=...` to be provided explicitly as an
        argument in each test, this fixture sets a global option to dictate
        which engine should be used to write Excel files. After executing
        the test it rolls back said change to the global option.

        It also uses a context manager to open a temporary excel file for
        the function to write to, accessible via `self.path`

        Notes
        -----
        This fixture will run as part of each test method defined in the
        class and any subclasses, on account of the `autouse=True`
        argument
        """
        option_name = 'io.excel.{ext}.writer'.format(ext=ext.strip('.'))
        prev_engine = get_option(option_name)
        set_option(option_name, engine)
        with ensure_clean(ext) as path:
            self.path = path
            yield
        set_option(option_name, prev_engine)  # Roll back option change


@pytest.mark.parametrize("merge_cells", [True, False])
@pytest.mark.parametrize("engine,ext", [
    pytest.param('openpyxl', '.xlsx', marks=pytest.mark.skipif(
        not td.safe_import('openpyxl'), reason='No openpyxl')),
    pytest.param('openpyxl', '.xlsm', marks=pytest.mark.skipif(
        not td.safe_import('openpyxl'), reason='No openpyxl')),
    pytest.param('xlwt', '.xls', marks=pytest.mark.skipif(
        not td.safe_import('xlwt'), reason='No xlwt')),
    pytest.param('xlsxwriter', '.xlsx', marks=pytest.mark.skipif(
        not td.safe_import('xlsxwriter'), reason='No xlsxwriter'))
])
class XlrdTestExcelWriter(_XlrdWriterBase):
    # Base class for test cases to run with different Excel writers.

    def test_excel_sheet_by_name_raise(self, *_):
        import xlrd

        gt = DataFrame(np.random.randn(10, 2))
        gt.to_excel(self.path)

        xl = ExcelFile(self.path)
        df = read_excel(xl, 0, index_col=0)

        tm.assert_frame_equal(gt, df)

        with pytest.raises(xlrd.XLRDError):
            read_excel(xl, "0")
