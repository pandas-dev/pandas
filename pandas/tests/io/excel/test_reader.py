from collections import OrderedDict
import contextlib
from datetime import datetime, time
from functools import partial
import os
import warnings

import numpy as np
import pytest

from pandas.compat import iteritems, range
import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series
import pandas.util.testing as tm
from pandas.util.testing import ensure_clean, makeCustomDataframe as mkdf

from pandas.io.common import URLError
from pandas.io.excel import ExcelFile, ExcelWriter, read_excel
from pandas.io.parsers import read_csv

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
class SharedItems(object):

    @pytest.fixture(autouse=True)
    def setup_method(self, datapath):
        self.dirpath = datapath("io", "data")
        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()

    def get_csv_refdf(self, basename):
        """
        Obtain the reference data from read_csv with the Python engine.
        Parameters
        ----------
        basename : str
            File base name, excluding file extension.
        Returns
        -------
        dfref : DataFrame
        """
        pref = os.path.join(self.dirpath, basename + '.csv')
        dfref = read_csv(pref, index_col=0, parse_dates=True, engine='python')
        return dfref

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


class ReadingTestsBase(SharedItems):
    # This is based on ExcelWriterBase

    @pytest.fixture(autouse=True, params=['xlrd', None])
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

    @pytest.mark.parametrize("usecols", [
        [0, 1, 3], [0, 3, 1],
        [1, 0, 3], [1, 3, 0],
        [3, 0, 1], [3, 1, 0],
    ])
    def test_usecols_diff_positional_int_columns_order(self, ext, usecols):
        expected = self.get_csv_refdf("test1")[["A", "C"]]
        result = self.get_exceldf("test1", ext, "Sheet1",
                                  index_col=0, usecols=usecols)
        tm.assert_frame_equal(result, expected, check_names=False)

    @pytest.mark.parametrize("usecols", [
        ["B", "D"], ["D", "B"]
    ])
    def test_usecols_diff_positional_str_columns_order(self, ext, usecols):
        expected = self.get_csv_refdf("test1")[["B", "D"]]
        expected.index = range(len(expected))

        result = self.get_exceldf("test1", ext, "Sheet1", usecols=usecols)
        tm.assert_frame_equal(result, expected, check_names=False)

    def test_read_excel_without_slicing(self, ext):
        expected = self.get_csv_refdf("test1")
        result = self.get_exceldf("test1", ext, "Sheet1", index_col=0)
        tm.assert_frame_equal(result, expected, check_names=False)

    def test_usecols_excel_range_str(self, ext):
        expected = self.get_csv_refdf("test1")[["C", "D"]]
        result = self.get_exceldf("test1", ext, "Sheet1",
                                  index_col=0, usecols="A,D:E")
        tm.assert_frame_equal(result, expected, check_names=False)

    def test_usecols_excel_range_str_invalid(self, ext):
        msg = "Invalid column name: E1"

        with pytest.raises(ValueError, match=msg):
            self.get_exceldf("test1", ext, "Sheet1", usecols="D:E1")

    def test_index_col_label_error(self, ext):
        msg = "list indices must be integers.*, not str"

        with pytest.raises(TypeError, match=msg):
            self.get_exceldf("test1", ext, "Sheet1", index_col=["A"],
                             usecols=["A", "C"])

    def test_index_col_empty(self, ext):
        # see gh-9208
        result = self.get_exceldf("test1", ext, "Sheet3",
                                  index_col=["A", "B", "C"])
        expected = DataFrame(columns=["D", "E", "F"],
                             index=MultiIndex(levels=[[]] * 3,
                                              codes=[[]] * 3,
                                              names=["A", "B", "C"]))
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("index_col", [None, 2])
    def test_index_col_with_unnamed(self, ext, index_col):
        # see gh-18792
        result = self.get_exceldf("test1", ext, "Sheet4",
                                  index_col=index_col)
        expected = DataFrame([["i1", "a", "x"], ["i2", "b", "y"]],
                             columns=["Unnamed: 0", "col1", "col2"])
        if index_col:
            expected = expected.set_index(expected.columns[index_col])

        tm.assert_frame_equal(result, expected)

    def test_usecols_pass_non_existent_column(self, ext):
        msg = ("Usecols do not match columns, "
               "columns expected but not found: " + r"\['E'\]")

        with pytest.raises(ValueError, match=msg):
            self.get_exceldf("test1", ext, usecols=["E"])

    def test_usecols_wrong_type(self, ext):
        msg = ("'usecols' must either be list-like of "
               "all strings, all unicode, all integers or a callable.")

        with pytest.raises(ValueError, match=msg):
            self.get_exceldf("test1", ext, usecols=["E1", 0])

    def test_excel_stop_iterator(self, ext):

        parsed = self.get_exceldf('test2', ext, 'Sheet1')
        expected = DataFrame([['aaaa', 'bbbbb']], columns=['Test', 'Test1'])
        tm.assert_frame_equal(parsed, expected)

    def test_excel_cell_error_na(self, ext):

        parsed = self.get_exceldf('test3', ext, 'Sheet1')
        expected = DataFrame([[np.nan]], columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

    def test_excel_passes_na(self, ext):

        excel = self.get_excelfile('test4', ext)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=False,
                            na_values=['apple'])
        expected = DataFrame([['NA'], [1], ['NA'], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=True,
                            na_values=['apple'])
        expected = DataFrame([[np.nan], [1], [np.nan], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

        # 13967
        excel = self.get_excelfile('test5', ext)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=False,
                            na_values=['apple'])
        expected = DataFrame([['1.#QNAN'], [1], ['nan'], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=True,
                            na_values=['apple'])
        expected = DataFrame([[np.nan], [1], [np.nan], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

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

    def test_excel_table(self, ext):

        dfref = self.get_csv_refdf('test1')

        df1 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0)
        df2 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0)
        # TODO add index to file
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)

        df3 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               skipfooter=1)
        tm.assert_frame_equal(df3, df1.iloc[:-1])

    def test_reader_special_dtypes(self, ext):

        expected = DataFrame.from_dict(OrderedDict([
            ("IntCol", [1, 2, -3, 4, 0]),
            ("FloatCol", [1.25, 2.25, 1.83, 1.92, 0.0000000005]),
            ("BoolCol", [True, False, True, True, False]),
            ("StrCol", [1, 2, 3, 4, 5]),
            # GH5394 - this is why convert_float isn't vectorized
            ("Str2Col", ["a", 3, "c", "d", "e"]),
            ("DateCol", [datetime(2013, 10, 30), datetime(2013, 10, 31),
                         datetime(1905, 1, 1), datetime(2013, 12, 14),
                         datetime(2015, 3, 14)])
        ]))
        basename = 'test_types'

        # should read in correctly and infer types
        actual = self.get_exceldf(basename, ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

        # if not coercing number, then int comes in as float
        float_expected = expected.copy()
        float_expected["IntCol"] = float_expected["IntCol"].astype(float)
        float_expected.loc[float_expected.index[1], "Str2Col"] = 3.0
        actual = self.get_exceldf(basename, ext, 'Sheet1', convert_float=False)
        tm.assert_frame_equal(actual, float_expected)

        # check setting Index (assuming xls and xlsx are the same here)
        for icol, name in enumerate(expected.columns):
            actual = self.get_exceldf(basename, ext, 'Sheet1', index_col=icol)
            exp = expected.set_index(name)
            tm.assert_frame_equal(actual, exp)

        # convert_float and converters should be different but both accepted
        expected["StrCol"] = expected["StrCol"].apply(str)
        actual = self.get_exceldf(
            basename, ext, 'Sheet1', converters={"StrCol": str})
        tm.assert_frame_equal(actual, expected)

        no_convert_float = float_expected.copy()
        no_convert_float["StrCol"] = no_convert_float["StrCol"].apply(str)
        actual = self.get_exceldf(basename, ext, 'Sheet1', convert_float=False,
                                  converters={"StrCol": str})
        tm.assert_frame_equal(actual, no_convert_float)

    # GH8212 - support for converters and missing values
    def test_reader_converters(self, ext):

        basename = 'test_converters'

        expected = DataFrame.from_dict(OrderedDict([
            ("IntCol", [1, 2, -3, -1000, 0]),
            ("FloatCol", [12.5, np.nan, 18.3, 19.2, 0.000000005]),
            ("BoolCol", ['Found', 'Found', 'Found', 'Not found', 'Found']),
            ("StrCol", ['1', np.nan, '3', '4', '5']),
        ]))

        converters = {'IntCol': lambda x: int(x) if x != '' else -1000,
                      'FloatCol': lambda x: 10 * x if x else np.nan,
                      2: lambda x: 'Found' if x != '' else 'Not found',
                      3: lambda x: str(x) if x else '',
                      }

        # should read in correctly and set types of single cells (not array
        # dtypes)
        actual = self.get_exceldf(basename, ext, 'Sheet1',
                                  converters=converters)
        tm.assert_frame_equal(actual, expected)

    def test_reader_dtype(self, ext):
        # GH 8212
        basename = 'testdtype'
        actual = self.get_exceldf(basename, ext)

        expected = DataFrame({
            'a': [1, 2, 3, 4],
            'b': [2.5, 3.5, 4.5, 5.5],
            'c': [1, 2, 3, 4],
            'd': [1.0, 2.0, np.nan, 4.0]}).reindex(
                columns=['a', 'b', 'c', 'd'])

        tm.assert_frame_equal(actual, expected)

        actual = self.get_exceldf(basename, ext,
                                  dtype={'a': 'float64',
                                         'b': 'float32',
                                         'c': str})

        expected['a'] = expected['a'].astype('float64')
        expected['b'] = expected['b'].astype('float32')
        expected['c'] = ['001', '002', '003', '004']
        tm.assert_frame_equal(actual, expected)

        with pytest.raises(ValueError):
            self.get_exceldf(basename, ext, dtype={'d': 'int64'})

    @pytest.mark.parametrize("dtype,expected", [
        (None,
         DataFrame({
             "a": [1, 2, 3, 4],
             "b": [2.5, 3.5, 4.5, 5.5],
             "c": [1, 2, 3, 4],
             "d": [1.0, 2.0, np.nan, 4.0]
         })),
        ({"a": "float64",
          "b": "float32",
          "c": str,
          "d": str
          },
         DataFrame({
             "a": Series([1, 2, 3, 4], dtype="float64"),
             "b": Series([2.5, 3.5, 4.5, 5.5], dtype="float32"),
             "c": ["001", "002", "003", "004"],
             "d": ["1", "2", np.nan, "4"]
         })),
    ])
    def test_reader_dtype_str(self, ext, dtype, expected):
        # see gh-20377
        basename = "testdtype"

        actual = self.get_exceldf(basename, ext, dtype=dtype)
        tm.assert_frame_equal(actual, expected)

    def test_reading_all_sheets(self, ext):
        # Test reading all sheetnames by setting sheetname to None,
        # Ensure a dict is returned.
        # See PR #9450
        basename = 'test_multisheet'
        dfs = self.get_exceldf(basename, ext, sheet_name=None)
        # ensure this is not alphabetical to test order preservation
        expected_keys = ['Charlie', 'Alpha', 'Beta']
        tm.assert_contains_all(expected_keys, dfs.keys())
        # Issue 9930
        # Ensure sheet order is preserved
        assert expected_keys == list(dfs.keys())

    def test_reading_multiple_specific_sheets(self, ext):
        # Test reading specific sheetnames by specifying a mixed list
        # of integers and strings, and confirm that duplicated sheet
        # references (positions/names) are removed properly.
        # Ensure a dict is returned
        # See PR #9450
        basename = 'test_multisheet'
        # Explicitly request duplicates. Only the set should be returned.
        expected_keys = [2, 'Charlie', 'Charlie']
        dfs = self.get_exceldf(basename, ext, sheet_name=expected_keys)
        expected_keys = list(set(expected_keys))
        tm.assert_contains_all(expected_keys, dfs.keys())
        assert len(expected_keys) == len(dfs.keys())

    def test_reading_all_sheets_with_blank(self, ext):
        # Test reading all sheetnames by setting sheetname to None,
        # In the case where some sheets are blank.
        # Issue #11711
        basename = 'blank_with_header'
        dfs = self.get_exceldf(basename, ext, sheet_name=None)
        expected_keys = ['Sheet1', 'Sheet2', 'Sheet3']
        tm.assert_contains_all(expected_keys, dfs.keys())

    # GH6403
    def test_read_excel_blank(self, ext):
        actual = self.get_exceldf('blank', ext, 'Sheet1')
        tm.assert_frame_equal(actual, DataFrame())

    def test_read_excel_blank_with_header(self, ext):
        expected = DataFrame(columns=['col_1', 'col_2'])
        actual = self.get_exceldf('blank_with_header', ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

    @td.skip_if_no("xlwt")
    @td.skip_if_no("openpyxl")
    @pytest.mark.parametrize("header,expected", [
        (None, DataFrame([np.nan] * 4)),
        (0, DataFrame({"Unnamed: 0": [np.nan] * 3}))
    ])
    def test_read_one_empty_col_no_header(self, ext, header, expected):
        # xref gh-12292
        filename = "no_header"
        df = pd.DataFrame(
            [["", 1, 100],
             ["", 2, 200],
             ["", 3, 300],
             ["", 4, 400]]
        )

        with ensure_clean(ext) as path:
            df.to_excel(path, filename, index=False, header=False)
            result = read_excel(path, filename, usecols=[0], header=header)

        tm.assert_frame_equal(result, expected)

    @td.skip_if_no("xlwt")
    @td.skip_if_no("openpyxl")
    @pytest.mark.parametrize("header,expected", [
        (None, DataFrame([0] + [np.nan] * 4)),
        (0, DataFrame([np.nan] * 4))
    ])
    def test_read_one_empty_col_with_header(self, ext, header, expected):
        filename = "with_header"
        df = pd.DataFrame(
            [["", 1, 100],
             ["", 2, 200],
             ["", 3, 300],
             ["", 4, 400]]
        )

        with ensure_clean(ext) as path:
            df.to_excel(path, 'with_header', index=False, header=True)
            result = read_excel(path, filename, usecols=[0], header=header)

        tm.assert_frame_equal(result, expected)

    @td.skip_if_no('openpyxl')
    @td.skip_if_no('xlwt')
    def test_set_column_names_in_parameter(self, ext):
        # GH 12870 : pass down column names associated with
        # keyword argument names
        refdf = pd.DataFrame([[1, 'foo'], [2, 'bar'],
                              [3, 'baz']], columns=['a', 'b'])

        with ensure_clean(ext) as pth:
            with ExcelWriter(pth) as writer:
                refdf.to_excel(writer, 'Data_no_head',
                               header=False, index=False)
                refdf.to_excel(writer, 'Data_with_head', index=False)

            refdf.columns = ['A', 'B']

            with ExcelFile(pth) as reader:
                xlsdf_no_head = read_excel(reader, 'Data_no_head',
                                           header=None, names=['A', 'B'])
                xlsdf_with_head = read_excel(reader, 'Data_with_head',
                                             index_col=None, names=['A', 'B'])

            tm.assert_frame_equal(xlsdf_no_head, refdf)
            tm.assert_frame_equal(xlsdf_with_head, refdf)

    def test_date_conversion_overflow(self, ext):
        # GH 10001 : pandas.ExcelFile ignore parse_dates=False
        expected = pd.DataFrame([[pd.Timestamp('2016-03-12'), 'Marc Johnson'],
                                 [pd.Timestamp('2016-03-16'), 'Jack Black'],
                                 [1e+20, 'Timothy Brown']],
                                columns=['DateColWithBigInt', 'StringCol'])

        result = self.get_exceldf('testdateoverflow', ext)
        tm.assert_frame_equal(result, expected)

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

    def test_sheet_name_both_raises(self, ext):
        with pytest.raises(TypeError, match="Cannot specify both"):
            self.get_exceldf('test1', ext, sheetname='Sheet1',
                             sheet_name='Sheet1')

        excel = self.get_excelfile('test1', ext)
        with pytest.raises(TypeError, match="Cannot specify both"):
            excel.parse(sheetname='Sheet1',
                        sheet_name='Sheet1')

    def test_excel_read_buffer(self, ext):

        pth = os.path.join(self.dirpath, 'test1' + ext)
        expected = read_excel(pth, 'Sheet1', index_col=0)
        with open(pth, 'rb') as f:
            actual = read_excel(f, 'Sheet1', index_col=0)
            tm.assert_frame_equal(expected, actual)

        with open(pth, 'rb') as f:
            xls = ExcelFile(f)
            actual = read_excel(xls, 'Sheet1', index_col=0)
            tm.assert_frame_equal(expected, actual)

    def test_bad_engine_raises(self, ext):
        bad_engine = 'foo'
        with pytest.raises(ValueError, match="Unknown engine: foo"):
            read_excel('', engine=bad_engine)

    @tm.network
    def test_read_from_http_url(self, ext):
        url = ('https://raw.github.com/pandas-dev/pandas/master/'
               'pandas/tests/io/data/test1' + ext)
        url_table = read_excel(url)
        local_table = self.get_exceldf('test1', ext)
        tm.assert_frame_equal(url_table, local_table)

    @td.skip_if_not_us_locale
    def test_read_from_s3_url(self, ext, s3_resource):
        # Bucket "pandas-test" created in tests/io/conftest.py
        file_name = os.path.join(self.dirpath, 'test1' + ext)

        with open(file_name, "rb") as f:
            s3_resource.Bucket("pandas-test").put_object(Key="test1" + ext,
                                                         Body=f)

        url = ('s3://pandas-test/test1' + ext)
        url_table = read_excel(url)
        local_table = self.get_exceldf('test1', ext)
        tm.assert_frame_equal(url_table, local_table)

    @pytest.mark.slow
    # ignore warning from old xlrd
    @pytest.mark.filterwarnings("ignore:This metho:PendingDeprecationWarning")
    def test_read_from_file_url(self, ext):

        # FILE
        localtable = os.path.join(self.dirpath, 'test1' + ext)
        local_table = read_excel(localtable)

        try:
            url_table = read_excel('file://localhost/' + localtable)
        except URLError:
            # fails on some systems
            import platform
            pytest.skip("failing on %s" %
                        ' '.join(platform.uname()).strip())

        tm.assert_frame_equal(url_table, local_table)

    @td.skip_if_no('pathlib')
    def test_read_from_pathlib_path(self, ext):

        # GH12655
        from pathlib import Path

        str_path = os.path.join(self.dirpath, 'test1' + ext)
        expected = read_excel(str_path, 'Sheet1', index_col=0)

        path_obj = Path(self.dirpath, 'test1' + ext)
        actual = read_excel(path_obj, 'Sheet1', index_col=0)

        tm.assert_frame_equal(expected, actual)

    @td.skip_if_no('py.path')
    def test_read_from_py_localpath(self, ext):

        # GH12655
        from py.path import local as LocalPath

        str_path = os.path.join(self.dirpath, 'test1' + ext)
        expected = read_excel(str_path, 'Sheet1', index_col=0)

        abs_dir = os.path.abspath(self.dirpath)
        path_obj = LocalPath(abs_dir).join('test1' + ext)
        actual = read_excel(path_obj, 'Sheet1', index_col=0)

        tm.assert_frame_equal(expected, actual)

    def test_reader_closes_file(self, ext):

        pth = os.path.join(self.dirpath, 'test1' + ext)
        f = open(pth, 'rb')
        with ExcelFile(f) as xlsx:
            # parses okay
            read_excel(xlsx, 'Sheet1', index_col=0)

        assert f.closed

    @td.skip_if_no("xlwt")
    @td.skip_if_no("openpyxl")
    def test_creating_and_reading_multiple_sheets(self, ext):
        # see gh-9450
        #
        # Test reading multiple sheets, from a runtime
        # created Excel file with multiple sheets.
        def tdf(col_sheet_name):
            d, i = [11, 22, 33], [1, 2, 3]
            return DataFrame(d, i, columns=[col_sheet_name])

        sheets = ["AAA", "BBB", "CCC"]

        dfs = [tdf(s) for s in sheets]
        dfs = dict(zip(sheets, dfs))

        with ensure_clean(ext) as pth:
            with ExcelWriter(pth) as ew:
                for sheetname, df in iteritems(dfs):
                    df.to_excel(ew, sheetname)

            dfs_returned = read_excel(pth, sheet_name=sheets, index_col=0)

            for s in sheets:
                tm.assert_frame_equal(dfs[s], dfs_returned[s])

    def test_reader_seconds(self, ext):

        # Test reading times with and without milliseconds. GH5945.
        expected = DataFrame.from_dict({"Time": [time(1, 2, 3),
                                                 time(2, 45, 56, 100000),
                                                 time(4, 29, 49, 200000),
                                                 time(6, 13, 42, 300000),
                                                 time(7, 57, 35, 400000),
                                                 time(9, 41, 28, 500000),
                                                 time(11, 25, 21, 600000),
                                                 time(13, 9, 14, 700000),
                                                 time(14, 53, 7, 800000),
                                                 time(16, 37, 0, 900000),
                                                 time(18, 20, 54)]})

        actual = self.get_exceldf('times_1900', ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

        actual = self.get_exceldf('times_1904', ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_multiindex(self, ext):
        # see gh-4679
        mi = MultiIndex.from_product([["foo", "bar"], ["a", "b"]])
        mi_file = os.path.join(self.dirpath, "testmultiindex" + ext)

        # "mi_column" sheet
        expected = DataFrame([[1, 2.5, pd.Timestamp("2015-01-01"), True],
                              [2, 3.5, pd.Timestamp("2015-01-02"), False],
                              [3, 4.5, pd.Timestamp("2015-01-03"), False],
                              [4, 5.5, pd.Timestamp("2015-01-04"), True]],
                             columns=mi)

        actual = read_excel(mi_file, "mi_column", header=[0, 1], index_col=0)
        tm.assert_frame_equal(actual, expected)

        # "mi_index" sheet
        expected.index = mi
        expected.columns = ["a", "b", "c", "d"]

        actual = read_excel(mi_file, "mi_index", index_col=[0, 1])
        tm.assert_frame_equal(actual, expected, check_names=False)

        # "both" sheet
        expected.columns = mi

        actual = read_excel(mi_file, "both", index_col=[0, 1], header=[0, 1])
        tm.assert_frame_equal(actual, expected, check_names=False)

        # "mi_index_name" sheet
        expected.columns = ["a", "b", "c", "d"]
        expected.index = mi.set_names(["ilvl1", "ilvl2"])

        actual = read_excel(mi_file, "mi_index_name", index_col=[0, 1])
        tm.assert_frame_equal(actual, expected)

        # "mi_column_name" sheet
        expected.index = list(range(4))
        expected.columns = mi.set_names(["c1", "c2"])
        actual = read_excel(mi_file, "mi_column_name",
                            header=[0, 1], index_col=0)
        tm.assert_frame_equal(actual, expected)

        # see gh-11317
        # "name_with_int" sheet
        expected.columns = mi.set_levels(
            [1, 2], level=1).set_names(["c1", "c2"])

        actual = read_excel(mi_file, "name_with_int",
                            index_col=0, header=[0, 1])
        tm.assert_frame_equal(actual, expected)

        # "both_name" sheet
        expected.columns = mi.set_names(["c1", "c2"])
        expected.index = mi.set_names(["ilvl1", "ilvl2"])

        actual = read_excel(mi_file, "both_name",
                            index_col=[0, 1], header=[0, 1])
        tm.assert_frame_equal(actual, expected)

        # "both_skiprows" sheet
        actual = read_excel(mi_file, "both_name_skiprows", index_col=[0, 1],
                            header=[0, 1], skiprows=2)
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_multiindex_header_only(self, ext):
        # see gh-11733.
        #
        # Don't try to parse a header name if there isn't one.
        mi_file = os.path.join(self.dirpath, "testmultiindex" + ext)
        result = read_excel(mi_file, "index_col_none", header=[0, 1])

        exp_columns = MultiIndex.from_product([("A", "B"), ("key", "val")])
        expected = DataFrame([[1, 2, 3, 4]] * 2, columns=exp_columns)
        tm.assert_frame_equal(result, expected)

    @td.skip_if_no("xlsxwriter")
    def test_read_excel_multiindex_empty_level(self, ext):
        # see gh-12453
        with ensure_clean(ext) as path:
            df = DataFrame({
                ("One", "x"): {0: 1},
                ("Two", "X"): {0: 3},
                ("Two", "Y"): {0: 7},
                ("Zero", ""): {0: 0}
            })

            expected = DataFrame({
                ("One", "x"): {0: 1},
                ("Two", "X"): {0: 3},
                ("Two", "Y"): {0: 7},
                ("Zero", "Unnamed: 4_level_1"): {0: 0}
            })

            df.to_excel(path)
            actual = pd.read_excel(path, header=[0, 1], index_col=0)
            tm.assert_frame_equal(actual, expected)

            df = pd.DataFrame({
                ("Beg", ""): {0: 0},
                ("Middle", "x"): {0: 1},
                ("Tail", "X"): {0: 3},
                ("Tail", "Y"): {0: 7}
            })

            expected = pd.DataFrame({
                ("Beg", "Unnamed: 1_level_1"): {0: 0},
                ("Middle", "x"): {0: 1},
                ("Tail", "X"): {0: 3},
                ("Tail", "Y"): {0: 7}
            })

            df.to_excel(path)
            actual = pd.read_excel(path, header=[0, 1], index_col=0)
            tm.assert_frame_equal(actual, expected)

    @td.skip_if_no("xlsxwriter")
    @pytest.mark.parametrize("c_idx_names", [True, False])
    @pytest.mark.parametrize("r_idx_names", [True, False])
    @pytest.mark.parametrize("c_idx_levels", [1, 3])
    @pytest.mark.parametrize("r_idx_levels", [1, 3])
    def test_excel_multindex_roundtrip(self, ext, c_idx_names, r_idx_names,
                                       c_idx_levels, r_idx_levels):
        # see gh-4679
        with ensure_clean(ext) as pth:
            if c_idx_levels == 1 and c_idx_names:
                pytest.skip("Column index name cannot be "
                            "serialized unless it's a MultiIndex")

            # Empty name case current read in as
            # unnamed levels, not Nones.
            check_names = r_idx_names or r_idx_levels <= 1

            df = mkdf(5, 5, c_idx_names, r_idx_names,
                      c_idx_levels, r_idx_levels)
            df.to_excel(pth)

            act = pd.read_excel(pth, index_col=list(range(r_idx_levels)),
                                header=list(range(c_idx_levels)))
            tm.assert_frame_equal(df, act, check_names=check_names)

            df.iloc[0, :] = np.nan
            df.to_excel(pth)

            act = pd.read_excel(pth, index_col=list(range(r_idx_levels)),
                                header=list(range(c_idx_levels)))
            tm.assert_frame_equal(df, act, check_names=check_names)

            df.iloc[-1, :] = np.nan
            df.to_excel(pth)
            act = pd.read_excel(pth, index_col=list(range(r_idx_levels)),
                                header=list(range(c_idx_levels)))
            tm.assert_frame_equal(df, act, check_names=check_names)

    def test_excel_old_index_format(self, ext):
        # see gh-4679
        filename = "test_index_name_pre17" + ext
        in_file = os.path.join(self.dirpath, filename)

        # We detect headers to determine if index names exist, so
        # that "index" name in the "names" version of the data will
        # now be interpreted as rows that include null data.
        data = np.array([[None, None, None, None, None],
                         ["R0C0", "R0C1", "R0C2", "R0C3", "R0C4"],
                         ["R1C0", "R1C1", "R1C2", "R1C3", "R1C4"],
                         ["R2C0", "R2C1", "R2C2", "R2C3", "R2C4"],
                         ["R3C0", "R3C1", "R3C2", "R3C3", "R3C4"],
                         ["R4C0", "R4C1", "R4C2", "R4C3", "R4C4"]])
        columns = ["C_l0_g0", "C_l0_g1", "C_l0_g2", "C_l0_g3", "C_l0_g4"]
        mi = MultiIndex(levels=[["R0", "R_l0_g0", "R_l0_g1",
                                 "R_l0_g2", "R_l0_g3", "R_l0_g4"],
                                ["R1", "R_l1_g0", "R_l1_g1",
                                 "R_l1_g2", "R_l1_g3", "R_l1_g4"]],
                        codes=[[0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]],
                        names=[None, None])
        si = Index(["R0", "R_l0_g0", "R_l0_g1", "R_l0_g2",
                    "R_l0_g3", "R_l0_g4"], name=None)

        expected = pd.DataFrame(data, index=si, columns=columns)

        actual = pd.read_excel(in_file, "single_names", index_col=0)
        tm.assert_frame_equal(actual, expected)

        expected.index = mi

        actual = pd.read_excel(in_file, "multi_names", index_col=[0, 1])
        tm.assert_frame_equal(actual, expected)

        # The analogous versions of the "names" version data
        # where there are explicitly no names for the indices.
        data = np.array([["R0C0", "R0C1", "R0C2", "R0C3", "R0C4"],
                         ["R1C0", "R1C1", "R1C2", "R1C3", "R1C4"],
                         ["R2C0", "R2C1", "R2C2", "R2C3", "R2C4"],
                         ["R3C0", "R3C1", "R3C2", "R3C3", "R3C4"],
                         ["R4C0", "R4C1", "R4C2", "R4C3", "R4C4"]])
        columns = ["C_l0_g0", "C_l0_g1", "C_l0_g2", "C_l0_g3", "C_l0_g4"]
        mi = MultiIndex(levels=[["R_l0_g0", "R_l0_g1", "R_l0_g2",
                                 "R_l0_g3", "R_l0_g4"],
                                ["R_l1_g0", "R_l1_g1", "R_l1_g2",
                                 "R_l1_g3", "R_l1_g4"]],
                        codes=[[0, 1, 2, 3, 4], [0, 1, 2, 3, 4]],
                        names=[None, None])
        si = Index(["R_l0_g0", "R_l0_g1", "R_l0_g2",
                    "R_l0_g3", "R_l0_g4"], name=None)

        expected = pd.DataFrame(data, index=si, columns=columns)

        actual = pd.read_excel(in_file, "single_no_names", index_col=0)
        tm.assert_frame_equal(actual, expected)

        expected.index = mi

        actual = pd.read_excel(in_file, "multi_no_names", index_col=[0, 1])
        tm.assert_frame_equal(actual, expected, check_names=False)

    def test_read_excel_bool_header_arg(self, ext):
        # GH 6114
        for arg in [True, False]:
            with pytest.raises(TypeError):
                pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                              header=arg)

    def test_read_excel_chunksize(self, ext):
        # GH 8011
        with pytest.raises(NotImplementedError):
            pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                          chunksize=100)

    @td.skip_if_no("xlwt")
    @td.skip_if_no("openpyxl")
    def test_read_excel_parse_dates(self, ext):
        # see gh-11544, gh-12051
        df = DataFrame(
            {"col": [1, 2, 3],
             "date_strings": pd.date_range("2012-01-01", periods=3)})
        df2 = df.copy()
        df2["date_strings"] = df2["date_strings"].dt.strftime("%m/%d/%Y")

        with ensure_clean(ext) as pth:
            df2.to_excel(pth)

            res = read_excel(pth, index_col=0)
            tm.assert_frame_equal(df2, res)

            res = read_excel(pth, parse_dates=["date_strings"], index_col=0)
            tm.assert_frame_equal(df, res)

            date_parser = lambda x: pd.datetime.strptime(x, "%m/%d/%Y")
            res = read_excel(pth, parse_dates=["date_strings"],
                             date_parser=date_parser, index_col=0)
            tm.assert_frame_equal(df, res)

    def test_read_excel_skiprows_list(self, ext):
        # GH 4903
        actual = pd.read_excel(os.path.join(self.dirpath,
                                            'testskiprows' + ext),
                               'skiprows_list', skiprows=[0, 2])
        expected = DataFrame([[1, 2.5, pd.Timestamp('2015-01-01'), True],
                              [2, 3.5, pd.Timestamp('2015-01-02'), False],
                              [3, 4.5, pd.Timestamp('2015-01-03'), False],
                              [4, 5.5, pd.Timestamp('2015-01-04'), True]],
                             columns=['a', 'b', 'c', 'd'])
        tm.assert_frame_equal(actual, expected)

        actual = pd.read_excel(os.path.join(self.dirpath,
                                            'testskiprows' + ext),
                               'skiprows_list', skiprows=np.array([0, 2]))
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_nrows(self, ext):
        # GH 16645
        num_rows_to_pull = 5
        actual = pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                               nrows=num_rows_to_pull)
        expected = pd.read_excel(os.path.join(self.dirpath,
                                              'test1' + ext))
        expected = expected[:num_rows_to_pull]
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_nrows_greater_than_nrows_in_file(self, ext):
        # GH 16645
        expected = pd.read_excel(os.path.join(self.dirpath,
                                              'test1' + ext))
        num_records_in_file = len(expected)
        num_rows_to_pull = num_records_in_file + 10
        actual = pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                               nrows=num_rows_to_pull)
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_nrows_non_integer_parameter(self, ext):
        # GH 16645
        msg = "'nrows' must be an integer >=0"
        with pytest.raises(ValueError, match=msg):
            pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                          nrows='5')

    def test_read_excel_squeeze(self, ext):
        # GH 12157
        f = os.path.join(self.dirpath, 'test_squeeze' + ext)

        actual = pd.read_excel(f, 'two_columns', index_col=0, squeeze=True)
        expected = pd.Series([2, 3, 4], [4, 5, 6], name='b')
        expected.index.name = 'a'
        tm.assert_series_equal(actual, expected)

        actual = pd.read_excel(f, 'two_columns', squeeze=True)
        expected = pd.DataFrame({'a': [4, 5, 6],
                                 'b': [2, 3, 4]})
        tm.assert_frame_equal(actual, expected)

        actual = pd.read_excel(f, 'one_column', squeeze=True)
        expected = pd.Series([1, 2, 3], name='a')
        tm.assert_series_equal(actual, expected)


@pytest.mark.parametrize("ext", ['.xls', '.xlsx', '.xlsm'])
class TestXlrdReader(ReadingTestsBase):
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
