# pylint: disable=E1101

from pandas.compat import u, range, map
import os
import unittest

import nose

from numpy import nan
import numpy as np

from pandas import DataFrame, Index, MultiIndex
from pandas.io.parsers import read_csv
from pandas.io.excel import (
    ExcelFile, ExcelWriter, read_excel, _XlwtWriter, _OpenpyxlWriter,
    register_writer
)
from pandas.util.testing import ensure_clean
from pandas.core.config import set_option, get_option
import pandas.util.testing as tm
import pandas as pd


def _skip_if_no_xlrd():
    try:
        import xlrd
        ver = tuple(map(int, xlrd.__VERSION__.split(".")[:2]))
        if ver < (0, 9):
            raise nose.SkipTest('xlrd < 0.9, skipping')
    except ImportError:
        raise nose.SkipTest('xlrd not installed, skipping')


def _skip_if_no_xlwt():
    try:
        import xlwt  # NOQA
    except ImportError:
        raise nose.SkipTest('xlwt not installed, skipping')


def _skip_if_no_openpyxl():
    try:
        import openpyxl  # NOQA
    except ImportError:
        raise nose.SkipTest('openpyxl not installed, skipping')


def _skip_if_no_xlsxwriter():
    try:
        import xlsxwriter  # NOQA
    except ImportError:
        raise nose.SkipTest('xlsxwriter not installed, skipping')


def _skip_if_no_excelsuite():
    _skip_if_no_xlrd()
    _skip_if_no_xlwt()
    _skip_if_no_openpyxl()


_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()
_frame = DataFrame(_seriesd)[:10]
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])[:10]
_tsframe = tm.makeTimeDataFrame()[:5]
_mixed_frame = _frame.copy()
_mixed_frame['foo'] = 'bar'


class SharedItems(object):
    def setUp(self):
        self.dirpath = tm.get_data_path()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')
        self.xlsx1 = os.path.join(self.dirpath, 'test.xlsx')
        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = 'python'
        return read_csv(*args, **kwds)


class ExcelReaderTests(SharedItems, unittest.TestCase):
    def test_parse_cols_int(self):
        _skip_if_no_openpyxl()
        _skip_if_no_xlrd()

        suffix = ['', 'x']

        for s in suffix:
            pth = os.path.join(self.dirpath, 'test.xls%s' % s)
            xls = ExcelFile(pth)
            df = xls.parse('Sheet1', index_col=0, parse_dates=True,
                           parse_cols=3)
            df2 = self.read_csv(self.csv1, index_col=0, parse_dates=True)
            df2 = df2.reindex(columns=['A', 'B', 'C'])
            df3 = xls.parse('Sheet2', skiprows=[1], index_col=0,
                            parse_dates=True, parse_cols=3)
            tm.assert_frame_equal(df, df2, check_names=False)  # TODO add index to xls file)
            tm.assert_frame_equal(df3, df2, check_names=False)

    def test_parse_cols_list(self):
        _skip_if_no_openpyxl()
        _skip_if_no_xlrd()

        suffix = ['', 'x']

        for s in suffix:
            pth = os.path.join(self.dirpath, 'test.xls%s' % s)
            xls = ExcelFile(pth)
            df = xls.parse('Sheet1', index_col=0, parse_dates=True,
                           parse_cols=[0, 2, 3])
            df2 = self.read_csv(self.csv1, index_col=0, parse_dates=True)
            df2 = df2.reindex(columns=['B', 'C'])
            df3 = xls.parse('Sheet2', skiprows=[1], index_col=0,
                            parse_dates=True,
                            parse_cols=[0, 2, 3])
            tm.assert_frame_equal(df, df2, check_names=False)  # TODO add index to xls file
            tm.assert_frame_equal(df3, df2, check_names=False)

    def test_parse_cols_str(self):
        _skip_if_no_openpyxl()
        _skip_if_no_xlrd()

        suffix = ['', 'x']

        for s in suffix:

            pth = os.path.join(self.dirpath, 'test.xls%s' % s)
            xls = ExcelFile(pth)

            df = xls.parse('Sheet1', index_col=0, parse_dates=True,
                           parse_cols='A:D')
            df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
            df2 = df2.reindex(columns=['A', 'B', 'C'])
            df3 = xls.parse('Sheet2', skiprows=[1], index_col=0,
                            parse_dates=True, parse_cols='A:D')
            tm.assert_frame_equal(df, df2, check_names=False)  # TODO add index to xls, read xls ignores index name ?
            tm.assert_frame_equal(df3, df2, check_names=False)
            del df, df2, df3

            df = xls.parse('Sheet1', index_col=0, parse_dates=True,
                           parse_cols='A,C,D')
            df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
            df2 = df2.reindex(columns=['B', 'C'])
            df3 = xls.parse('Sheet2', skiprows=[1], index_col=0,
                            parse_dates=True,
                            parse_cols='A,C,D')
            tm.assert_frame_equal(df, df2, check_names=False)  # TODO add index to xls file
            tm.assert_frame_equal(df3, df2, check_names=False)
            del df, df2, df3

            df = xls.parse('Sheet1', index_col=0, parse_dates=True,
                           parse_cols='A,C:D')
            df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
            df2 = df2.reindex(columns=['B', 'C'])
            df3 = xls.parse('Sheet2', skiprows=[1], index_col=0,
                            parse_dates=True,
                            parse_cols='A,C:D')
            tm.assert_frame_equal(df, df2, check_names=False)
            tm.assert_frame_equal(df3, df2, check_names=False)

    def test_excel_stop_iterator(self):
        _skip_if_no_xlrd()

        excel_data = ExcelFile(os.path.join(self.dirpath, 'test2.xls'))
        parsed = excel_data.parse('Sheet1')
        expected = DataFrame([['aaaa', 'bbbbb']], columns=['Test', 'Test1'])
        tm.assert_frame_equal(parsed, expected)

    def test_excel_cell_error_na(self):
        _skip_if_no_xlrd()

        excel_data = ExcelFile(os.path.join(self.dirpath, 'test3.xls'))
        parsed = excel_data.parse('Sheet1')
        expected = DataFrame([[np.nan]], columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

    def test_excel_passes_na(self):
        _skip_if_no_xlrd()

        excel_data = ExcelFile(os.path.join(self.dirpath, 'test2.xlsx'))
        parsed = excel_data.parse('Sheet1', keep_default_na=False,
                                  na_values=['apple'])
        expected = DataFrame([['NA'], [1], ['NA'], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

        parsed = excel_data.parse('Sheet1', keep_default_na=True,
                                  na_values=['apple'])
        expected = DataFrame([[np.nan], [1], [np.nan], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

    def check_excel_table_sheet_by_index(self, filename, csvfile):
        import xlrd

        pth = os.path.join(self.dirpath, filename)
        xls = ExcelFile(pth)
        df = xls.parse(0, index_col=0, parse_dates=True)
        df2 = self.read_csv(csvfile, index_col=0, parse_dates=True)
        df3 = xls.parse(1, skiprows=[1], index_col=0, parse_dates=True)
        tm.assert_frame_equal(df, df2, check_names=False)
        tm.assert_frame_equal(df3, df2, check_names=False)

        df4 = xls.parse(0, index_col=0, parse_dates=True, skipfooter=1)
        df5 = xls.parse(0, index_col=0, parse_dates=True, skip_footer=1)
        tm.assert_frame_equal(df4, df.ix[:-1])
        tm.assert_frame_equal(df4, df5)

        self.assertRaises(xlrd.XLRDError, xls.parse, 'asdf')

    def test_excel_table_sheet_by_index(self):
        _skip_if_no_xlrd()
        for filename, csvfile in [(self.xls1, self.csv1),
                                  (self.xlsx1, self.csv1)]:
            self.check_excel_table_sheet_by_index(filename, csvfile)


    def test_excel_table(self):
        _skip_if_no_xlrd()

        pth = os.path.join(self.dirpath, 'test.xls')
        xls = ExcelFile(pth)
        df = xls.parse('Sheet1', index_col=0, parse_dates=True)
        df2 = self.read_csv(self.csv1, index_col=0, parse_dates=True)
        df3 = xls.parse('Sheet2', skiprows=[1], index_col=0, parse_dates=True)
        tm.assert_frame_equal(df, df2, check_names=False)
        tm.assert_frame_equal(df3, df2, check_names=False)

        df4 = xls.parse('Sheet1', index_col=0, parse_dates=True,
                        skipfooter=1)
        df5 = xls.parse('Sheet1', index_col=0, parse_dates=True,
                        skip_footer=1)
        tm.assert_frame_equal(df4, df.ix[:-1])
        tm.assert_frame_equal(df4, df5)

    def test_excel_read_buffer(self):
        _skip_if_no_xlrd()
        _skip_if_no_openpyxl()

        pth = os.path.join(self.dirpath, 'test.xls')
        f = open(pth, 'rb')
        xls = ExcelFile(f)
        # it works
        xls.parse('Sheet1', index_col=0, parse_dates=True)

        pth = os.path.join(self.dirpath, 'test.xlsx')
        f = open(pth, 'rb')
        xl = ExcelFile(f)
        xl.parse('Sheet1', index_col=0, parse_dates=True)

    def test_read_xlrd_Book(self):
        _skip_if_no_xlrd()
        _skip_if_no_xlwt()

        import xlrd

        pth = '__tmp_excel_read_worksheet__.xls'
        df = self.frame

        with ensure_clean(pth) as pth:
            df.to_excel(pth, "SheetA")
            book = xlrd.open_workbook(pth)

            with ExcelFile(book, engine="xlrd") as xl:
                result = xl.parse("SheetA")
                tm.assert_frame_equal(df, result)

            result = read_excel(book, sheetname="SheetA", engine="xlrd")
            tm.assert_frame_equal(df, result)

    def test_xlsx_table(self):
        _skip_if_no_xlrd()
        _skip_if_no_openpyxl()

        pth = os.path.join(self.dirpath, 'test.xlsx')
        xlsx = ExcelFile(pth)
        df = xlsx.parse('Sheet1', index_col=0, parse_dates=True)
        df2 = self.read_csv(self.csv1, index_col=0, parse_dates=True)
        df3 = xlsx.parse('Sheet2', skiprows=[1], index_col=0, parse_dates=True)

        tm.assert_frame_equal(df, df2, check_names=False)  # TODO add index to xlsx file
        tm.assert_frame_equal(df3, df2, check_names=False)

        df4 = xlsx.parse('Sheet1', index_col=0, parse_dates=True,
                         skipfooter=1)
        df5 = xlsx.parse('Sheet1', index_col=0, parse_dates=True,
                         skip_footer=1)
        tm.assert_frame_equal(df4, df.ix[:-1])
        tm.assert_frame_equal(df4, df5)

    def test_reader_closes_file(self):
        _skip_if_no_xlrd()
        _skip_if_no_openpyxl()

        pth = os.path.join(self.dirpath, 'test.xlsx')
        f = open(pth, 'rb')
        with ExcelFile(f) as xlsx:
            # parses okay
            df = xlsx.parse('Sheet1', index_col=0)

        self.assertTrue(f.closed)


class ExcelWriterBase(SharedItems):
    # Base class for test cases to run with different Excel writers.
    # To add a writer test, define the following:
    # 1. A check_skip function that skips your tests if your writer isn't
    #    installed.
    # 2. Add a property ext, which is the file extension that your writer
    #    writes to.
    # 3. Add a property engine_name, which is the name of the writer class.
    def setUp(self):
        self.check_skip()
        super(ExcelWriterBase, self).setUp()
        self.option_name = 'io.excel.%s.writer' % self.ext
        self.prev_engine = get_option(self.option_name)
        set_option(self.option_name, self.engine_name)

    def tearDown(self):
        set_option(self.option_name, self.prev_engine)

    def test_excel_sheet_by_name_raise(self):
        _skip_if_no_xlrd()
        import xlrd

        ext = self.ext
        pth = os.path.join(self.dirpath, 'testit.{0}'.format(ext))

        with ensure_clean(pth) as pth:
            gt = DataFrame(np.random.randn(10, 2))
            gt.to_excel(pth)
            xl = ExcelFile(pth)
            df = xl.parse(0)
            tm.assert_frame_equal(gt, df)

            self.assertRaises(xlrd.XLRDError, xl.parse, '0')

    def test_excelwriter_contextmanager(self):
        ext = self.ext
        pth = os.path.join(self.dirpath, 'testit.{0}'.format(ext))

        with ensure_clean(pth) as pth:
            with ExcelWriter(pth) as writer:
                self.frame.to_excel(writer, 'Data1')
                self.frame2.to_excel(writer, 'Data2')

            with ExcelFile(pth) as reader:
                found_df = reader.parse('Data1')
                found_df2 = reader.parse('Data2')
                tm.assert_frame_equal(found_df, self.frame)
                tm.assert_frame_equal(found_df2, self.frame2)

    def test_roundtrip(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel__.' + ext

        with ensure_clean(path) as path:
            self.frame['A'][:5] = nan

            self.frame.to_excel(path, 'test1')
            self.frame.to_excel(path, 'test1', cols=['A', 'B'])
            self.frame.to_excel(path, 'test1', header=False)
            self.frame.to_excel(path, 'test1', index=False)

            # test roundtrip
            self.frame.to_excel(path, 'test1')
            recons = read_excel(path, 'test1', index_col=0)
            tm.assert_frame_equal(self.frame, recons)

            self.frame.to_excel(path, 'test1', index=False)
            recons = read_excel(path, 'test1', index_col=None)
            recons.index = self.frame.index
            tm.assert_frame_equal(self.frame, recons)

            self.frame.to_excel(path, 'test1', na_rep='NA')
            recons = read_excel(path, 'test1', index_col=0, na_values=['NA'])
            tm.assert_frame_equal(self.frame, recons)

            # GH 3611
            self.frame.to_excel(path, 'test1', na_rep='88')
            recons = read_excel(path, 'test1', index_col=0, na_values=['88'])
            tm.assert_frame_equal(self.frame, recons)

            self.frame.to_excel(path, 'test1', na_rep='88')
            recons = read_excel(path, 'test1', index_col=0, na_values=[88,88.0])
            tm.assert_frame_equal(self.frame, recons)

    def test_mixed(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_mixed__.' + ext

        with ensure_clean(path) as path:
            self.mixed_frame.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0)
            tm.assert_frame_equal(self.mixed_frame, recons)


    def test_tsframe(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_tsframe__.' + ext

        df = tm.makeTimeDataFrame()[:5]

        with ensure_clean(path) as path:
            df.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1')
            tm.assert_frame_equal(df, recons)

    def test_int64(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_int64__.' + ext

        with ensure_clean(path) as path:
            self.frame['A'][:5] = nan

            self.frame.to_excel(path, 'test1')
            self.frame.to_excel(path, 'test1', cols=['A', 'B'])
            self.frame.to_excel(path, 'test1', header=False)
            self.frame.to_excel(path, 'test1', index=False)

            # Test np.int64, values read come back as float
            frame = DataFrame(np.random.randint(-10, 10, size=(10, 2)), dtype=np.int64)
            frame.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1').astype(np.int64)
            tm.assert_frame_equal(frame, recons, check_dtype=False)

    def test_bool(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_bool__.' + ext

        with ensure_clean(path) as path:
            self.frame['A'][:5] = nan

            self.frame.to_excel(path, 'test1')
            self.frame.to_excel(path, 'test1', cols=['A', 'B'])
            self.frame.to_excel(path, 'test1', header=False)
            self.frame.to_excel(path, 'test1', index=False)

            # Test reading/writing np.bool8, roundtrip only works for xlsx
            frame = (DataFrame(np.random.randn(10, 2)) >= 0)
            frame.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1').astype(np.bool8)
            tm.assert_frame_equal(frame, recons)

    def test_sheets(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_sheets__.' + ext

        with ensure_clean(path) as path:
            self.frame['A'][:5] = nan

            self.frame.to_excel(path, 'test1')
            self.frame.to_excel(path, 'test1', cols=['A', 'B'])
            self.frame.to_excel(path, 'test1', header=False)
            self.frame.to_excel(path, 'test1', index=False)

            # Test writing to separate sheets
            writer = ExcelWriter(path)
            self.frame.to_excel(writer, 'test1')
            self.tsframe.to_excel(writer, 'test2')
            writer.save()
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0)
            tm.assert_frame_equal(self.frame, recons)
            recons = reader.parse('test2', index_col=0)
            tm.assert_frame_equal(self.tsframe, recons)
            np.testing.assert_equal(2, len(reader.sheet_names))
            np.testing.assert_equal('test1', reader.sheet_names[0])
            np.testing.assert_equal('test2', reader.sheet_names[1])

    def test_colaliases(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_aliases__.' + ext

        with ensure_clean(path) as path:
            self.frame['A'][:5] = nan

            self.frame.to_excel(path, 'test1')
            self.frame.to_excel(path, 'test1', cols=['A', 'B'])
            self.frame.to_excel(path, 'test1', header=False)
            self.frame.to_excel(path, 'test1', index=False)

            # column aliases
            col_aliases = Index(['AA', 'X', 'Y', 'Z'])
            self.frame2.to_excel(path, 'test1', header=col_aliases)
            reader = ExcelFile(path)
            rs = reader.parse('test1', index_col=0)
            xp = self.frame2.copy()
            xp.columns = col_aliases
            tm.assert_frame_equal(xp, rs)

    def test_roundtrip_indexlabels(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_indexlabels__.' + ext

        with ensure_clean(path) as path:

            self.frame['A'][:5] = nan

            self.frame.to_excel(path, 'test1')
            self.frame.to_excel(path, 'test1', cols=['A', 'B'])
            self.frame.to_excel(path, 'test1', header=False)
            self.frame.to_excel(path, 'test1', index=False)

            # test index_label
            frame = (DataFrame(np.random.randn(10, 2)) >= 0)
            frame.to_excel(path, 'test1', index_label=['test'])
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0).astype(np.int64)
            frame.index.names = ['test']
            self.assertEqual(frame.index.names, recons.index.names)

            frame = (DataFrame(np.random.randn(10, 2)) >= 0)
            frame.to_excel(
                path, 'test1', index_label=['test', 'dummy', 'dummy2'])
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0).astype(np.int64)
            frame.index.names = ['test']
            self.assertEqual(frame.index.names, recons.index.names)

            frame = (DataFrame(np.random.randn(10, 2)) >= 0)
            frame.to_excel(path, 'test1', index_label='test')
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0).astype(np.int64)
            frame.index.names = ['test']
            self.assertEqual(frame.index.names, recons.index.names)

        # test index_labels in same row as column names
        path = '%s.%s' % (tm.rands(10), ext)

        with ensure_clean(path) as path:

            self.frame.to_excel(path, 'test1',
                                cols=['A', 'B', 'C', 'D'], index=False)
            # take 'A' and 'B' as indexes (they are in same row as cols 'C',
            # 'D')
            df = self.frame.copy()
            df = df.set_index(['A', 'B'])

            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=[0, 1])
            tm.assert_frame_equal(df, recons)

    def test_excel_roundtrip_indexname(self):
        _skip_if_no_xlrd()

        path = '%s.%s' % (tm.rands(10), self.ext)

        df = DataFrame(np.random.randn(10, 4))
        df.index.name = 'foo'

        with ensure_clean(path) as path:
            df.to_excel(path)

            xf = ExcelFile(path)
            result = xf.parse(xf.sheet_names[0], index_col=0)

            tm.assert_frame_equal(result, df)
            self.assertEqual(result.index.name, 'foo')

    def test_excel_roundtrip_datetime(self):
        _skip_if_no_xlrd()

        # datetime.date, not sure what to test here exactly
        path = '__tmp_excel_roundtrip_datetime__.' + self.ext
        tsf = self.tsframe.copy()
        with ensure_clean(path) as path:

            tsf.index = [x.date() for x in self.tsframe.index]
            tsf.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1')
            tm.assert_frame_equal(self.tsframe, recons)

    def test_to_excel_periodindex(self):
        _skip_if_no_xlrd()
        path = '__tmp_to_excel_periodindex__.' + self.ext
        frame = self.tsframe
        xp = frame.resample('M', kind='period')

        with ensure_clean(path) as path:
            xp.to_excel(path, 'sht1')

            reader = ExcelFile(path)
            rs = reader.parse('sht1', index_col=0, parse_dates=True)
            tm.assert_frame_equal(xp, rs.to_period('M'))

    def test_to_excel_multiindex(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_multiindex__' + ext + '__.' + ext

        frame = self.frame
        old_index = frame.index
        arrays = np.arange(len(old_index) * 2).reshape(2, -1)
        new_index = MultiIndex.from_arrays(arrays,
                                           names=['first', 'second'])
        frame.index = new_index

        with ensure_clean(path) as path:
            frame.to_excel(path, 'test1', header=False)
            frame.to_excel(path, 'test1', cols=['A', 'B'])

            # round trip
            frame.to_excel(path, 'test1')
            reader = ExcelFile(path)
            df = reader.parse('test1', index_col=[0, 1], parse_dates=False)
            tm.assert_frame_equal(frame, df)
            self.assertEqual(frame.index.names, df.index.names)
            self.frame.index = old_index  # needed if setUP becomes a classmethod

    def test_to_excel_multiindex_dates(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_multiindex_dates__' + ext + '__.' + ext

        # try multiindex with dates
        tsframe = self.tsframe
        old_index = tsframe.index
        new_index = [old_index, np.arange(len(old_index))]
        tsframe.index = MultiIndex.from_arrays(new_index)

        with ensure_clean(path) as path:
            tsframe.to_excel(path, 'test1', index_label=['time', 'foo'])
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=[0, 1])

            tm.assert_frame_equal(tsframe, recons, check_names=False)
            self.assertEquals(recons.index.names, ('time', 'foo'))

            # infer index
            tsframe.to_excel(path, 'test1')
            reader = ExcelFile(path)
            recons = reader.parse('test1')
            tm.assert_frame_equal(tsframe, recons)

            self.tsframe.index = old_index  # needed if setUP becomes classmethod

    def test_to_excel_float_format(self):
        _skip_if_no_xlrd()
        ext = self.ext
        filename = '__tmp_to_excel_float_format__.' + ext
        df = DataFrame([[0.123456, 0.234567, 0.567567],
                        [12.32112, 123123.2, 321321.2]],
                        index=['A', 'B'], columns=['X', 'Y', 'Z'])

        with ensure_clean(filename) as filename:
            df.to_excel(filename, 'test1', float_format='%.2f')

            reader = ExcelFile(filename)
            rs = reader.parse('test1', index_col=None)
            xp = DataFrame([[0.12, 0.23, 0.57],
                            [12.32, 123123.20, 321321.20]],
                            index=['A', 'B'], columns=['X', 'Y', 'Z'])
            tm.assert_frame_equal(rs, xp)

    def test_to_excel_output_encoding(self):
        _skip_if_no_xlrd()
        ext = self.ext
        filename = '__tmp_to_excel_float_format__.' + ext
        df = DataFrame([[u'\u0192', u'\u0193', u'\u0194'],
                        [u'\u0195', u'\u0196', u'\u0197']],
                        index=[u'A\u0192', 'B'], columns=[u'X\u0193', 'Y', 'Z'])

        with ensure_clean(filename) as filename:
            df.to_excel(filename, sheet_name = 'TestSheet', encoding='utf8')
	    result = read_excel(filename, 'TestSheet', encoding = 'utf8')
            tm.assert_frame_equal(result,df)


    def test_to_excel_unicode_filename(self):
        _skip_if_no_xlrd()
        ext = self.ext
        filename = u('\u0192u.') + ext

        try:
            f = open(filename, 'wb')
        except UnicodeEncodeError:
            raise nose.SkipTest('no unicode file names on this system')
        else:
            f.close()

        df = DataFrame([[0.123456, 0.234567, 0.567567],
                        [12.32112, 123123.2, 321321.2]],
                        index=['A', 'B'], columns=['X', 'Y', 'Z'])

        with ensure_clean(filename) as filename:
            df.to_excel(filename, 'test1', float_format='%.2f')

            reader = ExcelFile(filename)
            rs = reader.parse('test1', index_col=None)
            xp = DataFrame([[0.12, 0.23, 0.57],
                            [12.32, 123123.20, 321321.20]],
                            index=['A', 'B'], columns=['X', 'Y', 'Z'])
            tm.assert_frame_equal(rs, xp)

    # def test_to_excel_header_styling_xls(self):

    #     import StringIO
    #     s = StringIO(
    #     """Date,ticker,type,value
    #     2001-01-01,x,close,12.2
    #     2001-01-01,x,open ,12.1
    #     2001-01-01,y,close,12.2
    #     2001-01-01,y,open ,12.1
    #     2001-02-01,x,close,12.2
    #     2001-02-01,x,open ,12.1
    #     2001-02-01,y,close,12.2
    #     2001-02-01,y,open ,12.1
    #     2001-03-01,x,close,12.2
    #     2001-03-01,x,open ,12.1
    #     2001-03-01,y,close,12.2
    #     2001-03-01,y,open ,12.1""")
    #     df = read_csv(s, parse_dates=["Date"])
    #     pdf = df.pivot_table(values="value", rows=["ticker"],
    #                                          cols=["Date", "type"])

    #     try:
    #         import xlwt
    #         import xlrd
    #     except ImportError:
    #         raise nose.SkipTest

    #     filename = '__tmp_to_excel_header_styling_xls__.xls'
    #     pdf.to_excel(filename, 'test1')

    #     wbk = xlrd.open_workbook(filename,
    #                              formatting_info=True)
    #     self.assertEquals(["test1"], wbk.sheet_names())
    #     ws = wbk.sheet_by_name('test1')
    #     self.assertEquals([(0, 1, 5, 7), (0, 1, 3, 5), (0, 1, 1, 3)],
    #                       ws.merged_cells)
    #     for i in range(0, 2):
    #         for j in range(0, 7):
    #             xfx = ws.cell_xf_index(0, 0)
    #             cell_xf = wbk.xf_list[xfx]
    #             font = wbk.font_list
    #             self.assertEquals(1, font[cell_xf.font_index].bold)
    #             self.assertEquals(1, cell_xf.border.top_line_style)
    #             self.assertEquals(1, cell_xf.border.right_line_style)
    #             self.assertEquals(1, cell_xf.border.bottom_line_style)
    #             self.assertEquals(1, cell_xf.border.left_line_style)
    #             self.assertEquals(2, cell_xf.alignment.hor_align)
    #     os.remove(filename)
    # def test_to_excel_header_styling_xlsx(self):
    #     import StringIO
    #     s = StringIO(
    #     """Date,ticker,type,value
    #     2001-01-01,x,close,12.2
    #     2001-01-01,x,open ,12.1
    #     2001-01-01,y,close,12.2
    #     2001-01-01,y,open ,12.1
    #     2001-02-01,x,close,12.2
    #     2001-02-01,x,open ,12.1
    #     2001-02-01,y,close,12.2
    #     2001-02-01,y,open ,12.1
    #     2001-03-01,x,close,12.2
    #     2001-03-01,x,open ,12.1
    #     2001-03-01,y,close,12.2
    #     2001-03-01,y,open ,12.1""")
    #     df = read_csv(s, parse_dates=["Date"])
    #     pdf = df.pivot_table(values="value", rows=["ticker"],
    #                                          cols=["Date", "type"])
    #     try:
    #         import openpyxl
    #         from openpyxl.cell import get_column_letter
    #     except ImportError:
    #         raise nose.SkipTest
    #     if openpyxl.__version__ < '1.6.1':
    #         raise nose.SkipTest
    #     # test xlsx_styling
    #     filename = '__tmp_to_excel_header_styling_xlsx__.xlsx'
    #     pdf.to_excel(filename, 'test1')
    #     wbk = openpyxl.load_workbook(filename)
    #     self.assertEquals(["test1"], wbk.get_sheet_names())
    #     ws = wbk.get_sheet_by_name('test1')
    #     xlsaddrs = ["%s2" % chr(i) for i in range(ord('A'), ord('H'))]
    #     xlsaddrs += ["A%s" % i for i in range(1, 6)]
    #     xlsaddrs += ["B1", "D1", "F1"]
    #     for xlsaddr in xlsaddrs:
    #         cell = ws.cell(xlsaddr)
    #         self.assertTrue(cell.style.font.bold)
    #         self.assertEquals(openpyxl.style.Border.BORDER_THIN,
    #                           cell.style.borders.top.border_style)
    #         self.assertEquals(openpyxl.style.Border.BORDER_THIN,
    #                           cell.style.borders.right.border_style)
    #         self.assertEquals(openpyxl.style.Border.BORDER_THIN,
    #                           cell.style.borders.bottom.border_style)
    #         self.assertEquals(openpyxl.style.Border.BORDER_THIN,
    #                           cell.style.borders.left.border_style)
    #         self.assertEquals(openpyxl.style.Alignment.HORIZONTAL_CENTER,
    #                           cell.style.alignment.horizontal)
    #     mergedcells_addrs = ["C1", "E1", "G1"]
    #     for maddr in mergedcells_addrs:
    #         self.assertTrue(ws.cell(maddr).merged)
    #     os.remove(filename)
    def test_excel_010_hemstring(self):
        _skip_if_no_xlrd()
        from pandas.util.testing import makeCustomDataframe as mkdf
        # ensure limited functionality in 0.10
        # override of #2370 until sorted out in 0.11

        def roundtrip(df, header=True, parser_hdr=0):
            path = '__tmp__test_xl_010_%s__.%s' % (np.random.randint(1, 10000), self.ext)
            df.to_excel(path, header=header)

            with ensure_clean(path) as path:
                xf = pd.ExcelFile(path)
                res = xf.parse(xf.sheet_names[0], header=parser_hdr)
                return res

        nrows = 5
        ncols = 3

        for i in range(1, 4):  # row multindex upto nlevel=3
            for j in range(1, 4):  # col ""
                df = mkdf(nrows, ncols, r_idx_nlevels=i, c_idx_nlevels=j)
                res = roundtrip(df)
                # shape
                self.assertEqual(res.shape, (nrows, ncols + i))

                # no nans
                for r in range(len(res.index)):
                    for c in range(len(res.columns)):
                        self.assertTrue(res.ix[r, c] is not np.nan)

        for i in range(1, 4):  # row multindex upto nlevel=3
            for j in range(1, 4):  # col ""
                df = mkdf(nrows, ncols, r_idx_nlevels=i, c_idx_nlevels=j)
                res = roundtrip(df, False)
                # shape
                self.assertEqual(res.shape, (
                    nrows - 1, ncols + i))  # first row taken as columns

                # no nans
                for r in range(len(res.index)):
                    for c in range(len(res.columns)):
                        self.assertTrue(res.ix[r, c] is not np.nan)

        res = roundtrip(DataFrame([0]))
        self.assertEqual(res.shape, (1, 1))
        self.assertTrue(res.ix[0, 0] is not np.nan)

        res = roundtrip(DataFrame([0]), False, None)
        self.assertEqual(res.shape, (1, 2))
        self.assertTrue(res.ix[0, 0] is not np.nan)


class OpenpyxlTests(ExcelWriterBase, unittest.TestCase):
    ext = 'xlsx'
    engine_name = 'openpyxl'
    check_skip = staticmethod(_skip_if_no_openpyxl)

    def test_to_excel_styleconverter(self):
        _skip_if_no_openpyxl()

        import openpyxl

        hstyle = {"font": {"bold": True},
                  "borders": {"top": "thin",
                              "right": "thin",
                              "bottom": "thin",
                              "left": "thin"},
                  "alignment": {"horizontal": "center"}}

        xlsx_style = _OpenpyxlWriter._convert_to_style(hstyle)
        self.assertTrue(xlsx_style.font.bold)
        self.assertEquals(openpyxl.style.Border.BORDER_THIN,
                          xlsx_style.borders.top.border_style)
        self.assertEquals(openpyxl.style.Border.BORDER_THIN,
                          xlsx_style.borders.right.border_style)
        self.assertEquals(openpyxl.style.Border.BORDER_THIN,
                          xlsx_style.borders.bottom.border_style)
        self.assertEquals(openpyxl.style.Border.BORDER_THIN,
                          xlsx_style.borders.left.border_style)
        self.assertEquals(openpyxl.style.Alignment.HORIZONTAL_CENTER,
                          xlsx_style.alignment.horizontal)


class XlwtTests(ExcelWriterBase, unittest.TestCase):
    ext = 'xls'
    engine_name = 'xlwt'
    check_skip = staticmethod(_skip_if_no_xlwt)

    def test_to_excel_styleconverter(self):
        _skip_if_no_xlwt()

        import xlwt

        hstyle = {"font": {"bold": True},
                  "borders": {"top": "thin",
                              "right": "thin",
                              "bottom": "thin",
                              "left": "thin"},
                  "alignment": {"horizontal": "center"}}
        xls_style = _XlwtWriter._convert_to_style(hstyle)
        self.assertTrue(xls_style.font.bold)
        self.assertEquals(xlwt.Borders.THIN, xls_style.borders.top)
        self.assertEquals(xlwt.Borders.THIN, xls_style.borders.right)
        self.assertEquals(xlwt.Borders.THIN, xls_style.borders.bottom)
        self.assertEquals(xlwt.Borders.THIN, xls_style.borders.left)
        self.assertEquals(xlwt.Alignment.HORZ_CENTER, xls_style.alignment.horz)


class XlsxWriterTests(ExcelWriterBase, unittest.TestCase):
    ext = 'xlsx'
    engine_name = 'xlsxwriter'
    check_skip = staticmethod(_skip_if_no_xlsxwriter)

    # Override test from the Superclass to use assertAlmostEqual on the
    # floating point values read back in from the output XlsxWriter file.
    def test_roundtrip_indexlabels(self):
        _skip_if_no_xlrd()
        ext = self.ext
        path = '__tmp_to_excel_from_excel_indexlabels__.' + ext

        with ensure_clean(path) as path:

            self.frame['A'][:5] = nan

            self.frame.to_excel(path, 'test1')
            self.frame.to_excel(path, 'test1', cols=['A', 'B'])
            self.frame.to_excel(path, 'test1', header=False)
            self.frame.to_excel(path, 'test1', index=False)

            # test index_label
            frame = (DataFrame(np.random.randn(10, 2)) >= 0)
            frame.to_excel(path, 'test1', index_label=['test'])
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0).astype(np.int64)
            frame.index.names = ['test']
            self.assertEqual(frame.index.names, recons.index.names)

            frame = (DataFrame(np.random.randn(10, 2)) >= 0)
            frame.to_excel(
                path, 'test1', index_label=['test', 'dummy', 'dummy2'])
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0).astype(np.int64)
            frame.index.names = ['test']
            self.assertEqual(frame.index.names, recons.index.names)

            frame = (DataFrame(np.random.randn(10, 2)) >= 0)
            frame.to_excel(path, 'test1', index_label='test')
            reader = ExcelFile(path)
            recons = reader.parse('test1', index_col=0).astype(np.int64)
            frame.index.names = ['test']
            self.assertAlmostEqual(frame.index.names, recons.index.names)


class ExcelWriterEngineTests(unittest.TestCase):
    def test_ExcelWriter_dispatch(self):
        with tm.assertRaisesRegexp(ValueError, 'No engine'):
            writer = ExcelWriter('nothing')

        _skip_if_no_openpyxl()
        writer = ExcelWriter('apple.xlsx')
        tm.assert_isinstance(writer, _OpenpyxlWriter)

        _skip_if_no_xlwt()
        writer = ExcelWriter('apple.xls')
        tm.assert_isinstance(writer, _XlwtWriter)


    def test_register_writer(self):
        # some awkward mocking to test out dispatch and such actually works
        called_save = []
        called_write_cells = []
        class DummyClass(ExcelWriter):
            called_save = False
            called_write_cells = False
            supported_extensions = ['test', 'xlsx', 'xls']
            engine = 'dummy'

            def save(self):
                called_save.append(True)

            def write_cells(self, *args, **kwargs):
                called_write_cells.append(True)

        def check_called(func):
            func()
            self.assert_(len(called_save) >= 1)
            self.assert_(len(called_write_cells) >= 1)
            del called_save[:]
            del called_write_cells[:]

        register_writer(DummyClass)
        writer = ExcelWriter('something.test')
        tm.assert_isinstance(writer, DummyClass)
        df = tm.makeCustomDataframe(1, 1)
        panel = tm.makePanel()
        func = lambda: df.to_excel('something.test')
        check_called(func)
        check_called(lambda: panel.to_excel('something.test'))
        from pandas import set_option, get_option
        val = get_option('io.excel.xlsx.writer')
        set_option('io.excel.xlsx.writer', 'dummy')
        check_called(lambda: df.to_excel('something.xlsx'))
        check_called(lambda: df.to_excel('something.xls', engine='dummy'))
        set_option('io.excel.xlsx.writer', val)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
