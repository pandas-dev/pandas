from cStringIO import StringIO
from datetime import datetime
import csv
import os
import unittest

import nose

from numpy import nan
import numpy as np

from pandas import DataFrame, Index
from pandas.io.parsers import read_csv, read_table, ExcelFile, TextParser
from pandas.util.testing import assert_almost_equal, assert_frame_equal
import pandas._tseries as lib

class TestParsers(unittest.TestCase):
    data1 = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""

    def setUp(self):
        self.dirpath = curpath()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')

    def test_read_csv(self):
        pass

    def test_custom_na_values(self):
        data = """A,B,C
ignore,this,row
1,NA,3
-1.#IND,5,baz
7,8,NaN
"""
        expected = [[1., nan, 3],
                    [nan, 5, nan],
                    [7, 8, nan]]

        df = read_csv(StringIO(data), na_values=['baz'], skiprows=[1])
        assert_almost_equal(df.values, expected)

        df2 = read_table(StringIO(data), sep=',', na_values=['baz'],
                         skiprows=[1])
        assert_almost_equal(df2.values, expected)

    def test_unnamed_columns(self):
        data = """A,B,C,,
1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        expected = [[1,2,3,4,5.],
                    [6,7,8,9,10],
                    [11,12,13,14,15]]
        df = read_table(StringIO(data), sep=',')
        assert_almost_equal(df.values, expected)
        self.assert_(np.array_equal(df.columns,
                                    ['A', 'B', 'C', 'Unnamed: 3',
                                     'Unnamed: 4']))

    def test_string_nas(self):
        data = """A,B,C
a,b,c
d,,f
,g,h
"""
        result = read_csv(StringIO(data))
        expected = DataFrame([['a', 'b', 'c'],
                              ['d', np.nan, 'f'],
                              [np.nan, 'g', 'h']],
                             columns=['A', 'B', 'C'])

        assert_frame_equal(result, expected)

    def test_duplicate_columns(self):
        data = """A,A,B,B,B
1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        df = read_table(StringIO(data), sep=',')
        self.assert_(np.array_equal(df.columns,
                                    ['A', 'A.1', 'B', 'B.1', 'B.2']))

    def test_csv_mixed_type(self):
        data = """A,B,C
a,1,2
b,3,4
c,4,5
"""
        df = read_csv(StringIO(data))
        # TODO

    def test_csv_custom_parser(self):
        data = """A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        df = read_csv(StringIO(data),
                      date_parser=lambda x: datetime.strptime(x, '%Y%m%d'))
        expected = read_csv(StringIO(data), parse_dates=True)
        assert_frame_equal(df, expected)

    def test_parse_dates_implicit_first_col(self):
        data = """A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        df = read_csv(StringIO(data), parse_dates=True)
        expected = read_csv(StringIO(data), index_col=0, parse_dates=True)
        self.assert_(isinstance(df.index[0], datetime))
        assert_frame_equal(df, expected)

    def test_no_header(self):
        data = """1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        df = read_table(StringIO(data), sep=',', header=None)
        names = ['foo', 'bar', 'baz', 'quux', 'panda']
        df2 = read_table(StringIO(data), sep=',', header=None, names=names)
        expected = [[1,2,3,4,5.],
                    [6,7,8,9,10],
                    [11,12,13,14,15]]
        assert_almost_equal(df.values, expected)
        assert_almost_equal(df.values, df2.values)
        self.assert_(np.array_equal(df.columns,
                                    ['X.1', 'X.2', 'X.3', 'X.4', 'X.5']))
        self.assert_(np.array_equal(df2.columns, names))

    def test_header_with_index_col(self):
        data = """foo,1,2,3
bar,4,5,6
baz,7,8,9
"""
        names = ['A', 'B', 'C']
        df = read_csv(StringIO(data), names=names)

        self.assertEqual(names, ['A', 'B', 'C'])

        values = [[1,2,3],[4,5,6],[7,8,9]]
        expected = DataFrame(values, index=['foo','bar','baz'],
                             columns=['A','B','C'])
        assert_frame_equal(df, expected)

    def test_read_csv_dataframe(self):
        df = read_csv(self.csv1, index_col=0, parse_dates=True)
        df2 = read_table(self.csv1, sep=',', index_col=0, parse_dates=True)
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D']))
        self.assert_(df.index.name == 'index')
        self.assert_(isinstance(df.index[0], datetime))
        self.assert_(df.values.dtype == np.float64)
        assert_frame_equal(df, df2)

    def test_read_csv_no_index_name(self):
        df = read_csv(self.csv2, index_col=0, parse_dates=True)
        df2 = read_table(self.csv2, sep=',', index_col=0, parse_dates=True)
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D', 'E']))
        self.assert_(isinstance(df.index[0], datetime))
        self.assert_(df.ix[:, ['A', 'B', 'C', 'D']].values.dtype == np.float64)
        assert_frame_equal(df, df2)

    def test_excel_table(self):
        try:
            import xlrd
        except ImportError:
            raise nose.SkipTest('xlrd not installed, skipping')

        pth = os.path.join(self.dirpath, 'test.xls')
        xls = ExcelFile(pth)
        df = xls.parse('Sheet1', index_col=0, parse_dates=True)
        df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
        df3 = xls.parse('Sheet2', skiprows=[1], index_col=0, parse_dates=True)
        assert_frame_equal(df, df2)
        assert_frame_equal(df3, df2)

    def test_read_table_wrong_num_columns(self):
        data = """A,B,C,D,E,F
1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        self.assertRaises(Exception, read_csv, StringIO(data))

    def test_read_table_duplicate_index(self):
        data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo,12,13,14,15
bar,12,13,14,15
"""

        self.assertRaises(Exception, read_csv, StringIO(data),
                          index_col=0)

    def test_parse_bools(self):
        data = """A,B
True,1
False,2
True,3
"""
        data = read_csv(StringIO(data))
        self.assert_(data['A'].dtype == np.bool_)

    def test_int_conversion(self):
        data = """A,B
1.0,1
2.0,2
3.0,3
"""
        data = read_csv(StringIO(data))
        self.assert_(data['A'].dtype == np.float64)
        self.assert_(data['B'].dtype == np.int64)

    def test_infer_index_col(self):
        data = """A,B,C
foo,1,2,3
bar,4,5,6
baz,7,8,9
"""
        data = read_csv(StringIO(data))
        self.assert_(data.index.equals(Index(['foo', 'bar', 'baz'])))

    def test_sniff_delimiter(self):
        data = """index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
        data = read_csv(StringIO(data), index_col=0)
        self.assert_(data.index.equals(Index(['foo', 'bar', 'baz'])))

    def test_read_nrows(self):
        df = read_csv(StringIO(self.data1), nrows=3)
        expected = read_csv(StringIO(self.data1))[:3]
        assert_frame_equal(df, expected)

    def test_read_chunksize(self):
        reader = read_csv(StringIO(self.data1), index_col=0, chunksize=2)
        df = read_csv(StringIO(self.data1), index_col=0)

        chunks = list(reader)

        assert_frame_equal(chunks[0], df[:2])
        assert_frame_equal(chunks[1], df[2:4])
        assert_frame_equal(chunks[2], df[4:])

    def test_iterator(self):
        reader = read_csv(StringIO(self.data1), index_col=0, iterator=True)

        df = read_csv(StringIO(self.data1), index_col=0)

        chunk = reader.get_chunk(3)
        assert_frame_equal(chunk, df[:3])

        last_chunk = reader.get_chunk(5)
        assert_frame_equal(last_chunk, df[3:])

        # pass list
        lines = list(csv.reader(StringIO(self.data1)))
        parser = TextParser(lines, index_col=0, chunksize=2)

        df = read_csv(StringIO(self.data1), index_col=0)

        chunks = list(parser)
        assert_frame_equal(chunks[0], df[:2])
        assert_frame_equal(chunks[1], df[2:4])
        assert_frame_equal(chunks[2], df[4:])

        treader = read_table(StringIO(self.data1), sep=',', index_col=0,
                             iterator=True)
        self.assert_(isinstance(treader, TextParser))

    def test_header_not_first_line(self):
        data = """got,to,ignore,this,line
got,to,ignore,this,line
index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
"""
        data2 = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
"""

        df = read_csv(StringIO(data), header=2, index_col=0)
        expected = read_csv(StringIO(data2), header=0, index_col=0)
        assert_frame_equal(df, expected)

    def test_pass_names_with_index(self):
        lines = self.data1.split('\n')
        no_header = '\n'.join(lines[1:])

        # regular index
        names = ['index', 'A', 'B', 'C', 'D']
        df = read_csv(StringIO(no_header), index_col=0, names=names)
        expected = read_csv(StringIO(self.data1), index_col=0)
        assert_frame_equal(df, expected)

        # multi index
        data = """index1,index2,A,B,C,D
foo,one,2,3,4,5
foo,two,7,8,9,10
foo,three,12,13,14,15
bar,one,12,13,14,15
bar,two,12,13,14,15
"""
        lines = data.split('\n')
        no_header = '\n'.join(lines[1:])
        names = ['index1', 'index2', 'A', 'B', 'C', 'D']
        df = read_csv(StringIO(no_header), index_col=[0, 1], names=names)
        expected = read_csv(StringIO(data), index_col=[0, 1])
        assert_frame_equal(df, expected)

    def test_multi_index_no_level_names(self):
        data = """index1,index2,A,B,C,D
foo,one,2,3,4,5
foo,two,7,8,9,10
foo,three,12,13,14,15
bar,one,12,13,14,15
bar,two,12,13,14,15
"""

        data2 = """A,B,C,D
foo,one,2,3,4,5
foo,two,7,8,9,10
foo,three,12,13,14,15
bar,one,12,13,14,15
bar,two,12,13,14,15
"""

        lines = data.split('\n')
        no_header = '\n'.join(lines[1:])
        names = ['A', 'B', 'C', 'D']
        df = read_csv(StringIO(no_header), index_col=[0, 1], names=names)
        expected = read_csv(StringIO(data), index_col=[0, 1])
        assert_frame_equal(df, expected)

        # 2 implicit first cols
        df2 = read_csv(StringIO(data2))
        assert_frame_equal(df2, df)

    def test_multi_index_parse_dates(self):
        data = """index1,index2,A,B,C
20090101,one,a,1,2
20090101,two,b,3,4
20090101,three,c,4,5
20090102,one,a,1,2
20090102,two,b,3,4
20090102,three,c,4,5
20090103,one,a,1,2
20090103,two,b,3,4
20090103,three,c,4,5
"""
        df = read_csv(StringIO(data), index_col=[0, 1], parse_dates=True)
        self.assert_(isinstance(df.index.levels[0][0], datetime))

class TestParseSQL(unittest.TestCase):

    def test_convert_sql_column_floats(self):
        arr = np.array([1.5, None, 3, 4.2], dtype=object)
        result = lib.convert_sql_column(arr)
        expected = np.array([1.5, np.nan, 3, 4.2], dtype='f8')
        assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_strings(self):
        arr = np.array(['1.5', None, '3', '4.2'], dtype=object)
        result = lib.convert_sql_column(arr)
        expected = np.array(['1.5', np.nan, '3', '4.2'], dtype=object)
        assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_unicode(self):
        arr = np.array([u'1.5', None, u'3', u'4.2'], dtype=object)
        result = lib.convert_sql_column(arr)
        expected = np.array([u'1.5', np.nan, u'3', u'4.2'], dtype=object)
        assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_ints(self):
        arr = np.array([1, 2, 3, 4], dtype='O')
        arr2 = np.array([1, 2, 3, 4], dtype='i4').astype('O')
        result = lib.convert_sql_column(arr)
        result2 = lib.convert_sql_column(arr2)
        expected = np.array([1, 2, 3, 4], dtype='i8')
        assert_same_values_and_dtype(result, expected)
        assert_same_values_and_dtype(result2, expected)

        arr = np.array([1, 2, 3, None, 4], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([1, 2, 3, np.nan, 4], dtype='f8')
        assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_longs(self):
        arr = np.array([1L, 2L, 3L, 4L], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([1, 2, 3, 4], dtype='i8')
        assert_same_values_and_dtype(result, expected)

        arr = np.array([1L, 2L, 3L, None, 4L], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([1, 2, 3, np.nan, 4], dtype='f8')
        assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_bools(self):
        arr = np.array([True, False, True, False], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([True, False, True, False], dtype=bool)
        assert_same_values_and_dtype(result, expected)

        arr = np.array([True, False, None, False], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([True, False, np.nan, False], dtype=object)
        assert_same_values_and_dtype(result, expected)

    def test_convert_sql_column_decimals(self):
        from decimal import Decimal
        arr = np.array([Decimal('1.5'), None, Decimal('3'), Decimal('4.2')])
        result = lib.convert_sql_column(arr)
        expected = np.array([1.5, np.nan, 3, 4.2], dtype='f8')
        assert_same_values_and_dtype(result, expected)

def assert_same_values_and_dtype(res, exp):
    assert(res.dtype == exp.dtype)
    assert_almost_equal(res, exp)

def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

