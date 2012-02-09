try:
    from io import BytesIO
except ImportError:  # pragma: no cover
    # Python < 2.6
    from cStringIO import StringIO as BytesIO

from cStringIO import StringIO
from datetime import datetime
import csv
import os
import sys
import re
import unittest

import nose

from numpy import nan
import numpy as np

from pandas import DataFrame, Index, isnull
from pandas.io.parsers import read_csv, read_table, ExcelFile, TextParser
from pandas.util.testing import assert_almost_equal, assert_frame_equal
import pandas._tseries as lib
from pandas.util import py3compat

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


    def test_skiprows_bug(self):
        # GH #505
        text = """#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
#foo,a,b,c
1/1/2000,1.,2.,3.
1/2/2000,4,5,6
1/3/2000,7,8,9
"""
        data = read_csv(StringIO(text), skiprows=range(6), header=None,
                        index_col=0, parse_dates=True)

        data2 = read_csv(StringIO(text), skiprows=6, header=None,
                         index_col=0, parse_dates=True)

        expected = DataFrame(np.arange(1., 10.).reshape((3,3)),
                             columns=['X.2', 'X.3', 'X.4'],
                             index=[datetime(2000, 1, 1), datetime(2000, 1, 2),
                                    datetime(2000, 1, 3)])
        assert_frame_equal(data, expected)
        assert_frame_equal(data, data2)


    def test_detect_string_na(self):
        data = """A,B
foo,bar
NA,baz
NaN,nan
"""
        expected = [['foo', 'bar'],
                    [nan, 'baz'],
                    [nan, nan]]

        df = read_csv(StringIO(data))
        assert_almost_equal(df.values, expected)

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
        text = """index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
        data = read_csv(StringIO(text), index_col=0, sep=None)
        self.assert_(data.index.equals(Index(['foo', 'bar', 'baz'])))

        data2 = read_csv(StringIO(text), index_col=0, delimiter='|')
        assert_frame_equal(data, data2)

        text = """ignore this
ignore this too
index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
        data3 = read_csv(StringIO(text), index_col=0, sep=None, skiprows=2)
        assert_frame_equal(data, data3)

        # can't get this to work on Python 3
        if not py3compat.PY3:
            text = u"""ignore this
ignore this too
index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
""".encode('utf-8')
            data4 = read_csv(BytesIO(text), index_col=0, sep=None, skiprows=2,
                             encoding='utf-8')
            assert_frame_equal(data, data4)

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

    def test_read_text_list(self):
        data = """A,B,C\nfoo,1,2,3\nbar,4,5,6"""
        as_list = [['A','B','C'],['foo','1','2','3'],['bar','4','5','6']]
        df = read_csv(StringIO(data), index_col=0)

        parser = TextParser(as_list, index_col=0, chunksize=2)
        chunk  = parser.get_chunk(None)

        assert_frame_equal(chunk, df)

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

        # pass skiprows
        parser = TextParser(lines, index_col=0, chunksize=2, skiprows=[1])
        chunks = list(parser)
        assert_frame_equal(chunks[0], df[1:3])

        # test bad parameter (skip_footer)
        reader = read_csv(StringIO(self.data1), index_col=0, iterator=True,
                          skip_footer=True)
        self.assertRaises(ValueError, reader.get_chunk, 3)

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

        # specify columns out of order!
        df2 = read_csv(StringIO(data), index_col=[1, 0], parse_dates=True)
        self.assert_(isinstance(df2.index.levels[1][0], datetime))

    def test_skip_footer(self):
        data = """A,B,C
1,2,3
4,5,6
7,8,9
want to skip this
also also skip this
and this
"""
        result = read_csv(StringIO(data), skip_footer=3)
        no_footer = '\n'.join(data.split('\n')[:-4])
        expected = read_csv(StringIO(no_footer))

        assert_frame_equal(result, expected)

    def test_no_unnamed_index(self):
        data = """ id c0 c1 c2
0 1 0 a b
1 2 0 c d
2 2 2 e f
"""
        df = read_table(StringIO(data), sep=' ')
        self.assert_(df.index.name is None)

    def test_converters(self):
        data = """A,B,C,D
a,1,2,01/01/2009
b,3,4,01/02/2009
c,4,5,01/03/2009
"""
        from dateutil import parser

        result = read_csv(StringIO(data), converters={'D' : parser.parse})
        result2 = read_csv(StringIO(data), converters={3 : parser.parse})

        expected = read_csv(StringIO(data))
        expected['D'] = expected['D'].map(parser.parse)

        self.assert_(isinstance(result['D'][0], datetime))
        assert_frame_equal(result, expected)
        assert_frame_equal(result2, expected)

        # produce integer
        converter = lambda x: int(x.split('/')[2])
        result = read_csv(StringIO(data), converters={'D' : converter})
        expected = read_csv(StringIO(data))
        expected['D'] = expected['D'].map(converter)
        assert_frame_equal(result, expected)

    def test_converters_euro_decimal_format(self):
        data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""
        f = lambda x : float(x.replace(",", "."))
        converter = {'Number1':f,'Number2':f, 'Number3':f}
        df2 = read_csv(StringIO(data), sep=';',converters=converter)
        self.assert_(df2['Number1'].dtype == float)
        self.assert_(df2['Number2'].dtype == float)
        self.assert_(df2['Number3'].dtype == float)

    def test_converter_return_string_bug(self):
        # GH #583
        data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""
        f = lambda x : x.replace(",", ".")
        converter = {'Number1':f,'Number2':f, 'Number3':f}
        df2 = read_csv(StringIO(data), sep=';',converters=converter)
        self.assert_(df2['Number1'].dtype == float)

    def test_regex_separator(self):
        data = """   A   B   C   D
a   1   2   3   4
b   1   2   3   4
c   1   2   3   4
"""
        df = read_table(StringIO(data), sep='\s+')
        expected = read_csv(StringIO(re.sub('[ ]+', ',', data)),
                            index_col=0)
        self.assert_(expected.index.name is None)
        assert_frame_equal(df, expected)

    def test_verbose_import(self):
        text = """a,b,c,d
one,1,2,3
one,1,2,3
,1,2,3
one,1,2,3
,1,2,3
,1,2,3
one,1,2,3
two,1,2,3"""

        buf = StringIO()
        sys.stdout = buf

        try:
            # it works!
            df = read_csv(StringIO(text), verbose=True)
            self.assert_(buf.getvalue() == 'Filled 3 NA values in column a\n')
        finally:
            sys.stdout = sys.__stdout__

        buf = StringIO()
        sys.stdout = buf

        text = """a,b,c,d
one,1,2,3
two,1,2,3
three,1,2,3
four,1,2,3
five,1,2,3
,1,2,3
seven,1,2,3
eight,1,2,3"""

        try:
            # it works!
            df = read_csv(StringIO(text), verbose=True, index_col=0)
            self.assert_(buf.getvalue() == 'Found 1 NA values in the index\n')
        finally:
            sys.stdout = sys.__stdout__

    def test_read_table_buglet_4x_multiindex(self):
        text = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""

        # it works!
        df = read_table(StringIO(text), sep='\s+')
        self.assertEquals(df.index.names, ['one', 'two', 'three', 'four'])

    def test_read_csv_parse_simple_list(self):
        text = """foo
bar baz
qux foo
foo
bar"""
        df = read_csv(StringIO(text), header=None)
        expected = DataFrame({'X.1' : ['foo', 'bar baz', 'qux foo',
                                       'foo', 'bar']})
        assert_frame_equal(df, expected)

    def test_converters_corner_with_nas(self):
        import StringIO
        import numpy as np
        import pandas
        csv = """id,score,days
1,2,12
2,2-5,
3,,14+
4,6-12,2"""

        def convert_days(x):
           x = x.strip()
           if not x: return np.nan

           is_plus = x.endswith('+')
           if is_plus:
               x = int(x[:-1]) + 1
           else:
               x = int(x)
           return x

        def convert_days_sentinel(x):
           x = x.strip()
           if not x: return -1

           is_plus = x.endswith('+')
           if is_plus:
               x = int(x[:-1]) + 1
           else:
               x = int(x)
           return x

        def convert_score(x):
           x = x.strip()
           if not x: return np.nan
           if x.find('-')>0:
               valmin, valmax = map(int, x.split('-'))
               val = 0.5*(valmin + valmax)
           else:
               val = float(x)

           return val

        fh = StringIO.StringIO(csv)
        result = pandas.read_csv(fh, converters={'score':convert_score,
                                                 'days':convert_days},
                                 na_values=[-1,'',None])
        self.assert_(isnull(result['days'][1]))

        fh = StringIO.StringIO(csv)
        result2 = pandas.read_csv(fh, converters={'score':convert_score,
                                                  'days':convert_days_sentinel},
                                  na_values=[-1,'',None])
        assert_frame_equal(result, result2)

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

