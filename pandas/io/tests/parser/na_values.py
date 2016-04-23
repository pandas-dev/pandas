# -*- coding: utf-8 -*-

"""
Tests that NA values are properly handled during
parsing for all of the parsers defined in parsers.py
"""

import numpy as np
from numpy import nan

import pandas.io.parsers as parsers
import pandas.util.testing as tm

from pandas import DataFrame, MultiIndex, read_csv
from pandas.compat import StringIO, range


class NAvaluesTests(object):

    def test_string_nas(self):
        data = """A,B,C
a,b,c
d,,f
,g,h
"""
        result = self.read_csv(StringIO(data))
        expected = DataFrame([['a', 'b', 'c'],
                              ['d', np.nan, 'f'],
                              [np.nan, 'g', 'h']],
                             columns=['A', 'B', 'C'])

        tm.assert_frame_equal(result, expected)

    def test_detect_string_na(self):
        data = """A,B
foo,bar
NA,baz
NaN,nan
"""
        expected = [['foo', 'bar'], [nan, 'baz'], [nan, nan]]
        df = self.read_csv(StringIO(data))
        tm.assert_almost_equal(df.values, expected)

    def test_non_string_na_values(self):
        # see gh-3611, na_values that are not a string are an issue
        with tm.ensure_clean('__non_string_na_values__.csv') as path:
            df = DataFrame({'A': [-999, 2, 3], 'B': [1.2, -999, 4.5]})
            df.to_csv(path, sep=' ', index=False)
            result1 = self.read_csv(path, sep=' ', header=0,
                                    na_values=['-999.0', '-999'])
            result2 = self.read_csv(path, sep=' ', header=0,
                                    na_values=[-999, -999.0])
            result3 = self.read_csv(path, sep=' ', header=0,
                                    na_values=[-999.0, -999])
            tm.assert_frame_equal(result1, result2)
            tm.assert_frame_equal(result2, result3)

            result4 = self.read_csv(
                path, sep=' ', header=0, na_values=['-999.0'])
            result5 = self.read_csv(
                path, sep=' ', header=0, na_values=['-999'])
            result6 = self.read_csv(
                path, sep=' ', header=0, na_values=[-999.0])
            result7 = self.read_csv(
                path, sep=' ', header=0, na_values=[-999])
            tm.assert_frame_equal(result4, result3)
            tm.assert_frame_equal(result5, result3)
            tm.assert_frame_equal(result6, result3)
            tm.assert_frame_equal(result7, result3)

            good_compare = result3

            # with an odd float format, so we can't match the string 999.0
            # exactly, but need float matching
            # TODO: change these to self.read_csv when Python bug is squashed
            df.to_csv(path, sep=' ', index=False, float_format='%.3f')
            result1 = read_csv(path, sep=' ', header=0,
                               na_values=['-999.0', '-999'])
            result2 = read_csv(path, sep=' ', header=0,
                               na_values=[-999.0, -999])
            tm.assert_frame_equal(result1, good_compare)
            tm.assert_frame_equal(result2, good_compare)

            result3 = read_csv(path, sep=' ',
                               header=0, na_values=['-999.0'])
            result4 = read_csv(path, sep=' ',
                               header=0, na_values=['-999'])
            result5 = read_csv(path, sep=' ',
                               header=0, na_values=[-999.0])
            result6 = read_csv(path, sep=' ',
                               header=0, na_values=[-999])
            tm.assert_frame_equal(result3, good_compare)
            tm.assert_frame_equal(result4, good_compare)
            tm.assert_frame_equal(result5, good_compare)
            tm.assert_frame_equal(result6, good_compare)

    def test_default_na_values(self):
        _NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN',
                          '#N/A', 'N/A', 'NA', '#NA', 'NULL', 'NaN',
                          'nan', '-NaN', '-nan', '#N/A N/A', ''])
        self.assertEqual(_NA_VALUES, parsers._NA_VALUES)
        nv = len(_NA_VALUES)

        def f(i, v):
            if i == 0:
                buf = ''
            elif i > 0:
                buf = ''.join([','] * i)

            buf = "{0}{1}".format(buf, v)

            if i < nv - 1:
                buf = "{0}{1}".format(buf, ''.join([','] * (nv - i - 1)))

            return buf

        data = StringIO('\n'.join([f(i, v) for i, v in enumerate(_NA_VALUES)]))
        expected = DataFrame(np.nan, columns=range(nv), index=range(nv))
        df = self.read_csv(data, header=None)
        tm.assert_frame_equal(df, expected)

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

        df = self.read_csv(StringIO(data), na_values=['baz'], skiprows=[1])
        tm.assert_almost_equal(df.values, expected)

        df2 = self.read_table(StringIO(data), sep=',', na_values=['baz'],
                              skiprows=[1])
        tm.assert_almost_equal(df2.values, expected)

        df3 = self.read_table(StringIO(data), sep=',', na_values='baz',
                              skiprows=[1])
        tm.assert_almost_equal(df3.values, expected)

    def test_bool_na_values(self):
        data = """A,B,C
True,False,True
NA,True,False
False,NA,True"""

        result = self.read_csv(StringIO(data))
        expected = DataFrame({'A': np.array([True, nan, False], dtype=object),
                              'B': np.array([False, True, nan], dtype=object),
                              'C': [True, False, True]})

        tm.assert_frame_equal(result, expected)

    def test_na_value_dict(self):
        data = """A,B,C
foo,bar,NA
bar,foo,foo
foo,bar,NA
bar,foo,foo"""

        df = self.read_csv(StringIO(data),
                           na_values={'A': ['foo'], 'B': ['bar']})
        expected = DataFrame({'A': [np.nan, 'bar', np.nan, 'bar'],
                              'B': [np.nan, 'foo', np.nan, 'foo'],
                              'C': [np.nan, 'foo', np.nan, 'foo']})
        tm.assert_frame_equal(df, expected)

        data = """\
a,b,c,d
0,NA,1,5
"""
        xp = DataFrame({'b': [np.nan], 'c': [1], 'd': [5]}, index=[0])
        xp.index.name = 'a'
        df = self.read_csv(StringIO(data), na_values={}, index_col=0)
        tm.assert_frame_equal(df, xp)

        xp = DataFrame({'b': [np.nan], 'd': [5]},
                       MultiIndex.from_tuples([(0, 1)]))
        xp.index.names = ['a', 'c']
        df = self.read_csv(StringIO(data), na_values={}, index_col=[0, 2])
        tm.assert_frame_equal(df, xp)

        xp = DataFrame({'b': [np.nan], 'd': [5]},
                       MultiIndex.from_tuples([(0, 1)]))
        xp.index.names = ['a', 'c']
        df = self.read_csv(StringIO(data), na_values={}, index_col=['a', 'c'])
        tm.assert_frame_equal(df, xp)

    def test_na_values_keep_default(self):
        data = """\
One,Two,Three
a,1,one
b,2,two
,3,three
d,4,nan
e,5,five
nan,6,
g,7,seven
"""
        df = self.read_csv(StringIO(data))
        xp = DataFrame({'One': ['a', 'b', np.nan, 'd', 'e', np.nan, 'g'],
                        'Two': [1, 2, 3, 4, 5, 6, 7],
                        'Three': ['one', 'two', 'three', np.nan, 'five',
                                  np.nan, 'seven']})
        tm.assert_frame_equal(xp.reindex(columns=df.columns), df)

        df = self.read_csv(StringIO(data), na_values={'One': [], 'Three': []},
                           keep_default_na=False)
        xp = DataFrame({'One': ['a', 'b', '', 'd', 'e', 'nan', 'g'],
                        'Two': [1, 2, 3, 4, 5, 6, 7],
                        'Three': ['one', 'two', 'three', 'nan', 'five',
                                  '', 'seven']})
        tm.assert_frame_equal(xp.reindex(columns=df.columns), df)

        df = self.read_csv(
            StringIO(data), na_values=['a'], keep_default_na=False)
        xp = DataFrame({'One': [np.nan, 'b', '', 'd', 'e', 'nan', 'g'],
                        'Two': [1, 2, 3, 4, 5, 6, 7],
                        'Three': ['one', 'two', 'three', 'nan', 'five', '',
                                  'seven']})
        tm.assert_frame_equal(xp.reindex(columns=df.columns), df)

        df = self.read_csv(StringIO(data), na_values={'One': [], 'Three': []})
        xp = DataFrame({'One': ['a', 'b', np.nan, 'd', 'e', np.nan, 'g'],
                        'Two': [1, 2, 3, 4, 5, 6, 7],
                        'Three': ['one', 'two', 'three', np.nan, 'five',
                                  np.nan, 'seven']})
        tm.assert_frame_equal(xp.reindex(columns=df.columns), df)

        # see gh-4318: passing na_values=None and
        # keep_default_na=False yields 'None' as a na_value
        data = """\
One,Two,Three
a,1,None
b,2,two
,3,None
d,4,nan
e,5,five
nan,6,
g,7,seven
"""
        df = self.read_csv(
            StringIO(data), keep_default_na=False)
        xp = DataFrame({'One': ['a', 'b', '', 'd', 'e', 'nan', 'g'],
                        'Two': [1, 2, 3, 4, 5, 6, 7],
                        'Three': ['None', 'two', 'None', 'nan', 'five', '',
                                  'seven']})
        tm.assert_frame_equal(xp.reindex(columns=df.columns), df)

    def test_skiprow_with_newline(self):
        # see gh-12775 and gh-10911
        data = """id,text,num_lines
1,"line 11
line 12",2
2,"line 21
line 22",2
3,"line 31",1"""
        expected = [[2, 'line 21\nline 22', 2],
                    [3, 'line 31', 1]]
        expected = DataFrame(expected, columns=[
            'id', 'text', 'num_lines'])
        df = self.read_csv(StringIO(data), skiprows=[1])
        tm.assert_frame_equal(df, expected)

        data = ('a,b,c\n~a\n b~,~e\n d~,'
                '~f\n f~\n1,2,~12\n 13\n 14~')
        expected = [['a\n b', 'e\n d', 'f\n f']]
        expected = DataFrame(expected, columns=[
            'a', 'b', 'c'])
        df = self.read_csv(StringIO(data),
                           quotechar="~",
                           skiprows=[2])
        tm.assert_frame_equal(df, expected)

        data = ('Text,url\n~example\n '
                'sentence\n one~,url1\n~'
                'example\n sentence\n two~,url2\n~'
                'example\n sentence\n three~,url3')
        expected = [['example\n sentence\n two', 'url2']]
        expected = DataFrame(expected, columns=[
            'Text', 'url'])
        df = self.read_csv(StringIO(data),
                           quotechar="~",
                           skiprows=[1, 3])
        tm.assert_frame_equal(df, expected)

    def test_skiprow_with_quote(self):
        # see gh-12775 and gh-10911
        data = """id,text,num_lines
1,"line '11' line 12",2
2,"line '21' line 22",2
3,"line '31' line 32",1"""
        expected = [[2, "line '21' line 22", 2],
                    [3, "line '31' line 32", 1]]
        expected = DataFrame(expected, columns=[
            'id', 'text', 'num_lines'])
        df = self.read_csv(StringIO(data), skiprows=[1])
        tm.assert_frame_equal(df, expected)

    def test_skiprow_with_newline_and_quote(self):
        # see gh-12775 and gh-10911
        data = """id,text,num_lines
1,"line \n'11' line 12",2
2,"line \n'21' line 22",2
3,"line \n'31' line 32",1"""
        expected = [[2, "line \n'21' line 22", 2],
                    [3, "line \n'31' line 32", 1]]
        expected = DataFrame(expected, columns=[
            'id', 'text', 'num_lines'])
        df = self.read_csv(StringIO(data), skiprows=[1])
        tm.assert_frame_equal(df, expected)

        data = """id,text,num_lines
1,"line '11\n' line 12",2
2,"line '21\n' line 22",2
3,"line '31\n' line 32",1"""
        expected = [[2, "line '21\n' line 22", 2],
                    [3, "line '31\n' line 32", 1]]
        expected = DataFrame(expected, columns=[
            'id', 'text', 'num_lines'])
        df = self.read_csv(StringIO(data), skiprows=[1])
        tm.assert_frame_equal(df, expected)

        data = """id,text,num_lines
1,"line '11\n' \r\tline 12",2
2,"line '21\n' \r\tline 22",2
3,"line '31\n' \r\tline 32",1"""
        expected = [[2, "line '21\n' \r\tline 22", 2],
                    [3, "line '31\n' \r\tline 32", 1]]
        expected = DataFrame(expected, columns=[
            'id', 'text', 'num_lines'])
        df = self.read_csv(StringIO(data), skiprows=[1])
        tm.assert_frame_equal(df, expected)

    def test_skiprows_lineterminator(self):
        # see gh-9079
        data = '\n'.join(['SMOSMANIA ThetaProbe-ML2X ',
                          '2007/01/01 01:00   0.2140 U M ',
                          '2007/01/01 02:00   0.2141 M O ',
                          '2007/01/01 04:00   0.2142 D M '])
        expected = DataFrame([['2007/01/01', '01:00', 0.2140, 'U', 'M'],
                              ['2007/01/01', '02:00', 0.2141, 'M', 'O'],
                              ['2007/01/01', '04:00', 0.2142, 'D', 'M']],
                             columns=['date', 'time', 'var', 'flag',
                                      'oflag'])

        # test with default line terminators "LF" and "CRLF"
        df = self.read_csv(StringIO(data), skiprows=1, delim_whitespace=True,
                           names=['date', 'time', 'var', 'flag', 'oflag'])
        tm.assert_frame_equal(df, expected)

        df = self.read_csv(StringIO(data.replace('\n', '\r\n')),
                           skiprows=1, delim_whitespace=True,
                           names=['date', 'time', 'var', 'flag', 'oflag'])
        tm.assert_frame_equal(df, expected)

        # "CR" is not respected with the Python parser yet
        if self.engine == 'c':
            df = self.read_csv(StringIO(data.replace('\n', '\r')),
                               skiprows=1, delim_whitespace=True,
                               names=['date', 'time', 'var', 'flag', 'oflag'])
            tm.assert_frame_equal(df, expected)
