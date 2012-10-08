from pandas.util.py3compat import StringIO, BytesIO
from datetime import datetime
import csv
import os
import sys
import re
import unittest

import nose

from numpy import nan
import numpy as np

from pandas import DataFrame, Series, Index, isnull, MultiIndex, DatetimeIndex
import pandas.io.parsers as parsers
from pandas.io.parsers import (read_csv, read_table, read_fwf,
                               ExcelFile, TextParser)
from pandas.util.testing import (assert_almost_equal, assert_frame_equal,
                                 assert_series_equal, network)
import pandas.lib as lib
from pandas.util import py3compat
from pandas.lib import Timestamp
from pandas.tseries.index import date_range
import pandas.tseries.tools as tools

from numpy.testing.decorators import slow
from pandas.io.date_converters import (
    parse_date_time, parse_date_fields, parse_all_fields
)

def _skip_if_no_xlrd():
    try:
        import xlrd
    except ImportError:
        raise nose.SkipTest('xlrd not installed, skipping')

def _skip_if_no_openpyxl():
    try:
        import openpyxl
    except ImportError:
        raise nose.SkipTest('openpyxl not installed, skipping')


class TestParsers(unittest.TestCase):
    data1 = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""
    ts_data = """\
ID,date,nominalTime,actualTime,A,B,C,D,E
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000
"""

    def setUp(self):
        self.dirpath = curpath()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')

    def test_empty_string(self):
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
        df = read_csv(StringIO(data))
        xp = DataFrame({'One' : ['a', 'b', np.nan, 'd', 'e', np.nan, 'g'],
                        'Two' : [1,2,3,4,5,6,7],
                        'Three' : ['one', 'two', 'three', np.nan, 'five',
                                   np.nan, 'seven']})
        assert_frame_equal(xp.reindex(columns=df.columns), df)

        df = read_csv(StringIO(data), na_values={'One': [], 'Three': []},
                      keep_default_na=False)
        xp = DataFrame({'One' : ['a', 'b', '', 'd', 'e', 'nan', 'g'],
                        'Two' : [1,2,3,4,5,6,7],
                        'Three' : ['one', 'two', 'three', 'nan', 'five',
                                   '', 'seven']})
        assert_frame_equal(xp.reindex(columns=df.columns), df)

        df = read_csv(StringIO(data), na_values=['a'], keep_default_na=False)
        xp = DataFrame({'One' : [np.nan, 'b', '', 'd', 'e', 'nan', 'g'],
                        'Two' : [1, 2, 3, 4, 5, 6, 7],
                        'Three' : ['one', 'two', 'three', 'nan', 'five', '',
                                   'seven']})
        assert_frame_equal(xp.reindex(columns=df.columns), df)

        df = read_csv(StringIO(data), na_values={'One': [], 'Three': []})
        xp = DataFrame({'One' : ['a', 'b', np.nan, 'd', 'e', np.nan, 'g'],
                        'Two' : [1,2,3,4,5,6,7],
                        'Three' : ['one', 'two', 'three', np.nan, 'five',
                                   np.nan, 'seven']})
        assert_frame_equal(xp.reindex(columns=df.columns), df)


    def test_read_csv(self):
        if not py3compat.PY3:
            fname=u"file:///"+unicode(self.csv1)
            try:
                df1 = read_csv(fname, index_col=0, parse_dates=True)
            except IOError:
                assert(False), "read_csv should accept unicode objects as urls"

    def test_dialect(self):
        data = """\
label1,label2,label3
index1,"a,c,e
index2,b,d,f
"""

        dia = csv.excel()
        dia.quoting = csv.QUOTE_NONE
        df = read_csv(StringIO(data), dialect=dia)

        data = '''\
label1,label2,label3
index1,a,c,e
index2,b,d,f
'''
        exp = read_csv(StringIO(data))
        exp.replace('a', '"a', inplace=True)
        assert_frame_equal(df, exp)

    def test_1000_sep(self):
        data = """A|B|C
1|2,334.0|5
10|13|10.
"""
        expected = [[1, 2334., 5],
                    [10, 13, 10]]

        df = read_csv(StringIO(data), sep='|', thousands=',')
        assert_almost_equal(df.values, expected)

        df = read_table(StringIO(data), sep='|', thousands=',')
        assert_almost_equal(df.values, expected)

    def test_1000_fwf(self):
        data = """
 1 2,334.0    5
10   13     10.
"""
        expected = [[1, 2334., 5],
                    [10, 13, 10]]
        df = read_fwf(StringIO(data), colspecs=[(0,3),(3,11),(12,16)],
                      thousands=',')
        assert_almost_equal(df.values, expected)

    def test_comment(self):
        data = """A,B,C
1,2.,4.#hello world
5.,NaN,10.0
"""
        expected = [[1., 2., 4.],
                    [5., np.nan, 10.]]
        df = read_csv(StringIO(data), comment='#')
        assert_almost_equal(df.values, expected)

        df = read_table(StringIO(data), sep=',', comment='#', na_values=['NaN'])
        assert_almost_equal(df.values, expected)

    def test_comment_fwf(self):
        data = """
  1   2.   4  #hello world
  5  NaN  10.0
"""
        expected = [[1, 2., 4],
                    [5, np.nan, 10.]]
        df = read_fwf(StringIO(data), colspecs=[(0,3),(4,9),(9,25)],
                      comment='#')
        assert_almost_equal(df.values, expected)

    def test_squeeze(self):
        data = """\
a,1
b,2
c,3
"""
        expected = Series([1,2,3], ['a', 'b', 'c'])
        result = read_table(StringIO(data), sep=',', index_col=0,
                            header=None, squeeze=True)
        self.assert_(isinstance(result, Series))
        assert_series_equal(result, expected)

    def test_multiple_date_col(self):
        # Can use multiple date parsers
        data = """\
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000
"""
        def func(*date_cols):
            return lib.try_parse_dates(parsers._concat_date_cols(date_cols))

        df = read_csv(StringIO(data), header=None,
                      date_parser=func,
                      parse_dates={'nominal' : [1, 2],
                                   'actual' : [1,3]})
        self.assert_('nominal' in df)
        self.assert_('actual' in df)
        self.assert_('X1' not in df)
        self.assert_('X2' not in df)
        self.assert_('X3' not in df)

        d = datetime(1999, 1, 27, 19, 0)
        self.assert_(df.ix[0, 'nominal'] == d)

        df = read_csv(StringIO(data), header=None,
                      date_parser=func,
                      parse_dates={'nominal' : [1, 2],
                                     'actual' : [1,3]},
                      keep_date_col=True)
        self.assert_('nominal' in df)
        self.assert_('actual' in df)

        self.assert_('X1' in df)
        self.assert_('X2' in df)
        self.assert_('X3' in df)

        data = """\
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000
"""
        df = read_csv(StringIO(data), header=None,
                      parse_dates=[[1, 2], [1,3]])

        self.assert_('X1_X2' in df)
        self.assert_('X1_X3' in df)
        self.assert_('X1' not in df)
        self.assert_('X2' not in df)
        self.assert_('X3' not in df)

        d = datetime(1999, 1, 27, 19, 0)
        self.assert_(df.ix[0, 'X1_X2'] == d)

        df = read_csv(StringIO(data), header=None,
                      parse_dates=[[1, 2], [1,3]], keep_date_col=True)

        self.assert_('X1_X2' in df)
        self.assert_('X1_X3' in df)
        self.assert_('X1' in df)
        self.assert_('X2' in df)
        self.assert_('X3' in df)

        data = '''\
KORD,19990127 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
'''
        df = read_csv(StringIO(data), sep=',', header=None,
                      parse_dates=[1], index_col=1)
        d = datetime(1999, 1, 27, 19, 0)
        self.assert_(df.index[0] == d)

    def test_multiple_date_cols_int_cast(self):
        data =  ("KORD,19990127, 19:00:00, 18:56:00, 0.8100\n"
                 "KORD,19990127, 20:00:00, 19:56:00, 0.0100\n"
                 "KORD,19990127, 21:00:00, 20:56:00, -0.5900\n"
                 "KORD,19990127, 21:00:00, 21:18:00, -0.9900\n"
                 "KORD,19990127, 22:00:00, 21:56:00, -0.5900\n"
                 "KORD,19990127, 23:00:00, 22:56:00, -0.5900")
        date_spec = {'nominal': [1, 2], 'actual': [1, 3]}
        import pandas.io.date_converters as conv

        # it works!
        df = read_csv(StringIO(data), header=None, parse_dates=date_spec,
                      date_parser=conv.parse_date_time)
        self.assert_('nominal' in df)

    def test_single_line(self):
        df = read_csv(StringIO('1,2'), names=['a', 'b'], sep=None)
        assert_frame_equal(DataFrame({'a': [1], 'b': [2]}), df)

    def test_multiple_date_cols_with_header(self):
        data = """\
ID,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""

        df = read_csv(StringIO(data), parse_dates={'nominal': [1, 2]})
        self.assert_(not isinstance(df.nominal[0], basestring))

    def test_multiple_date_cols_index(self):
        data = """\
ID,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir
KORD1,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD2,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD3,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD4,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD5,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD6,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""

        xp = read_csv(StringIO(data), parse_dates={'nominal': [1, 2]})
        df = read_csv(StringIO(data), parse_dates={'nominal': [1, 2]},
                      index_col='nominal')
        assert_frame_equal(xp.set_index('nominal'), df)
        df2 = read_csv(StringIO(data), parse_dates={'nominal': [1, 2]},
                      index_col=0)
        assert_frame_equal(df2, df)

        df3 = read_csv(StringIO(data), parse_dates=[[1, 2]], index_col=0)
        assert_frame_equal(df3, df)

    def test_multiple_date_cols_chunked(self):
        df = read_csv(StringIO(self.ts_data), parse_dates={'nominal': [1,2]},
                      index_col='nominal')
        reader = read_csv(StringIO(self.ts_data), parse_dates={'nominal': [1,2]},
                          index_col='nominal', chunksize=2)

        chunks = list(reader)

        assert_frame_equal(chunks[0], df[:2])
        assert_frame_equal(chunks[1], df[2:4])
        assert_frame_equal(chunks[2], df[4:])

    def test_multiple_date_col_multiple_index(self):
        df = read_csv(StringIO(self.ts_data), parse_dates={'nominal' : [1, 2]},
                      index_col=['nominal', 'ID'])
        xp = read_csv(StringIO(self.ts_data), parse_dates={'nominal' : [1, 2]})
        assert_frame_equal(xp.set_index(['nominal', 'ID']), df)

    def test_multiple_date_col_name_collision(self):
        self.assertRaises(ValueError, read_csv, StringIO(self.ts_data),
                          parse_dates={'ID' : [1, 2]})

        data = """\
date_NominalTime,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir
KORD1,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD2,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD3,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD4,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD5,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD6,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""

        self.assertRaises(ValueError, read_csv, StringIO(data),
                          parse_dates=[[1, 2]])

    def test_multiple_date_col_named_components(self):
        xp = read_csv(StringIO(self.ts_data), parse_dates={'nominal': [1,2]},
                      index_col='nominal')
        colspec = {'nominal' : ['date', 'nominalTime']}
        df = read_csv(StringIO(self.ts_data), parse_dates=colspec,
                      index_col='nominal')
        assert_frame_equal(df, xp)

    def test_index_col_named(self):
        no_header = """\
KORD1,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD2,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD3,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD4,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD5,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD6,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""

        h = "ID,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir\n"
        data = h + no_header
        #import pdb; pdb.set_trace()
        rs = read_csv(StringIO(data), index_col='ID')
        xp = read_csv(StringIO(data), header=0).set_index('ID')
        assert_frame_equal(rs, xp)

        self.assertRaises(ValueError, read_csv, StringIO(no_header),
                          index_col='ID')

        data = """\
1,2,3,4,hello
5,6,7,8,world
9,10,11,12,foo
"""
        names = ['a', 'b', 'c', 'd', 'message']
        xp = DataFrame({'a' : [1, 5, 9], 'b' : [2, 6, 10], 'c' : [3, 7, 11],
                        'd' : [4, 8, 12]},
                       index=Index(['hello', 'world', 'foo'], name='message'))
        rs = read_csv(StringIO(data), names=names, index_col=['message'])
        assert_frame_equal(xp, rs)
        self.assert_(xp.index.name == rs.index.name)

        rs = read_csv(StringIO(data), names=names, index_col='message')
        assert_frame_equal(xp, rs)
        self.assert_(xp.index.name == rs.index.name)

    def test_converter_index_col_bug(self):
        #1835
        data = "A;B\n1;2\n3;4"
        rs = read_csv(StringIO(data), sep=';', index_col='A',
                      converters={'A' : lambda x: x})
        xp = DataFrame({'B' : [2, 4]}, index=Index([1, 3], name='A'))
        assert_frame_equal(rs, xp)
        self.assert_(rs.index.name == xp.index.name)

    def test_multiple_skts_example(self):
        data = "year, month, a, b\n 2001, 01, 0.0, 10.\n 2001, 02, 1.1, 11."
        pass

    def test_malformed(self):
        # all
        data = """ignore
A,B,C
1,2,3 # comment
1,2,3,4,5
2,3,4
"""

        try:
            df = read_table(StringIO(data), sep=',', header=1, comment='#')
            self.assert_(False)
        except ValueError, inst:
            self.assert_('Expecting 3 columns, got 5 in row 3' in str(inst))

        #skip_footer
        data = """ignore
A,B,C
1,2,3 # comment
1,2,3,4,5
2,3,4
footer
"""

        try:
            df = read_table(StringIO(data), sep=',', header=1, comment='#',
                            skip_footer=1)
            self.assert_(False)
        except ValueError, inst:
            self.assert_('Expecting 3 columns, got 5 in row 3' in str(inst))

        # first chunk
        data = """ignore
A,B,C
skip
1,2,3
3,5,10 # comment
1,2,3,4,5
2,3,4
"""
        try:
            it = read_table(StringIO(data), sep=',',
                            header=1, comment='#', iterator=True, chunksize=1,
                            skiprows=[2])
            df = it.get_chunk(5)
            self.assert_(False)
        except ValueError, inst:
            self.assert_('Expecting 3 columns, got 5 in row 5' in str(inst))


        # middle chunk
        data = """ignore
A,B,C
skip
1,2,3
3,5,10 # comment
1,2,3,4,5
2,3,4
"""
        try:
            it = read_table(StringIO(data), sep=',',
                            header=1, comment='#', iterator=True, chunksize=1,
                            skiprows=[2])
            df = it.get_chunk(1)
            it.get_chunk(2)
            self.assert_(False)
        except ValueError, inst:
            self.assert_('Expecting 3 columns, got 5 in row 5' in str(inst))


        # last chunk
        data = """ignore
A,B,C
skip
1,2,3
3,5,10 # comment
1,2,3,4,5
2,3,4
"""
        try:
            it = read_table(StringIO(data), sep=',',
                            header=1, comment='#', iterator=True, chunksize=1,
                            skiprows=[2])
            df = it.get_chunk(1)
            it.get_chunk()
            self.assert_(False)
        except ValueError, inst:
            self.assert_('Expecting 3 columns, got 5 in row 5' in str(inst))

    def test_quoting(self):
        bad_line_small = """printer\tresult\tvariant_name
Klosterdruckerei\tKlosterdruckerei <Salem> (1611-1804)\tMuller, Jacob
Klosterdruckerei\tKlosterdruckerei <Salem> (1611-1804)\tMuller, Jakob
Klosterdruckerei\tKlosterdruckerei <Kempten> (1609-1805)\t"Furststiftische Hofdruckerei,  <Kempten""
Klosterdruckerei\tKlosterdruckerei <Kempten> (1609-1805)\tGaller, Alois
Klosterdruckerei\tKlosterdruckerei <Kempten> (1609-1805)\tHochfurstliche Buchhandlung <Kempten>"""
        self.assertRaises(Exception, read_table, StringIO(bad_line_small),
                          sep='\t')

        good_line_small = bad_line_small + '"'
        df = read_table(StringIO(good_line_small), sep='\t')
        self.assert_(len(df) == 3)

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

        df3 = read_table(StringIO(data), sep=',', na_values='baz',
                         skiprows=[1])
        assert_almost_equal(df3.values, expected)

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
                             columns=['X1', 'X2', 'X3'],
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
        self.assert_(isinstance(df.index[0], (datetime, np.datetime64, Timestamp)))
        assert_frame_equal(df, expected)

    def test_parse_dates_string(self):
        data = """date,A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        rs = read_csv(StringIO(data), index_col='date', parse_dates='date')
        idx = date_range('1/1/2009', periods=3).asobject
        idx.name = 'date'
        xp = DataFrame({'A': ['a', 'b', 'c'],
                        'B': [1, 3, 4],
                        'C': [2, 4, 5]}, idx)
        assert_frame_equal(rs, xp)

    def test_yy_format(self):
        data = """date,time,B,C
090131,0010,1,2
090228,1020,3,4
090331,0830,5,6
"""
        rs = read_csv(StringIO(data), index_col=0,
                      parse_dates=[['date', 'time']])
        idx = DatetimeIndex([datetime(2009,1,31,0,10,0),
                             datetime(2009,2,28,10,20,0),
                             datetime(2009,3,31,8,30,0)]).asobject
        idx.name = 'date'
        xp = DataFrame({'B': [1, 3, 5], 'C': [2, 4, 6]}, idx)
        assert_frame_equal(rs, xp)

    def test_parse_dates_column_list(self):
        from pandas.core.datetools import to_datetime

        data = '''date;destination;ventilationcode;unitcode;units;aux_date
01/01/2010;P;P;50;1;12/1/2011
01/01/2010;P;R;50;1;13/1/2011
15/01/2010;P;P;50;1;14/1/2011
01/05/2010;P;P;50;1;15/1/2011'''

        expected = read_csv(StringIO(data), sep=";", index_col=range(4))

        lev = expected.index.levels[0]
        expected.index.levels[0] = lev.to_datetime(dayfirst=True)
        expected['aux_date'] = to_datetime(expected['aux_date'],
                                           dayfirst=True)
        expected['aux_date'] = map(Timestamp, expected['aux_date'])
        self.assert_(isinstance(expected['aux_date'][0], datetime))

        df = read_csv(StringIO(data), sep=";", index_col = range(4),
                      parse_dates=[0, 5], dayfirst=True)
        assert_frame_equal(df, expected)

        df = read_csv(StringIO(data), sep=";", index_col = range(4),
                      parse_dates=['date', 'aux_date'], dayfirst=True)
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
                                    ['X0', 'X1', 'X2', 'X3', 'X4']))
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
        self.assert_(isinstance(df.index[0], (datetime, np.datetime64, Timestamp)))
        self.assert_(df.values.dtype == np.float64)
        assert_frame_equal(df, df2)

    def test_read_csv_no_index_name(self):
        df = read_csv(self.csv2, index_col=0, parse_dates=True)
        df2 = read_table(self.csv2, sep=',', index_col=0, parse_dates=True)
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D', 'E']))
        self.assert_(isinstance(df.index[0], (datetime, np.datetime64, Timestamp)))
        self.assert_(df.ix[:, ['A', 'B', 'C', 'D']].values.dtype == np.float64)
        assert_frame_equal(df, df2)

    def test_excel_stop_iterator(self):
        _skip_if_no_xlrd()

        excel_data = ExcelFile(os.path.join(self.dirpath, 'test2.xls'))
        parsed = excel_data.parse('Sheet1')
        expected = DataFrame([['aaaa','bbbbb']], columns=['Test', 'Test1'])
        assert_frame_equal(parsed, expected)

    def test_excel_cell_error_na(self):
        _skip_if_no_xlrd()

        excel_data = ExcelFile(os.path.join(self.dirpath, 'test3.xls'))
        parsed = excel_data.parse('Sheet1')
        expected = DataFrame([[np.nan]], columns=['Test'])
        assert_frame_equal(parsed, expected)

    def test_excel_table(self):
        _skip_if_no_xlrd()

        pth = os.path.join(self.dirpath, 'test.xls')
        xls = ExcelFile(pth)
        df = xls.parse('Sheet1', index_col=0, parse_dates=True)
        df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
        df3 = xls.parse('Sheet2', skiprows=[1], index_col=0, parse_dates=True)
        assert_frame_equal(df, df2)
        assert_frame_equal(df3, df2)

        df4 = xls.parse('Sheet1', index_col=0, parse_dates=True,
                        skipfooter=1)
        df5 = xls.parse('Sheet1', index_col=0, parse_dates=True,
                        skip_footer=1)
        assert_frame_equal(df4, df.ix[:-1])
        assert_frame_equal(df4, df5)

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
        df = xl.parse('Sheet1', index_col=0, parse_dates=True)

    def test_xlsx_table(self):
        _skip_if_no_openpyxl()

        pth = os.path.join(self.dirpath, 'test.xlsx')
        xlsx = ExcelFile(pth)
        df = xlsx.parse('Sheet1', index_col=0, parse_dates=True)
        df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
        df3 = xlsx.parse('Sheet2', skiprows=[1], index_col=0, parse_dates=True)
        assert_frame_equal(df, df2)
        assert_frame_equal(df3, df2)

        df4 = xlsx.parse('Sheet1', index_col=0, parse_dates=True,
                         skipfooter=1)
        df5 = xlsx.parse('Sheet1', index_col=0, parse_dates=True,
                         skip_footer=1)
        assert_frame_equal(df4, df.ix[:-1])
        assert_frame_equal(df4, df5)

    def test_parse_cols_int(self):
        _skip_if_no_openpyxl()
        _skip_if_no_xlrd()

        suffix = ['', 'x']

        for s in suffix:
            pth = os.path.join(self.dirpath, 'test.xls%s' % s)
            xls = ExcelFile(pth)
            df = xls.parse('Sheet1', index_col=0, parse_dates=True,
                            parse_cols=3)
            df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
            df2 = df2.reindex(columns=['A', 'B', 'C'])
            df3 = xls.parse('Sheet2', skiprows=[1], index_col=0,
                            parse_dates=True, parse_cols=3)
            assert_frame_equal(df, df2)
            assert_frame_equal(df3, df2)

    def test_parse_cols_list(self):
        _skip_if_no_openpyxl()
        _skip_if_no_xlrd()

        suffix = ['', 'x']

        for s in suffix:
            pth = os.path.join(self.dirpath, 'test.xls%s' % s)
            xlsx = ExcelFile(pth)
            df = xlsx.parse('Sheet1', index_col=0, parse_dates=True,
                            parse_cols=[0, 2, 3])
            df2 = read_csv(self.csv1, index_col=0, parse_dates=True)
            df2 = df2.reindex(columns=['B', 'C'])
            df3 = xlsx.parse('Sheet2', skiprows=[1], index_col=0,
                             parse_dates=True,
                             parse_cols=[0, 2, 3])
            assert_frame_equal(df, df2)
            assert_frame_equal(df3, df2)

    def test_read_table_unicode(self):
        fin = StringIO('\u0141aski, Jan;1')
        df1 = read_table(fin, sep=";", encoding="utf-8", header=None)
        self.assert_(isinstance(df1['X0'].values[0], unicode))

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

        result = read_csv(StringIO(data), index_col=0)
        expected = read_csv(StringIO(data)).set_index('index',
                                                      verify_integrity=False)
        assert_frame_equal(result, expected)

    def test_read_table_duplicate_index_implicit(self):
        data = """A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo,12,13,14,15
bar,12,13,14,15
"""

        # it works!
        result = read_csv(StringIO(data))

    def test_parse_bools(self):
        data = """A,B
True,1
False,2
True,3
"""
        data = read_csv(StringIO(data))
        self.assert_(data['A'].dtype == np.bool_)

        data = """A,B
YES,1
no,2
yes,3
No,3
Yes,3
"""
        data = read_csv(StringIO(data))
        self.assert_(data['A'].dtype == np.bool_)

        data = """A,B
TRUE,1
FALSE,2
TRUE,3
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

    def test_read_chunksize_named(self):
        reader = read_csv(StringIO(self.data1), index_col='index', chunksize=2)
        df = read_csv(StringIO(self.data1), index_col='index')

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

        df = read_csv(StringIO(data), index_col=['index1', 'index2'])
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
        self.assert_(isinstance(df.index.levels[0][0],
                     (datetime, np.datetime64, Timestamp)))

        # specify columns out of order!
        df2 = read_csv(StringIO(data), index_col=[1, 0], parse_dates=True)
        self.assert_(isinstance(df2.index.levels[1][0],
                     (datetime, np.datetime64, Timestamp)))

    def test_skip_footer(self):
        data = """A,B,C
1,2,3
4,5,6
7,8,9
want to skip this
also also skip this
"""
        result = read_csv(StringIO(data), skip_footer=2)
        no_footer = '\n'.join(data.split('\n')[:-3])
        expected = read_csv(StringIO(no_footer))

        assert_frame_equal(result, expected)

        # equivalent to nrows
        result = read_csv(StringIO(data), nrows=3)
        assert_frame_equal(result, expected)

        # skipfooter alias
        result = read_csv(StringIO(data), skipfooter=2)
        no_footer = '\n'.join(data.split('\n')[:-3])
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

        self.assert_(isinstance(result['D'][0], (datetime, Timestamp)))
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
            self.assert_(buf.getvalue() == 'Filled 1 NA values in column a\n')
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
        expected = DataFrame({'X0' : ['foo', 'bar baz', 'qux foo',
                                      'foo', 'bar']})
        assert_frame_equal(df, expected)

    def test_parse_dates_custom_euroformat(self):
        from dateutil.parser import parse
        text = """foo,bar,baz
31/01/2010,1,2
01/02/2010,1,NA
02/02/2010,1,2
"""
        parser = lambda d: parse(d, dayfirst=True)
        df = read_csv(StringIO(text), skiprows=[0],
                      names=['time', 'Q', 'NTU'], index_col=0,
                      parse_dates=True, date_parser=parser,
                      na_values=['NA'])

        exp_index = Index([datetime(2010, 1, 31), datetime(2010, 2, 1),
                           datetime(2010, 2, 2)], name='time')
        expected = DataFrame({'Q' : [1, 1, 1], 'NTU' : [2, np.nan, 2]},
                             index=exp_index, columns=['Q', 'NTU'])
        assert_frame_equal(df, expected)

        parser = lambda d: parse(d, day_first=True)
        self.assertRaises(Exception, read_csv,
                          StringIO(text), skiprows=[0],
                          names=['time', 'Q', 'NTU'], index_col=0,
                          parse_dates=True, date_parser=parser,
                          na_values=['NA'])

    def test_converters_corner_with_nas(self):
        # skip aberration observed on Win64 Python 3.2.2
        if hash(np.int64(-1)) != -2:
            raise nose.SkipTest

        import StringIO
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

    def test_fwf(self):
        data_expected = """\
2011,58,360.242940,149.910199,11950.7
2011,59,444.953632,166.985655,11788.4
2011,60,364.136849,183.628767,11806.2
2011,61,413.836124,184.375703,11916.8
2011,62,502.953953,173.237159,12468.3
"""
        expected = read_csv(StringIO(data_expected), header=None)

        data1 = """\
201158    360.242940   149.910199   11950.7
201159    444.953632   166.985655   11788.4
201160    364.136849   183.628767   11806.2
201161    413.836124   184.375703   11916.8
201162    502.953953   173.237159   12468.3
"""
        colspecs = [(0, 4), (4, 8), (8, 20), (21, 33), (34, 43)]
        df = read_fwf(StringIO(data1), colspecs=colspecs, header=None)
        assert_frame_equal(df, expected)

        data2 = """\
2011 58   360.242940   149.910199   11950.7
2011 59   444.953632   166.985655   11788.4
2011 60   364.136849   183.628767   11806.2
2011 61   413.836124   184.375703   11916.8
2011 62   502.953953   173.237159   12468.3
"""
        df = read_fwf(StringIO(data2), widths=[5, 5, 13, 13, 7], header=None)
        assert_frame_equal(df, expected)

        # From Thomas Kluyver: apparently some non-space filler characters can
        # be seen, this is supported by specifying the 'delimiter' character:
        # http://publib.boulder.ibm.com/infocenter/dmndhelp/v6r1mx/index.jsp?topic=/com.ibm.wbit.612.help.config.doc/topics/rfixwidth.html
        data3 = """\
201158~~~~360.242940~~~149.910199~~~11950.7
201159~~~~444.953632~~~166.985655~~~11788.4
201160~~~~364.136849~~~183.628767~~~11806.2
201161~~~~413.836124~~~184.375703~~~11916.8
201162~~~~502.953953~~~173.237159~~~12468.3
"""
        df = read_fwf(StringIO(data3), colspecs=colspecs, delimiter='~', header=None)
        assert_frame_equal(df, expected)

        self.assertRaises(ValueError, read_fwf, StringIO(data3),
                          colspecs=colspecs, widths=[6, 10, 10, 7])
    def test_na_value_dict(self):
        data = """A,B,C
foo,bar,NA
bar,foo,foo
foo,bar,NA
bar,foo,foo"""

        df = read_csv(StringIO(data),
                      na_values={'A': ['foo'], 'B': ['bar']})
        expected = DataFrame({'A': [np.nan, 'bar', np.nan, 'bar'],
                              'B': [np.nan, 'foo', np.nan, 'foo'],
                              'C': [np.nan, 'foo', np.nan, 'foo']})
        assert_frame_equal(df, expected)

        data = """\
a,b,c,d
0,NA,1,5
"""
        xp = DataFrame({'b': [np.nan], 'c': [1], 'd': [5]}, index=[0])
        xp.index.name = 'a'
        df = read_csv(StringIO(data), na_values={}, index_col=0)
        assert_frame_equal(df, xp)

        xp = DataFrame({'b': [np.nan], 'd': [5]},
                       MultiIndex.from_tuples([(0, 1)]))
        df = read_csv(StringIO(data), na_values={}, index_col=[0, 2])
        assert_frame_equal(df, xp)

        xp = DataFrame({'b': [np.nan], 'd': [5]},
                       MultiIndex.from_tuples([(0, 1)]))
        df = read_csv(StringIO(data), na_values={}, index_col=['a', 'c'])
        assert_frame_equal(df, xp)

    @slow
    @network
    def test_url(self):
        import urllib2
        try:
            # HTTP(S)
            url = ('https://raw.github.com/pydata/pandas/master/'
                   'pandas/io/tests/salary.table')
            url_table = read_table(url)
            dirpath = curpath()
            localtable = os.path.join(dirpath, 'salary.table')
            local_table = read_table(localtable)
            assert_frame_equal(url_table, local_table)
            #TODO: ftp testing

        except urllib2.URLError:
            try:
                urllib2.urlopen('http://www.google.com')
            except urllib2.URLError:
                raise nose.SkipTest
            else:
                raise

    @slow
    def test_file(self):
        import urllib2

        # FILE
        if sys.version_info[:2] < (2, 6):
            raise nose.SkipTest("file:// not supported with Python < 2.6")
        dirpath = curpath()
        localtable = os.path.join(dirpath, 'salary.table')
        local_table = read_table(localtable)

        try:
            url_table = read_table('file://localhost/'+localtable)
        except urllib2.URLError:
            # fails on some systems
            raise nose.SkipTest

        assert_frame_equal(url_table, local_table)

    def test_parse_tz_aware(self):
        import pytz
        # #1693
        data = StringIO("Date,x\n2012-06-13T01:39:00Z,0.5")

        # it works
        result = read_csv(data, index_col=0, parse_dates=True)
        stamp = result.index[0]
        self.assert_(stamp.minute == 39)
        try:
            self.assert_(result.index.tz is pytz.utc)
        except AssertionError: # hello Yaroslav
            arr = result.index.to_pydatetime()
            result = tools.to_datetime(arr, utc=True)[0]
            self.assert_(stamp.minute == result.minute)
            self.assert_(stamp.hour == result.hour)
            self.assert_(stamp.day == result.day)

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
