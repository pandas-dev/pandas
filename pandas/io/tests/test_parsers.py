# -*- coding: utf-8 -*-
# pylint: disable=E1101

from datetime import datetime
import csv
import os
import sys
import re
import nose
import platform

from numpy import nan
import numpy as np
from pandas.io.common import DtypeWarning

from pandas import DataFrame, Series, Index, MultiIndex, DatetimeIndex
from pandas.compat import(
    StringIO, BytesIO, PY3, range, long, lrange, lmap, u
)
from pandas.io.common import URLError
import pandas.io.parsers as parsers
from pandas.io.parsers import (read_csv, read_table, read_fwf,
                               TextFileReader, TextParser)

import pandas.util.testing as tm
import pandas as pd

from pandas.compat import parse_date
import pandas.lib as lib
from pandas import compat
from pandas.lib import Timestamp
from pandas.tseries.index import date_range
import pandas.tseries.tools as tools

from numpy.testing.decorators import slow
from numpy.testing import assert_array_equal

from pandas.parser import OverflowError


class ParserTests(object):
    """
    Want to be able to test either C+Cython or Python+Cython parsers
    """
    data1 = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""

    def read_csv(self, *args, **kwargs):
        raise NotImplementedError

    def read_table(self, *args, **kwargs):
        raise NotImplementedError

    def setUp(self):
        import warnings
        warnings.filterwarnings(action='ignore', category=FutureWarning)

        self.dirpath = tm.get_data_path()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')

    def test_converters_type_must_be_dict(self):
        with tm.assertRaisesRegexp(TypeError, 'Type converters.+'):
            self.read_csv(StringIO(self.data1), converters=0)

    def test_multi_character_decimal_marker(self):
        data = """A|B|C
1|2,334|5
10|13|10.
"""
        self.assertRaises(ValueError, read_csv, StringIO(data), decimal=',,')

    def test_empty_decimal_marker(self):
        data = """A|B|C
1|2,334|5
10|13|10.
"""
        self.assertRaises(ValueError, read_csv, StringIO(data), decimal='')

    def test_empty_thousands_marker(self):
        data = """A|B|C
1|2,334|5
10|13|10.
"""
        self.assertRaises(ValueError, read_csv, StringIO(data), thousands='')


    def test_multi_character_decimal_marker(self):
        data = """A|B|C
1|2,334|5
10|13|10.
"""
        self.assertRaises(ValueError, read_csv, StringIO(data), thousands=',,')

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


        # GH4318, passing na_values=None and keep_default_na=False yields 'None' as a na_value
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


    def test_read_csv(self):
        if not compat.PY3:
            if 'win' in sys.platform:
                prefix = u("file:///")
            else:
                prefix = u("file://")
            fname = prefix + compat.text_type(self.csv1)
            # it works!
            df1 = read_csv(fname, index_col=0, parse_dates=True)

    def test_dialect(self):
        data = """\
label1,label2,label3
index1,"a,c,e
index2,b,d,f
"""

        dia = csv.excel()
        dia.quoting = csv.QUOTE_NONE
        df = self.read_csv(StringIO(data), dialect=dia)

        data = '''\
label1,label2,label3
index1,a,c,e
index2,b,d,f
'''
        exp = self.read_csv(StringIO(data))
        exp.replace('a', '"a', inplace=True)
        tm.assert_frame_equal(df, exp)

    def test_1000_sep(self):
        data = """A|B|C
1|2,334|5
10|13|10.
"""
        expected = DataFrame({
            'A': [1, 10],
            'B': [2334, 13],
            'C': [5, 10.]
        })

        df = self.read_csv(StringIO(data), sep='|', thousands=',')
        tm.assert_frame_equal(df, expected)

        df = self.read_table(StringIO(data), sep='|', thousands=',')
        tm.assert_frame_equal(df, expected)

    def test_1000_sep_with_decimal(self):
        data = """A|B|C
1|2,334.01|5
10|13|10.
"""
        expected = DataFrame({
            'A': [1, 10],
            'B': [2334.01, 13],
            'C': [5, 10.]
        })

        tm.assert_equal(expected.A.dtype, 'int64')
        tm.assert_equal(expected.B.dtype, 'float')
        tm.assert_equal(expected.C.dtype, 'float')

        df = self.read_csv(StringIO(data), sep='|', thousands=',', decimal='.')
        tm.assert_frame_equal(df, expected)

        df = self.read_table(StringIO(data), sep='|', thousands=',', decimal='.')
        tm.assert_frame_equal(df, expected)

        data_with_odd_sep = """A|B|C
1|2.334,01|5
10|13|10,
"""
        df = self.read_csv(StringIO(data_with_odd_sep), sep='|', thousands='.', decimal=',')
        tm.assert_frame_equal(df, expected)

        df = self.read_table(StringIO(data_with_odd_sep), sep='|', thousands='.', decimal=',')
        tm.assert_frame_equal(df, expected)

    def test_separator_date_conflict(self):
        # Regression test for issue #4678: make sure thousands separator and
        # date parsing do not conflict.
        data = '06-02-2013;13:00;1-000.215'
        expected = DataFrame(
            [[datetime(2013, 6, 2, 13, 0, 0), 1000.215]],
            columns=['Date', 2]
        )

        df = self.read_csv(StringIO(data), sep=';', thousands='-', parse_dates={'Date': [0, 1]}, header=None)
        tm.assert_frame_equal(df, expected)

    def test_squeeze(self):
        data = """\
a,1
b,2
c,3
"""
        expected = Series([1, 2, 3], ['a', 'b', 'c'])
        result = self.read_table(StringIO(data), sep=',', index_col=0,
                                 header=None, squeeze=True)
        tm.assert_isinstance(result, Series)
        tm.assert_series_equal(result, expected)

    def test_inf_parsing(self):
        data = """\
,A
a,inf
b,-inf
c,Inf
d,-Inf
e,INF
f,-INF
g,INf
h,-INf
i,inF
j,-inF"""
        inf = float('inf')
        expected = Series([inf, -inf] * 5)
        df = read_csv(StringIO(data), index_col=0)
        tm.assert_almost_equal(df['A'].values, expected.values)
        df = read_csv(StringIO(data), index_col=0, na_filter=False)
        tm.assert_almost_equal(df['A'].values, expected.values)

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

        df = self.read_csv(StringIO(data), header=None,
                           date_parser=func,
                           prefix='X',
                           parse_dates={'nominal': [1, 2],
                                        'actual': [1, 3]})
        self.assert_('nominal' in df)
        self.assert_('actual' in df)
        self.assert_('X1' not in df)
        self.assert_('X2' not in df)
        self.assert_('X3' not in df)

        d = datetime(1999, 1, 27, 19, 0)
        self.assert_(df.ix[0, 'nominal'] == d)

        df = self.read_csv(StringIO(data), header=None,
                           date_parser=func,
                           parse_dates={'nominal': [1, 2],
                                        'actual': [1, 3]},
                           keep_date_col=True)
        self.assert_('nominal' in df)
        self.assert_('actual' in df)

        self.assert_(1 in df)
        self.assert_(2 in df)
        self.assert_(3 in df)

        data = """\
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000
"""
        df = read_csv(StringIO(data), header=None,
                      prefix='X',
                      parse_dates=[[1, 2], [1, 3]])

        self.assert_('X1_X2' in df)
        self.assert_('X1_X3' in df)
        self.assert_('X1' not in df)
        self.assert_('X2' not in df)
        self.assert_('X3' not in df)

        d = datetime(1999, 1, 27, 19, 0)
        self.assert_(df.ix[0, 'X1_X2'] == d)

        df = read_csv(StringIO(data), header=None,
                      parse_dates=[[1, 2], [1, 3]], keep_date_col=True)

        self.assert_('1_2' in df)
        self.assert_('1_3' in df)
        self.assert_(1 in df)
        self.assert_(2 in df)
        self.assert_(3 in df)

        data = '''\
KORD,19990127 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
'''
        df = self.read_csv(StringIO(data), sep=',', header=None,
                           parse_dates=[1], index_col=1)
        d = datetime(1999, 1, 27, 19, 0)
        self.assert_(df.index[0] == d)

    def test_multiple_date_cols_int_cast(self):
        data = ("KORD,19990127, 19:00:00, 18:56:00, 0.8100\n"
                "KORD,19990127, 20:00:00, 19:56:00, 0.0100\n"
                "KORD,19990127, 21:00:00, 20:56:00, -0.5900\n"
                "KORD,19990127, 21:00:00, 21:18:00, -0.9900\n"
                "KORD,19990127, 22:00:00, 21:56:00, -0.5900\n"
                "KORD,19990127, 23:00:00, 22:56:00, -0.5900")
        date_spec = {'nominal': [1, 2], 'actual': [1, 3]}
        import pandas.io.date_converters as conv

        # it works!
        df = self.read_csv(StringIO(data), header=None, parse_dates=date_spec,
                           date_parser=conv.parse_date_time)
        self.assert_('nominal' in df)

    def test_multiple_date_col_timestamp_parse(self):
        data = """05/31/2012,15:30:00.029,1306.25,1,E,0,,1306.25
05/31/2012,15:30:00.029,1306.25,8,E,0,,1306.25"""
        result = self.read_csv(StringIO(data), sep=',', header=None,
                               parse_dates=[[0,1]], date_parser=Timestamp)

        ex_val = Timestamp('05/31/2012 15:30:00.029')
        self.assertEqual(result['0_1'][0], ex_val)

    def test_single_line(self):
        # sniff separator
        buf = StringIO()
        sys.stdout = buf

        # printing warning message when engine == 'c' for now

        try:
            # it works!
            df = self.read_csv(StringIO('1,2'), names=['a', 'b'],
                               header=None, sep=None)
            tm.assert_frame_equal(DataFrame({'a': [1], 'b': [2]}), df)
        finally:
            sys.stdout = sys.__stdout__

    def test_multiple_date_cols_with_header(self):
        data = """\
ID,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""

        df = self.read_csv(StringIO(data), parse_dates={'nominal': [1, 2]})
        self.assert_(not isinstance(df.nominal[0], compat.string_types))

    ts_data = """\
ID,date,nominalTime,actualTime,A,B,C,D,E
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000
"""

    def test_multiple_date_col_name_collision(self):
        self.assertRaises(ValueError, self.read_csv, StringIO(self.ts_data),
                          parse_dates={'ID': [1, 2]})

        data = """\
date_NominalTime,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir
KORD1,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD2,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD3,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD4,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD5,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD6,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""

        self.assertRaises(ValueError, self.read_csv, StringIO(data),
                          parse_dates=[[1, 2]])

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
        # import pdb; pdb.set_trace()
        rs = self.read_csv(StringIO(data), index_col='ID')
        xp = self.read_csv(StringIO(data), header=0).set_index('ID')
        tm.assert_frame_equal(rs, xp)

        self.assertRaises(ValueError, self.read_csv, StringIO(no_header),
                          index_col='ID')

        data = """\
1,2,3,4,hello
5,6,7,8,world
9,10,11,12,foo
"""
        names = ['a', 'b', 'c', 'd', 'message']
        xp = DataFrame({'a': [1, 5, 9], 'b': [2, 6, 10], 'c': [3, 7, 11],
                        'd': [4, 8, 12]},
                       index=Index(['hello', 'world', 'foo'], name='message'))
        rs = self.read_csv(StringIO(data), names=names, index_col=['message'])
        tm.assert_frame_equal(xp, rs)
        self.assert_(xp.index.name == rs.index.name)

        rs = self.read_csv(StringIO(data), names=names, index_col='message')
        tm.assert_frame_equal(xp, rs)
        self.assert_(xp.index.name == rs.index.name)

    def test_converter_index_col_bug(self):
        # 1835
        data = "A;B\n1;2\n3;4"

        rs = self.read_csv(StringIO(data), sep=';', index_col='A',
                           converters={'A': lambda x: x})

        xp = DataFrame({'B': [2, 4]}, index=Index([1, 3], name='A'))
        tm.assert_frame_equal(rs, xp)
        self.assert_(rs.index.name == xp.index.name)

    def test_date_parser_int_bug(self):
        # #3071
        log_file = StringIO(
            'posix_timestamp,elapsed,sys,user,queries,query_time,rows,'
                'accountid,userid,contactid,level,silo,method\n'
            '1343103150,0.062353,0,4,6,0.01690,3,'
                '12345,1,-1,3,invoice_InvoiceResource,search\n'
        )

        def f(posix_string):
            return datetime.utcfromtimestamp(int(posix_string))

        # it works!
        read_csv(log_file, index_col=0, parse_dates=0, date_parser=f)

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
            df = self.read_table(
                StringIO(data), sep=',', header=1, comment='#')
            self.assert_(False)
        except Exception as inst:
            self.assert_('Expected 3 fields in line 4, saw 5' in str(inst))

        # skip_footer
        data = """ignore
A,B,C
1,2,3 # comment
1,2,3,4,5
2,3,4
footer
"""

        try:
            df = self.read_table(
                StringIO(data), sep=',', header=1, comment='#',
                skip_footer=1)
            self.assert_(False)
        except Exception as inst:
            self.assert_('Expected 3 fields in line 4, saw 5' in str(inst))

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
            it = self.read_table(StringIO(data), sep=',',
                                 header=1, comment='#', iterator=True, chunksize=1,
                                 skiprows=[2])
            df = it.read(5)
            self.assert_(False)
        except Exception as inst:
            self.assert_('Expected 3 fields in line 6, saw 5' in str(inst))

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
            it = self.read_table(StringIO(data), sep=',', header=1,
                                 comment='#', iterator=True, chunksize=1,
                                 skiprows=[2])
            df = it.read(1)
            it.read(2)
            self.assert_(False)
        except Exception as inst:
            self.assert_('Expected 3 fields in line 6, saw 5' in str(inst))

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
            it = self.read_table(StringIO(data), sep=',',
                                 header=1, comment='#', iterator=True, chunksize=1,
                                 skiprows=[2])
            df = it.read(1)
            it.read()
            self.assert_(False)
        except Exception as inst:
            self.assert_('Expected 3 fields in line 6, saw 5' in str(inst))

    def test_passing_dtype(self):

        df = DataFrame(np.random.rand(5,2),columns=list('AB'),index=['1A','1B','1C','1D','1E'])

        with tm.ensure_clean('__passing_str_as_dtype__.csv') as path:
            df.to_csv(path)

            # GH 3795
            # passing 'str' as the dtype
            result = pd.read_csv(path, dtype=str, index_col=0)
            tm.assert_series_equal(result.dtypes,Series({ 'A' : 'object', 'B' : 'object' }))

            # we expect all object columns, so need to convert to test for equivalence
            result = result.astype(float)
            tm.assert_frame_equal(result,df)

            # invalid dtype
            self.assertRaises(TypeError, pd.read_csv, path, dtype={'A' : 'foo', 'B' : 'float64' },
                              index_col=0)

            # valid but we don't support it (date)
            self.assertRaises(TypeError, pd.read_csv, path, dtype={'A' : 'datetime64', 'B' : 'float64' },
                              index_col=0)
            self.assertRaises(TypeError, pd.read_csv, path, dtype={'A' : 'datetime64', 'B' : 'float64' },
                              index_col=0, parse_dates=['B'])

            # valid but we don't support it
            self.assertRaises(TypeError, pd.read_csv, path, dtype={'A' : 'timedelta64', 'B' : 'float64' },
                              index_col=0)

    def test_quoting(self):
        bad_line_small = """printer\tresult\tvariant_name
Klosterdruckerei\tKlosterdruckerei <Salem> (1611-1804)\tMuller, Jacob
Klosterdruckerei\tKlosterdruckerei <Salem> (1611-1804)\tMuller, Jakob
Klosterdruckerei\tKlosterdruckerei <Kempten> (1609-1805)\t"Furststiftische Hofdruckerei,  <Kempten""
Klosterdruckerei\tKlosterdruckerei <Kempten> (1609-1805)\tGaller, Alois
Klosterdruckerei\tKlosterdruckerei <Kempten> (1609-1805)\tHochfurstliche Buchhandlung <Kempten>"""
        self.assertRaises(Exception, self.read_table, StringIO(bad_line_small),
                          sep='\t')

        good_line_small = bad_line_small + '"'
        df = self.read_table(StringIO(good_line_small), sep='\t')
        self.assert_(len(df) == 3)

    def test_non_string_na_values(self):
        # GH3611, na_values that are not a string are an issue
        with tm.ensure_clean('__non_string_na_values__.csv') as path:
            df = DataFrame({'A' : [-999, 2, 3], 'B' : [1.2, -999, 4.5]})
            df.to_csv(path, sep=' ', index=False)
            result1 = read_csv(path, sep= ' ', header=0, na_values=['-999.0','-999'])
            result2 = read_csv(path, sep= ' ', header=0, na_values=[-999,-999.0])
            result3 = read_csv(path, sep= ' ', header=0, na_values=[-999.0,-999])
            tm.assert_frame_equal(result1,result2)
            tm.assert_frame_equal(result2,result3)

            result4 = read_csv(path, sep= ' ', header=0, na_values=['-999.0'])
            result5 = read_csv(path, sep= ' ', header=0, na_values=['-999'])
            result6 = read_csv(path, sep= ' ', header=0, na_values=[-999.0])
            result7 = read_csv(path, sep= ' ', header=0, na_values=[-999])
            tm.assert_frame_equal(result4,result3)
            tm.assert_frame_equal(result5,result3)
            tm.assert_frame_equal(result6,result3)
            tm.assert_frame_equal(result7,result3)

            good_compare = result3

            # with an odd float format, so we can't match the string 999.0 exactly,
            # but need float matching
            df.to_csv(path, sep=' ', index=False, float_format = '%.3f')
            result1 = read_csv(path, sep= ' ', header=0, na_values=['-999.0','-999'])
            result2 = read_csv(path, sep= ' ', header=0, na_values=[-999,-999.0])
            result3 = read_csv(path, sep= ' ', header=0, na_values=[-999.0,-999])
            tm.assert_frame_equal(result1,good_compare)
            tm.assert_frame_equal(result2,good_compare)
            tm.assert_frame_equal(result3,good_compare)

            result4 = read_csv(path, sep= ' ', header=0, na_values=['-999.0'])
            result5 = read_csv(path, sep= ' ', header=0, na_values=['-999'])
            result6 = read_csv(path, sep= ' ', header=0, na_values=[-999.0])
            result7 = read_csv(path, sep= ' ', header=0, na_values=[-999])
            tm.assert_frame_equal(result4,good_compare)
            tm.assert_frame_equal(result5,good_compare)
            tm.assert_frame_equal(result6,good_compare)
            tm.assert_frame_equal(result7,good_compare)

    def test_default_na_values(self):
        _NA_VALUES = set(['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN',
                          '#N/A','N/A', 'NA', '#NA', 'NULL', 'NaN',
                          'nan', '-NaN', '-nan', ''])
        assert_array_equal (_NA_VALUES, parsers._NA_VALUES)
        nv = len(_NA_VALUES)
        def f(i, v):
            if i == 0:
                buf = ''
            elif i > 0:
                buf = ''.join([','] * i)

            buf = "{0}{1}".format(buf,v)

            if i < nv-1:
                buf = "{0}{1}".format(buf,''.join([','] * (nv-i-1)))

            return buf

        data = StringIO('\n'.join([ f(i, v) for i, v in enumerate(_NA_VALUES) ]))

        expected = DataFrame(np.nan,columns=range(nv),index=range(nv))
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

    def test_nat_parse(self):

        # GH 3062
        df = DataFrame(dict({
                    'A' : np.asarray(lrange(10),dtype='float64'),
                    'B' : pd.Timestamp('20010101') }))
        df.iloc[3:6,:] = np.nan

        with tm.ensure_clean('__nat_parse_.csv') as path:
            df.to_csv(path)
            result = read_csv(path,index_col=0,parse_dates=['B'])
            tm.assert_frame_equal(result,df)

            expected = Series(dict( A = 'float64',B = 'datetime64[ns]'))
            tm.assert_series_equal(expected,result.dtypes)

            # test with NaT for the nan_rep
            # we don't have a method to specif the Datetime na_rep (it defaults to '')
            df.to_csv(path)
            result = read_csv(path,index_col=0,parse_dates=['B'])
            tm.assert_frame_equal(result,df)

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
        data = self.read_csv(StringIO(text), skiprows=lrange(6), header=None,
                             index_col=0, parse_dates=True)

        data2 = self.read_csv(StringIO(text), skiprows=6, header=None,
                              index_col=0, parse_dates=True)

        expected = DataFrame(np.arange(1., 10.).reshape((3, 3)),
                             columns=[1, 2, 3],
                             index=[datetime(2000, 1, 1), datetime(2000, 1, 2),
                                    datetime(2000, 1, 3)])
        expected.index.name = 0
        tm.assert_frame_equal(data, expected)
        tm.assert_frame_equal(data, data2)

    def test_deep_skiprows(self):
        # GH #4382
        text = "a,b,c\n" + "\n".join([",".join([str(i), str(i+1), str(i+2)]) for i in range(10)])
        condensed_text = "a,b,c\n" + "\n".join([",".join([str(i), str(i+1), str(i+2)]) for i in [0, 1, 2, 3, 4, 6, 8, 9]])
        data = self.read_csv(StringIO(text), skiprows=[6, 8])
        condensed_data = self.read_csv(StringIO(condensed_text))
        tm.assert_frame_equal(data, condensed_data)

    def test_detect_string_na(self):
        data = """A,B
foo,bar
NA,baz
NaN,nan
"""
        expected = [['foo', 'bar'],
                    [nan, 'baz'],
                    [nan, nan]]

        df = self.read_csv(StringIO(data))
        tm.assert_almost_equal(df.values, expected)

    def test_unnamed_columns(self):
        data = """A,B,C,,
1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        expected = [[1, 2, 3, 4, 5.],
                    [6, 7, 8, 9, 10],
                    [11, 12, 13, 14, 15]]
        df = self.read_table(StringIO(data), sep=',')
        tm.assert_almost_equal(df.values, expected)
        self.assert_(np.array_equal(df.columns,
                                    ['A', 'B', 'C', 'Unnamed: 3',
                                     'Unnamed: 4']))

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

    def test_duplicate_columns(self):
        for engine in ['python', 'c']:
            data = """A,A,B,B,B
    1,2,3,4,5
    6,7,8,9,10
    11,12,13,14,15
    """
            # check default beahviour
            df = self.read_table(StringIO(data), sep=',',engine=engine)
            self.assertEqual(list(df.columns), ['A', 'A.1', 'B', 'B.1', 'B.2'])

            df = self.read_table(StringIO(data), sep=',',engine=engine,mangle_dupe_cols=False)
            self.assertEqual(list(df.columns), ['A', 'A', 'B', 'B', 'B'])

            df = self.read_table(StringIO(data), sep=',',engine=engine,mangle_dupe_cols=True)
            self.assertEqual(list(df.columns), ['A', 'A.1', 'B', 'B.1', 'B.2'])

    def test_csv_mixed_type(self):
        data = """A,B,C
a,1,2
b,3,4
c,4,5
"""
        df = self.read_csv(StringIO(data))
        # TODO

    def test_csv_custom_parser(self):
        data = """A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        f = lambda x: datetime.strptime(x, '%Y%m%d')
        df = self.read_csv(StringIO(data), date_parser=f)
        expected = self.read_csv(StringIO(data), parse_dates=True)
        tm.assert_frame_equal(df, expected)

    def test_parse_dates_implicit_first_col(self):
        data = """A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        df = self.read_csv(StringIO(data), parse_dates=True)
        expected = self.read_csv(StringIO(data), index_col=0, parse_dates=True)
        self.assert_(
            isinstance(df.index[0], (datetime, np.datetime64, Timestamp)))
        tm.assert_frame_equal(df, expected)

    def test_parse_dates_string(self):
        data = """date,A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        rs = self.read_csv(
            StringIO(data), index_col='date', parse_dates='date')
        idx = date_range('1/1/2009', periods=3)
        idx.name = 'date'
        xp = DataFrame({'A': ['a', 'b', 'c'],
                        'B': [1, 3, 4],
                        'C': [2, 4, 5]}, idx)
        tm.assert_frame_equal(rs, xp)

    def test_yy_format(self):
        data = """date,time,B,C
090131,0010,1,2
090228,1020,3,4
090331,0830,5,6
"""
        rs = self.read_csv(StringIO(data), index_col=0,
                           parse_dates=[['date', 'time']])
        idx = DatetimeIndex([datetime(2009, 1, 31, 0, 10, 0),
                             datetime(2009, 2, 28, 10, 20, 0),
                             datetime(2009, 3, 31, 8, 30, 0)]).asobject
        idx.name = 'date_time'
        xp = DataFrame({'B': [1, 3, 5], 'C': [2, 4, 6]}, idx)
        tm.assert_frame_equal(rs, xp)

        rs = self.read_csv(StringIO(data), index_col=0,
                           parse_dates=[[0, 1]])
        idx = DatetimeIndex([datetime(2009, 1, 31, 0, 10, 0),
                             datetime(2009, 2, 28, 10, 20, 0),
                             datetime(2009, 3, 31, 8, 30, 0)]).asobject
        idx.name = 'date_time'
        xp = DataFrame({'B': [1, 3, 5], 'C': [2, 4, 6]}, idx)
        tm.assert_frame_equal(rs, xp)

    def test_parse_dates_column_list(self):
        from pandas.core.datetools import to_datetime

        data = '''date;destination;ventilationcode;unitcode;units;aux_date
01/01/2010;P;P;50;1;12/1/2011
01/01/2010;P;R;50;1;13/1/2011
15/01/2010;P;P;50;1;14/1/2011
01/05/2010;P;P;50;1;15/1/2011'''

        expected = self.read_csv(StringIO(data), sep=";", index_col=lrange(4))

        lev = expected.index.levels[0]
        levels = list(expected.index.levels)
        levels[0] = lev.to_datetime(dayfirst=True)
        # hack to get this to work - remove for final test
        levels[0].name = lev.name
        expected.index.set_levels(levels, inplace=True)
        expected['aux_date'] = to_datetime(expected['aux_date'],
                                           dayfirst=True)
        expected['aux_date'] = lmap(Timestamp, expected['aux_date'])
        tm.assert_isinstance(expected['aux_date'][0], datetime)

        df = self.read_csv(StringIO(data), sep=";", index_col=lrange(4),
                           parse_dates=[0, 5], dayfirst=True)
        tm.assert_frame_equal(df, expected)

        df = self.read_csv(StringIO(data), sep=";", index_col=lrange(4),
                           parse_dates=['date', 'aux_date'], dayfirst=True)
        tm.assert_frame_equal(df, expected)

    def test_no_header(self):
        data = """1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        df = self.read_table(StringIO(data), sep=',', header=None)
        df_pref = self.read_table(StringIO(data), sep=',', prefix='X',
                                  header=None)

        names = ['foo', 'bar', 'baz', 'quux', 'panda']
        df2 = self.read_table(StringIO(data), sep=',', names=names)
        expected = [[1, 2, 3, 4, 5.],
                    [6, 7, 8, 9, 10],
                    [11, 12, 13, 14, 15]]
        tm.assert_almost_equal(df.values, expected)
        tm.assert_almost_equal(df.values, df2.values)

        self.assert_(np.array_equal(df_pref.columns,
                                    ['X0', 'X1', 'X2', 'X3', 'X4']))
        self.assert_(np.array_equal(df.columns, lrange(5)))

        self.assert_(np.array_equal(df2.columns, names))

    def test_no_header_prefix(self):
        data = """1,2,3,4,5
6,7,8,9,10
11,12,13,14,15
"""
        df_pref = self.read_table(StringIO(data), sep=',', prefix='Field',
                                  header=None)

        expected = [[1, 2, 3, 4, 5.],
                    [6, 7, 8, 9, 10],
                    [11, 12, 13, 14, 15]]
        tm.assert_almost_equal(df_pref.values, expected)

        self.assert_(np.array_equal(df_pref.columns,
                                    ['Field0', 'Field1', 'Field2', 'Field3', 'Field4']))

    def test_header_with_index_col(self):
        data = """foo,1,2,3
bar,4,5,6
baz,7,8,9
"""
        names = ['A', 'B', 'C']
        df = self.read_csv(StringIO(data), names=names)

        self.assertEqual(names, ['A', 'B', 'C'])

        values = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        expected = DataFrame(values, index=['foo', 'bar', 'baz'],
                             columns=['A', 'B', 'C'])
        tm.assert_frame_equal(df, expected)

    def test_read_csv_dataframe(self):
        df = self.read_csv(self.csv1, index_col=0, parse_dates=True)
        df2 = self.read_table(self.csv1, sep=',', index_col=0,
                              parse_dates=True)
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D']))
        self.assert_(df.index.name == 'index')
        self.assert_(isinstance(df.index[0], (datetime, np.datetime64,
                                              Timestamp)))
        self.assert_(df.values.dtype == np.float64)
        tm.assert_frame_equal(df, df2)

    def test_read_csv_no_index_name(self):
        df = self.read_csv(self.csv2, index_col=0, parse_dates=True)
        df2 = self.read_table(self.csv2, sep=',', index_col=0,
                              parse_dates=True)
        self.assert_(np.array_equal(df.columns, ['A', 'B', 'C', 'D', 'E']))
        self.assert_(isinstance(df.index[0], (datetime, np.datetime64,
                                              Timestamp)))
        self.assert_(df.ix[:, ['A', 'B', 'C', 'D']].values.dtype == np.float64)
        tm.assert_frame_equal(df, df2)

    def test_read_table_unicode(self):
        fin = BytesIO(u('\u0141aski, Jan;1').encode('utf-8'))
        df1 = read_table(fin, sep=";", encoding="utf-8", header=None)
        tm.assert_isinstance(df1[0].values[0], compat.text_type)

    def test_read_table_wrong_num_columns(self):
        # too few!
        data = """A,B,C,D,E,F
1,2,3,4,5,6
6,7,8,9,10,11,12
11,12,13,14,15,16
"""
        self.assertRaises(Exception, self.read_csv, StringIO(data))

    def test_read_table_duplicate_index(self):
        data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo,12,13,14,15
bar,12,13,14,15
"""

        result = self.read_csv(StringIO(data), index_col=0)
        expected = self.read_csv(StringIO(data)).set_index('index',
                                                           verify_integrity=False)
        tm.assert_frame_equal(result, expected)

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
        result = self.read_csv(StringIO(data))

    def test_parse_bools(self):
        data = """A,B
True,1
False,2
True,3
"""
        data = self.read_csv(StringIO(data))
        self.assert_(data['A'].dtype == np.bool_)

        data = """A,B
YES,1
no,2
yes,3
No,3
Yes,3
"""
        data = self.read_csv(StringIO(data),
                             true_values=['yes', 'Yes', 'YES'],
                             false_values=['no', 'NO', 'No'])
        self.assert_(data['A'].dtype == np.bool_)

        data = """A,B
TRUE,1
FALSE,2
TRUE,3
"""
        data = self.read_csv(StringIO(data))
        self.assert_(data['A'].dtype == np.bool_)

        data = """A,B
foo,bar
bar,foo"""
        result = self.read_csv(StringIO(data), true_values=['foo'],
                               false_values=['bar'])
        expected = DataFrame({'A': [True, False], 'B': [False, True]})
        tm.assert_frame_equal(result, expected)

    def test_int_conversion(self):
        data = """A,B
1.0,1
2.0,2
3.0,3
"""
        data = self.read_csv(StringIO(data))
        self.assert_(data['A'].dtype == np.float64)
        self.assert_(data['B'].dtype == np.int64)

    def test_infer_index_col(self):
        data = """A,B,C
foo,1,2,3
bar,4,5,6
baz,7,8,9
"""
        data = self.read_csv(StringIO(data))
        self.assert_(data.index.equals(Index(['foo', 'bar', 'baz'])))

    def test_read_nrows(self):
        df = self.read_csv(StringIO(self.data1), nrows=3)
        expected = self.read_csv(StringIO(self.data1))[:3]
        tm.assert_frame_equal(df, expected)

    def test_read_chunksize(self):
        reader = self.read_csv(StringIO(self.data1), index_col=0, chunksize=2)
        df = self.read_csv(StringIO(self.data1), index_col=0)

        chunks = list(reader)

        tm.assert_frame_equal(chunks[0], df[:2])
        tm.assert_frame_equal(chunks[1], df[2:4])
        tm.assert_frame_equal(chunks[2], df[4:])

    def test_read_chunksize_named(self):
        reader = self.read_csv(
            StringIO(self.data1), index_col='index', chunksize=2)
        df = self.read_csv(StringIO(self.data1), index_col='index')

        chunks = list(reader)

        tm.assert_frame_equal(chunks[0], df[:2])
        tm.assert_frame_equal(chunks[1], df[2:4])
        tm.assert_frame_equal(chunks[2], df[4:])

    def test_get_chunk_passed_chunksize(self):
        data = """A,B,C
1,2,3
4,5,6
7,8,9
1,2,3"""
        result = self.read_csv(StringIO(data), chunksize=2)

        piece = result.get_chunk()
        self.assertEqual(len(piece), 2)

    def test_read_text_list(self):
        data = """A,B,C\nfoo,1,2,3\nbar,4,5,6"""
        as_list = [['A', 'B', 'C'], ['foo', '1', '2', '3'], ['bar',
                                                             '4', '5', '6']]
        df = self.read_csv(StringIO(data), index_col=0)

        parser = TextParser(as_list, index_col=0, chunksize=2)
        chunk = parser.read(None)

        tm.assert_frame_equal(chunk, df)

    def test_iterator(self):
        reader = self.read_csv(StringIO(self.data1), index_col=0,
                               iterator=True)
        df = self.read_csv(StringIO(self.data1), index_col=0)

        chunk = reader.read(3)
        tm.assert_frame_equal(chunk, df[:3])

        last_chunk = reader.read(5)
        tm.assert_frame_equal(last_chunk, df[3:])

        # pass list
        lines = list(csv.reader(StringIO(self.data1)))
        parser = TextParser(lines, index_col=0, chunksize=2)

        df = self.read_csv(StringIO(self.data1), index_col=0)

        chunks = list(parser)
        tm.assert_frame_equal(chunks[0], df[:2])
        tm.assert_frame_equal(chunks[1], df[2:4])
        tm.assert_frame_equal(chunks[2], df[4:])

        # pass skiprows
        parser = TextParser(lines, index_col=0, chunksize=2, skiprows=[1])
        chunks = list(parser)
        tm.assert_frame_equal(chunks[0], df[1:3])

        # test bad parameter (skip_footer)
        reader = self.read_csv(StringIO(self.data1), index_col=0,
                               iterator=True, skip_footer=True)
        self.assertRaises(ValueError, reader.read, 3)

        treader = self.read_table(StringIO(self.data1), sep=',', index_col=0,
                                  iterator=True)
        tm.assert_isinstance(treader, TextFileReader)

        # stopping iteration when on chunksize is specified, GH 3967
        data = """A,B,C
foo,1,2,3
bar,4,5,6
baz,7,8,9
"""
        reader = self.read_csv(StringIO(data), iterator=True)
        result = list(reader)
        expected = DataFrame(dict(A = [1,4,7], B = [2,5,8], C = [3,6,9]), index=['foo','bar','baz'])
        tm.assert_frame_equal(result[0], expected)

        # chunksize = 1
        reader = self.read_csv(StringIO(data), chunksize=1)
        result = list(reader)
        expected = DataFrame(dict(A = [1,4,7], B = [2,5,8], C = [3,6,9]), index=['foo','bar','baz'])
        self.assert_(len(result) == 3)
        tm.assert_frame_equal(pd.concat(result), expected)

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

        df = self.read_csv(StringIO(data), header=2, index_col=0)
        expected = self.read_csv(StringIO(data2), header=0, index_col=0)
        tm.assert_frame_equal(df, expected)

    def test_header_multi_index(self):
        expected = tm.makeCustomDataframe(5,3,r_idx_nlevels=2,c_idx_nlevels=4)

        data = """\
C0,,C_l0_g0,C_l0_g1,C_l0_g2

C1,,C_l1_g0,C_l1_g1,C_l1_g2
C2,,C_l2_g0,C_l2_g1,C_l2_g2
C3,,C_l3_g0,C_l3_g1,C_l3_g2
R0,R1,,,
R_l0_g0,R_l1_g0,R0C0,R0C1,R0C2
R_l0_g1,R_l1_g1,R1C0,R1C1,R1C2
R_l0_g2,R_l1_g2,R2C0,R2C1,R2C2
R_l0_g3,R_l1_g3,R3C0,R3C1,R3C2
R_l0_g4,R_l1_g4,R4C0,R4C1,R4C2
"""

        df = self.read_csv(StringIO(data), header=[0, 2, 3, 4], index_col=[0, 1], tupleize_cols=False)
        tm.assert_frame_equal(df, expected)

        # skipping lines in the header
        df = self.read_csv(StringIO(data), header=[0, 2, 3, 4], index_col=[0, 1], tupleize_cols=False)
        tm.assert_frame_equal(df, expected)

        #### invalid options ####

        # no as_recarray
        self.assertRaises(ValueError, self.read_csv, StringIO(data), header=[0,1,2,3],
                          index_col=[0,1], as_recarray=True, tupleize_cols=False)

        # names
        self.assertRaises(ValueError, self.read_csv, StringIO(data), header=[0,1,2,3],
                          index_col=[0,1], names=['foo','bar'], tupleize_cols=False)
        # usecols
        self.assertRaises(ValueError, self.read_csv, StringIO(data), header=[0,1,2,3],
                          index_col=[0,1], usecols=['foo','bar'], tupleize_cols=False)
        # non-numeric index_col
        self.assertRaises(ValueError, self.read_csv, StringIO(data), header=[0,1,2,3],
                          index_col=['foo','bar'], tupleize_cols=False)

    def test_header_multiindex_common_format(self):

        df = DataFrame([[1,2,3,4,5,6],[7,8,9,10,11,12]],
                       index=['one','two'],
                       columns=MultiIndex.from_tuples([('a','q'),('a','r'),('a','s'),
                                                       ('b','t'),('c','u'),('c','v')]))

        # to_csv
        data = """,a,a,a,b,c,c
,q,r,s,t,u,v
,,,,,,
one,1,2,3,4,5,6
two,7,8,9,10,11,12"""

        result = self.read_csv(StringIO(data),header=[0,1],index_col=0)
        tm.assert_frame_equal(df,result)

        # common
        data = """,a,a,a,b,c,c
,q,r,s,t,u,v
one,1,2,3,4,5,6
two,7,8,9,10,11,12"""

        result = self.read_csv(StringIO(data),header=[0,1],index_col=0)
        tm.assert_frame_equal(df,result)

        # common, no index_col
        data = """a,a,a,b,c,c
q,r,s,t,u,v
1,2,3,4,5,6
7,8,9,10,11,12"""

        result = self.read_csv(StringIO(data),header=[0,1],index_col=None)
        tm.assert_frame_equal(df.reset_index(drop=True),result)

        # malformed case 1
        expected = DataFrame(np.array([[2,  3,  4,  5,  6],
                                       [8,  9, 10, 11, 12]], dtype='int64'),
                             index=Index([1, 7]),
                             columns=MultiIndex(levels=[[u('a'), u('b'), u('c')], [u('r'), u('s'), u('t'), u('u'), u('v')]],
                                                labels=[[0, 0, 1, 2, 2], [0, 1, 2, 3, 4]],
                                                names=[u('a'), u('q')]))

        data = """a,a,a,b,c,c
q,r,s,t,u,v
1,2,3,4,5,6
7,8,9,10,11,12"""

        result = self.read_csv(StringIO(data),header=[0,1],index_col=0)
        tm.assert_frame_equal(expected,result)

        # malformed case 2
        expected = DataFrame(np.array([[2,  3,  4,  5,  6],
                                       [8,  9, 10, 11, 12]], dtype='int64'),
                             index=Index([1, 7]),
                             columns=MultiIndex(levels=[[u('a'), u('b'), u('c')], [u('r'), u('s'), u('t'), u('u'), u('v')]],
                                                labels=[[0, 0, 1, 2, 2], [0, 1, 2, 3, 4]],
                                                names=[None, u('q')]))

        data = """,a,a,b,c,c
q,r,s,t,u,v
1,2,3,4,5,6
7,8,9,10,11,12"""

        result = self.read_csv(StringIO(data),header=[0,1],index_col=0)
        tm.assert_frame_equal(expected,result)

        # mi on columns and index (malformed)
        expected = DataFrame(np.array([[ 3,  4,  5,  6],
                                       [ 9, 10, 11, 12]], dtype='int64'),
                             index=MultiIndex(levels=[[1, 7], [2, 8]],
                                              labels=[[0, 1], [0, 1]]),
                             columns=MultiIndex(levels=[[u('a'), u('b'), u('c')], [u('s'), u('t'), u('u'), u('v')]],
                                                labels=[[0, 1, 2, 2], [0, 1, 2, 3]],
                                                names=[None, u('q')]))

        data = """,a,a,b,c,c
q,r,s,t,u,v
1,2,3,4,5,6
7,8,9,10,11,12"""

        result = self.read_csv(StringIO(data),header=[0,1],index_col=[0, 1])
        tm.assert_frame_equal(expected,result)

    def test_pass_names_with_index(self):
        lines = self.data1.split('\n')
        no_header = '\n'.join(lines[1:])

        # regular index
        names = ['index', 'A', 'B', 'C', 'D']
        df = self.read_csv(StringIO(no_header), index_col=0, names=names)
        expected = self.read_csv(StringIO(self.data1), index_col=0)
        tm.assert_frame_equal(df, expected)

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
        df = self.read_csv(StringIO(no_header), index_col=[0, 1],
                           names=names)
        expected = self.read_csv(StringIO(data), index_col=[0, 1])
        tm.assert_frame_equal(df, expected)

        df = self.read_csv(StringIO(data), index_col=['index1', 'index2'])
        tm.assert_frame_equal(df, expected)

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

        df = self.read_csv(StringIO(no_header), index_col=[0, 1],
                           header=None, names=names)
        expected = self.read_csv(StringIO(data), index_col=[0, 1])
        tm.assert_frame_equal(df, expected, check_names=False)

        # 2 implicit first cols
        df2 = self.read_csv(StringIO(data2))
        tm.assert_frame_equal(df2, df)

        # reverse order of index
        df = self.read_csv(StringIO(no_header), index_col=[1, 0], names=names,
                           header=None)
        expected = self.read_csv(StringIO(data), index_col=[1, 0])
        tm.assert_frame_equal(df, expected, check_names=False)

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
        df = self.read_csv(StringIO(data), index_col=[0, 1], parse_dates=True)
        self.assert_(isinstance(df.index.levels[0][0],
                     (datetime, np.datetime64, Timestamp)))

        # specify columns out of order!
        df2 = self.read_csv(StringIO(data), index_col=[1, 0], parse_dates=True)
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
        result = self.read_csv(StringIO(data), skip_footer=2)
        no_footer = '\n'.join(data.split('\n')[:-3])
        expected = self.read_csv(StringIO(no_footer))

        tm.assert_frame_equal(result, expected)

        result = self.read_csv(StringIO(data), nrows=3)
        tm.assert_frame_equal(result, expected)

        # skipfooter alias
        result = read_csv(StringIO(data), skipfooter=2)
        no_footer = '\n'.join(data.split('\n')[:-3])
        expected = read_csv(StringIO(no_footer))

        tm.assert_frame_equal(result, expected)

    def test_no_unnamed_index(self):
        data = """ id c0 c1 c2
0 1 0 a b
1 2 0 c d
2 2 2 e f
"""
        df = self.read_table(StringIO(data), sep=' ')
        self.assert_(df.index.name is None)

    def test_converters(self):
        data = """A,B,C,D
a,1,2,01/01/2009
b,3,4,01/02/2009
c,4,5,01/03/2009
"""
        from pandas.compat import parse_date

        result = self.read_csv(StringIO(data), converters={'D': parse_date})
        result2 = self.read_csv(StringIO(data), converters={3: parse_date})

        expected = self.read_csv(StringIO(data))
        expected['D'] = expected['D'].map(parse_date)

        tm.assert_isinstance(result['D'][0], (datetime, Timestamp))
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(result2, expected)

        # produce integer
        converter = lambda x: int(x.split('/')[2])
        result = self.read_csv(StringIO(data), converters={'D': converter})
        expected = self.read_csv(StringIO(data))
        expected['D'] = expected['D'].map(converter)
        tm.assert_frame_equal(result, expected)

    def test_converters_no_implicit_conv(self):
        # GH2184
        data = """000102,1.2,A\n001245,2,B"""
        f = lambda x: x.strip()
        converter = {0: f}
        df = self.read_csv(StringIO(data), header=None, converters=converter)
        self.assert_(df[0].dtype == object)

    def test_converters_euro_decimal_format(self):
        data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""
        f = lambda x: float(x.replace(",", "."))
        converter = {'Number1': f, 'Number2': f, 'Number3': f}
        df2 = self.read_csv(StringIO(data), sep=';', converters=converter)
        self.assert_(df2['Number1'].dtype == float)
        self.assert_(df2['Number2'].dtype == float)
        self.assert_(df2['Number3'].dtype == float)

    def test_converter_return_string_bug(self):
        # GH #583
        data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""
        f = lambda x: float(x.replace(",", "."))
        converter = {'Number1': f, 'Number2': f, 'Number3': f}
        df2 = self.read_csv(StringIO(data), sep=';', converters=converter)
        self.assert_(df2['Number1'].dtype == float)

    def test_read_table_buglet_4x_multiindex(self):
        text = """                      A       B       C       D        E
one two three   four
a   b   10.0032 5    -0.5109 -2.3358 -0.4645  0.05076  0.3640
a   q   20      4     0.4473  1.4152  0.2834  1.00661  0.1744
x   q   30      3    -0.6662 -0.5243 -0.3580  0.89145  2.5838"""

        # it works!
        df = self.read_table(StringIO(text), sep='\s+')
        self.assertEquals(df.index.names, ('one', 'two', 'three', 'four'))

    def test_read_csv_parse_simple_list(self):
        text = """foo
bar baz
qux foo
foo
bar"""
        df = read_csv(StringIO(text), header=None)
        expected = DataFrame({0: ['foo', 'bar baz', 'qux foo',
                                  'foo', 'bar']})
        tm.assert_frame_equal(df, expected)

    def test_parse_dates_custom_euroformat(self):
        text = """foo,bar,baz
31/01/2010,1,2
01/02/2010,1,NA
02/02/2010,1,2
"""
        parser = lambda d: parse_date(d, dayfirst=True)
        df = self.read_csv(StringIO(text),
                           names=['time', 'Q', 'NTU'], header=0,
                           index_col=0, parse_dates=True,
                           date_parser=parser, na_values=['NA'])

        exp_index = Index([datetime(2010, 1, 31), datetime(2010, 2, 1),
                           datetime(2010, 2, 2)], name='time')
        expected = DataFrame({'Q': [1, 1, 1], 'NTU': [2, np.nan, 2]},
                             index=exp_index, columns=['Q', 'NTU'])
        tm.assert_frame_equal(df, expected)

        parser = lambda d: parse_date(d, day_first=True)
        self.assertRaises(Exception, self.read_csv,
                          StringIO(text), skiprows=[0],
                          names=['time', 'Q', 'NTU'], index_col=0,
                          parse_dates=True, date_parser=parser,
                          na_values=['NA'])

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

    @tm.network
    def test_url(self):
        # HTTP(S)
        url = ('https://raw.github.com/pydata/pandas/master/'
                'pandas/io/tests/data/salary.table')
        url_table = self.read_table(url)
        dirpath = tm.get_data_path()
        localtable = os.path.join(dirpath, 'salary.table')
        local_table = self.read_table(localtable)
        tm.assert_frame_equal(url_table, local_table)
        # TODO: ftp testing

    @slow
    def test_file(self):

        # FILE
        if sys.version_info[:2] < (2, 6):
            raise nose.SkipTest("file:// not supported with Python < 2.6")
        dirpath = tm.get_data_path()
        localtable = os.path.join(dirpath, 'salary.table')
        local_table = self.read_table(localtable)

        try:
            url_table = self.read_table('file://localhost/' + localtable)
        except URLError:
            # fails on some systems
            raise nose.SkipTest("failing on %s" %
                                ' '.join(platform.uname()).strip())

        tm.assert_frame_equal(url_table, local_table)

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
        except AssertionError:  # hello Yaroslav
            arr = result.index.to_pydatetime()
            result = tools.to_datetime(arr, utc=True)[0]
            self.assert_(stamp.minute == result.minute)
            self.assert_(stamp.hour == result.hour)
            self.assert_(stamp.day == result.day)

    def test_multiple_date_cols_index(self):
        data = """\
ID,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir
KORD1,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD2,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD3,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD4,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD5,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD6,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""

        xp = self.read_csv(StringIO(data), parse_dates={'nominal': [1, 2]})
        df = self.read_csv(StringIO(data), parse_dates={'nominal': [1, 2]},
                           index_col='nominal')
        tm.assert_frame_equal(xp.set_index('nominal'), df)
        df2 = self.read_csv(StringIO(data), parse_dates={'nominal': [1, 2]},
                            index_col=0)
        tm.assert_frame_equal(df2, df)

        df3 = self.read_csv(StringIO(data), parse_dates=[[1, 2]], index_col=0)
        tm.assert_frame_equal(df3, df, check_names=False)

    def test_multiple_date_cols_chunked(self):
        df = self.read_csv(StringIO(self.ts_data), parse_dates={
                           'nominal': [1, 2]}, index_col='nominal')
        reader = self.read_csv(StringIO(self.ts_data), parse_dates={'nominal':
                               [1, 2]}, index_col='nominal', chunksize=2)

        chunks = list(reader)

        self.assert_('nominalTime' not in df)

        tm.assert_frame_equal(chunks[0], df[:2])
        tm.assert_frame_equal(chunks[1], df[2:4])
        tm.assert_frame_equal(chunks[2], df[4:])

    def test_multiple_date_col_named_components(self):
        xp = self.read_csv(StringIO(self.ts_data),
                           parse_dates={'nominal': [1, 2]},
                           index_col='nominal')
        colspec = {'nominal': ['date', 'nominalTime']}
        df = self.read_csv(StringIO(self.ts_data), parse_dates=colspec,
                           index_col='nominal')
        tm.assert_frame_equal(df, xp)

    def test_multiple_date_col_multiple_index(self):
        df = self.read_csv(StringIO(self.ts_data),
                           parse_dates={'nominal': [1, 2]},
                           index_col=['nominal', 'ID'])

        xp = self.read_csv(StringIO(self.ts_data),
                           parse_dates={'nominal': [1, 2]})

        tm.assert_frame_equal(xp.set_index(['nominal', 'ID']), df)

    def test_comment(self):
        data = """A,B,C
1,2.,4.#hello world
5.,NaN,10.0
"""
        expected = [[1., 2., 4.],
                    [5., np.nan, 10.]]
        df = self.read_csv(StringIO(data), comment='#')
        tm.assert_almost_equal(df.values, expected)

        df = self.read_table(StringIO(data), sep=',', comment='#',
                             na_values=['NaN'])
        tm.assert_almost_equal(df.values, expected)

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

    def test_nonexistent_path(self):
        # don't segfault pls #2428
        path = '%s.csv' % tm.rands(10)
        self.assertRaises(Exception, self.read_csv, path)

    def test_missing_trailing_delimiters(self):
        data = """A,B,C,D
1,2,3,4
1,3,3,
1,4,5"""
        result = self.read_csv(StringIO(data))
        self.assertTrue(result['D'].isnull()[1:].all())

    def test_skipinitialspace(self):
        s = ('"09-Apr-2012", "01:10:18.300", 2456026.548822908, 12849, '
             '1.00361,  1.12551, 330.65659, 0355626618.16711,  73.48821, '
             '314.11625,  1917.09447,   179.71425,  80.000, 240.000, -350,  '
             '70.06056, 344.98370, 1,   1, -0.689265, -0.692787,  '
             '0.212036,    14.7674,   41.605,   -9999.0,   -9999.0,   '
             '-9999.0,   -9999.0,   -9999.0,  -9999.0, 000, 012, 128')

        sfile = StringIO(s)
        # it's 33 columns
        result = self.read_csv(sfile, names=lrange(33), na_values=['-9999.0'],
                               header=None, skipinitialspace=True)
        self.assertTrue(pd.isnull(result.ix[0, 29]))

    def test_utf16_bom_skiprows(self):
        # #2298
        data = u("""skip this
skip this too
A\tB\tC
1\t2\t3
4\t5\t6""")

        data2 = u("""skip this
skip this too
A,B,C
1,2,3
4,5,6""")

        path = '__%s__.csv' % tm.rands(10)

        with tm.ensure_clean(path) as path:
            for sep, dat in [('\t', data), (',', data2)]:
                for enc in ['utf-16', 'utf-16le', 'utf-16be']:
                    bytes = dat.encode(enc)
                    with open(path, 'wb') as f:
                        f.write(bytes)

                    s = BytesIO(dat.encode('utf-8'))
                    if compat.PY3:
                        # somewhat False since the code never sees bytes
                        from io import TextIOWrapper
                        s = TextIOWrapper(s, encoding='utf-8')

                    result = self.read_csv(path, encoding=enc, skiprows=2,
                                           sep=sep)
                    expected = self.read_csv(s, encoding='utf-8', skiprows=2,
                                             sep=sep)

                    tm.assert_frame_equal(result, expected)

    def test_utf16_example(self):
        path = tm.get_data_path('utf16_ex.txt')

        # it works! and is the right length
        result = self.read_table(path, encoding='utf-16')
        self.assertEquals(len(result), 50)

        if not compat.PY3:
            buf = BytesIO(open(path, 'rb').read())
            result = self.read_table(buf, encoding='utf-16')
            self.assertEquals(len(result), 50)

    def test_converters_corner_with_nas(self):
        # skip aberration observed on Win64 Python 3.2.2
        if hash(np.int64(-1)) != -2:
            raise nose.SkipTest("skipping because of windows hash on Python"
                                " 3.2.2")

        csv = """id,score,days
1,2,12
2,2-5,
3,,14+
4,6-12,2"""

        def convert_days(x):
            x = x.strip()
            if not x:
                return np.nan

            is_plus = x.endswith('+')
            if is_plus:
                x = int(x[:-1]) + 1
            else:
                x = int(x)
            return x

        def convert_days_sentinel(x):
            x = x.strip()
            if not x:
                return np.nan

            is_plus = x.endswith('+')
            if is_plus:
                x = int(x[:-1]) + 1
            else:
                x = int(x)
            return x

        def convert_score(x):
            x = x.strip()
            if not x:
                return np.nan
            if x.find('-') > 0:
                valmin, valmax = lmap(int, x.split('-'))
                val = 0.5 * (valmin + valmax)
            else:
                val = float(x)

            return val

        fh = StringIO(csv)
        result = self.read_csv(fh, converters={'score': convert_score,
                                               'days': convert_days},
                               na_values=['', None])
        self.assert_(pd.isnull(result['days'][1]))

        fh = StringIO(csv)
        result2 = self.read_csv(fh, converters={'score': convert_score,
                                                'days': convert_days_sentinel},
                                na_values=['', None])
        tm.assert_frame_equal(result, result2)

    def test_unicode_encoding(self):
        pth = tm.get_data_path('unicode_series.csv')

        result = self.read_csv(pth, header=None, encoding='latin-1')
        result = result.set_index(0)

        got = result[1][1632]
        expected = u('\xc1 k\xf6ldum klaka (Cold Fever) (1994)')

        self.assertEquals(got, expected)

    def test_trailing_delimiters(self):
        # #2442. grumble grumble
        data = """A,B,C
1,2,3,
4,5,6,
7,8,9,"""
        result = self.read_csv(StringIO(data), index_col=False)

        expected = DataFrame({'A': [1, 4, 7], 'B': [2, 5, 8],
                              'C': [3, 6, 9]})

        tm.assert_frame_equal(result, expected)

    def test_escapechar(self):
        # http://stackoverflow.com/questions/13824840/feature-request-for-
        # pandas-read-csv
        data = '''SEARCH_TERM,ACTUAL_URL
"bra tv bord","http://www.ikea.com/se/sv/catalog/categories/departments/living_room/10475/?se%7cps%7cnonbranded%7cvardagsrum%7cgoogle%7ctv_bord"
"tv p\xc3\xa5 hjul","http://www.ikea.com/se/sv/catalog/categories/departments/living_room/10475/?se%7cps%7cnonbranded%7cvardagsrum%7cgoogle%7ctv_bord"
"SLAGBORD, \\"Bergslagen\\", IKEA:s 1700-tals serie","http://www.ikea.com/se/sv/catalog/categories/departments/living_room/10475/?se%7cps%7cnonbranded%7cvardagsrum%7cgoogle%7ctv_bord"'''

        result = self.read_csv(StringIO(data), escapechar='\\',
                               quotechar='"', encoding='utf-8')
        self.assertEqual(result['SEARCH_TERM'][2],
                         'SLAGBORD, "Bergslagen", IKEA:s 1700-tals serie')
        self.assertTrue(np.array_equal(result.columns,
                                       ['SEARCH_TERM', 'ACTUAL_URL']))

    def test_header_names_backward_compat(self):
        # #2539
        data = '1,2,3\n4,5,6'

        result = self.read_csv(StringIO(data), names=['a', 'b', 'c'])
        expected = self.read_csv(StringIO(data), names=['a', 'b', 'c'],
                                 header=None)
        tm.assert_frame_equal(result, expected)

        data2 = 'foo,bar,baz\n' + data
        result = self.read_csv(StringIO(data2), names=['a', 'b', 'c'],
                               header=0)
        tm.assert_frame_equal(result, expected)

    def test_integer_overflow_bug(self):
        # #2601
        data = "65248E10 11\n55555E55 22\n"

        result = self.read_csv(StringIO(data), header=None, sep=' ')
        self.assertTrue(result[0].dtype == np.float64)

    def test_int64_min_issues(self):
        # #2599
        data = 'A,B\n0,0\n0,'

        result = self.read_csv(StringIO(data))
        expected = DataFrame({'A': [0, 0], 'B': [0, np.nan]})

        tm.assert_frame_equal(result, expected)

    def test_parse_integers_above_fp_precision(self):
        data = """Numbers
17007000002000191
17007000002000191
17007000002000191
17007000002000191
17007000002000192
17007000002000192
17007000002000192
17007000002000192
17007000002000192
17007000002000194"""

        result = self.read_csv(StringIO(data))
        expected = DataFrame({'Numbers': [17007000002000191,
                                          17007000002000191,
                                          17007000002000191,
                                          17007000002000191,
                                          17007000002000192,
                                          17007000002000192,
                                          17007000002000192,
                                          17007000002000192,
                                          17007000002000192,
                                          17007000002000194]})

        self.assertTrue(np.array_equal(result['Numbers'], expected['Numbers']))

    def test_usecols_index_col_conflict(self):
        # Issue 4201  Test that index_col as integer reflects usecols
        data = """SecId,Time,Price,P2,P3
10000,2013-5-11,100,10,1
500,2013-5-12,101,11,1
"""
        expected = DataFrame({'Price': [100, 101]}, index=[datetime(2013, 5, 11), datetime(2013, 5, 12)])
        expected.index.name = 'Time'

        df = pd.read_csv(StringIO(data), usecols=['Time', 'Price'], parse_dates=True, index_col=0)
        tm.assert_frame_equal(expected, df)

        df = pd.read_csv(StringIO(data), usecols=['Time', 'Price'], parse_dates=True, index_col='Time')
        tm.assert_frame_equal(expected, df)

        df = pd.read_csv(StringIO(data), usecols=[1, 2], parse_dates=True, index_col='Time')
        tm.assert_frame_equal(expected, df)

        df = pd.read_csv(StringIO(data), usecols=[1, 2], parse_dates=True, index_col=0)
        tm.assert_frame_equal(expected, df)

        expected = DataFrame({'P3': [1, 1], 'Price': (100, 101), 'P2': (10, 11)})
        expected = expected.set_index(['Price', 'P2'])
        df = pd.read_csv(StringIO(data), usecols=['Price', 'P2', 'P3'], parse_dates=True, index_col=['Price', 'P2'])
        tm.assert_frame_equal(expected, df)

    def test_chunks_have_consistent_numerical_type(self):
        integers = [str(i) for i in range(499999)]
        data = "a\n" + "\n".join(integers + ["1.0", "2.0"] + integers)

        with tm.assert_produces_warning(False):
            df = self.read_csv(StringIO(data))
        self.assertTrue(type(df.a[0]) is np.float64)  # Assert that types were coerced.
        self.assertEqual(df.a.dtype, np.float)

    def test_warn_if_chunks_have_mismatched_type(self):
        # See test in TestCParserLowMemory.
        integers = [str(i) for i in range(499999)]
        data = "a\n" + "\n".join(integers + ['a', 'b'] + integers)

        with tm.assert_produces_warning(False):
            df = self.read_csv(StringIO(data))
        self.assertEqual(df.a.dtype, np.object)

    def test_usecols(self):
        data = """\
a,b,c
1,2,3
4,5,6
7,8,9
10,11,12"""

        result = self.read_csv(StringIO(data), usecols=(1, 2))
        result2 = self.read_csv(StringIO(data), usecols=('b', 'c'))
        exp = self.read_csv(StringIO(data))

        self.assertEquals(len(result.columns), 2)
        self.assertTrue((result['b'] == exp['b']).all())
        self.assertTrue((result['c'] == exp['c']).all())

        tm.assert_frame_equal(result, result2)

        result = self.read_csv(StringIO(data), usecols=[1, 2], header=0,
                               names=['foo', 'bar'])
        expected = self.read_csv(StringIO(data), usecols=[1, 2])
        expected.columns = ['foo', 'bar']
        tm.assert_frame_equal(result, expected)

        data = """\
1,2,3
4,5,6
7,8,9
10,11,12"""
        result = self.read_csv(StringIO(data), names=['b', 'c'],
                               header=None, usecols=[1, 2])

        expected = self.read_csv(StringIO(data), names=['a', 'b', 'c'],
                                 header=None)
        expected = expected[['b', 'c']]
        tm.assert_frame_equal(result, expected)

        result2 = self.read_csv(StringIO(data), names=['a', 'b', 'c'],
                                header=None, usecols=['b', 'c'])
        tm.assert_frame_equal(result2, result)


        # 5766
        result = self.read_csv(StringIO(data), names=['a', 'b'],
                               header=None, usecols=[0, 1])

        expected = self.read_csv(StringIO(data), names=['a', 'b', 'c'],
                                 header=None)
        expected = expected[['a', 'b']]
        tm.assert_frame_equal(result, expected)

        # length conflict, passed names and usecols disagree
        self.assertRaises(ValueError, self.read_csv, StringIO(data),
                          names=['a', 'b'], usecols=[1], header=None)

    def test_integer_overflow_bug(self):
        # #2601
        data = "65248E10 11\n55555E55 22\n"

        result = self.read_csv(StringIO(data), header=None, sep=' ')
        self.assertTrue(result[0].dtype == np.float64)

        result = self.read_csv(StringIO(data), header=None, sep='\s+')
        self.assertTrue(result[0].dtype == np.float64)

    def test_catch_too_many_names(self):
        # Issue 5156
        data = """\
1,2,3
4,,6
7,8,9
10,11,12\n"""
        tm.assertRaises(Exception, read_csv, StringIO(data), header=0, names=['a', 'b', 'c', 'd'])


class TestPythonParser(ParserTests, tm.TestCase):
    def test_negative_skipfooter_raises(self):
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

        with tm.assertRaisesRegexp(ValueError,
                                   'skip footer cannot be negative'):
            df = self.read_csv(StringIO(text), skipfooter=-1)

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = 'python'
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = 'python'
        return read_table(*args, **kwds)

    def test_sniff_delimiter(self):
        text = """index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
        data = self.read_csv(StringIO(text), index_col=0, sep=None)
        self.assert_(data.index.equals(Index(['foo', 'bar', 'baz'])))

        data2 = self.read_csv(StringIO(text), index_col=0, delimiter='|')
        tm.assert_frame_equal(data, data2)

        text = """ignore this
ignore this too
index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
"""
        data3 = self.read_csv(StringIO(text), index_col=0,
                              sep=None, skiprows=2)
        tm.assert_frame_equal(data, data3)

        text = u("""ignore this
ignore this too
index|A|B|C
foo|1|2|3
bar|4|5|6
baz|7|8|9
""").encode('utf-8')

        s = BytesIO(text)
        if compat.PY3:
            # somewhat False since the code never sees bytes
            from io import TextIOWrapper
            s = TextIOWrapper(s, encoding='utf-8')

        data4 = self.read_csv(s, index_col=0, sep=None, skiprows=2,
                              encoding='utf-8')
        tm.assert_frame_equal(data, data4)

    def test_regex_separator(self):
        data = """   A   B   C   D
a   1   2   3   4
b   1   2   3   4
c   1   2   3   4
"""
        df = self.read_table(StringIO(data), sep='\s+')
        expected = self.read_csv(StringIO(re.sub('[ ]+', ',', data)),
                                 index_col=0)
        self.assert_(expected.index.name is None)
        tm.assert_frame_equal(df, expected)

    def test_1000_fwf(self):
        data = """
 1 2,334.0    5
10   13     10.
"""
        expected = [[1, 2334., 5],
                    [10, 13, 10]]
        df = read_fwf(StringIO(data), colspecs=[(0, 3), (3, 11), (12, 16)],
                      thousands=',')
        tm.assert_almost_equal(df.values, expected)

    def test_1000_sep_with_decimal(self):
        data = """A|B|C
1|2,334.01|5
10|13|10.
"""

        expected = DataFrame({
            'A': [1, 10],
            'B': [2334.01, 13],
            'C': [5, 10.]
        })

        df = self.read_csv(StringIO(data), sep='|', thousands=',')
        tm.assert_frame_equal(df, expected)

        df = self.read_table(StringIO(data), sep='|', thousands=',')
        tm.assert_frame_equal(df, expected)

    def test_comment_fwf(self):
        data = """
  1   2.   4  #hello world
  5  NaN  10.0
"""
        expected = [[1, 2., 4],
                    [5, np.nan, 10.]]
        df = read_fwf(StringIO(data), colspecs=[(0, 3), (4, 9), (9, 25)],
                      comment='#')
        tm.assert_almost_equal(df.values, expected)

    def test_fwf(self):
        data_expected = """\
2011,58,360.242940,149.910199,11950.7
2011,59,444.953632,166.985655,11788.4
2011,60,364.136849,183.628767,11806.2
2011,61,413.836124,184.375703,11916.8
2011,62,502.953953,173.237159,12468.3
"""
        expected = self.read_csv(StringIO(data_expected), header=None)

        data1 = """\
201158    360.242940   149.910199   11950.7
201159    444.953632   166.985655   11788.4
201160    364.136849   183.628767   11806.2
201161    413.836124   184.375703   11916.8
201162    502.953953   173.237159   12468.3
"""
        colspecs = [(0, 4), (4, 8), (8, 20), (21, 33), (34, 43)]
        df = read_fwf(StringIO(data1), colspecs=colspecs, header=None)
        tm.assert_frame_equal(df, expected)

        data2 = """\
2011 58   360.242940   149.910199   11950.7
2011 59   444.953632   166.985655   11788.4
2011 60   364.136849   183.628767   11806.2
2011 61   413.836124   184.375703   11916.8
2011 62   502.953953   173.237159   12468.3
"""
        df = read_fwf(StringIO(data2), widths=[5, 5, 13, 13, 7], header=None)
        tm.assert_frame_equal(df, expected)

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
        df = read_fwf(
            StringIO(data3), colspecs=colspecs, delimiter='~', header=None)
        tm.assert_frame_equal(df, expected)

        with tm.assertRaisesRegexp(ValueError, "must specify only one of"):
            read_fwf(StringIO(data3), colspecs=colspecs, widths=[6, 10, 10, 7])

        with tm.assertRaisesRegexp(ValueError, "Must specify either"):
            read_fwf(StringIO(data3), colspecs=None, widths=None)

    def test_fwf_colspecs_is_list_or_tuple(self):
        with tm.assertRaisesRegexp(TypeError,
                                   'column specifications must be a list or '
                                   'tuple.+'):
            pd.io.parsers.FixedWidthReader(StringIO(self.data1),
                                           {'a': 1}, ',', '#')

    def test_fwf_colspecs_is_list_or_tuple_of_two_element_tuples(self):
        with tm.assertRaisesRegexp(TypeError,
                                   'Each column specification must be.+'):
            read_fwf(StringIO(self.data1), [('a', 1)])

    def test_fwf_regression(self):
        # GH 3594
        #### turns out 'T060' is parsable as a datetime slice!

        tzlist = [1,10,20,30,60,80,100]
        ntz = len(tzlist)
        tcolspecs = [16]+[8]*ntz
        tcolnames = ['SST'] + ["T%03d" % z for z in tzlist[1:]]
        data = """  2009164202000   9.5403  9.4105  8.6571  7.8372  6.0612  5.8843  5.5192
  2009164203000   9.5435  9.2010  8.6167  7.8176  6.0804  5.8728  5.4869
  2009164204000   9.5873  9.1326  8.4694  7.5889  6.0422  5.8526  5.4657
  2009164205000   9.5810  9.0896  8.4009  7.4652  6.0322  5.8189  5.4379
  2009164210000   9.6034  9.0897  8.3822  7.4905  6.0908  5.7904  5.4039
"""

        df = read_fwf(StringIO(data),
                      index_col=0,
                      header=None,
                      names=tcolnames,
                      widths=tcolspecs,
                      parse_dates=True,
                      date_parser=lambda s: datetime.strptime(s,'%Y%j%H%M%S'))

        for c in df.columns:
            res = df.loc[:,c]
            self.assert_(len(res))

    def test_fwf_compression(self):
        try:
            import gzip
            import bz2
        except ImportError:
            raise nose.SkipTest("Need gzip and bz2 to run this test")

        data = """1111111111
        2222222222
        3333333333""".strip()
        widths = [5, 5]
        names = ['one', 'two']
        expected = read_fwf(StringIO(data), widths=widths, names=names)
        if compat.PY3:
            data = bytes(data, encoding='utf-8')
        comps = [('gzip', gzip.GzipFile), ('bz2', bz2.BZ2File)]
        for comp_name, compresser in comps:
            with tm.ensure_clean() as path:
                tmp = compresser(path, mode='wb')
                tmp.write(data)
                tmp.close()
                result = read_fwf(path, widths=widths, names=names,
                                  compression=comp_name)
                tm.assert_frame_equal(result, expected)

    def test_BytesIO_input(self):
        if not compat.PY3:
            raise nose.SkipTest("Bytes-related test - only needs to work on Python 3")
        result = pd.read_fwf(BytesIO("שלום\nשלום".encode('utf8')), widths=[2,2], encoding='utf8')
        expected = pd.DataFrame([["של", "ום"]], columns=["של", "ום"])
        tm.assert_frame_equal(result, expected)
        data = BytesIO("שלום::1234\n562::123".encode('cp1255'))
        result = pd.read_table(data, sep="::", engine='python',
                               encoding='cp1255')
        expected = pd.DataFrame([[562, 123]], columns=["שלום","1234"])
        tm.assert_frame_equal(result, expected)

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
            df = self.read_csv(StringIO(text), verbose=True)
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
            df = self.read_csv(StringIO(text), verbose=True, index_col=0)
            self.assert_(buf.getvalue() == 'Filled 1 NA values in column a\n')
        finally:
            sys.stdout = sys.__stdout__

    def test_iteration_open_handle(self):
        if PY3:
            raise nose.SkipTest("won't work in Python 3 {0}".format(sys.version_info))

        with tm.ensure_clean() as path:
            with open(path, 'wb') as f:
                f.write('AAA\nBBB\nCCC\nDDD\nEEE\nFFF\nGGG')

            with open(path, 'rb') as f:
                for line in f:
                    if 'CCC' in line:
                        break

                try:
                    read_table(f, squeeze=True, header=None, engine='c')
                except Exception:
                    pass
                else:
                    raise ValueError('this should not happen')

                result = read_table(f, squeeze=True, header=None,
                                    engine='python')

                expected = Series(['DDD', 'EEE', 'FFF', 'GGG'])
                tm.assert_series_equal(result, expected)


class TestFwfColspaceSniffing(tm.TestCase):
    def test_full_file(self):
        # File with all values
        test = '''index                             A    B    C
2000-01-03T00:00:00  0.980268513777    3  foo
2000-01-04T00:00:00  1.04791624281    -4  bar
2000-01-05T00:00:00  0.498580885705   73  baz
2000-01-06T00:00:00  1.12020151869     1  foo
2000-01-07T00:00:00  0.487094399463    0  bar
2000-01-10T00:00:00  0.836648671666    2  baz
2000-01-11T00:00:00  0.157160753327   34  foo'''
        colspecs = ((0, 19), (21, 35), (38, 40), (42, 45))
        expected = read_fwf(StringIO(test), colspecs=colspecs)
        tm.assert_frame_equal(expected, read_fwf(StringIO(test)))

    def test_full_file_with_missing(self):
        # File with missing values
        test = '''index                             A    B    C
2000-01-03T00:00:00  0.980268513777    3  foo
2000-01-04T00:00:00  1.04791624281    -4  bar
                     0.498580885705   73  baz
2000-01-06T00:00:00  1.12020151869     1  foo
2000-01-07T00:00:00                    0  bar
2000-01-10T00:00:00  0.836648671666    2  baz
                                      34'''
        colspecs = ((0, 19), (21, 35), (38, 40), (42, 45))
        expected = read_fwf(StringIO(test), colspecs=colspecs)
        tm.assert_frame_equal(expected, read_fwf(StringIO(test)))

    def test_full_file_with_spaces(self):
        # File with spaces in columns
        test = '''
Account                 Name  Balance     CreditLimit   AccountCreated
101     Keanu Reeves          9315.45     10000.00           1/17/1998
312     Gerard Butler         90.00       1000.00             8/6/2003
868     Jennifer Love Hewitt  0           17000.00           5/25/1985
761     Jada Pinkett-Smith    49654.87    100000.00          12/5/2006
317     Bill Murray           789.65      5000.00             2/5/2007
'''.strip('\r\n')
        colspecs = ((0, 7), (8, 28), (30, 38), (42, 53), (56, 70))
        expected = read_fwf(StringIO(test), colspecs=colspecs)
        tm.assert_frame_equal(expected, read_fwf(StringIO(test)))

    def test_full_file_with_spaces_and_missing(self):
        # File with spaces and missing values in columsn
        test = '''
Account               Name    Balance     CreditLimit   AccountCreated
101                           10000.00                       1/17/1998
312     Gerard Butler         90.00       1000.00             8/6/2003
868                                                          5/25/1985
761     Jada Pinkett-Smith    49654.87    100000.00          12/5/2006
317     Bill Murray           789.65
'''.strip('\r\n')
        colspecs = ((0, 7), (8, 28), (30, 38), (42, 53), (56, 70))
        expected = read_fwf(StringIO(test), colspecs=colspecs)
        tm.assert_frame_equal(expected, read_fwf(StringIO(test)))

    def test_messed_up_data(self):
        # Completely messed up file
        test = '''
   Account          Name             Balance     Credit Limit   Account Created
       101                           10000.00                       1/17/1998
       312     Gerard Butler         90.00       1000.00

       761     Jada Pinkett-Smith    49654.87    100000.00          12/5/2006
  317          Bill Murray           789.65
'''.strip('\r\n')
        colspecs = ((2, 10), (15, 33), (37, 45), (49, 61), (64, 79))
        expected = read_fwf(StringIO(test), colspecs=colspecs)
        tm.assert_frame_equal(expected, read_fwf(StringIO(test)))

    def test_multiple_delimiters(self):
        test = r'''
col1~~~~~col2  col3++++++++++++++++++col4
~~22.....11.0+++foo~~~~~~~~~~Keanu Reeves
  33+++122.33\\\bar.........Gerard Butler
++44~~~~12.01   baz~~Jennifer Love Hewitt
~~55       11+++foo++++Jada Pinkett-Smith
..66++++++.03~~~bar           Bill Murray
'''.strip('\r\n')
        colspecs = ((0, 4), (7, 13), (15, 19), (21, 41))
        expected = read_fwf(StringIO(test), colspecs=colspecs,
                            delimiter=' +~.\\')
        tm.assert_frame_equal(expected, read_fwf(StringIO(test),
                                                 delimiter=' +~.\\'))

    def test_variable_width_unicode(self):
        if not compat.PY3:
            raise nose.SkipTest('Bytes-related test - only needs to work on Python 3')
        test = '''
שלום שלום
ום   שלל
של   ום
'''.strip('\r\n')
        expected = pd.read_fwf(BytesIO(test.encode('utf8')),
                               colspecs=[(0, 4), (5, 9)], header=None, encoding='utf8')
        tm.assert_frame_equal(expected, read_fwf(BytesIO(test.encode('utf8')),
                                                 header=None, encoding='utf8'))


class TestCParserHighMemory(ParserTests, tm.TestCase):

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = 'c'
        kwds['low_memory'] = False
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = 'c'
        kwds['low_memory'] = False
        return read_table(*args, **kwds)

    def test_compact_ints(self):
        data = ('0,1,0,0\n'
                '1,1,0,0\n'
                '0,1,0,1')

        result = read_csv(StringIO(data), delimiter=',', header=None,
                          compact_ints=True, as_recarray=True)
        ex_dtype = np.dtype([(str(i), 'i1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

        result = read_csv(StringIO(data), delimiter=',', header=None,
                          as_recarray=True, compact_ints=True,
                          use_unsigned=True)
        ex_dtype = np.dtype([(str(i), 'u1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

    def test_parse_dates_empty_string(self):
        # #2263
        s = StringIO("Date, test\n2012-01-01, 1\n,2")
        result = pd.read_csv(s, parse_dates=["Date"], na_filter=False)
        self.assertTrue(result['Date'].isnull()[1])

    def test_usecols(self):
        raise nose.SkipTest("Usecols is not supported in C High Memory engine.")


class TestCParserLowMemory(ParserTests, tm.TestCase):

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = 'c'
        kwds['low_memory'] = True
        kwds['buffer_lines'] = 2
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = 'c'
        kwds['low_memory'] = True
        kwds['buffer_lines'] = 2
        return read_table(*args, **kwds)

    def test_compact_ints(self):
        data = ('0,1,0,0\n'
                '1,1,0,0\n'
                '0,1,0,1')

        result = read_csv(StringIO(data), delimiter=',', header=None,
                          compact_ints=True, as_recarray=True)
        ex_dtype = np.dtype([(str(i), 'i1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

        result = read_csv(StringIO(data), delimiter=',', header=None,
                          as_recarray=True, compact_ints=True,
                          use_unsigned=True)
        ex_dtype = np.dtype([(str(i), 'u1') for i in range(4)])
        self.assertEqual(result.dtype, ex_dtype)

    def test_pass_dtype(self):
        data = """\
one,two
1,2.5
2,3.5
3,4.5
4,5.5"""

        result = self.read_csv(StringIO(data), dtype={'one': 'u1', 1: 'S1'},
                               as_recarray=True)
        self.assert_(result['one'].dtype == 'u1')
        self.assert_(result['two'].dtype == 'S1')

    def test_usecols_dtypes(self):
        data = """\
1,2,3
4,5,6
7,8,9
10,11,12"""
        result = self.read_csv(StringIO(data), usecols=(0, 1, 2),
                               names=('a', 'b', 'c'),
                               header=None,
                               converters={'a': str},
                               dtype={'b': int, 'c': float},
                              )
        result2 = self.read_csv(StringIO(data), usecols=(0, 2),
                               names=('a', 'b', 'c'),
                               header=None,
                               converters={'a': str},
                               dtype={'b': int, 'c': float},
                              )
        self.assertTrue((result.dtypes == [object, np.int, np.float]).all())
        self.assertTrue((result2.dtypes == [object, np.float]).all())

    def test_usecols_implicit_index_col(self):
        # #2654
        data = 'a,b,c\n4,apple,bat,5.7\n8,orange,cow,10'

        result = self.read_csv(StringIO(data), usecols=['a', 'b'])
        expected = DataFrame({'a': ['apple', 'orange'],
                              'b': ['bat', 'cow']}, index=[4, 8])

        tm.assert_frame_equal(result, expected)

    def test_usecols_with_whitespace(self):
        data = 'a  b  c\n4  apple  bat  5.7\n8  orange  cow  10'

        result = self.read_csv(StringIO(data), delim_whitespace=True,
                               usecols=('a', 'b'))
        expected = DataFrame({'a': ['apple', 'orange'],
                              'b': ['bat', 'cow']}, index=[4, 8])

        tm.assert_frame_equal(result, expected)

    def test_usecols_regex_sep(self):
        # #2733
        data = 'a  b  c\n4  apple  bat  5.7\n8  orange  cow  10'

        df = self.read_csv(StringIO(data), sep='\s+', usecols=('a', 'b'))

        expected = DataFrame({'a': ['apple', 'orange'],
                              'b': ['bat', 'cow']}, index=[4, 8])
        tm.assert_frame_equal(df, expected)

    def test_pure_python_failover(self):
        data = "a,b,c\n1,2,3#ignore this!\n4,5,6#ignorethistoo"

        result = self.read_csv(StringIO(data), comment='#')
        expected = DataFrame({'a': [1, 4], 'b': [2, 5], 'c': [3, 6]})
        tm.assert_frame_equal(result, expected)

    def test_decompression(self):
        try:
            import gzip
            import bz2
        except ImportError:
            raise nose.SkipTest('need gzip and bz2 to run')

        data = open(self.csv1, 'rb').read()
        expected = self.read_csv(self.csv1)

        with tm.ensure_clean() as path:
            tmp = gzip.GzipFile(path, mode='wb')
            tmp.write(data)
            tmp.close()

            result = self.read_csv(path, compression='gzip')
            tm.assert_frame_equal(result, expected)

            result = self.read_csv(open(path, 'rb'), compression='gzip')
            tm.assert_frame_equal(result, expected)

        with tm.ensure_clean() as path:
            tmp = bz2.BZ2File(path, mode='wb')
            tmp.write(data)
            tmp.close()

            result = self.read_csv(path, compression='bz2')
            tm.assert_frame_equal(result, expected)

            # result = self.read_csv(open(path, 'rb'), compression='bz2')
            # tm.assert_frame_equal(result, expected)

            self.assertRaises(ValueError, self.read_csv,
                              path, compression='bz3')

    def test_decompression_regex_sep(self):
        try:
            import gzip
            import bz2
        except ImportError:
            raise nose.SkipTest('need gzip and bz2 to run')

        data = open(self.csv1, 'rb').read()
        data = data.replace(b',', b'::')
        expected = self.read_csv(self.csv1)

        with tm.ensure_clean() as path:
            tmp = gzip.GzipFile(path, mode='wb')
            tmp.write(data)
            tmp.close()

            result = self.read_csv(path, sep='::', compression='gzip')
            tm.assert_frame_equal(result, expected)

        with tm.ensure_clean() as path:
            tmp = bz2.BZ2File(path, mode='wb')
            tmp.write(data)
            tmp.close()

            result = self.read_csv(path, sep='::', compression='bz2')
            tm.assert_frame_equal(result, expected)

            self.assertRaises(ValueError, self.read_csv,
                              path, compression='bz3')

    def test_memory_map(self):
        # it works!
        result = self.read_csv(self.csv1, memory_map=True)

    def test_disable_bool_parsing(self):
        # #2090

        data = """A,B,C
Yes,No,Yes
No,Yes,Yes
Yes,,Yes
No,No,No"""

        result = read_csv(StringIO(data), dtype=object)
        self.assertTrue((result.dtypes == object).all())

        result = read_csv(StringIO(data), dtype=object, na_filter=False)
        self.assertEquals(result['B'][2], '')

    def test_int64_overflow(self):
        data = """ID
00013007854817840016671868
00013007854817840016749251
00013007854817840016754630
00013007854817840016781876
00013007854817840017028824
00013007854817840017963235
00013007854817840018860166"""

        result = read_csv(StringIO(data))
        self.assertTrue(result['ID'].dtype == object)

        self.assertRaises(OverflowError, read_csv, StringIO(data),
                          dtype='i8')

    def test_euro_decimal_format(self):
        data = """Id;Number1;Number2;Text1;Text2;Number3
1;1521,1541;187101,9543;ABC;poi;4,738797819
2;121,12;14897,76;DEF;uyt;0,377320872
3;878,158;108013,434;GHI;rez;2,735694704"""

        df2 = self.read_csv(StringIO(data), sep=';', decimal=',')
        self.assert_(df2['Number1'].dtype == float)
        self.assert_(df2['Number2'].dtype == float)
        self.assert_(df2['Number3'].dtype == float)

    def test_custom_lineterminator(self):
        data = 'a,b,c~1,2,3~4,5,6'

        result = self.read_csv(StringIO(data), lineterminator='~')
        expected = self.read_csv(StringIO(data.replace('~', '\n')))

        tm.assert_frame_equal(result, expected)

        data2 = data.replace('~', '~~')
        result = self.assertRaises(ValueError, read_csv, StringIO(data2),
                                   lineterminator='~~')

    def test_raise_on_passed_int_dtype_with_nas(self):
        # #2631
        data = """YEAR, DOY, a
2001,106380451,10
2001,,11
2001,106380451,67"""
        self.assertRaises(Exception, read_csv, StringIO(data), sep=",",
                          skipinitialspace=True,
                          dtype={'DOY': np.int64})

    def test_na_trailing_columns(self):
        data = """Date,Currenncy,Symbol,Type,Units,UnitPrice,Cost,Tax
2012-03-14,USD,AAPL,BUY,1000
2012-05-12,USD,SBUX,SELL,500"""

        result = self.read_csv(StringIO(data))
        self.assertEquals(result['Date'][1], '2012-05-12')
        self.assertTrue(result['UnitPrice'].isnull().all())

    def test_parse_ragged_csv(self):
        data = """1,2,3
1,2,3,4
1,2,3,4,5
1,2
1,2,3,4"""

        nice_data = """1,2,3,,
1,2,3,4,
1,2,3,4,5
1,2,,,
1,2,3,4,"""
        result = self.read_csv(StringIO(data), header=None,
                               names=['a', 'b', 'c', 'd', 'e'])

        expected = self.read_csv(StringIO(nice_data), header=None,
                                 names=['a', 'b', 'c', 'd', 'e'])

        tm.assert_frame_equal(result, expected)

        # too many columns, cause segfault if not careful
        data = "1,2\n3,4,5"

        result = self.read_csv(StringIO(data), header=None,
                               names=lrange(50))
        expected = self.read_csv(StringIO(data), header=None,
                                 names=lrange(3)).reindex(columns=lrange(50))

        tm.assert_frame_equal(result, expected)

    def test_tokenize_CR_with_quoting(self):
        # #3453, this doesn't work with Python parser for some reason

        data = ' a,b,c\r"a,b","e,d","f,f"'

        result = self.read_csv(StringIO(data), header=None)
        expected = self.read_csv(StringIO(data.replace('\r', '\n')),
                                 header=None)
        tm.assert_frame_equal(result, expected)

        result = self.read_csv(StringIO(data))
        expected = self.read_csv(StringIO(data.replace('\r', '\n')))
        tm.assert_frame_equal(result, expected)

    def test_raise_on_no_columns(self):
        # single newline
        data = "\n"
        self.assertRaises(ValueError, self.read_csv, StringIO(data))

        # test with more than a single newline
        data = "\n\n\n"
        self.assertRaises(ValueError, self.read_csv, StringIO(data))

    def test_warn_if_chunks_have_mismatched_type(self):
        # Issue #3866 If chunks are different types and can't
        # be coerced using numerical types, then issue warning.
        integers = [str(i) for i in range(499999)]
        data = "a\n" + "\n".join(integers + ['a', 'b'] + integers)

        with tm.assert_produces_warning(DtypeWarning):
            df = self.read_csv(StringIO(data))
        self.assertEqual(df.a.dtype, np.object)

    def test_invalid_c_parser_opts_with_not_c_parser(self):
        from pandas.io.parsers import _c_parser_defaults as c_defaults

        data = """1,2,3,,
1,2,3,4,
1,2,3,4,5
1,2,,,
1,2,3,4,"""

        engines = 'python', 'python-fwf'
        for default in c_defaults:
            for engine in engines:
                kwargs = {default: object()}
                with tm.assertRaisesRegexp(ValueError,
                                           'The %r option is not supported '
                                           'with the %r engine' % (default,
                                                                   engine)):
                    read_csv(StringIO(data), engine=engine, **kwargs)

class TestParseSQL(tm.TestCase):

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
        arr = np.array([u('1.5'), None, u('3'), u('4.2')],
                       dtype=object)
        result = lib.convert_sql_column(arr)
        expected = np.array([u('1.5'), np.nan, u('3'), u('4.2')],
                            dtype=object)
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
        arr = np.array([long(1), long(2), long(3), long(4)], dtype='O')
        result = lib.convert_sql_column(arr)
        expected = np.array([1, 2, 3, 4], dtype='i8')
        assert_same_values_and_dtype(result, expected)

        arr = np.array([long(1), long(2), long(3), None, long(4)], dtype='O')
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
    tm.assert_equal(res.dtype, exp.dtype)
    tm.assert_almost_equal(res, exp)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
