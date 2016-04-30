# -*- coding: utf-8 -*-

"""
Tests date parsing functionality for all of the
parsers defined in parsers.py
"""

from distutils.version import LooseVersion
from datetime import datetime

import nose
import numpy as np
import pandas.lib as lib
from pandas.lib import Timestamp

import pandas as pd
import pandas.io.parsers as parsers
import pandas.tseries.tools as tools
import pandas.util.testing as tm

from pandas import DataFrame, Series, Index, DatetimeIndex
from pandas import compat
from pandas.compat import(parse_date, StringIO,
                          lrange, lmap)
from pandas.tseries.index import date_range


class ParseDatesTests(object):
    def test_separator_date_conflict(self):
        # Regression test for gh-4678: make sure thousands separator and
        # date parsing do not conflict.
        data = '06-02-2013;13:00;1-000.215'
        expected = DataFrame(
            [[datetime(2013, 6, 2, 13, 0, 0), 1000.215]],
            columns=['Date', 2]
        )

        df = self.read_csv(StringIO(data), sep=';', thousands='-',
                           parse_dates={'Date': [0, 1]}, header=None)
        tm.assert_frame_equal(df, expected)

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
        self.assertIn('nominal', df)
        self.assertIn('actual', df)
        self.assertNotIn('X1', df)
        self.assertNotIn('X2', df)
        self.assertNotIn('X3', df)

        d = datetime(1999, 1, 27, 19, 0)
        self.assertEqual(df.ix[0, 'nominal'], d)

        df = self.read_csv(StringIO(data), header=None,
                           date_parser=func,
                           parse_dates={'nominal': [1, 2],
                                        'actual': [1, 3]},
                           keep_date_col=True)
        self.assertIn('nominal', df)
        self.assertIn('actual', df)

        self.assertIn(1, df)
        self.assertIn(2, df)
        self.assertIn(3, df)

        data = """\
KORD,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000
"""
        df = self.read_csv(StringIO(data), header=None,
                           prefix='X', parse_dates=[[1, 2], [1, 3]])

        self.assertIn('X1_X2', df)
        self.assertIn('X1_X3', df)
        self.assertNotIn('X1', df)
        self.assertNotIn('X2', df)
        self.assertNotIn('X3', df)

        d = datetime(1999, 1, 27, 19, 0)
        self.assertEqual(df.ix[0, 'X1_X2'], d)

        df = self.read_csv(StringIO(data), header=None,
                           parse_dates=[[1, 2], [1, 3]], keep_date_col=True)

        self.assertIn('1_2', df)
        self.assertIn('1_3', df)
        self.assertIn(1, df)
        self.assertIn(2, df)
        self.assertIn(3, df)

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
        self.assertEqual(df.index[0], d)

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
        self.assertIn('nominal', df)

    def test_multiple_date_col_timestamp_parse(self):
        data = """05/31/2012,15:30:00.029,1306.25,1,E,0,,1306.25
05/31/2012,15:30:00.029,1306.25,8,E,0,,1306.25"""
        result = self.read_csv(StringIO(data), sep=',', header=None,
                               parse_dates=[[0, 1]], date_parser=Timestamp)

        ex_val = Timestamp('05/31/2012 15:30:00.029')
        self.assertEqual(result['0_1'][0], ex_val)

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
        self.assertNotIsInstance(df.nominal[0], compat.string_types)

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
KORD6,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000"""  # noqa

        self.assertRaises(ValueError, self.read_csv, StringIO(data),
                          parse_dates=[[1, 2]])

    def test_date_parser_int_bug(self):
        # See gh-3071
        log_file = StringIO(
            'posix_timestamp,elapsed,sys,user,queries,query_time,rows,'
            'accountid,userid,contactid,level,silo,method\n'
            '1343103150,0.062353,0,4,6,0.01690,3,'
            '12345,1,-1,3,invoice_InvoiceResource,search\n'
        )

        def f(posix_string):
            return datetime.utcfromtimestamp(int(posix_string))

        # it works!
        self.read_csv(log_file, index_col=0, parse_dates=[0], date_parser=f)

    def test_nat_parse(self):
        # See gh-3062
        df = DataFrame(dict({
            'A': np.asarray(lrange(10), dtype='float64'),
            'B': pd.Timestamp('20010101')}))
        df.iloc[3:6, :] = np.nan

        with tm.ensure_clean('__nat_parse_.csv') as path:
            df.to_csv(path)
            result = self.read_csv(path, index_col=0, parse_dates=['B'])
            tm.assert_frame_equal(result, df)

            expected = Series(dict(A='float64', B='datetime64[ns]'))
            tm.assert_series_equal(expected, result.dtypes)

            # test with NaT for the nan_rep
            # we don't have a method to specif the Datetime na_rep (it defaults
            # to '')
            df.to_csv(path)
            result = self.read_csv(path, index_col=0, parse_dates=['B'])
            tm.assert_frame_equal(result, df)

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
        self.assertIsInstance(
            df.index[0], (datetime, np.datetime64, Timestamp))
        tm.assert_frame_equal(df, expected)

    def test_parse_dates_string(self):
        data = """date,A,B,C
20090101,a,1,2
20090102,b,3,4
20090103,c,4,5
"""
        rs = self.read_csv(
            StringIO(data), index_col='date', parse_dates=['date'])
        idx = date_range('1/1/2009', periods=3)
        idx.name = 'date'
        xp = DataFrame({'A': ['a', 'b', 'c'],
                        'B': [1, 3, 4],
                        'C': [2, 4, 5]}, idx)
        tm.assert_frame_equal(rs, xp)

    def test_yy_format_with_yearfirst(self):
        data = """date,time,B,C
090131,0010,1,2
090228,1020,3,4
090331,0830,5,6
"""

        # See gh-217
        import dateutil
        if dateutil.__version__ >= LooseVersion('2.5.0'):
            raise nose.SkipTest("testing yearfirst=True not-support"
                                "on datetutil < 2.5.0 this works but"
                                "is wrong")

        rs = self.read_csv(StringIO(data), index_col=0,
                           parse_dates=[['date', 'time']])
        idx = DatetimeIndex([datetime(2009, 1, 31, 0, 10, 0),
                             datetime(2009, 2, 28, 10, 20, 0),
                             datetime(2009, 3, 31, 8, 30, 0)],
                            dtype=object, name='date_time')
        xp = DataFrame({'B': [1, 3, 5], 'C': [2, 4, 6]}, idx)
        tm.assert_frame_equal(rs, xp)

        rs = self.read_csv(StringIO(data), index_col=0,
                           parse_dates=[[0, 1]])
        idx = DatetimeIndex([datetime(2009, 1, 31, 0, 10, 0),
                             datetime(2009, 2, 28, 10, 20, 0),
                             datetime(2009, 3, 31, 8, 30, 0)],
                            dtype=object, name='date_time')
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
        tm.assertIsInstance(expected['aux_date'][0], datetime)

        df = self.read_csv(StringIO(data), sep=";", index_col=lrange(4),
                           parse_dates=[0, 5], dayfirst=True)
        tm.assert_frame_equal(df, expected)

        df = self.read_csv(StringIO(data), sep=";", index_col=lrange(4),
                           parse_dates=['date', 'aux_date'], dayfirst=True)
        tm.assert_frame_equal(df, expected)

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
        self.assertIsInstance(df.index.levels[0][0],
                              (datetime, np.datetime64, Timestamp))

        # specify columns out of order!
        df2 = self.read_csv(StringIO(data), index_col=[1, 0], parse_dates=True)
        self.assertIsInstance(df2.index.levels[1][0],
                              (datetime, np.datetime64, Timestamp))

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
        self.assertRaises(TypeError, self.read_csv,
                          StringIO(text), skiprows=[0],
                          names=['time', 'Q', 'NTU'], index_col=0,
                          parse_dates=True, date_parser=parser,
                          na_values=['NA'])

    def test_parse_tz_aware(self):
        # See gh-1693
        import pytz
        data = StringIO("Date,x\n2012-06-13T01:39:00Z,0.5")

        # it works
        result = self.read_csv(data, index_col=0, parse_dates=True)
        stamp = result.index[0]
        self.assertEqual(stamp.minute, 39)
        try:
            self.assertIs(result.index.tz, pytz.utc)
        except AssertionError:  # hello Yaroslav
            arr = result.index.to_pydatetime()
            result = tools.to_datetime(arr, utc=True)[0]
            self.assertEqual(stamp.minute, result.minute)
            self.assertEqual(stamp.hour, result.hour)
            self.assertEqual(stamp.day, result.day)

    def test_multiple_date_cols_index(self):
        data = """
ID,date,NominalTime,ActualTime,TDew,TAir,Windspeed,Precip,WindDir
KORD1,19990127, 19:00:00, 18:56:00, 0.8100, 2.8100, 7.2000, 0.0000, 280.0000
KORD2,19990127, 20:00:00, 19:56:00, 0.0100, 2.2100, 7.2000, 0.0000, 260.0000
KORD3,19990127, 21:00:00, 20:56:00, -0.5900, 2.2100, 5.7000, 0.0000, 280.0000
KORD4,19990127, 21:00:00, 21:18:00, -0.9900, 2.0100, 3.6000, 0.0000, 270.0000
KORD5,19990127, 22:00:00, 21:56:00, -0.5900, 1.7100, 5.1000, 0.0000, 290.0000
KORD6,19990127, 23:00:00, 22:56:00, -0.5900, 1.7100, 4.6000, 0.0000, 280.0000
"""

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
        reader = self.read_csv(StringIO(self.ts_data),
                               parse_dates={'nominal': [1, 2]},
                               index_col='nominal', chunksize=2)

        chunks = list(reader)

        self.assertNotIn('nominalTime', df)

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

    def test_read_with_parse_dates_scalar_non_bool(self):
        # See gh-5636
        errmsg = ("Only booleans, lists, and "
                  "dictionaries are accepted "
                  "for the 'parse_dates' parameter")
        data = """A,B,C
        1,2,2003-11-1"""

        tm.assertRaisesRegexp(TypeError, errmsg, self.read_csv,
                              StringIO(data), parse_dates="C")
        tm.assertRaisesRegexp(TypeError, errmsg, self.read_csv,
                              StringIO(data), parse_dates="C",
                              index_col="C")

    def test_read_with_parse_dates_invalid_type(self):
        errmsg = ("Only booleans, lists, and "
                  "dictionaries are accepted "
                  "for the 'parse_dates' parameter")
        data = """A,B,C
        1,2,2003-11-1"""

        tm.assertRaisesRegexp(TypeError, errmsg, self.read_csv,
                              StringIO(data), parse_dates=(1,))
        tm.assertRaisesRegexp(TypeError, errmsg, self.read_csv,
                              StringIO(data), parse_dates=np.array([4, 5]))
        tm.assertRaisesRegexp(TypeError, errmsg, self.read_csv,
                              StringIO(data), parse_dates=set([1, 3, 3]))
