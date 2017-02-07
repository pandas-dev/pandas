""" test to_datetime """

import nose

import sys
import calendar
import locale
from datetime import datetime

import numpy as np
from pandas.types.common import is_datetime64_ns_dtype
from pandas import (isnull, to_datetime, Timestamp, Series, DataFrame,
                    Index, DatetimeIndex, NaT, date_range, bdate_range)
from pandas import tslib
from pandas.compat import lmap
import pandas as pd
from pandas.tseries import tools
from pandas.util import testing as tm
from pandas.util.testing import assert_series_equal


class TimeConversionFormats(tm.TestCase):

    def test_to_datetime_format(self):
        values = ['1/1/2000', '1/2/2000', '1/3/2000']

        results1 = [Timestamp('20000101'), Timestamp('20000201'),
                    Timestamp('20000301')]
        results2 = [Timestamp('20000101'), Timestamp('20000102'),
                    Timestamp('20000103')]
        for vals, expecteds in [(values, (Index(results1), Index(results2))),
                                (Series(values),
                                 (Series(results1), Series(results2))),
                                (values[0], (results1[0], results2[0])),
                                (values[1], (results1[1], results2[1])),
                                (values[2], (results1[2], results2[2]))]:

            for i, fmt in enumerate(['%d/%m/%Y', '%m/%d/%Y']):
                result = to_datetime(vals, format=fmt)
                expected = expecteds[i]

                if isinstance(expected, Series):
                    assert_series_equal(result, Series(expected))
                elif isinstance(expected, Timestamp):
                    self.assertEqual(result, expected)
                else:
                    tm.assert_index_equal(result, expected)

    def test_to_datetime_format_YYYYMMDD(self):
        s = Series([19801222, 19801222] + [19810105] * 5)
        expected = Series([Timestamp(x) for x in s.apply(str)])

        result = to_datetime(s, format='%Y%m%d')
        assert_series_equal(result, expected)

        result = to_datetime(s.apply(str), format='%Y%m%d')
        assert_series_equal(result, expected)

        # with NaT
        expected = Series([Timestamp("19801222"), Timestamp("19801222")] +
                          [Timestamp("19810105")] * 5)
        expected[2] = np.nan
        s[2] = np.nan

        result = to_datetime(s, format='%Y%m%d')
        assert_series_equal(result, expected)

        # string with NaT
        s = s.apply(str)
        s[2] = 'nat'
        result = to_datetime(s, format='%Y%m%d')
        assert_series_equal(result, expected)

        # coercion
        # GH 7930
        s = Series([20121231, 20141231, 99991231])
        result = pd.to_datetime(s, format='%Y%m%d', errors='ignore')
        expected = Series([datetime(2012, 12, 31),
                           datetime(2014, 12, 31), datetime(9999, 12, 31)],
                          dtype=object)
        self.assert_series_equal(result, expected)

        result = pd.to_datetime(s, format='%Y%m%d', errors='coerce')
        expected = Series(['20121231', '20141231', 'NaT'], dtype='M8[ns]')
        assert_series_equal(result, expected)

    # GH 10178
    def test_to_datetime_format_integer(self):
        s = Series([2000, 2001, 2002])
        expected = Series([Timestamp(x) for x in s.apply(str)])

        result = to_datetime(s, format='%Y')
        assert_series_equal(result, expected)

        s = Series([200001, 200105, 200206])
        expected = Series([Timestamp(x[:4] + '-' + x[4:]) for x in s.apply(str)
                           ])

        result = to_datetime(s, format='%Y%m')
        assert_series_equal(result, expected)

    def test_to_datetime_format_microsecond(self):

        # these are locale dependent
        lang, _ = locale.getlocale()
        month_abbr = calendar.month_abbr[4]
        val = '01-{}-2011 00:00:01.978'.format(month_abbr)

        format = '%d-%b-%Y %H:%M:%S.%f'
        result = to_datetime(val, format=format)
        exp = datetime.strptime(val, format)
        self.assertEqual(result, exp)

    def test_to_datetime_format_time(self):
        data = [
            ['01/10/2010 15:20', '%m/%d/%Y %H:%M',
             Timestamp('2010-01-10 15:20')],
            ['01/10/2010 05:43', '%m/%d/%Y %I:%M',
             Timestamp('2010-01-10 05:43')],
            ['01/10/2010 13:56:01', '%m/%d/%Y %H:%M:%S',
             Timestamp('2010-01-10 13:56:01')]  # ,
            # ['01/10/2010 08:14 PM', '%m/%d/%Y %I:%M %p',
            #  Timestamp('2010-01-10 20:14')],
            # ['01/10/2010 07:40 AM', '%m/%d/%Y %I:%M %p',
            #  Timestamp('2010-01-10 07:40')],
            # ['01/10/2010 09:12:56 AM', '%m/%d/%Y %I:%M:%S %p',
            #  Timestamp('2010-01-10 09:12:56')]
        ]
        for s, format, dt in data:
            self.assertEqual(to_datetime(s, format=format), dt)

    def test_to_datetime_with_non_exact(self):
        # GH 10834
        tm._skip_if_has_locale()

        # 8904
        # exact kw
        if sys.version_info < (2, 7):
            raise nose.SkipTest('on python version < 2.7')

        s = Series(['19MAY11', 'foobar19MAY11', '19MAY11:00:00:00',
                    '19MAY11 00:00:00Z'])
        result = to_datetime(s, format='%d%b%y', exact=False)
        expected = to_datetime(s.str.extract(r'(\d+\w+\d+)', expand=False),
                               format='%d%b%y')
        assert_series_equal(result, expected)

    def test_parse_nanoseconds_with_formula(self):

        # GH8989
        # trunctaing the nanoseconds when a format was provided
        for v in ["2012-01-01 09:00:00.000000001",
                  "2012-01-01 09:00:00.000001",
                  "2012-01-01 09:00:00.001",
                  "2012-01-01 09:00:00.001000",
                  "2012-01-01 09:00:00.001000000", ]:
            expected = pd.to_datetime(v)
            result = pd.to_datetime(v, format="%Y-%m-%d %H:%M:%S.%f")
            self.assertEqual(result, expected)

    def test_to_datetime_format_weeks(self):
        data = [
            ['2009324', '%Y%W%w', Timestamp('2009-08-13')],
            ['2013020', '%Y%U%w', Timestamp('2013-01-13')]
        ]
        for s, format, dt in data:
            self.assertEqual(to_datetime(s, format=format), dt)


class TestToDatetime(tm.TestCase):

    def test_to_datetime_dt64s(self):
        in_bound_dts = [
            np.datetime64('2000-01-01'),
            np.datetime64('2000-01-02'),
        ]

        for dt in in_bound_dts:
            self.assertEqual(pd.to_datetime(dt), Timestamp(dt))

        oob_dts = [np.datetime64('1000-01-01'), np.datetime64('5000-01-02'), ]

        for dt in oob_dts:
            self.assertRaises(ValueError, pd.to_datetime, dt, errors='raise')
            self.assertRaises(ValueError, Timestamp, dt)
            self.assertIs(pd.to_datetime(dt, errors='coerce'), NaT)

    def test_to_datetime_array_of_dt64s(self):
        dts = [np.datetime64('2000-01-01'), np.datetime64('2000-01-02'), ]

        # Assuming all datetimes are in bounds, to_datetime() returns
        # an array that is equal to Timestamp() parsing
        self.assert_numpy_array_equal(
            pd.to_datetime(dts, box=False),
            np.array([Timestamp(x).asm8 for x in dts])
        )

        # A list of datetimes where the last one is out of bounds
        dts_with_oob = dts + [np.datetime64('9999-01-01')]

        self.assertRaises(ValueError, pd.to_datetime, dts_with_oob,
                          errors='raise')

        self.assert_numpy_array_equal(
            pd.to_datetime(dts_with_oob, box=False, errors='coerce'),
            np.array(
                [
                    Timestamp(dts_with_oob[0]).asm8,
                    Timestamp(dts_with_oob[1]).asm8,
                    tslib.iNaT,
                ],
                dtype='M8'
            )
        )

        # With errors='ignore', out of bounds datetime64s
        # are converted to their .item(), which depending on the version of
        # numpy is either a python datetime.datetime or datetime.date
        self.assert_numpy_array_equal(
            pd.to_datetime(dts_with_oob, box=False, errors='ignore'),
            np.array(
                [dt.item() for dt in dts_with_oob],
                dtype='O'
            )
        )

    def test_to_datetime_tz(self):

        # xref 8260
        # uniform returns a DatetimeIndex
        arr = [pd.Timestamp('2013-01-01 13:00:00-0800', tz='US/Pacific'),
               pd.Timestamp('2013-01-02 14:00:00-0800', tz='US/Pacific')]
        result = pd.to_datetime(arr)
        expected = DatetimeIndex(
            ['2013-01-01 13:00:00', '2013-01-02 14:00:00'], tz='US/Pacific')
        tm.assert_index_equal(result, expected)

        # mixed tzs will raise
        arr = [pd.Timestamp('2013-01-01 13:00:00', tz='US/Pacific'),
               pd.Timestamp('2013-01-02 14:00:00', tz='US/Eastern')]
        self.assertRaises(ValueError, lambda: pd.to_datetime(arr))

    def test_to_datetime_tz_pytz(self):

        # xref 8260
        tm._skip_if_no_pytz()
        import pytz

        us_eastern = pytz.timezone('US/Eastern')
        arr = np.array([us_eastern.localize(datetime(year=2000, month=1, day=1,
                                                     hour=3, minute=0)),
                        us_eastern.localize(datetime(year=2000, month=6, day=1,
                                                     hour=3, minute=0))],
                       dtype=object)
        result = pd.to_datetime(arr, utc=True)
        expected = DatetimeIndex(['2000-01-01 08:00:00+00:00',
                                  '2000-06-01 07:00:00+00:00'],
                                 dtype='datetime64[ns, UTC]', freq=None)
        tm.assert_index_equal(result, expected)

    def test_to_datetime_utc_is_true(self):
        # See gh-11934
        start = pd.Timestamp('2014-01-01', tz='utc')
        end = pd.Timestamp('2014-01-03', tz='utc')
        date_range = pd.bdate_range(start, end)

        result = pd.to_datetime(date_range, utc=True)
        expected = pd.DatetimeIndex(data=date_range)
        tm.assert_index_equal(result, expected)

    def test_to_datetime_tz_psycopg2(self):

        # xref 8260
        try:
            import psycopg2
        except ImportError:
            raise nose.SkipTest("no psycopg2 installed")

        # misc cases
        tz1 = psycopg2.tz.FixedOffsetTimezone(offset=-300, name=None)
        tz2 = psycopg2.tz.FixedOffsetTimezone(offset=-240, name=None)
        arr = np.array([datetime(2000, 1, 1, 3, 0, tzinfo=tz1),
                        datetime(2000, 6, 1, 3, 0, tzinfo=tz2)],
                       dtype=object)

        result = pd.to_datetime(arr, errors='coerce', utc=True)
        expected = DatetimeIndex(['2000-01-01 08:00:00+00:00',
                                  '2000-06-01 07:00:00+00:00'],
                                 dtype='datetime64[ns, UTC]', freq=None)
        tm.assert_index_equal(result, expected)

        # dtype coercion
        i = pd.DatetimeIndex([
            '2000-01-01 08:00:00+00:00'
        ], tz=psycopg2.tz.FixedOffsetTimezone(offset=-300, name=None))
        self.assertTrue(is_datetime64_ns_dtype(i))

        # tz coerceion
        result = pd.to_datetime(i, errors='coerce')
        tm.assert_index_equal(result, i)

        result = pd.to_datetime(i, errors='coerce', utc=True)
        expected = pd.DatetimeIndex(['2000-01-01 13:00:00'],
                                    dtype='datetime64[ns, UTC]')
        tm.assert_index_equal(result, expected)

    def test_datetime_bool(self):
        # GH13176
        with self.assertRaises(TypeError):
            to_datetime(False)
        self.assertTrue(to_datetime(False, errors="coerce") is NaT)
        self.assertEqual(to_datetime(False, errors="ignore"), False)
        with self.assertRaises(TypeError):
            to_datetime(True)
        self.assertTrue(to_datetime(True, errors="coerce") is NaT)
        self.assertEqual(to_datetime(True, errors="ignore"), True)
        with self.assertRaises(TypeError):
            to_datetime([False, datetime.today()])
        with self.assertRaises(TypeError):
            to_datetime(['20130101', True])
        tm.assert_index_equal(to_datetime([0, False, NaT, 0.0],
                                          errors="coerce"),
                              DatetimeIndex([to_datetime(0), NaT,
                                             NaT, to_datetime(0)]))

    def test_datetime_invalid_datatype(self):
        # GH13176

        with self.assertRaises(TypeError):
            pd.to_datetime(bool)
        with self.assertRaises(TypeError):
            pd.to_datetime(pd.to_datetime)


class ToDatetimeUnit(tm.TestCase):

    def test_unit(self):
        # GH 11758
        # test proper behavior with erros

        with self.assertRaises(ValueError):
            to_datetime([1], unit='D', format='%Y%m%d')

        values = [11111111, 1, 1.0, tslib.iNaT, NaT, np.nan,
                  'NaT', '']
        result = to_datetime(values, unit='D', errors='ignore')
        expected = Index([11111111, Timestamp('1970-01-02'),
                          Timestamp('1970-01-02'), NaT,
                          NaT, NaT, NaT, NaT],
                         dtype=object)
        tm.assert_index_equal(result, expected)

        result = to_datetime(values, unit='D', errors='coerce')
        expected = DatetimeIndex(['NaT', '1970-01-02', '1970-01-02',
                                  'NaT', 'NaT', 'NaT', 'NaT', 'NaT'])
        tm.assert_index_equal(result, expected)

        with self.assertRaises(tslib.OutOfBoundsDatetime):
            to_datetime(values, unit='D', errors='raise')

        values = [1420043460000, tslib.iNaT, NaT, np.nan, 'NaT']

        result = to_datetime(values, errors='ignore', unit='s')
        expected = Index([1420043460000, NaT, NaT,
                          NaT, NaT], dtype=object)
        tm.assert_index_equal(result, expected)

        result = to_datetime(values, errors='coerce', unit='s')
        expected = DatetimeIndex(['NaT', 'NaT', 'NaT', 'NaT', 'NaT'])
        tm.assert_index_equal(result, expected)

        with self.assertRaises(tslib.OutOfBoundsDatetime):
            to_datetime(values, errors='raise', unit='s')

        # if we have a string, then we raise a ValueError
        # and NOT an OutOfBoundsDatetime
        for val in ['foo', Timestamp('20130101')]:
            try:
                to_datetime(val, errors='raise', unit='s')
            except tslib.OutOfBoundsDatetime:
                raise AssertionError("incorrect exception raised")
            except ValueError:
                pass

    def test_unit_consistency(self):

        # consistency of conversions
        expected = Timestamp('1970-05-09 14:25:11')
        result = pd.to_datetime(11111111, unit='s', errors='raise')
        self.assertEqual(result, expected)
        self.assertIsInstance(result, Timestamp)

        result = pd.to_datetime(11111111, unit='s', errors='coerce')
        self.assertEqual(result, expected)
        self.assertIsInstance(result, Timestamp)

        result = pd.to_datetime(11111111, unit='s', errors='ignore')
        self.assertEqual(result, expected)
        self.assertIsInstance(result, Timestamp)

    def test_unit_with_numeric(self):

        # GH 13180
        # coercions from floats/ints are ok
        expected = DatetimeIndex(['2015-06-19 05:33:20',
                                  '2015-05-27 22:33:20'])
        arr1 = [1.434692e+18, 1.432766e+18]
        arr2 = np.array(arr1).astype('int64')
        for errors in ['ignore', 'raise', 'coerce']:
            result = pd.to_datetime(arr1, errors=errors)
            tm.assert_index_equal(result, expected)

            result = pd.to_datetime(arr2, errors=errors)
            tm.assert_index_equal(result, expected)

        # but we want to make sure that we are coercing
        # if we have ints/strings
        expected = DatetimeIndex(['NaT',
                                  '2015-06-19 05:33:20',
                                  '2015-05-27 22:33:20'])
        arr = ['foo', 1.434692e+18, 1.432766e+18]
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

        expected = DatetimeIndex(['2015-06-19 05:33:20',
                                  '2015-05-27 22:33:20',
                                  'NaT',
                                  'NaT'])
        arr = [1.434692e+18, 1.432766e+18, 'foo', 'NaT']
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

    def test_unit_mixed(self):

        # mixed integers/datetimes
        expected = DatetimeIndex(['2013-01-01', 'NaT', 'NaT'])
        arr = [pd.Timestamp('20130101'), 1.434692e+18, 1.432766e+18]
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

        with self.assertRaises(ValueError):
            pd.to_datetime(arr, errors='raise')

        expected = DatetimeIndex(['NaT',
                                  'NaT',
                                  '2013-01-01'])
        arr = [1.434692e+18, 1.432766e+18, pd.Timestamp('20130101')]
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

        with self.assertRaises(ValueError):
            pd.to_datetime(arr, errors='raise')

    def test_dataframe(self):

        df = DataFrame({'year': [2015, 2016],
                        'month': [2, 3],
                        'day': [4, 5],
                        'hour': [6, 7],
                        'minute': [58, 59],
                        'second': [10, 11],
                        'ms': [1, 1],
                        'us': [2, 2],
                        'ns': [3, 3]})

        result = to_datetime({'year': df['year'],
                              'month': df['month'],
                              'day': df['day']})
        expected = Series([Timestamp('20150204 00:00:00'),
                           Timestamp('20160305 00:0:00')])
        assert_series_equal(result, expected)

        # dict-like
        result = to_datetime(df[['year', 'month', 'day']].to_dict())
        assert_series_equal(result, expected)

        # dict but with constructable
        df2 = df[['year', 'month', 'day']].to_dict()
        df2['month'] = 2
        result = to_datetime(df2)
        expected2 = Series([Timestamp('20150204 00:00:00'),
                            Timestamp('20160205 00:0:00')])
        assert_series_equal(result, expected2)

        # unit mappings
        units = [{'year': 'years',
                  'month': 'months',
                  'day': 'days',
                  'hour': 'hours',
                  'minute': 'minutes',
                  'second': 'seconds'},
                 {'year': 'year',
                  'month': 'month',
                  'day': 'day',
                  'hour': 'hour',
                  'minute': 'minute',
                  'second': 'second'},
                 ]

        for d in units:
            result = to_datetime(df[list(d.keys())].rename(columns=d))
            expected = Series([Timestamp('20150204 06:58:10'),
                               Timestamp('20160305 07:59:11')])
            assert_series_equal(result, expected)

        d = {'year': 'year',
             'month': 'month',
             'day': 'day',
             'hour': 'hour',
             'minute': 'minute',
             'second': 'second',
             'ms': 'ms',
             'us': 'us',
             'ns': 'ns'}

        result = to_datetime(df.rename(columns=d))
        expected = Series([Timestamp('20150204 06:58:10.001002003'),
                           Timestamp('20160305 07:59:11.001002003')])
        assert_series_equal(result, expected)

        # coerce back to int
        result = to_datetime(df.astype(str))
        assert_series_equal(result, expected)

        # passing coerce
        df2 = DataFrame({'year': [2015, 2016],
                         'month': [2, 20],
                         'day': [4, 5]})
        with self.assertRaises(ValueError):
            to_datetime(df2)
        result = to_datetime(df2, errors='coerce')
        expected = Series([Timestamp('20150204 00:00:00'),
                           NaT])
        assert_series_equal(result, expected)

        # extra columns
        with self.assertRaises(ValueError):
            df2 = df.copy()
            df2['foo'] = 1
            to_datetime(df2)

        # not enough
        for c in [['year'],
                  ['year', 'month'],
                  ['year', 'month', 'second'],
                  ['month', 'day'],
                  ['year', 'day', 'second']]:
            with self.assertRaises(ValueError):
                to_datetime(df[c])

        # duplicates
        df2 = DataFrame({'year': [2015, 2016],
                         'month': [2, 20],
                         'day': [4, 5]})
        df2.columns = ['year', 'year', 'day']
        with self.assertRaises(ValueError):
            to_datetime(df2)

        df2 = DataFrame({'year': [2015, 2016],
                         'month': [2, 20],
                         'day': [4, 5],
                         'hour': [4, 5]})
        df2.columns = ['year', 'month', 'day', 'day']
        with self.assertRaises(ValueError):
            to_datetime(df2)

    def test_dataframe_dtypes(self):
        # #13451
        df = DataFrame({'year': [2015, 2016],
                        'month': [2, 3],
                        'day': [4, 5]})

        # int16
        result = to_datetime(df.astype('int16'))
        expected = Series([Timestamp('20150204 00:00:00'),
                           Timestamp('20160305 00:00:00')])
        assert_series_equal(result, expected)

        # mixed dtypes
        df['month'] = df['month'].astype('int8')
        df['day'] = df['day'].astype('int8')
        result = to_datetime(df)
        expected = Series([Timestamp('20150204 00:00:00'),
                           Timestamp('20160305 00:00:00')])
        assert_series_equal(result, expected)

        # float
        df = DataFrame({'year': [2000, 2001],
                        'month': [1.5, 1],
                        'day': [1, 1]})
        with self.assertRaises(ValueError):
            to_datetime(df)


class ToDatetimeMisc(tm.TestCase):

    def test_index_to_datetime(self):
        idx = Index(['1/1/2000', '1/2/2000', '1/3/2000'])

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = idx.to_datetime()
            expected = DatetimeIndex(pd.to_datetime(idx.values))
            tm.assert_index_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            today = datetime.today()
            idx = Index([today], dtype=object)
            result = idx.to_datetime()
            expected = DatetimeIndex([today])
            tm.assert_index_equal(result, expected)

    def test_to_datetime_iso8601(self):
        result = to_datetime(["2012-01-01 00:00:00"])
        exp = Timestamp("2012-01-01 00:00:00")
        self.assertEqual(result[0], exp)

        result = to_datetime(['20121001'])  # bad iso 8601
        exp = Timestamp('2012-10-01')
        self.assertEqual(result[0], exp)

    def test_to_datetime_default(self):
        rs = to_datetime('2001')
        xp = datetime(2001, 1, 1)
        self.assertTrue(rs, xp)

        # dayfirst is essentially broken

        # to_datetime('01-13-2012', dayfirst=True)
        # self.assertRaises(ValueError, to_datetime('01-13-2012',
        #                   dayfirst=True))

    def test_to_datetime_on_datetime64_series(self):
        # #2699
        s = Series(date_range('1/1/2000', periods=10))

        result = to_datetime(s)
        self.assertEqual(result[0], s[0])

    def test_to_datetime_with_space_in_series(self):
        # GH 6428
        s = Series(['10/18/2006', '10/18/2008', ' '])
        tm.assertRaises(ValueError, lambda: to_datetime(s, errors='raise'))
        result_coerce = to_datetime(s, errors='coerce')
        expected_coerce = Series([datetime(2006, 10, 18),
                                  datetime(2008, 10, 18),
                                  NaT])
        tm.assert_series_equal(result_coerce, expected_coerce)
        result_ignore = to_datetime(s, errors='ignore')
        tm.assert_series_equal(result_ignore, s)

    def test_to_datetime_with_apply(self):
        # this is only locale tested with US/None locales
        tm._skip_if_has_locale()

        # GH 5195
        # with a format and coerce a single item to_datetime fails
        td = Series(['May 04', 'Jun 02', 'Dec 11'], index=[1, 2, 3])
        expected = pd.to_datetime(td, format='%b %y')
        result = td.apply(pd.to_datetime, format='%b %y')
        assert_series_equal(result, expected)

        td = pd.Series(['May 04', 'Jun 02', ''], index=[1, 2, 3])
        self.assertRaises(ValueError,
                          lambda: pd.to_datetime(td, format='%b %y',
                                                 errors='raise'))
        self.assertRaises(ValueError,
                          lambda: td.apply(pd.to_datetime, format='%b %y',
                                           errors='raise'))
        expected = pd.to_datetime(td, format='%b %y', errors='coerce')

        result = td.apply(
            lambda x: pd.to_datetime(x, format='%b %y', errors='coerce'))
        assert_series_equal(result, expected)

    def test_to_datetime_types(self):

        # empty string
        result = to_datetime('')
        self.assertIs(result, NaT)

        result = to_datetime(['', ''])
        self.assertTrue(isnull(result).all())

        # ints
        result = Timestamp(0)
        expected = to_datetime(0)
        self.assertEqual(result, expected)

        # GH 3888 (strings)
        expected = to_datetime(['2012'])[0]
        result = to_datetime('2012')
        self.assertEqual(result, expected)

        # array = ['2012','20120101','20120101 12:01:01']
        array = ['20120101', '20120101 12:01:01']
        expected = list(to_datetime(array))
        result = lmap(Timestamp, array)
        tm.assert_almost_equal(result, expected)

        # currently fails ###
        # result = Timestamp('2012')
        # expected = to_datetime('2012')
        # self.assertEqual(result, expected)

    def test_to_datetime_unprocessable_input(self):
        # GH 4928
        self.assert_numpy_array_equal(
            to_datetime([1, '1'], errors='ignore'),
            np.array([1, '1'], dtype='O')
        )
        self.assertRaises(TypeError, to_datetime, [1, '1'], errors='raise')

    def test_to_datetime_other_datetime64_units(self):
        # 5/25/2012
        scalar = np.int64(1337904000000000).view('M8[us]')
        as_obj = scalar.astype('O')

        index = DatetimeIndex([scalar])
        self.assertEqual(index[0], scalar.astype('O'))

        value = Timestamp(scalar)
        self.assertEqual(value, as_obj)

    def test_to_datetime_list_of_integers(self):
        rng = date_range('1/1/2000', periods=20)
        rng = DatetimeIndex(rng.values)

        ints = list(rng.asi8)

        result = DatetimeIndex(ints)

        tm.assert_index_equal(rng, result)

    def test_to_datetime_freq(self):
        xp = bdate_range('2000-1-1', periods=10, tz='UTC')
        rs = xp.to_datetime()
        self.assertEqual(xp.freq, rs.freq)
        self.assertEqual(xp.tzinfo, rs.tzinfo)

    def test_string_na_nat_conversion(self):
        # GH #999, #858

        from pandas.compat import parse_date

        strings = np.array(['1/1/2000', '1/2/2000', np.nan,
                            '1/4/2000, 12:34:56'], dtype=object)

        expected = np.empty(4, dtype='M8[ns]')
        for i, val in enumerate(strings):
            if isnull(val):
                expected[i] = tslib.iNaT
            else:
                expected[i] = parse_date(val)

        result = tslib.array_to_datetime(strings)
        tm.assert_almost_equal(result, expected)

        result2 = to_datetime(strings)
        tm.assertIsInstance(result2, DatetimeIndex)
        tm.assert_numpy_array_equal(result, result2.values)

        malformed = np.array(['1/100/2000', np.nan], dtype=object)

        # GH 10636, default is now 'raise'
        self.assertRaises(ValueError,
                          lambda: to_datetime(malformed, errors='raise'))

        result = to_datetime(malformed, errors='ignore')
        tm.assert_numpy_array_equal(result, malformed)

        self.assertRaises(ValueError, to_datetime, malformed, errors='raise')

        idx = ['a', 'b', 'c', 'd', 'e']
        series = Series(['1/1/2000', np.nan, '1/3/2000', np.nan,
                         '1/5/2000'], index=idx, name='foo')
        dseries = Series([to_datetime('1/1/2000'), np.nan,
                          to_datetime('1/3/2000'), np.nan,
                          to_datetime('1/5/2000')], index=idx, name='foo')

        result = to_datetime(series)
        dresult = to_datetime(dseries)

        expected = Series(np.empty(5, dtype='M8[ns]'), index=idx)
        for i in range(5):
            x = series[i]
            if isnull(x):
                expected[i] = tslib.iNaT
            else:
                expected[i] = to_datetime(x)

        assert_series_equal(result, expected, check_names=False)
        self.assertEqual(result.name, 'foo')

        assert_series_equal(dresult, expected, check_names=False)
        self.assertEqual(dresult.name, 'foo')

    def test_dti_constructor_numpy_timeunits(self):
        # GH 9114
        base = pd.to_datetime(['2000-01-01T00:00', '2000-01-02T00:00', 'NaT'])

        for dtype in ['datetime64[h]', 'datetime64[m]', 'datetime64[s]',
                      'datetime64[ms]', 'datetime64[us]', 'datetime64[ns]']:
            values = base.values.astype(dtype)

            tm.assert_index_equal(DatetimeIndex(values), base)
            tm.assert_index_equal(to_datetime(values), base)

    def test_dayfirst(self):
        # GH 5917
        arr = ['10/02/2014', '11/02/2014', '12/02/2014']
        expected = DatetimeIndex([datetime(2014, 2, 10), datetime(2014, 2, 11),
                                  datetime(2014, 2, 12)])
        idx1 = DatetimeIndex(arr, dayfirst=True)
        idx2 = DatetimeIndex(np.array(arr), dayfirst=True)
        idx3 = to_datetime(arr, dayfirst=True)
        idx4 = to_datetime(np.array(arr), dayfirst=True)
        idx5 = DatetimeIndex(Index(arr), dayfirst=True)
        idx6 = DatetimeIndex(Series(arr), dayfirst=True)
        tm.assert_index_equal(expected, idx1)
        tm.assert_index_equal(expected, idx2)
        tm.assert_index_equal(expected, idx3)
        tm.assert_index_equal(expected, idx4)
        tm.assert_index_equal(expected, idx5)
        tm.assert_index_equal(expected, idx6)


class TestGuessDatetimeFormat(tm.TestCase):

    def test_guess_datetime_format_with_parseable_formats(self):
        tm._skip_if_not_us_locale()
        dt_string_to_format = (('20111230', '%Y%m%d'),
                               ('2011-12-30', '%Y-%m-%d'),
                               ('30-12-2011', '%d-%m-%Y'),
                               ('2011-12-30 00:00:00', '%Y-%m-%d %H:%M:%S'),
                               ('2011-12-30T00:00:00', '%Y-%m-%dT%H:%M:%S'),
                               ('2011-12-30 00:00:00.000000',
                                '%Y-%m-%d %H:%M:%S.%f'), )

        for dt_string, dt_format in dt_string_to_format:
            self.assertEqual(
                tools._guess_datetime_format(dt_string),
                dt_format
            )

    def test_guess_datetime_format_with_dayfirst(self):
        ambiguous_string = '01/01/2011'
        self.assertEqual(
            tools._guess_datetime_format(ambiguous_string, dayfirst=True),
            '%d/%m/%Y'
        )
        self.assertEqual(
            tools._guess_datetime_format(ambiguous_string, dayfirst=False),
            '%m/%d/%Y'
        )

    def test_guess_datetime_format_with_locale_specific_formats(self):
        # The month names will vary depending on the locale, in which
        # case these wont be parsed properly (dateutil can't parse them)
        tm._skip_if_has_locale()

        dt_string_to_format = (('30/Dec/2011', '%d/%b/%Y'),
                               ('30/December/2011', '%d/%B/%Y'),
                               ('30/Dec/2011 00:00:00', '%d/%b/%Y %H:%M:%S'), )

        for dt_string, dt_format in dt_string_to_format:
            self.assertEqual(
                tools._guess_datetime_format(dt_string),
                dt_format
            )

    def test_guess_datetime_format_invalid_inputs(self):
        # A datetime string must include a year, month and a day for it
        # to be guessable, in addition to being a string that looks like
        # a datetime
        invalid_dts = [
            '2013',
            '01/2013',
            '12:00:00',
            '1/1/1/1',
            'this_is_not_a_datetime',
            '51a',
            9,
            datetime(2011, 1, 1),
        ]

        for invalid_dt in invalid_dts:
            self.assertTrue(tools._guess_datetime_format(invalid_dt) is None)

    def test_guess_datetime_format_nopadding(self):
        # GH 11142
        dt_string_to_format = (('2011-1-1', '%Y-%m-%d'),
                               ('30-1-2011', '%d-%m-%Y'),
                               ('1/1/2011', '%m/%d/%Y'),
                               ('2011-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'),
                               ('2011-1-1 0:0:0', '%Y-%m-%d %H:%M:%S'),
                               ('2011-1-3T00:00:0', '%Y-%m-%dT%H:%M:%S'))

        for dt_string, dt_format in dt_string_to_format:
            self.assertEqual(
                tools._guess_datetime_format(dt_string),
                dt_format
            )

    def test_guess_datetime_format_for_array(self):
        tm._skip_if_not_us_locale()
        expected_format = '%Y-%m-%d %H:%M:%S.%f'
        dt_string = datetime(2011, 12, 30, 0, 0, 0).strftime(expected_format)

        test_arrays = [
            np.array([dt_string, dt_string, dt_string], dtype='O'),
            np.array([np.nan, np.nan, dt_string], dtype='O'),
            np.array([dt_string, 'random_string'], dtype='O'),
        ]

        for test_array in test_arrays:
            self.assertEqual(
                tools._guess_datetime_format_for_array(test_array),
                expected_format
            )

        format_for_string_of_nans = tools._guess_datetime_format_for_array(
            np.array(
                [np.nan, np.nan, np.nan], dtype='O'))
        self.assertTrue(format_for_string_of_nans is None)


class TestToDatetimeInferFormat(tm.TestCase):

    def test_to_datetime_infer_datetime_format_consistent_format(self):
        s = pd.Series(pd.date_range('20000101', periods=50, freq='H'))

        test_formats = ['%m-%d-%Y', '%m/%d/%Y %H:%M:%S.%f',
                        '%Y-%m-%dT%H:%M:%S.%f']

        for test_format in test_formats:
            s_as_dt_strings = s.apply(lambda x: x.strftime(test_format))

            with_format = pd.to_datetime(s_as_dt_strings, format=test_format)
            no_infer = pd.to_datetime(s_as_dt_strings,
                                      infer_datetime_format=False)
            yes_infer = pd.to_datetime(s_as_dt_strings,
                                       infer_datetime_format=True)

            # Whether the format is explicitly passed, it is inferred, or
            # it is not inferred, the results should all be the same
            self.assert_series_equal(with_format, no_infer)
            self.assert_series_equal(no_infer, yes_infer)

    def test_to_datetime_infer_datetime_format_inconsistent_format(self):
        s = pd.Series(np.array(['01/01/2011 00:00:00',
                                '01-02-2011 00:00:00',
                                '2011-01-03T00:00:00']))

        # When the format is inconsistent, infer_datetime_format should just
        # fallback to the default parsing
        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

        s = pd.Series(np.array(['Jan/01/2011', 'Feb/01/2011', 'Mar/01/2011']))

        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

    def test_to_datetime_infer_datetime_format_series_with_nans(self):
        s = pd.Series(np.array(['01/01/2011 00:00:00', np.nan,
                                '01/03/2011 00:00:00', np.nan]))
        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

    def test_to_datetime_infer_datetime_format_series_starting_with_nans(self):
        s = pd.Series(np.array([np.nan, np.nan, '01/01/2011 00:00:00',
                                '01/02/2011 00:00:00', '01/03/2011 00:00:00']))

        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

    def test_to_datetime_iso8601_noleading_0s(self):
        # GH 11871
        s = pd.Series(['2014-1-1', '2014-2-2', '2015-3-3'])
        expected = pd.Series([pd.Timestamp('2014-01-01'),
                              pd.Timestamp('2014-02-02'),
                              pd.Timestamp('2015-03-03')])
        tm.assert_series_equal(pd.to_datetime(s), expected)
        tm.assert_series_equal(pd.to_datetime(s, format='%Y-%m-%d'), expected)


class TestDaysInMonth(tm.TestCase):
    # tests for issue #10154

    def test_day_not_in_month_coerce(self):
        self.assertTrue(isnull(to_datetime('2015-02-29', errors='coerce')))
        self.assertTrue(isnull(to_datetime('2015-02-29', format="%Y-%m-%d",
                                           errors='coerce')))
        self.assertTrue(isnull(to_datetime('2015-02-32', format="%Y-%m-%d",
                                           errors='coerce')))
        self.assertTrue(isnull(to_datetime('2015-04-31', format="%Y-%m-%d",
                                           errors='coerce')))

    def test_day_not_in_month_raise(self):
        self.assertRaises(ValueError, to_datetime, '2015-02-29',
                          errors='raise')
        self.assertRaises(ValueError, to_datetime, '2015-02-29',
                          errors='raise', format="%Y-%m-%d")
        self.assertRaises(ValueError, to_datetime, '2015-02-32',
                          errors='raise', format="%Y-%m-%d")
        self.assertRaises(ValueError, to_datetime, '2015-04-31',
                          errors='raise', format="%Y-%m-%d")

    def test_day_not_in_month_ignore(self):
        self.assertEqual(to_datetime(
            '2015-02-29', errors='ignore'), '2015-02-29')
        self.assertEqual(to_datetime(
            '2015-02-29', errors='ignore', format="%Y-%m-%d"), '2015-02-29')
        self.assertEqual(to_datetime(
            '2015-02-32', errors='ignore', format="%Y-%m-%d"), '2015-02-32')
        self.assertEqual(to_datetime(
            '2015-04-31', errors='ignore', format="%Y-%m-%d"), '2015-04-31')
