import nose
from distutils.version import LooseVersion
import numpy as np

from pandas import tslib
import pandas._period as period
import datetime

import pandas as pd
from pandas.core.api import Timestamp, Series, Timedelta, Period, to_datetime
from pandas.tslib import get_timezone
from pandas._period import period_asfreq, period_ordinal
from pandas.tseries.index import date_range, DatetimeIndex
from pandas.tseries.frequencies import (
    get_freq,
    US_RESO, MS_RESO, S_RESO, H_RESO, D_RESO, T_RESO
)
import pandas.tseries.tools as tools
import pandas.tseries.offsets as offsets
import pandas.util.testing as tm
import pandas.compat as compat
from pandas.compat.numpy import (np_datetime64_compat,
                                 np_array_datetime64_compat)

from pandas.util.testing import assert_series_equal, _skip_if_has_locale


class TestTimestamp(tm.TestCase):

    def test_constructor(self):
        base_str = '2014-07-01 09:00'
        base_dt = datetime.datetime(2014, 7, 1, 9)
        base_expected = 1404205200000000000

        # confirm base representation is correct
        import calendar
        self.assertEqual(calendar.timegm(base_dt.timetuple()) * 1000000000,
                         base_expected)

        tests = [(base_str, base_dt, base_expected),
                 ('2014-07-01 10:00', datetime.datetime(2014, 7, 1, 10),
                  base_expected + 3600 * 1000000000),
                 ('2014-07-01 09:00:00.000008000',
                  datetime.datetime(2014, 7, 1, 9, 0, 0, 8),
                  base_expected + 8000),
                 ('2014-07-01 09:00:00.000000005',
                  Timestamp('2014-07-01 09:00:00.000000005'),
                  base_expected + 5)]

        tm._skip_if_no_pytz()
        tm._skip_if_no_dateutil()
        import pytz
        import dateutil
        timezones = [(None, 0), ('UTC', 0), (pytz.utc, 0), ('Asia/Tokyo', 9),
                     ('US/Eastern', -4), ('dateutil/US/Pacific', -7),
                     (pytz.FixedOffset(-180), -3),
                     (dateutil.tz.tzoffset(None, 18000), 5)]

        for date_str, date, expected in tests:
            for result in [Timestamp(date_str), Timestamp(date)]:
                # only with timestring
                self.assertEqual(result.value, expected)
                self.assertEqual(tslib.pydt_to_i8(result), expected)

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                self.assertEqual(result.value, expected)
                self.assertEqual(tslib.pydt_to_i8(result), expected)

            # with timezone
            for tz, offset in timezones:
                for result in [Timestamp(date_str, tz=tz), Timestamp(date,
                                                                     tz=tz)]:
                    expected_tz = expected - offset * 3600 * 1000000000
                    self.assertEqual(result.value, expected_tz)
                    self.assertEqual(tslib.pydt_to_i8(result), expected_tz)

                    # should preserve tz
                    result = Timestamp(result)
                    self.assertEqual(result.value, expected_tz)
                    self.assertEqual(tslib.pydt_to_i8(result), expected_tz)

                    # should convert to UTC
                    result = Timestamp(result, tz='UTC')
                    expected_utc = expected - offset * 3600 * 1000000000
                    self.assertEqual(result.value, expected_utc)
                    self.assertEqual(tslib.pydt_to_i8(result), expected_utc)

    def test_constructor_with_stringoffset(self):
        # GH 7833
        base_str = '2014-07-01 11:00:00+02:00'
        base_dt = datetime.datetime(2014, 7, 1, 9)
        base_expected = 1404205200000000000

        # confirm base representation is correct
        import calendar
        self.assertEqual(calendar.timegm(base_dt.timetuple()) * 1000000000,
                         base_expected)

        tests = [(base_str, base_expected),
                 ('2014-07-01 12:00:00+02:00',
                  base_expected + 3600 * 1000000000),
                 ('2014-07-01 11:00:00.000008000+02:00', base_expected + 8000),
                 ('2014-07-01 11:00:00.000000005+02:00', base_expected + 5)]

        tm._skip_if_no_pytz()
        tm._skip_if_no_dateutil()
        import pytz
        import dateutil
        timezones = [(None, 0), ('UTC', 0), (pytz.utc, 0), ('Asia/Tokyo', 9),
                     ('US/Eastern', -4), ('dateutil/US/Pacific', -7),
                     (pytz.FixedOffset(-180), -3),
                     (dateutil.tz.tzoffset(None, 18000), 5)]

        for date_str, expected in tests:
            for result in [Timestamp(date_str)]:
                # only with timestring
                self.assertEqual(result.value, expected)
                self.assertEqual(tslib.pydt_to_i8(result), expected)

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                self.assertEqual(result.value, expected)
                self.assertEqual(tslib.pydt_to_i8(result), expected)

            # with timezone
            for tz, offset in timezones:
                result = Timestamp(date_str, tz=tz)
                expected_tz = expected
                self.assertEqual(result.value, expected_tz)
                self.assertEqual(tslib.pydt_to_i8(result), expected_tz)

                # should preserve tz
                result = Timestamp(result)
                self.assertEqual(result.value, expected_tz)
                self.assertEqual(tslib.pydt_to_i8(result), expected_tz)

                # should convert to UTC
                result = Timestamp(result, tz='UTC')
                expected_utc = expected
                self.assertEqual(result.value, expected_utc)
                self.assertEqual(tslib.pydt_to_i8(result), expected_utc)

        # This should be 2013-11-01 05:00 in UTC
        # converted to Chicago tz
        result = Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')
        self.assertEqual(result.value, Timestamp('2013-11-01 05:00').value)
        expected = "Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')"  # noqa
        self.assertEqual(repr(result), expected)
        self.assertEqual(result, eval(repr(result)))

        # This should be 2013-11-01 05:00 in UTC
        # converted to Tokyo tz (+09:00)
        result = Timestamp('2013-11-01 00:00:00-0500', tz='Asia/Tokyo')
        self.assertEqual(result.value, Timestamp('2013-11-01 05:00').value)
        expected = "Timestamp('2013-11-01 14:00:00+0900', tz='Asia/Tokyo')"
        self.assertEqual(repr(result), expected)
        self.assertEqual(result, eval(repr(result)))

        # GH11708
        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Katmandu
        result = Timestamp("2015-11-18 15:45:00+05:45", tz="Asia/Katmandu")
        self.assertEqual(result.value, Timestamp("2015-11-18 10:00").value)
        expected = "Timestamp('2015-11-18 15:45:00+0545', tz='Asia/Katmandu')"
        self.assertEqual(repr(result), expected)
        self.assertEqual(result, eval(repr(result)))

        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Kolkata
        result = Timestamp("2015-11-18 15:30:00+05:30", tz="Asia/Kolkata")
        self.assertEqual(result.value, Timestamp("2015-11-18 10:00").value)
        expected = "Timestamp('2015-11-18 15:30:00+0530', tz='Asia/Kolkata')"
        self.assertEqual(repr(result), expected)
        self.assertEqual(result, eval(repr(result)))

    def test_constructor_invalid(self):
        with tm.assertRaisesRegexp(TypeError, 'Cannot convert input'):
            Timestamp(slice(2))
        with tm.assertRaisesRegexp(ValueError, 'Cannot convert Period'):
            Timestamp(Period('1000-01-01'))

    def test_conversion(self):
        # GH 9255
        ts = Timestamp('2000-01-01')

        result = ts.to_pydatetime()
        expected = datetime.datetime(2000, 1, 1)
        self.assertEqual(result, expected)
        self.assertEqual(type(result), type(expected))

        result = ts.to_datetime64()
        expected = np.datetime64(ts.value, 'ns')
        self.assertEqual(result, expected)
        self.assertEqual(type(result), type(expected))
        self.assertEqual(result.dtype, expected.dtype)

    def test_repr(self):
        tm._skip_if_no_pytz()
        tm._skip_if_no_dateutil()

        dates = ['2014-03-07', '2014-01-01 09:00',
                 '2014-01-01 00:00:00.000000001']

        # dateutil zone change (only matters for repr)
        import dateutil
        if dateutil.__version__ >= LooseVersion(
                '2.3') and dateutil.__version__ <= LooseVersion('2.4.0'):
            timezones = ['UTC', 'Asia/Tokyo', 'US/Eastern',
                         'dateutil/US/Pacific']
        else:
            timezones = ['UTC', 'Asia/Tokyo', 'US/Eastern',
                         'dateutil/America/Los_Angeles']

        freqs = ['D', 'M', 'S', 'N']

        for date in dates:
            for tz in timezones:
                for freq in freqs:

                    # avoid to match with timezone name
                    freq_repr = "'{0}'".format(freq)
                    if tz.startswith('dateutil'):
                        tz_repr = tz.replace('dateutil', '')
                    else:
                        tz_repr = tz

                    date_only = Timestamp(date)
                    self.assertIn(date, repr(date_only))
                    self.assertNotIn(tz_repr, repr(date_only))
                    self.assertNotIn(freq_repr, repr(date_only))
                    self.assertEqual(date_only, eval(repr(date_only)))

                    date_tz = Timestamp(date, tz=tz)
                    self.assertIn(date, repr(date_tz))
                    self.assertIn(tz_repr, repr(date_tz))
                    self.assertNotIn(freq_repr, repr(date_tz))
                    self.assertEqual(date_tz, eval(repr(date_tz)))

                    date_freq = Timestamp(date, offset=freq)
                    self.assertIn(date, repr(date_freq))
                    self.assertNotIn(tz_repr, repr(date_freq))
                    self.assertIn(freq_repr, repr(date_freq))
                    self.assertEqual(date_freq, eval(repr(date_freq)))

                    date_tz_freq = Timestamp(date, tz=tz, offset=freq)
                    self.assertIn(date, repr(date_tz_freq))
                    self.assertIn(tz_repr, repr(date_tz_freq))
                    self.assertIn(freq_repr, repr(date_tz_freq))
                    self.assertEqual(date_tz_freq, eval(repr(date_tz_freq)))

        # this can cause the tz field to be populated, but it's redundant to
        # information in the datestring
        tm._skip_if_no_pytz()
        import pytz  # noqa
        date_with_utc_offset = Timestamp('2014-03-13 00:00:00-0400', tz=None)
        self.assertIn('2014-03-13 00:00:00-0400', repr(date_with_utc_offset))
        self.assertNotIn('tzoffset', repr(date_with_utc_offset))
        self.assertIn('pytz.FixedOffset(-240)', repr(date_with_utc_offset))
        expr = repr(date_with_utc_offset).replace("'pytz.FixedOffset(-240)'",
                                                  'pytz.FixedOffset(-240)')
        self.assertEqual(date_with_utc_offset, eval(expr))

    def test_bounds_with_different_units(self):
        out_of_bounds_dates = ('1677-09-21', '2262-04-12', )

        time_units = ('D', 'h', 'm', 's', 'ms', 'us')

        for date_string in out_of_bounds_dates:
            for unit in time_units:
                self.assertRaises(ValueError, Timestamp, np.datetime64(
                    date_string, dtype='M8[%s]' % unit))

        in_bounds_dates = ('1677-09-23', '2262-04-11', )

        for date_string in in_bounds_dates:
            for unit in time_units:
                Timestamp(np.datetime64(date_string, dtype='M8[%s]' % unit))

    def test_tz(self):
        t = '2014-02-01 09:00'
        ts = Timestamp(t)
        local = ts.tz_localize('Asia/Tokyo')
        self.assertEqual(local.hour, 9)
        self.assertEqual(local, Timestamp(t, tz='Asia/Tokyo'))
        conv = local.tz_convert('US/Eastern')
        self.assertEqual(conv, Timestamp('2014-01-31 19:00', tz='US/Eastern'))
        self.assertEqual(conv.hour, 19)

        # preserves nanosecond
        ts = Timestamp(t) + offsets.Nano(5)
        local = ts.tz_localize('Asia/Tokyo')
        self.assertEqual(local.hour, 9)
        self.assertEqual(local.nanosecond, 5)
        conv = local.tz_convert('US/Eastern')
        self.assertEqual(conv.nanosecond, 5)
        self.assertEqual(conv.hour, 19)

    def test_tz_localize_ambiguous(self):

        ts = Timestamp('2014-11-02 01:00')
        ts_dst = ts.tz_localize('US/Eastern', ambiguous=True)
        ts_no_dst = ts.tz_localize('US/Eastern', ambiguous=False)

        rng = date_range('2014-11-02', periods=3, freq='H', tz='US/Eastern')
        self.assertEqual(rng[1], ts_dst)
        self.assertEqual(rng[2], ts_no_dst)
        self.assertRaises(ValueError, ts.tz_localize, 'US/Eastern',
                          ambiguous='infer')

        # GH 8025
        with tm.assertRaisesRegexp(TypeError,
                                   'Cannot localize tz-aware Timestamp, use '
                                   'tz_convert for conversions'):
            Timestamp('2011-01-01', tz='US/Eastern').tz_localize('Asia/Tokyo')

        with tm.assertRaisesRegexp(TypeError,
                                   'Cannot convert tz-naive Timestamp, use '
                                   'tz_localize to localize'):
            Timestamp('2011-01-01').tz_convert('Asia/Tokyo')

    def test_tz_localize_roundtrip(self):
        for tz in ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']:
            for t in ['2014-02-01 09:00', '2014-07-08 09:00',
                      '2014-11-01 17:00', '2014-11-05 00:00']:
                ts = Timestamp(t)
                localized = ts.tz_localize(tz)
                self.assertEqual(localized, Timestamp(t, tz=tz))

                with tm.assertRaises(TypeError):
                    localized.tz_localize(tz)

                reset = localized.tz_localize(None)
                self.assertEqual(reset, ts)
                self.assertTrue(reset.tzinfo is None)

    def test_tz_convert_roundtrip(self):
        for tz in ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']:
            for t in ['2014-02-01 09:00', '2014-07-08 09:00',
                      '2014-11-01 17:00', '2014-11-05 00:00']:
                ts = Timestamp(t, tz='UTC')
                converted = ts.tz_convert(tz)

                reset = converted.tz_convert(None)
                self.assertEqual(reset, Timestamp(t))
                self.assertTrue(reset.tzinfo is None)
                self.assertEqual(reset,
                                 converted.tz_convert('UTC').tz_localize(None))

    def test_barely_oob_dts(self):
        one_us = np.timedelta64(1).astype('timedelta64[us]')

        # By definition we can't go out of bounds in [ns], so we
        # convert the datetime64s to [us] so we can go out of bounds
        min_ts_us = np.datetime64(Timestamp.min).astype('M8[us]')
        max_ts_us = np.datetime64(Timestamp.max).astype('M8[us]')

        # No error for the min/max datetimes
        Timestamp(min_ts_us)
        Timestamp(max_ts_us)

        # One us less than the minimum is an error
        self.assertRaises(ValueError, Timestamp, min_ts_us - one_us)

        # One us more than the maximum is an error
        self.assertRaises(ValueError, Timestamp, max_ts_us + one_us)

    def test_utc_z_designator(self):
        self.assertEqual(get_timezone(
            Timestamp('2014-11-02 01:00Z').tzinfo), 'UTC')

    def test_now(self):
        # #9000
        ts_from_string = Timestamp('now')
        ts_from_method = Timestamp.now()
        ts_datetime = datetime.datetime.now()

        ts_from_string_tz = Timestamp('now', tz='US/Eastern')
        ts_from_method_tz = Timestamp.now(tz='US/Eastern')

        # Check that the delta between the times is less than 1s (arbitrarily
        # small)
        delta = Timedelta(seconds=1)
        self.assertTrue(abs(ts_from_method - ts_from_string) < delta)
        self.assertTrue(abs(ts_datetime - ts_from_method) < delta)
        self.assertTrue(abs(ts_from_method_tz - ts_from_string_tz) < delta)
        self.assertTrue(abs(ts_from_string_tz.tz_localize(None) -
                            ts_from_method_tz.tz_localize(None)) < delta)

    def test_today(self):

        ts_from_string = Timestamp('today')
        ts_from_method = Timestamp.today()
        ts_datetime = datetime.datetime.today()

        ts_from_string_tz = Timestamp('today', tz='US/Eastern')
        ts_from_method_tz = Timestamp.today(tz='US/Eastern')

        # Check that the delta between the times is less than 1s (arbitrarily
        # small)
        delta = Timedelta(seconds=1)
        self.assertTrue(abs(ts_from_method - ts_from_string) < delta)
        self.assertTrue(abs(ts_datetime - ts_from_method) < delta)
        self.assertTrue(abs(ts_from_method_tz - ts_from_string_tz) < delta)
        self.assertTrue(abs(ts_from_string_tz.tz_localize(None) -
                            ts_from_method_tz.tz_localize(None)) < delta)

    def test_asm8(self):
        np.random.seed(7960929)
        ns = [Timestamp.min.value, Timestamp.max.value, 1000, ]
        for n in ns:
            self.assertEqual(Timestamp(n).asm8.view('i8'),
                             np.datetime64(n, 'ns').view('i8'), n)
        self.assertEqual(Timestamp('nat').asm8.view('i8'),
                         np.datetime64('nat', 'ns').view('i8'))

    def test_fields(self):
        def check(value, equal):
            # that we are int/long like
            self.assertTrue(isinstance(value, (int, compat.long)))
            self.assertEqual(value, equal)

        # GH 10050
        ts = Timestamp('2015-05-10 09:06:03.000100001')
        check(ts.year, 2015)
        check(ts.month, 5)
        check(ts.day, 10)
        check(ts.hour, 9)
        check(ts.minute, 6)
        check(ts.second, 3)
        self.assertRaises(AttributeError, lambda: ts.millisecond)
        check(ts.microsecond, 100)
        check(ts.nanosecond, 1)
        check(ts.dayofweek, 6)
        check(ts.quarter, 2)
        check(ts.dayofyear, 130)
        check(ts.week, 19)
        check(ts.daysinmonth, 31)
        check(ts.daysinmonth, 31)

    def test_nat_fields(self):
        # GH 10050
        ts = Timestamp('NaT')
        self.assertTrue(np.isnan(ts.year))
        self.assertTrue(np.isnan(ts.month))
        self.assertTrue(np.isnan(ts.day))
        self.assertTrue(np.isnan(ts.hour))
        self.assertTrue(np.isnan(ts.minute))
        self.assertTrue(np.isnan(ts.second))
        self.assertTrue(np.isnan(ts.microsecond))
        self.assertTrue(np.isnan(ts.nanosecond))
        self.assertTrue(np.isnan(ts.dayofweek))
        self.assertTrue(np.isnan(ts.quarter))
        self.assertTrue(np.isnan(ts.dayofyear))
        self.assertTrue(np.isnan(ts.week))
        self.assertTrue(np.isnan(ts.daysinmonth))
        self.assertTrue(np.isnan(ts.days_in_month))

    def test_pprint(self):
        # GH12622
        import pprint
        nested_obj = {'foo': 1,
                      'bar': [{'w': {'a': Timestamp('2011-01-01')}}] * 10}
        result = pprint.pformat(nested_obj, width=50)
        expected = r"""{'bar': [{'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}}],
 'foo': 1}"""
        self.assertEqual(result, expected)


class TestDatetimeParsingWrappers(tm.TestCase):
    def test_does_not_convert_mixed_integer(self):
        bad_date_strings = ('-50000', '999', '123.1234', 'm', 'T')

        for bad_date_string in bad_date_strings:
            self.assertFalse(tslib._does_string_look_like_datetime(
                bad_date_string))

        good_date_strings = ('2012-01-01',
                             '01/01/2012',
                             'Mon Sep 16, 2013',
                             '01012012',
                             '0101',
                             '1-1', )

        for good_date_string in good_date_strings:
            self.assertTrue(tslib._does_string_look_like_datetime(
                good_date_string))

    def test_parsers(self):

        # https://github.com/dateutil/dateutil/issues/217
        import dateutil
        yearfirst = dateutil.__version__ >= LooseVersion('2.5.0')

        cases = {'2011-01-01': datetime.datetime(2011, 1, 1),
                 '2Q2005': datetime.datetime(2005, 4, 1),
                 '2Q05': datetime.datetime(2005, 4, 1),
                 '2005Q1': datetime.datetime(2005, 1, 1),
                 '05Q1': datetime.datetime(2005, 1, 1),
                 '2011Q3': datetime.datetime(2011, 7, 1),
                 '11Q3': datetime.datetime(2011, 7, 1),
                 '3Q2011': datetime.datetime(2011, 7, 1),
                 '3Q11': datetime.datetime(2011, 7, 1),

                 # quarterly without space
                 '2000Q4': datetime.datetime(2000, 10, 1),
                 '00Q4': datetime.datetime(2000, 10, 1),
                 '4Q2000': datetime.datetime(2000, 10, 1),
                 '4Q00': datetime.datetime(2000, 10, 1),
                 '2000q4': datetime.datetime(2000, 10, 1),
                 '2000-Q4': datetime.datetime(2000, 10, 1),
                 '00-Q4': datetime.datetime(2000, 10, 1),
                 '4Q-2000': datetime.datetime(2000, 10, 1),
                 '4Q-00': datetime.datetime(2000, 10, 1),
                 '2000q4': datetime.datetime(2000, 10, 1),
                 '00q4': datetime.datetime(2000, 10, 1),
                 '2005': datetime.datetime(2005, 1, 1),
                 '2005-11': datetime.datetime(2005, 11, 1),
                 '2005 11': datetime.datetime(2005, 11, 1),
                 '11-2005': datetime.datetime(2005, 11, 1),
                 '11 2005': datetime.datetime(2005, 11, 1),
                 '200511': datetime.datetime(2020, 5, 11),
                 '20051109': datetime.datetime(2005, 11, 9),
                 '20051109 10:15': datetime.datetime(2005, 11, 9, 10, 15),
                 '20051109 08H': datetime.datetime(2005, 11, 9, 8, 0),
                 '2005-11-09 10:15': datetime.datetime(2005, 11, 9, 10, 15),
                 '2005-11-09 08H': datetime.datetime(2005, 11, 9, 8, 0),
                 '2005/11/09 10:15': datetime.datetime(2005, 11, 9, 10, 15),
                 '2005/11/09 08H': datetime.datetime(2005, 11, 9, 8, 0),
                 "Thu Sep 25 10:36:28 2003": datetime.datetime(2003, 9, 25, 10,
                                                               36, 28),
                 "Thu Sep 25 2003": datetime.datetime(2003, 9, 25),
                 "Sep 25 2003": datetime.datetime(2003, 9, 25),
                 "January 1 2014": datetime.datetime(2014, 1, 1),

                 # GH 10537
                 '2014-06': datetime.datetime(2014, 6, 1),
                 '06-2014': datetime.datetime(2014, 6, 1),
                 '2014-6': datetime.datetime(2014, 6, 1),
                 '6-2014': datetime.datetime(2014, 6, 1),

                 '20010101 12': datetime.datetime(2001, 1, 1, 12),
                 '20010101 1234': datetime.datetime(2001, 1, 1, 12, 34),
                 '20010101 123456': datetime.datetime(2001, 1, 1, 12, 34, 56),
                 }

        for date_str, expected in compat.iteritems(cases):
            result1, _, _ = tools.parse_time_string(date_str,
                                                    yearfirst=yearfirst)
            result2 = to_datetime(date_str, yearfirst=yearfirst)
            result3 = to_datetime([date_str], yearfirst=yearfirst)
            result4 = to_datetime(np.array([date_str], dtype=object),
                                  yearfirst=yearfirst)
            result6 = DatetimeIndex([date_str], yearfirst=yearfirst)[0]
            self.assertEqual(result1, expected)
            self.assertEqual(result2, expected)
            self.assertEqual(result3, expected)
            self.assertEqual(result4, expected)
            self.assertEqual(result6, expected)

            # these really need to have yearfist, but we don't support
            if not yearfirst:
                result5 = Timestamp(date_str)
                self.assertEqual(result5, expected)
                result7 = date_range(date_str, freq='S', periods=1,
                                     yearfirst=yearfirst)
                self.assertEqual(result7, expected)

        # NaT
        result1, _, _ = tools.parse_time_string('NaT')
        result2 = to_datetime('NaT')
        result3 = Timestamp('NaT')
        result4 = DatetimeIndex(['NaT'])[0]
        self.assertTrue(result1 is tslib.NaT)
        self.assertTrue(result1 is tslib.NaT)
        self.assertTrue(result1 is tslib.NaT)
        self.assertTrue(result1 is tslib.NaT)

    def test_parsers_quarter_invalid(self):

        cases = ['2Q 2005', '2Q-200A', '2Q-200', '22Q2005', '6Q-20', '2Q200.']
        for case in cases:
            self.assertRaises(ValueError, tools.parse_time_string, case)

    def test_parsers_dayfirst_yearfirst(self):
        tm._skip_if_no_dateutil()

        # OK
        # 2.5.1 10-11-12   [dayfirst=0, yearfirst=0] -> 2012-10-11 00:00:00
        # 2.5.2 10-11-12   [dayfirst=0, yearfirst=1] -> 2012-10-11 00:00:00
        # 2.5.3 10-11-12   [dayfirst=0, yearfirst=0] -> 2012-10-11 00:00:00

        # OK
        # 2.5.1 10-11-12   [dayfirst=0, yearfirst=1] -> 2010-11-12 00:00:00
        # 2.5.2 10-11-12   [dayfirst=0, yearfirst=1] -> 2010-11-12 00:00:00
        # 2.5.3 10-11-12   [dayfirst=0, yearfirst=1] -> 2010-11-12 00:00:00

        # bug fix in 2.5.2
        # 2.5.1 10-11-12   [dayfirst=1, yearfirst=1] -> 2010-11-12 00:00:00
        # 2.5.2 10-11-12   [dayfirst=1, yearfirst=1] -> 2010-12-11 00:00:00
        # 2.5.3 10-11-12   [dayfirst=1, yearfirst=1] -> 2010-12-11 00:00:00

        # OK
        # 2.5.1 10-11-12   [dayfirst=1, yearfirst=0] -> 2012-11-10 00:00:00
        # 2.5.2 10-11-12   [dayfirst=1, yearfirst=0] -> 2012-11-10 00:00:00
        # 2.5.3 10-11-12   [dayfirst=1, yearfirst=0] -> 2012-11-10 00:00:00

        # OK
        # 2.5.1 20/12/21   [dayfirst=0, yearfirst=0] -> 2021-12-20 00:00:00
        # 2.5.2 20/12/21   [dayfirst=0, yearfirst=0] -> 2021-12-20 00:00:00
        # 2.5.3 20/12/21   [dayfirst=0, yearfirst=0] -> 2021-12-20 00:00:00

        # OK
        # 2.5.1 20/12/21   [dayfirst=0, yearfirst=1] -> 2020-12-21 00:00:00
        # 2.5.2 20/12/21   [dayfirst=0, yearfirst=1] -> 2020-12-21 00:00:00
        # 2.5.3 20/12/21   [dayfirst=0, yearfirst=1] -> 2020-12-21 00:00:00

        # revert of bug in 2.5.2
        # 2.5.1 20/12/21   [dayfirst=1, yearfirst=1] -> 2020-12-21 00:00:00
        # 2.5.2 20/12/21   [dayfirst=1, yearfirst=1] -> month must be in 1..12
        # 2.5.3 20/12/21   [dayfirst=1, yearfirst=1] -> 2020-12-21 00:00:00

        # OK
        # 2.5.1 20/12/21   [dayfirst=1, yearfirst=0] -> 2021-12-20 00:00:00
        # 2.5.2 20/12/21   [dayfirst=1, yearfirst=0] -> 2021-12-20 00:00:00
        # 2.5.3 20/12/21   [dayfirst=1, yearfirst=0] -> 2021-12-20 00:00:00

        import dateutil
        is_lt_253 = dateutil.__version__ < LooseVersion('2.5.3')

        # str : dayfirst, yearfirst, expected
        cases = {'10-11-12': [(False, False,
                               datetime.datetime(2012, 10, 11)),
                              (True, False,
                               datetime.datetime(2012, 11, 10)),
                              (False, True,
                               datetime.datetime(2010, 11, 12)),
                              (True, True,
                               datetime.datetime(2010, 12, 11))],
                 '20/12/21': [(False, False,
                               datetime.datetime(2021, 12, 20)),
                              (True, False,
                               datetime.datetime(2021, 12, 20)),
                              (False, True,
                               datetime.datetime(2020, 12, 21)),
                              (True, True,
                               datetime.datetime(2020, 12, 21))]}

        from dateutil.parser import parse
        for date_str, values in compat.iteritems(cases):
            for dayfirst, yearfirst, expected in values:

                # odd comparisons across version
                # let's just skip
                if dayfirst and yearfirst and is_lt_253:
                    continue

                # compare with dateutil result
                dateutil_result = parse(date_str, dayfirst=dayfirst,
                                        yearfirst=yearfirst)
                self.assertEqual(dateutil_result, expected)

                result1, _, _ = tools.parse_time_string(date_str,
                                                        dayfirst=dayfirst,
                                                        yearfirst=yearfirst)

                # we don't support dayfirst/yearfirst here:
                if not dayfirst and not yearfirst:
                    result2 = Timestamp(date_str)
                    self.assertEqual(result2, expected)

                result3 = to_datetime(date_str, dayfirst=dayfirst,
                                      yearfirst=yearfirst)

                result4 = DatetimeIndex([date_str], dayfirst=dayfirst,
                                        yearfirst=yearfirst)[0]

                self.assertEqual(result1, expected)
                self.assertEqual(result3, expected)
                self.assertEqual(result4, expected)

    def test_parsers_timestring(self):
        tm._skip_if_no_dateutil()
        from dateutil.parser import parse

        # must be the same as dateutil result
        cases = {'10:15': (parse('10:15'), datetime.datetime(1, 1, 1, 10, 15)),
                 '9:05': (parse('9:05'), datetime.datetime(1, 1, 1, 9, 5))}

        for date_str, (exp_now, exp_def) in compat.iteritems(cases):
            result1, _, _ = tools.parse_time_string(date_str)
            result2 = to_datetime(date_str)
            result3 = to_datetime([date_str])
            result4 = Timestamp(date_str)
            result5 = DatetimeIndex([date_str])[0]
            # parse time string return time string based on default date
            # others are not, and can't be changed because it is used in
            # time series plot
            self.assertEqual(result1, exp_def)
            self.assertEqual(result2, exp_now)
            self.assertEqual(result3, exp_now)
            self.assertEqual(result4, exp_now)
            self.assertEqual(result5, exp_now)

    def test_parsers_time(self):
        # GH11818
        _skip_if_has_locale()
        strings = ["14:15", "1415", "2:15pm", "0215pm", "14:15:00", "141500",
                   "2:15:00pm", "021500pm", datetime.time(14, 15)]
        expected = datetime.time(14, 15)

        for time_string in strings:
            self.assertEqual(tools.to_time(time_string), expected)

        new_string = "14.15"
        self.assertRaises(ValueError, tools.to_time, new_string)
        self.assertEqual(tools.to_time(new_string, format="%H.%M"), expected)

        arg = ["14:15", "20:20"]
        expected_arr = [datetime.time(14, 15), datetime.time(20, 20)]
        self.assertEqual(tools.to_time(arg), expected_arr)
        self.assertEqual(tools.to_time(arg, format="%H:%M"), expected_arr)
        self.assertEqual(tools.to_time(arg, infer_time_format=True),
                         expected_arr)
        self.assertEqual(tools.to_time(arg, format="%I:%M%p", errors="coerce"),
                         [None, None])
        self.assert_numpy_array_equal(tools.to_time(arg, format="%I:%M%p",
                                                    errors="ignore"),
                                      np.array(arg))
        self.assertRaises(ValueError,
                          lambda: tools.to_time(arg, format="%I:%M%p",
                                                errors="raise"))
        self.assert_series_equal(tools.to_time(Series(arg, name="test")),
                                 Series(expected_arr, name="test"))
        self.assert_numpy_array_equal(tools.to_time(np.array(arg)),
                                      np.array(expected_arr))

    def test_parsers_monthfreq(self):
        cases = {'201101': datetime.datetime(2011, 1, 1, 0, 0),
                 '200005': datetime.datetime(2000, 5, 1, 0, 0)}

        for date_str, expected in compat.iteritems(cases):
            result1, _, _ = tools.parse_time_string(date_str, freq='M')
            result2 = tools._to_datetime(date_str, freq='M')
            self.assertEqual(result1, expected)
            self.assertEqual(result2, expected)

    def test_parsers_quarterly_with_freq(self):
        msg = ('Incorrect quarterly string is given, quarter '
               'must be between 1 and 4: 2013Q5')
        with tm.assertRaisesRegexp(tslib.DateParseError, msg):
            tools.parse_time_string('2013Q5')

        # GH 5418
        msg = ('Unable to retrieve month information from given freq: '
               'INVLD-L-DEC-SAT')
        with tm.assertRaisesRegexp(tslib.DateParseError, msg):
            tools.parse_time_string('2013Q1', freq='INVLD-L-DEC-SAT')

        cases = {('2013Q2', None): datetime.datetime(2013, 4, 1),
                 ('2013Q2', 'A-APR'): datetime.datetime(2012, 8, 1),
                 ('2013-Q2', 'A-DEC'): datetime.datetime(2013, 4, 1)}

        for (date_str, freq), exp in compat.iteritems(cases):
            result, _, _ = tools.parse_time_string(date_str, freq=freq)
            self.assertEqual(result, exp)

    def test_parsers_timezone_minute_offsets_roundtrip(self):
        # GH11708
        base = to_datetime("2013-01-01 00:00:00")
        dt_strings = [
            ('2013-01-01 05:45+0545',
             "Asia/Katmandu",
             "Timestamp('2013-01-01 05:45:00+0545', tz='Asia/Katmandu')"),
            ('2013-01-01 05:30+0530',
             "Asia/Kolkata",
             "Timestamp('2013-01-01 05:30:00+0530', tz='Asia/Kolkata')")
        ]

        for dt_string, tz, dt_string_repr in dt_strings:
            dt_time = to_datetime(dt_string)
            self.assertEqual(base, dt_time)
            converted_time = dt_time.tz_localize('UTC').tz_convert(tz)
            self.assertEqual(dt_string_repr, repr(converted_time))

    def test_parsers_iso8601(self):
        # GH 12060
        # test only the iso parser - flexibility to different
        # separators and leadings 0s
        # Timestamp construction falls back to dateutil
        cases = {'2011-01-02': datetime.datetime(2011, 1, 2),
                 '2011-1-2': datetime.datetime(2011, 1, 2),
                 '2011-01': datetime.datetime(2011, 1, 1),
                 '2011-1': datetime.datetime(2011, 1, 1),
                 '2011 01 02': datetime.datetime(2011, 1, 2),
                 '2011.01.02': datetime.datetime(2011, 1, 2),
                 '2011/01/02': datetime.datetime(2011, 1, 2),
                 '2011\\01\\02': datetime.datetime(2011, 1, 2),
                 '2013-01-01 05:30:00': datetime.datetime(2013, 1, 1, 5, 30),
                 '2013-1-1 5:30:00': datetime.datetime(2013, 1, 1, 5, 30)}
        for date_str, exp in compat.iteritems(cases):
            actual = tslib._test_parse_iso8601(date_str)
            self.assertEqual(actual, exp)

        # seperators must all match - YYYYMM not valid
        invalid_cases = ['2011-01/02', '2011^11^11',
                         '201401', '201111', '200101',
                         # mixed separated and unseparated
                         '2005-0101', '200501-01',
                         '20010101 12:3456', '20010101 1234:56',
                         # HHMMSS must have two digits in each component
                         # if unseparated
                         '20010101 1', '20010101 123', '20010101 12345',
                         '20010101 12345Z',
                         # wrong separator for HHMMSS
                         '2001-01-01 12-34-56']
        for date_str in invalid_cases:
            with tm.assertRaises(ValueError):
                tslib._test_parse_iso8601(date_str)
                # If no ValueError raised, let me know which case failed.
                raise Exception(date_str)


class TestArrayToDatetime(tm.TestCase):
    def test_parsing_valid_dates(self):
        arr = np.array(['01-01-2013', '01-02-2013'], dtype=object)
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr),
            np_array_datetime64_compat(
                [
                    '2013-01-01T00:00:00.000000000-0000',
                    '2013-01-02T00:00:00.000000000-0000'
                ],
                dtype='M8[ns]'
            )
        )

        arr = np.array(['Mon Sep 16 2013', 'Tue Sep 17 2013'], dtype=object)
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr),
            np_array_datetime64_compat(
                [
                    '2013-09-16T00:00:00.000000000-0000',
                    '2013-09-17T00:00:00.000000000-0000'
                ],
                dtype='M8[ns]'
            )
        )

    def test_number_looking_strings_not_into_datetime(self):
        # #4601
        # These strings don't look like datetimes so they shouldn't be
        # attempted to be converted
        arr = np.array(['-352.737091', '183.575577'], dtype=object)
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr, errors='ignore'), arr)

        arr = np.array(['1', '2', '3', '4', '5'], dtype=object)
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr, errors='ignore'), arr)

    def test_coercing_dates_outside_of_datetime64_ns_bounds(self):
        invalid_dates = [
            datetime.date(1000, 1, 1),
            datetime.datetime(1000, 1, 1),
            '1000-01-01',
            'Jan 1, 1000',
            np.datetime64('1000-01-01'),
        ]

        for invalid_date in invalid_dates:
            self.assertRaises(ValueError,
                              tslib.array_to_datetime,
                              np.array(
                                  [invalid_date], dtype='object'),
                              errors='raise', )
            self.assert_numpy_array_equal(
                tslib.array_to_datetime(
                    np.array([invalid_date], dtype='object'),
                    errors='coerce'),
                np.array([tslib.iNaT], dtype='M8[ns]')
            )

        arr = np.array(['1/1/1000', '1/1/2000'], dtype=object)
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr, errors='coerce'),
            np_array_datetime64_compat(
                [
                    tslib.iNaT,
                    '2000-01-01T00:00:00.000000000-0000'
                ],
                dtype='M8[ns]'
            )
        )

    def test_coerce_of_invalid_datetimes(self):
        arr = np.array(['01-01-2013', 'not_a_date', '1'], dtype=object)

        # Without coercing, the presence of any invalid dates prevents
        # any values from being converted
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr, errors='ignore'), arr)

        # With coercing, the invalid dates becomes iNaT
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr, errors='coerce'),
            np_array_datetime64_compat(
                [
                    '2013-01-01T00:00:00.000000000-0000',
                    tslib.iNaT,
                    tslib.iNaT
                ],
                dtype='M8[ns]'
            )
        )

    def test_parsing_timezone_offsets(self):
        # All of these datetime strings with offsets are equivalent
        # to the same datetime after the timezone offset is added
        dt_strings = [
            '01-01-2013 08:00:00+08:00',
            '2013-01-01T08:00:00.000000000+0800',
            '2012-12-31T16:00:00.000000000-0800',
            '12-31-2012 23:00:00-01:00'
        ]

        expected_output = tslib.array_to_datetime(np.array(
            ['01-01-2013 00:00:00'], dtype=object))

        for dt_string in dt_strings:
            self.assert_numpy_array_equal(
                tslib.array_to_datetime(
                    np.array([dt_string], dtype=object)
                ),
                expected_output
            )


class TestTimestampNsOperations(tm.TestCase):
    def setUp(self):
        self.timestamp = Timestamp(datetime.datetime.utcnow())

    def assert_ns_timedelta(self, modified_timestamp, expected_value):
        value = self.timestamp.value
        modified_value = modified_timestamp.value

        self.assertEqual(modified_value - value, expected_value)

    def test_timedelta_ns_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'ns'),
                                 -123)

    def test_timedelta_ns_based_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(
            1234567898, 'ns'), 1234567898)

    def test_timedelta_us_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'us'),
                                 -123000)

    def test_timedelta_ms_arithmetic(self):
        time = self.timestamp + np.timedelta64(-123, 'ms')
        self.assert_ns_timedelta(time, -123000000)

    def test_nanosecond_string_parsing(self):
        ts = Timestamp('2013-05-01 07:15:45.123456789')
        # GH 7878
        expected_repr = '2013-05-01 07:15:45.123456789'
        expected_value = 1367392545123456789
        self.assertEqual(ts.value, expected_value)
        self.assertIn(expected_repr, repr(ts))

        ts = Timestamp('2013-05-01 07:15:45.123456789+09:00', tz='Asia/Tokyo')
        self.assertEqual(ts.value, expected_value - 9 * 3600 * 1000000000)
        self.assertIn(expected_repr, repr(ts))

        ts = Timestamp('2013-05-01 07:15:45.123456789', tz='UTC')
        self.assertEqual(ts.value, expected_value)
        self.assertIn(expected_repr, repr(ts))

        ts = Timestamp('2013-05-01 07:15:45.123456789', tz='US/Eastern')
        self.assertEqual(ts.value, expected_value + 4 * 3600 * 1000000000)
        self.assertIn(expected_repr, repr(ts))

        # GH 10041
        ts = Timestamp('20130501T071545.123456789')
        self.assertEqual(ts.value, expected_value)
        self.assertIn(expected_repr, repr(ts))

    def test_nanosecond_timestamp(self):
        # GH 7610
        expected = 1293840000000000005
        t = Timestamp('2011-01-01') + offsets.Nano(5)
        self.assertEqual(repr(t), "Timestamp('2011-01-01 00:00:00.000000005')")
        self.assertEqual(t.value, expected)
        self.assertEqual(t.nanosecond, 5)

        t = Timestamp(t)
        self.assertEqual(repr(t), "Timestamp('2011-01-01 00:00:00.000000005')")
        self.assertEqual(t.value, expected)
        self.assertEqual(t.nanosecond, 5)

        t = Timestamp(np_datetime64_compat('2011-01-01 00:00:00.000000005Z'))
        self.assertEqual(repr(t), "Timestamp('2011-01-01 00:00:00.000000005')")
        self.assertEqual(t.value, expected)
        self.assertEqual(t.nanosecond, 5)

        expected = 1293840000000000010
        t = t + offsets.Nano(5)
        self.assertEqual(repr(t), "Timestamp('2011-01-01 00:00:00.000000010')")
        self.assertEqual(t.value, expected)
        self.assertEqual(t.nanosecond, 10)

        t = Timestamp(t)
        self.assertEqual(repr(t), "Timestamp('2011-01-01 00:00:00.000000010')")
        self.assertEqual(t.value, expected)
        self.assertEqual(t.nanosecond, 10)

        t = Timestamp(np_datetime64_compat('2011-01-01 00:00:00.000000010Z'))
        self.assertEqual(repr(t), "Timestamp('2011-01-01 00:00:00.000000010')")
        self.assertEqual(t.value, expected)
        self.assertEqual(t.nanosecond, 10)

    def test_nat_arithmetic(self):
        # GH 6873
        i = 2
        f = 1.5

        for (left, right) in [(pd.NaT, i), (pd.NaT, f), (pd.NaT, np.nan)]:
            self.assertIs(left / right, pd.NaT)
            self.assertIs(left * right, pd.NaT)
            self.assertIs(right * left, pd.NaT)
            with tm.assertRaises(TypeError):
                right / left

        # Timestamp / datetime
        t = Timestamp('2014-01-01')
        dt = datetime.datetime(2014, 1, 1)
        for (left, right) in [(pd.NaT, pd.NaT), (pd.NaT, t), (pd.NaT, dt)]:
            # NaT __add__ or __sub__ Timestamp-like (or inverse) returns NaT
            self.assertIs(right + left, pd.NaT)
            self.assertIs(left + right, pd.NaT)
            self.assertIs(left - right, pd.NaT)
            self.assertIs(right - left, pd.NaT)

        # timedelta-like
        # offsets are tested in test_offsets.py

        delta = datetime.timedelta(3600)
        td = Timedelta('5s')

        for (left, right) in [(pd.NaT, delta), (pd.NaT, td)]:
            # NaT + timedelta-like returns NaT
            self.assertIs(right + left, pd.NaT)
            self.assertIs(left + right, pd.NaT)
            self.assertIs(right - left, pd.NaT)
            self.assertIs(left - right, pd.NaT)

        # GH 11718
        tm._skip_if_no_pytz()
        import pytz

        t_utc = Timestamp('2014-01-01', tz='UTC')
        t_tz = Timestamp('2014-01-01', tz='US/Eastern')
        dt_tz = pytz.timezone('Asia/Tokyo').localize(dt)

        for (left, right) in [(pd.NaT, t_utc), (pd.NaT, t_tz),
                              (pd.NaT, dt_tz)]:
            # NaT __add__ or __sub__ Timestamp-like (or inverse) returns NaT
            self.assertIs(right + left, pd.NaT)
            self.assertIs(left + right, pd.NaT)
            self.assertIs(left - right, pd.NaT)
            self.assertIs(right - left, pd.NaT)

    def test_nat_arithmetic_index(self):
        # GH 11718

        # datetime
        tm._skip_if_no_pytz()

        dti = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], name='x')
        exp = pd.DatetimeIndex([pd.NaT, pd.NaT], name='x')
        self.assert_index_equal(dti + pd.NaT, exp)
        self.assert_index_equal(pd.NaT + dti, exp)

        dti_tz = pd.DatetimeIndex(['2011-01-01', '2011-01-02'],
                                  tz='US/Eastern', name='x')
        exp = pd.DatetimeIndex([pd.NaT, pd.NaT], name='x', tz='US/Eastern')
        self.assert_index_equal(dti_tz + pd.NaT, exp)
        self.assert_index_equal(pd.NaT + dti_tz, exp)

        exp = pd.TimedeltaIndex([pd.NaT, pd.NaT], name='x')
        for (left, right) in [(pd.NaT, dti), (pd.NaT, dti_tz)]:
            self.assert_index_equal(left - right, exp)
            self.assert_index_equal(right - left, exp)

        # timedelta
        tdi = pd.TimedeltaIndex(['1 day', '2 day'], name='x')
        exp = pd.DatetimeIndex([pd.NaT, pd.NaT], name='x')
        for (left, right) in [(pd.NaT, tdi)]:
            self.assert_index_equal(left + right, exp)
            self.assert_index_equal(right + left, exp)
            self.assert_index_equal(left - right, exp)
            self.assert_index_equal(right - left, exp)


class TestTslib(tm.TestCase):
    def test_intraday_conversion_factors(self):
        self.assertEqual(period_asfreq(
            1, get_freq('D'), get_freq('H'), False), 24)
        self.assertEqual(period_asfreq(
            1, get_freq('D'), get_freq('T'), False), 1440)
        self.assertEqual(period_asfreq(
            1, get_freq('D'), get_freq('S'), False), 86400)
        self.assertEqual(period_asfreq(1, get_freq(
            'D'), get_freq('L'), False), 86400000)
        self.assertEqual(period_asfreq(1, get_freq(
            'D'), get_freq('U'), False), 86400000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'D'), get_freq('N'), False), 86400000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('H'), get_freq('T'), False), 60)
        self.assertEqual(period_asfreq(
            1, get_freq('H'), get_freq('S'), False), 3600)
        self.assertEqual(period_asfreq(1, get_freq('H'),
                                       get_freq('L'), False), 3600000)
        self.assertEqual(period_asfreq(1, get_freq(
            'H'), get_freq('U'), False), 3600000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'H'), get_freq('N'), False), 3600000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('T'), get_freq('S'), False), 60)
        self.assertEqual(period_asfreq(
            1, get_freq('T'), get_freq('L'), False), 60000)
        self.assertEqual(period_asfreq(1, get_freq(
            'T'), get_freq('U'), False), 60000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'T'), get_freq('N'), False), 60000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('S'), get_freq('L'), False), 1000)
        self.assertEqual(period_asfreq(1, get_freq('S'),
                                       get_freq('U'), False), 1000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'S'), get_freq('N'), False), 1000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('L'), get_freq('U'), False), 1000)
        self.assertEqual(period_asfreq(1, get_freq('L'),
                                       get_freq('N'), False), 1000000)

        self.assertEqual(period_asfreq(
            1, get_freq('U'), get_freq('N'), False), 1000)

    def test_period_ordinal_start_values(self):
        # information for 1.1.1970
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('A')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('M')))
        self.assertEqual(1, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('W')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('D')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('B')))

    def test_period_ordinal_week(self):
        self.assertEqual(1, period_ordinal(1970, 1, 4, 0, 0, 0, 0, 0,
                                           get_freq('W')))
        self.assertEqual(2, period_ordinal(1970, 1, 5, 0, 0, 0, 0, 0,
                                           get_freq('W')))

        self.assertEqual(2284, period_ordinal(2013, 10, 6, 0, 0, 0, 0, 0,
                                              get_freq('W')))
        self.assertEqual(2285, period_ordinal(2013, 10, 7, 0, 0, 0, 0, 0,
                                              get_freq('W')))

    def test_period_ordinal_business_day(self):
        # Thursday
        self.assertEqual(11415, period_ordinal(2013, 10, 3, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Friday
        self.assertEqual(11416, period_ordinal(2013, 10, 4, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Saturday
        self.assertEqual(11417, period_ordinal(2013, 10, 5, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Sunday
        self.assertEqual(11417, period_ordinal(2013, 10, 6, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Monday
        self.assertEqual(11417, period_ordinal(2013, 10, 7, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Tuesday
        self.assertEqual(11418, period_ordinal(2013, 10, 8, 0, 0, 0, 0, 0,
                                               get_freq('B')))

    def test_tslib_tz_convert(self):
        def compare_utc_to_local(tz_didx, utc_didx):
            f = lambda x: tslib.tz_convert_single(x, 'UTC', tz_didx.tz)
            result = tslib.tz_convert(tz_didx.asi8, 'UTC', tz_didx.tz)
            result_single = np.vectorize(f)(tz_didx.asi8)
            self.assert_numpy_array_equal(result, result_single)

        def compare_local_to_utc(tz_didx, utc_didx):
            f = lambda x: tslib.tz_convert_single(x, tz_didx.tz, 'UTC')
            result = tslib.tz_convert(utc_didx.asi8, tz_didx.tz, 'UTC')
            result_single = np.vectorize(f)(utc_didx.asi8)
            self.assert_numpy_array_equal(result, result_single)

        for tz in ['UTC', 'Asia/Tokyo', 'US/Eastern', 'Europe/Moscow']:
            # US: 2014-03-09 - 2014-11-11
            # MOSCOW: 2014-10-26  /  2014-12-31
            tz_didx = date_range('2014-03-01', '2015-01-10', freq='H', tz=tz)
            utc_didx = date_range('2014-03-01', '2015-01-10', freq='H')
            compare_utc_to_local(tz_didx, utc_didx)
            # local tz to UTC can be differ in hourly (or higher) freqs because
            # of DST
            compare_local_to_utc(tz_didx, utc_didx)

            tz_didx = date_range('2000-01-01', '2020-01-01', freq='D', tz=tz)
            utc_didx = date_range('2000-01-01', '2020-01-01', freq='D')
            compare_utc_to_local(tz_didx, utc_didx)
            compare_local_to_utc(tz_didx, utc_didx)

            tz_didx = date_range('2000-01-01', '2100-01-01', freq='A', tz=tz)
            utc_didx = date_range('2000-01-01', '2100-01-01', freq='A')
            compare_utc_to_local(tz_didx, utc_didx)
            compare_local_to_utc(tz_didx, utc_didx)

        # Check empty array
        result = tslib.tz_convert(np.array([], dtype=np.int64),
                                  tslib.maybe_get_tz('US/Eastern'),
                                  tslib.maybe_get_tz('Asia/Tokyo'))
        self.assert_numpy_array_equal(result, np.array([], dtype=np.int64))

        # Check all-NaT array
        result = tslib.tz_convert(np.array([tslib.iNaT], dtype=np.int64),
                                  tslib.maybe_get_tz('US/Eastern'),
                                  tslib.maybe_get_tz('Asia/Tokyo'))
        self.assert_numpy_array_equal(result, np.array(
            [tslib.iNaT], dtype=np.int64))

    def test_shift_months(self):
        s = DatetimeIndex([Timestamp('2000-01-05 00:15:00'), Timestamp(
            '2000-01-31 00:23:00'), Timestamp('2000-01-01'), Timestamp(
                '2000-02-29'), Timestamp('2000-12-31')])
        for years in [-1, 0, 1]:
            for months in [-2, 0, 2]:
                actual = DatetimeIndex(tslib.shift_months(s.asi8, years * 12 +
                                                          months))
                expected = DatetimeIndex([x + offsets.DateOffset(
                    years=years, months=months) for x in s])
                tm.assert_index_equal(actual, expected)


class TestTimestampOps(tm.TestCase):
    def test_timestamp_and_datetime(self):
        self.assertEqual((Timestamp(datetime.datetime(
            2013, 10, 13)) - datetime.datetime(2013, 10, 12)).days, 1)
        self.assertEqual((datetime.datetime(2013, 10, 12) -
                          Timestamp(datetime.datetime(2013, 10, 13))).days, -1)

    def test_timestamp_and_series(self):
        timestamp_series = Series(date_range('2014-03-17', periods=2, freq='D',
                                             tz='US/Eastern'))
        first_timestamp = timestamp_series[0]

        delta_series = Series([np.timedelta64(0, 'D'), np.timedelta64(1, 'D')])
        assert_series_equal(timestamp_series - first_timestamp, delta_series)
        assert_series_equal(first_timestamp - timestamp_series, -delta_series)

    def test_addition_subtraction_types(self):
        # Assert on the types resulting from Timestamp +/- various date/time
        # objects
        datetime_instance = datetime.datetime(2014, 3, 4)
        timedelta_instance = datetime.timedelta(seconds=1)
        # build a timestamp with a frequency, since then it supports
        # addition/subtraction of integers
        timestamp_instance = date_range(datetime_instance, periods=1,
                                        freq='D')[0]

        self.assertEqual(type(timestamp_instance + 1), Timestamp)
        self.assertEqual(type(timestamp_instance - 1), Timestamp)

        # Timestamp + datetime not supported, though subtraction is supported
        # and yields timedelta more tests in tseries/base/tests/test_base.py
        self.assertEqual(
            type(timestamp_instance - datetime_instance), Timedelta)
        self.assertEqual(
            type(timestamp_instance + timedelta_instance), Timestamp)
        self.assertEqual(
            type(timestamp_instance - timedelta_instance), Timestamp)

        # Timestamp +/- datetime64 not supported, so not tested (could possibly
        # assert error raised?)
        timedelta64_instance = np.timedelta64(1, 'D')
        self.assertEqual(
            type(timestamp_instance + timedelta64_instance), Timestamp)
        self.assertEqual(
            type(timestamp_instance - timedelta64_instance), Timestamp)

    def test_addition_subtraction_preserve_frequency(self):
        timestamp_instance = date_range('2014-03-05', periods=1, freq='D')[0]
        timedelta_instance = datetime.timedelta(days=1)
        original_freq = timestamp_instance.freq
        self.assertEqual((timestamp_instance + 1).freq, original_freq)
        self.assertEqual((timestamp_instance - 1).freq, original_freq)
        self.assertEqual(
            (timestamp_instance + timedelta_instance).freq, original_freq)
        self.assertEqual(
            (timestamp_instance - timedelta_instance).freq, original_freq)

        timedelta64_instance = np.timedelta64(1, 'D')
        self.assertEqual(
            (timestamp_instance + timedelta64_instance).freq, original_freq)
        self.assertEqual(
            (timestamp_instance - timedelta64_instance).freq, original_freq)

    def test_resolution(self):

        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H', 'T',
                                   'S', 'L', 'U'],
                                  [D_RESO, D_RESO,
                                   D_RESO, D_RESO,
                                   H_RESO, T_RESO,
                                   S_RESO, MS_RESO,
                                   US_RESO]):
            for tz in [None, 'Asia/Tokyo', 'US/Eastern',
                       'dateutil/US/Eastern']:
                idx = date_range(start='2013-04-01', periods=30, freq=freq,
                                 tz=tz)
                result = period.resolution(idx.asi8, idx.tz)
                self.assertEqual(result, expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
