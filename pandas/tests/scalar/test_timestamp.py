""" test the scalar Timestamp """

import sys
import operator
import calendar
import numpy as np
from datetime import datetime, timedelta
from distutils.version import LooseVersion

import pandas as pd
import pandas.util.testing as tm
import pandas._period as period
from pandas.tseries import offsets, frequencies
from pandas.tslib import get_timezone, iNaT
from pandas.compat import lrange, long
from pandas.util.testing import assert_series_equal
from pandas.compat.numpy import np_datetime64_compat
from pandas import (Timestamp, date_range, Period, Timedelta, tslib, compat,
                    Series, NaT, isnull, DataFrame, DatetimeIndex)
from pandas.tseries.frequencies import (RESO_DAY, RESO_HR, RESO_MIN, RESO_US,
                                        RESO_MS, RESO_SEC)

randn = np.random.randn


class TestTimestamp(tm.TestCase):

    def test_constructor(self):
        base_str = '2014-07-01 09:00'
        base_dt = datetime(2014, 7, 1, 9)
        base_expected = 1404205200000000000

        # confirm base representation is correct
        import calendar
        self.assertEqual(calendar.timegm(base_dt.timetuple()) * 1000000000,
                         base_expected)

        tests = [(base_str, base_dt, base_expected),
                 ('2014-07-01 10:00', datetime(2014, 7, 1, 10),
                  base_expected + 3600 * 1000000000),
                 ('2014-07-01 09:00:00.000008000',
                  datetime(2014, 7, 1, 9, 0, 0, 8),
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
        base_dt = datetime(2014, 7, 1, 9)
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

    def test_constructor_positional(self):
        # GH 10758
        with tm.assertRaises(TypeError):
            Timestamp(2000, 1)
        with tm.assertRaises(ValueError):
            Timestamp(2000, 0, 1)
        with tm.assertRaises(ValueError):
            Timestamp(2000, 13, 1)
        with tm.assertRaises(ValueError):
            Timestamp(2000, 1, 0)
        with tm.assertRaises(ValueError):
            Timestamp(2000, 1, 32)

        # GH 11630
        self.assertEqual(
            repr(Timestamp(2015, 11, 12)),
            repr(Timestamp('20151112')))

        self.assertEqual(
            repr(Timestamp(2015, 11, 12, 1, 2, 3, 999999)),
            repr(Timestamp('2015-11-12 01:02:03.999999')))

        self.assertIs(Timestamp(None), pd.NaT)

    def test_constructor_keyword(self):
        # GH 10758
        with tm.assertRaises(TypeError):
            Timestamp(year=2000, month=1)
        with tm.assertRaises(ValueError):
            Timestamp(year=2000, month=0, day=1)
        with tm.assertRaises(ValueError):
            Timestamp(year=2000, month=13, day=1)
        with tm.assertRaises(ValueError):
            Timestamp(year=2000, month=1, day=0)
        with tm.assertRaises(ValueError):
            Timestamp(year=2000, month=1, day=32)

        self.assertEqual(
            repr(Timestamp(year=2015, month=11, day=12)),
            repr(Timestamp('20151112')))

        self.assertEqual(
            repr(Timestamp(year=2015, month=11, day=12,
                           hour=1, minute=2, second=3, microsecond=999999)),
            repr(Timestamp('2015-11-12 01:02:03.999999')))

    def test_constructor_fromordinal(self):
        base = datetime(2000, 1, 1)

        ts = Timestamp.fromordinal(base.toordinal(), freq='D')
        self.assertEqual(base, ts)
        self.assertEqual(ts.freq, 'D')
        self.assertEqual(base.toordinal(), ts.toordinal())

        ts = Timestamp.fromordinal(base.toordinal(), tz='US/Eastern')
        self.assertEqual(pd.Timestamp('2000-01-01', tz='US/Eastern'), ts)
        self.assertEqual(base.toordinal(), ts.toordinal())

    def test_constructor_offset_depr(self):
        # GH 12160
        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            ts = Timestamp('2011-01-01', offset='D')
        self.assertEqual(ts.freq, 'D')

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            self.assertEqual(ts.offset, 'D')

        msg = "Can only specify freq or offset, not both"
        with tm.assertRaisesRegexp(TypeError, msg):
            Timestamp('2011-01-01', offset='D', freq='D')

    def test_constructor_offset_depr_fromordinal(self):
        # GH 12160
        base = datetime(2000, 1, 1)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            ts = Timestamp.fromordinal(base.toordinal(), offset='D')
        self.assertEqual(pd.Timestamp('2000-01-01'), ts)
        self.assertEqual(ts.freq, 'D')
        self.assertEqual(base.toordinal(), ts.toordinal())

        msg = "Can only specify freq or offset, not both"
        with tm.assertRaisesRegexp(TypeError, msg):
            Timestamp.fromordinal(base.toordinal(), offset='D', freq='D')

    def test_conversion(self):
        # GH 9255
        ts = Timestamp('2000-01-01')

        result = ts.to_pydatetime()
        expected = datetime(2000, 1, 1)
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
        if (dateutil.__version__ >= LooseVersion('2.3') and
            (dateutil.__version__ <= LooseVersion('2.4.0') or
             dateutil.__version__ >= LooseVersion('2.6.0'))):
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

                    date_freq = Timestamp(date, freq=freq)
                    self.assertIn(date, repr(date_freq))
                    self.assertNotIn(tz_repr, repr(date_freq))
                    self.assertIn(freq_repr, repr(date_freq))
                    self.assertEqual(date_freq, eval(repr(date_freq)))

                    date_tz_freq = Timestamp(date, tz=tz, freq=freq)
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

    def test_tz_localize_nonexistent(self):
        # See issue 13057
        from pytz.exceptions import NonExistentTimeError
        times = ['2015-03-08 02:00', '2015-03-08 02:30',
                 '2015-03-29 02:00', '2015-03-29 02:30']
        timezones = ['US/Eastern', 'US/Pacific',
                     'Europe/Paris', 'Europe/Belgrade']
        for t, tz in zip(times, timezones):
            ts = Timestamp(t)
            self.assertRaises(NonExistentTimeError, ts.tz_localize,
                              tz)
            self.assertRaises(NonExistentTimeError, ts.tz_localize,
                              tz, errors='raise')
            self.assertIs(ts.tz_localize(tz, errors='coerce'),
                          pd.NaT)

    def test_tz_localize_errors_ambiguous(self):
        # See issue 13057
        from pytz.exceptions import AmbiguousTimeError
        ts = pd.Timestamp('2015-11-1 01:00')
        self.assertRaises(AmbiguousTimeError,
                          ts.tz_localize, 'US/Pacific', errors='coerce')

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
        ts_datetime = datetime.now()

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
        ts_datetime = datetime.today()

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

    def test_nat_vector_field_access(self):
        idx = DatetimeIndex(['1/1/2000', None, None, '1/4/2000'])

        fields = ['year', 'quarter', 'month', 'day', 'hour', 'minute',
                  'second', 'microsecond', 'nanosecond', 'week', 'dayofyear',
                  'days_in_month', 'is_leap_year']

        for field in fields:
            result = getattr(idx, field)
            expected = [getattr(x, field) for x in idx]
            self.assert_numpy_array_equal(result, np.array(expected))

        s = pd.Series(idx)

        for field in fields:
            result = getattr(s.dt, field)
            expected = [getattr(x, field) for x in idx]
            self.assert_series_equal(result, pd.Series(expected))

    def test_nat_scalar_field_access(self):
        fields = ['year', 'quarter', 'month', 'day', 'hour', 'minute',
                  'second', 'microsecond', 'nanosecond', 'week', 'dayofyear',
                  'days_in_month', 'daysinmonth', 'dayofweek', 'weekday_name']
        for field in fields:
            result = getattr(NaT, field)
            self.assertTrue(np.isnan(result))

    def test_NaT_methods(self):
        # GH 9513
        raise_methods = ['astimezone', 'combine', 'ctime', 'dst',
                         'fromordinal', 'fromtimestamp', 'isocalendar',
                         'strftime', 'strptime', 'time', 'timestamp',
                         'timetuple', 'timetz', 'toordinal', 'tzname',
                         'utcfromtimestamp', 'utcnow', 'utcoffset',
                         'utctimetuple']
        nat_methods = ['date', 'now', 'replace', 'to_datetime', 'today']
        nan_methods = ['weekday', 'isoweekday']

        for method in raise_methods:
            if hasattr(NaT, method):
                self.assertRaises(ValueError, getattr(NaT, method))

        for method in nan_methods:
            if hasattr(NaT, method):
                self.assertTrue(np.isnan(getattr(NaT, method)()))

        for method in nat_methods:
            if hasattr(NaT, method):
                # see gh-8254
                exp_warning = None
                if method == 'to_datetime':
                    exp_warning = FutureWarning
                with tm.assert_produces_warning(
                        exp_warning, check_stacklevel=False):
                    self.assertIs(getattr(NaT, method)(), NaT)

        # GH 12300
        self.assertEqual(NaT.isoformat(), 'NaT')

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

    def to_datetime_depr(self):
        # see gh-8254
        ts = Timestamp('2011-01-01')

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            expected = datetime(2011, 1, 1)
            result = ts.to_datetime()
            self.assertEqual(result, expected)

    def to_pydatetime_nonzero_nano(self):
        ts = Timestamp('2011-01-01 9:00:00.123456789')

        # Warn the user of data loss (nanoseconds).
        with tm.assert_produces_warning(UserWarning,
                                        check_stacklevel=False):
            expected = datetime(2011, 1, 1, 9, 0, 0, 123456)
            result = ts.to_pydatetime()
            self.assertEqual(result, expected)

    def test_round(self):

        # round
        dt = Timestamp('20130101 09:10:11')
        result = dt.round('D')
        expected = Timestamp('20130101')
        self.assertEqual(result, expected)

        dt = Timestamp('20130101 19:10:11')
        result = dt.round('D')
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

        dt = Timestamp('20130201 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130202')
        self.assertEqual(result, expected)

        dt = Timestamp('20130104 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130105')
        self.assertEqual(result, expected)

        dt = Timestamp('20130104 12:32:00')
        result = dt.round('30Min')
        expected = Timestamp('20130104 12:30:00')
        self.assertEqual(result, expected)

        dti = date_range('20130101 09:10:11', periods=5)
        result = dti.round('D')
        expected = date_range('20130101', periods=5)
        tm.assert_index_equal(result, expected)

        # floor
        dt = Timestamp('20130101 09:10:11')
        result = dt.floor('D')
        expected = Timestamp('20130101')
        self.assertEqual(result, expected)

        # ceil
        dt = Timestamp('20130101 09:10:11')
        result = dt.ceil('D')
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

        # round with tz
        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('D')
        expected = Timestamp('20130101', tz='US/Eastern')
        self.assertEqual(result, expected)

        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('s')
        self.assertEqual(result, dt)

        dti = date_range('20130101 09:10:11',
                         periods=5).tz_localize('UTC').tz_convert('US/Eastern')
        result = dti.round('D')
        expected = date_range('20130101', periods=5).tz_localize('US/Eastern')
        tm.assert_index_equal(result, expected)

        result = dti.round('s')
        tm.assert_index_equal(result, dti)

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            self.assertRaises(ValueError, lambda: dti.round(freq))

    def test_class_ops_pytz(self):
        tm._skip_if_no_pytz()
        from pytz import timezone

        def compare(x, y):
            self.assertEqual(int(Timestamp(x).value / 1e9),
                             int(Timestamp(y).value / 1e9))

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now('UTC'), datetime.now(timezone('UTC')))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(Timestamp.utcfromtimestamp(current_time),
                datetime.utcfromtimestamp(current_time))
        compare(Timestamp.fromtimestamp(current_time),
                datetime.fromtimestamp(current_time))

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(Timestamp.combine(date_component, time_component),
                datetime.combine(date_component, time_component))

    def test_class_ops_dateutil(self):
        tm._skip_if_no_dateutil()
        from dateutil.tz import tzutc

        def compare(x, y):
            self.assertEqual(int(np.round(Timestamp(x).value / 1e9)),
                             int(np.round(Timestamp(y).value / 1e9)))

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now('UTC'), datetime.now(tzutc()))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(Timestamp.utcfromtimestamp(current_time),
                datetime.utcfromtimestamp(current_time))
        compare(Timestamp.fromtimestamp(current_time),
                datetime.fromtimestamp(current_time))

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(Timestamp.combine(date_component, time_component),
                datetime.combine(date_component, time_component))

    def test_basics_nanos(self):
        val = np.int64(946684800000000000).view('M8[ns]')
        stamp = Timestamp(val.view('i8') + 500)
        self.assertEqual(stamp.year, 2000)
        self.assertEqual(stamp.month, 1)
        self.assertEqual(stamp.microsecond, 0)
        self.assertEqual(stamp.nanosecond, 500)

        # GH 14415
        val = np.iinfo(np.int64).min + 80000000000000
        stamp = Timestamp(val)
        self.assertEqual(stamp.year, 1677)
        self.assertEqual(stamp.month, 9)
        self.assertEqual(stamp.day, 21)
        self.assertEqual(stamp.microsecond, 145224)
        self.assertEqual(stamp.nanosecond, 192)

    def test_unit(self):

        def check(val, unit=None, h=1, s=1, us=0):
            stamp = Timestamp(val, unit=unit)
            self.assertEqual(stamp.year, 2000)
            self.assertEqual(stamp.month, 1)
            self.assertEqual(stamp.day, 1)
            self.assertEqual(stamp.hour, h)
            if unit != 'D':
                self.assertEqual(stamp.minute, 1)
                self.assertEqual(stamp.second, s)
                self.assertEqual(stamp.microsecond, us)
            else:
                self.assertEqual(stamp.minute, 0)
                self.assertEqual(stamp.second, 0)
                self.assertEqual(stamp.microsecond, 0)
            self.assertEqual(stamp.nanosecond, 0)

        ts = Timestamp('20000101 01:01:01')
        val = ts.value
        days = (ts - Timestamp('1970-01-01')).days

        check(val)
        check(val / long(1000), unit='us')
        check(val / long(1000000), unit='ms')
        check(val / long(1000000000), unit='s')
        check(days, unit='D', h=0)

        # using truediv, so these are like floats
        if compat.PY3:
            check((val + 500000) / long(1000000000), unit='s', us=500)
            check((val + 500000000) / long(1000000000), unit='s', us=500000)
            check((val + 500000) / long(1000000), unit='ms', us=500)

        # get chopped in py2
        else:
            check((val + 500000) / long(1000000000), unit='s')
            check((val + 500000000) / long(1000000000), unit='s')
            check((val + 500000) / long(1000000), unit='ms')

        # ok
        check((val + 500000) / long(1000), unit='us', us=500)
        check((val + 500000000) / long(1000000), unit='ms', us=500000)

        # floats
        check(val / 1000.0 + 5, unit='us', us=5)
        check(val / 1000.0 + 5000, unit='us', us=5000)
        check(val / 1000000.0 + 0.5, unit='ms', us=500)
        check(val / 1000000.0 + 0.005, unit='ms', us=5)
        check(val / 1000000000.0 + 0.5, unit='s', us=500000)
        check(days + 0.5, unit='D', h=12)

        # nan
        result = Timestamp(np.nan)
        self.assertIs(result, NaT)

        result = Timestamp(None)
        self.assertIs(result, NaT)

        result = Timestamp(iNaT)
        self.assertIs(result, NaT)

        result = Timestamp(NaT)
        self.assertIs(result, NaT)

        result = Timestamp('NaT')
        self.assertIs(result, NaT)

        self.assertTrue(isnull(Timestamp('nat')))

    def test_roundtrip(self):

        # test value to string and back conversions
        # further test accessors
        base = Timestamp('20140101 00:00:00')

        result = Timestamp(base.value + pd.Timedelta('5ms').value)
        self.assertEqual(result, Timestamp(str(base) + ".005000"))
        self.assertEqual(result.microsecond, 5000)

        result = Timestamp(base.value + pd.Timedelta('5us').value)
        self.assertEqual(result, Timestamp(str(base) + ".000005"))
        self.assertEqual(result.microsecond, 5)

        result = Timestamp(base.value + pd.Timedelta('5ns').value)
        self.assertEqual(result, Timestamp(str(base) + ".000000005"))
        self.assertEqual(result.nanosecond, 5)
        self.assertEqual(result.microsecond, 0)

        result = Timestamp(base.value + pd.Timedelta('6ms 5us').value)
        self.assertEqual(result, Timestamp(str(base) + ".006005"))
        self.assertEqual(result.microsecond, 5 + 6 * 1000)

        result = Timestamp(base.value + pd.Timedelta('200ms 5us').value)
        self.assertEqual(result, Timestamp(str(base) + ".200005"))
        self.assertEqual(result.microsecond, 5 + 200 * 1000)

    def test_comparison(self):
        # 5-18-2012 00:00:00.000
        stamp = long(1337299200000000000)

        val = Timestamp(stamp)

        self.assertEqual(val, val)
        self.assertFalse(val != val)
        self.assertFalse(val < val)
        self.assertTrue(val <= val)
        self.assertFalse(val > val)
        self.assertTrue(val >= val)

        other = datetime(2012, 5, 18)
        self.assertEqual(val, other)
        self.assertFalse(val != other)
        self.assertFalse(val < other)
        self.assertTrue(val <= other)
        self.assertFalse(val > other)
        self.assertTrue(val >= other)

        other = Timestamp(stamp + 100)

        self.assertNotEqual(val, other)
        self.assertNotEqual(val, other)
        self.assertTrue(val < other)
        self.assertTrue(val <= other)
        self.assertTrue(other > val)
        self.assertTrue(other >= val)

    def test_compare_invalid(self):

        # GH 8058
        val = Timestamp('20130101 12:01:02')
        self.assertFalse(val == 'foo')
        self.assertFalse(val == 10.0)
        self.assertFalse(val == 1)
        self.assertFalse(val == long(1))
        self.assertFalse(val == [])
        self.assertFalse(val == {'foo': 1})
        self.assertFalse(val == np.float64(1))
        self.assertFalse(val == np.int64(1))

        self.assertTrue(val != 'foo')
        self.assertTrue(val != 10.0)
        self.assertTrue(val != 1)
        self.assertTrue(val != long(1))
        self.assertTrue(val != [])
        self.assertTrue(val != {'foo': 1})
        self.assertTrue(val != np.float64(1))
        self.assertTrue(val != np.int64(1))

        # ops testing
        df = DataFrame(randn(5, 2))
        a = df[0]
        b = Series(randn(5))
        b.name = Timestamp('2000-01-01')
        tm.assert_series_equal(a / b, 1 / (b / a))

    def test_cant_compare_tz_naive_w_aware(self):
        tm._skip_if_no_pytz()
        # #1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz='utc')

        self.assertRaises(Exception, a.__eq__, b)
        self.assertRaises(Exception, a.__ne__, b)
        self.assertRaises(Exception, a.__lt__, b)
        self.assertRaises(Exception, a.__gt__, b)
        self.assertRaises(Exception, b.__eq__, a)
        self.assertRaises(Exception, b.__ne__, a)
        self.assertRaises(Exception, b.__lt__, a)
        self.assertRaises(Exception, b.__gt__, a)

        if sys.version_info < (3, 3):
            self.assertRaises(Exception, a.__eq__, b.to_pydatetime())
            self.assertRaises(Exception, a.to_pydatetime().__eq__, b)
        else:
            self.assertFalse(a == b.to_pydatetime())
            self.assertFalse(a.to_pydatetime() == b)

    def test_cant_compare_tz_naive_w_aware_explicit_pytz(self):
        tm._skip_if_no_pytz()
        from pytz import utc
        # #1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz=utc)

        self.assertRaises(Exception, a.__eq__, b)
        self.assertRaises(Exception, a.__ne__, b)
        self.assertRaises(Exception, a.__lt__, b)
        self.assertRaises(Exception, a.__gt__, b)
        self.assertRaises(Exception, b.__eq__, a)
        self.assertRaises(Exception, b.__ne__, a)
        self.assertRaises(Exception, b.__lt__, a)
        self.assertRaises(Exception, b.__gt__, a)

        if sys.version_info < (3, 3):
            self.assertRaises(Exception, a.__eq__, b.to_pydatetime())
            self.assertRaises(Exception, a.to_pydatetime().__eq__, b)
        else:
            self.assertFalse(a == b.to_pydatetime())
            self.assertFalse(a.to_pydatetime() == b)

    def test_cant_compare_tz_naive_w_aware_dateutil(self):
        tm._skip_if_no_dateutil()
        from dateutil.tz import tzutc
        utc = tzutc()
        # #1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz=utc)

        self.assertRaises(Exception, a.__eq__, b)
        self.assertRaises(Exception, a.__ne__, b)
        self.assertRaises(Exception, a.__lt__, b)
        self.assertRaises(Exception, a.__gt__, b)
        self.assertRaises(Exception, b.__eq__, a)
        self.assertRaises(Exception, b.__ne__, a)
        self.assertRaises(Exception, b.__lt__, a)
        self.assertRaises(Exception, b.__gt__, a)

        if sys.version_info < (3, 3):
            self.assertRaises(Exception, a.__eq__, b.to_pydatetime())
            self.assertRaises(Exception, a.to_pydatetime().__eq__, b)
        else:
            self.assertFalse(a == b.to_pydatetime())
            self.assertFalse(a.to_pydatetime() == b)

    def test_delta_preserve_nanos(self):
        val = Timestamp(long(1337299200000000123))
        result = val + timedelta(1)
        self.assertEqual(result.nanosecond, val.nanosecond)

    def test_frequency_misc(self):
        self.assertEqual(frequencies.get_freq_group('T'),
                         frequencies.FreqGroup.FR_MIN)

        code, stride = frequencies.get_freq_code(offsets.Hour())
        self.assertEqual(code, frequencies.FreqGroup.FR_HR)

        code, stride = frequencies.get_freq_code((5, 'T'))
        self.assertEqual(code, frequencies.FreqGroup.FR_MIN)
        self.assertEqual(stride, 5)

        offset = offsets.Hour()
        result = frequencies.to_offset(offset)
        self.assertEqual(result, offset)

        result = frequencies.to_offset((5, 'T'))
        expected = offsets.Minute(5)
        self.assertEqual(result, expected)

        self.assertRaises(ValueError, frequencies.get_freq_code, (5, 'baz'))

        self.assertRaises(ValueError, frequencies.to_offset, '100foo')

        self.assertRaises(ValueError, frequencies.to_offset, ('', ''))

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = frequencies.get_standard_freq(offsets.Hour())
        self.assertEqual(result, 'H')

    def test_hash_equivalent(self):
        d = {datetime(2011, 1, 1): 5}
        stamp = Timestamp(datetime(2011, 1, 1))
        self.assertEqual(d[stamp], 5)

    def test_timestamp_compare_scalars(self):
        # case where ndim == 0
        lhs = np.datetime64(datetime(2013, 12, 6))
        rhs = Timestamp('now')
        nat = Timestamp('nat')

        ops = {'gt': 'lt',
               'lt': 'gt',
               'ge': 'le',
               'le': 'ge',
               'eq': 'eq',
               'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)
            expected = left_f(lhs, rhs)

            result = right_f(rhs, lhs)
            self.assertEqual(result, expected)

            expected = left_f(rhs, nat)
            result = right_f(nat, rhs)
            self.assertEqual(result, expected)

    def test_timestamp_compare_series(self):
        # make sure we can compare Timestamps on the right AND left hand side
        # GH4982
        s = Series(date_range('20010101', periods=10), name='dates')
        s_nat = s.copy(deep=True)

        s[0] = pd.Timestamp('nat')
        s[3] = pd.Timestamp('nat')

        ops = {'lt': 'gt', 'le': 'ge', 'eq': 'eq', 'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)

            # no nats
            expected = left_f(s, Timestamp('20010109'))
            result = right_f(Timestamp('20010109'), s)
            tm.assert_series_equal(result, expected)

            # nats
            expected = left_f(s, Timestamp('nat'))
            result = right_f(Timestamp('nat'), s)
            tm.assert_series_equal(result, expected)

            # compare to timestamp with series containing nats
            expected = left_f(s_nat, Timestamp('20010109'))
            result = right_f(Timestamp('20010109'), s_nat)
            tm.assert_series_equal(result, expected)

            # compare to nat with series containing nats
            expected = left_f(s_nat, Timestamp('nat'))
            result = right_f(Timestamp('nat'), s_nat)
            tm.assert_series_equal(result, expected)

    def test_is_leap_year(self):
        # GH 13727
        for tz in [None, 'UTC', 'US/Eastern', 'Asia/Tokyo']:
            dt = Timestamp('2000-01-01 00:00:00', tz=tz)
            self.assertTrue(dt.is_leap_year)
            self.assertIsInstance(dt.is_leap_year, bool)

            dt = Timestamp('1999-01-01 00:00:00', tz=tz)
            self.assertFalse(dt.is_leap_year)

            dt = Timestamp('2004-01-01 00:00:00', tz=tz)
            self.assertTrue(dt.is_leap_year)

            dt = Timestamp('2100-01-01 00:00:00', tz=tz)
            self.assertFalse(dt.is_leap_year)

        self.assertFalse(pd.NaT.is_leap_year)
        self.assertIsInstance(pd.NaT.is_leap_year, bool)

    def test_round_nat(self):
        # GH14940
        ts = Timestamp('nat')
        print(dir(ts))
        for method in ["round", "floor", "ceil"]:
            round_method = getattr(ts, method)
            for freq in ["s", "5s", "min", "5min", "h", "5h"]:
                self.assertIs(round_method(freq), ts)


class TestTimestampNsOperations(tm.TestCase):

    def setUp(self):
        self.timestamp = Timestamp(datetime.utcnow())

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
        dt = datetime(2014, 1, 1)
        for (left, right) in [(pd.NaT, pd.NaT), (pd.NaT, t), (pd.NaT, dt)]:
            # NaT __add__ or __sub__ Timestamp-like (or inverse) returns NaT
            self.assertIs(right + left, pd.NaT)
            self.assertIs(left + right, pd.NaT)
            self.assertIs(left - right, pd.NaT)
            self.assertIs(right - left, pd.NaT)

        # timedelta-like
        # offsets are tested in test_offsets.py

        delta = timedelta(3600)
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

        # int addition / subtraction
        for (left, right) in [(pd.NaT, 2), (pd.NaT, 0), (pd.NaT, -3)]:
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


class TestTimestampOps(tm.TestCase):

    def test_timestamp_and_datetime(self):
        self.assertEqual((Timestamp(datetime(
            2013, 10, 13)) - datetime(2013, 10, 12)).days, 1)
        self.assertEqual((datetime(2013, 10, 12) -
                          Timestamp(datetime(2013, 10, 13))).days, -1)

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
        datetime_instance = datetime(2014, 3, 4)
        timedelta_instance = timedelta(seconds=1)
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
        timedelta_instance = timedelta(days=1)
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
                                  [RESO_DAY, RESO_DAY,
                                   RESO_DAY, RESO_DAY,
                                   RESO_HR, RESO_MIN,
                                   RESO_SEC, RESO_MS,
                                   RESO_US]):
            for tz in [None, 'Asia/Tokyo', 'US/Eastern',
                       'dateutil/US/Eastern']:
                idx = date_range(start='2013-04-01', periods=30, freq=freq,
                                 tz=tz)
                result = period.resolution(idx.asi8, idx.tz)
                self.assertEqual(result, expected)


class TestTimestampToJulianDate(tm.TestCase):

    def test_compare_1700(self):
        r = Timestamp('1700-06-23').to_julian_date()
        self.assertEqual(r, 2342145.5)

    def test_compare_2000(self):
        r = Timestamp('2000-04-12').to_julian_date()
        self.assertEqual(r, 2451646.5)

    def test_compare_2100(self):
        r = Timestamp('2100-08-12').to_julian_date()
        self.assertEqual(r, 2488292.5)

    def test_compare_hour01(self):
        r = Timestamp('2000-08-12T01:00:00').to_julian_date()
        self.assertEqual(r, 2451768.5416666666666666)

    def test_compare_hour13(self):
        r = Timestamp('2000-08-12T13:00:00').to_julian_date()
        self.assertEqual(r, 2451769.0416666666666666)


class TestTimeSeries(tm.TestCase):

    def test_timestamp_to_datetime(self):
        tm._skip_if_no_pytz()
        rng = date_range('20090415', '20090519', tz='US/Eastern')

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_timestamp_to_datetime_dateutil(self):
        tm._skip_if_no_pytz()
        rng = date_range('20090415', '20090519', tz='dateutil/US/Eastern')

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_timestamp_to_datetime_explicit_pytz(self):
        tm._skip_if_no_pytz()
        import pytz
        rng = date_range('20090415', '20090519',
                         tz=pytz.timezone('US/Eastern'))

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_timestamp_to_datetime_explicit_dateutil(self):
        tm._skip_if_windows_python_3()
        tm._skip_if_no_dateutil()
        from pandas.tslib import _dateutil_gettz as gettz
        rng = date_range('20090415', '20090519', tz=gettz('US/Eastern'))

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_timestamp_fields(self):
        # extra fields from DatetimeIndex like quarter and week
        idx = tm.makeDateIndex(100)

        fields = ['dayofweek', 'dayofyear', 'week', 'weekofyear', 'quarter',
                  'days_in_month', 'is_month_start', 'is_month_end',
                  'is_quarter_start', 'is_quarter_end', 'is_year_start',
                  'is_year_end', 'weekday_name']
        for f in fields:
            expected = getattr(idx, f)[-1]
            result = getattr(Timestamp(idx[-1]), f)
            self.assertEqual(result, expected)

        self.assertEqual(idx.freq, Timestamp(idx[-1], idx.freq).freq)
        self.assertEqual(idx.freqstr, Timestamp(idx[-1], idx.freq).freqstr)

    def test_timestamp_date_out_of_range(self):
        self.assertRaises(ValueError, Timestamp, '1676-01-01')
        self.assertRaises(ValueError, Timestamp, '2263-01-01')

        # 1475
        self.assertRaises(ValueError, DatetimeIndex, ['1400-01-01'])
        self.assertRaises(ValueError, DatetimeIndex, [datetime(1400, 1, 1)])

    def test_timestamp_repr(self):
        # pre-1900
        stamp = Timestamp('1850-01-01', tz='US/Eastern')
        repr(stamp)

        iso8601 = '1850-01-01 01:23:45.012345'
        stamp = Timestamp(iso8601, tz='US/Eastern')
        result = repr(stamp)
        self.assertIn(iso8601, result)

    def test_timestamp_from_ordinal(self):

        # GH 3042
        dt = datetime(2011, 4, 16, 0, 0)
        ts = Timestamp.fromordinal(dt.toordinal())
        self.assertEqual(ts.to_pydatetime(), dt)

        # with a tzinfo
        stamp = Timestamp('2011-4-16', tz='US/Eastern')
        dt_tz = stamp.to_pydatetime()
        ts = Timestamp.fromordinal(dt_tz.toordinal(), tz='US/Eastern')
        self.assertEqual(ts.to_pydatetime(), dt_tz)

    def test_timestamp_compare_with_early_datetime(self):
        # e.g. datetime.min
        stamp = Timestamp('2012-01-01')

        self.assertFalse(stamp == datetime.min)
        self.assertFalse(stamp == datetime(1600, 1, 1))
        self.assertFalse(stamp == datetime(2700, 1, 1))
        self.assertNotEqual(stamp, datetime.min)
        self.assertNotEqual(stamp, datetime(1600, 1, 1))
        self.assertNotEqual(stamp, datetime(2700, 1, 1))
        self.assertTrue(stamp > datetime(1600, 1, 1))
        self.assertTrue(stamp >= datetime(1600, 1, 1))
        self.assertTrue(stamp < datetime(2700, 1, 1))
        self.assertTrue(stamp <= datetime(2700, 1, 1))

    def test_timestamp_equality(self):

        # GH 11034
        s = Series([Timestamp('2000-01-29 01:59:00'), 'NaT'])
        result = s != s
        assert_series_equal(result, Series([False, True]))
        result = s != s[0]
        assert_series_equal(result, Series([False, True]))
        result = s != s[1]
        assert_series_equal(result, Series([True, True]))

        result = s == s
        assert_series_equal(result, Series([True, False]))
        result = s == s[0]
        assert_series_equal(result, Series([True, False]))
        result = s == s[1]
        assert_series_equal(result, Series([False, False]))

    def test_series_box_timestamp(self):
        rng = date_range('20090415', '20090519', freq='B')
        s = Series(rng)

        tm.assertIsInstance(s[5], Timestamp)

        rng = date_range('20090415', '20090519', freq='B')
        s = Series(rng, index=rng)
        tm.assertIsInstance(s[5], Timestamp)

        tm.assertIsInstance(s.iat[5], Timestamp)

    def test_frame_setitem_timestamp(self):
        # 2155
        columns = DatetimeIndex(start='1/1/2012', end='2/1/2012',
                                freq=offsets.BDay())
        index = lrange(10)
        data = DataFrame(columns=columns, index=index)
        t = datetime(2012, 11, 1)
        ts = Timestamp(t)
        data[ts] = np.nan  # works

    def test_to_html_timestamp(self):
        rng = date_range('2000-01-01', periods=10)
        df = DataFrame(np.random.randn(10, 4), index=rng)

        result = df.to_html()
        self.assertIn('2000-01-01', result)

    def test_series_map_box_timestamps(self):
        # #2689, #2627
        s = Series(date_range('1/1/2000', periods=10))

        def f(x):
            return (x.hour, x.day, x.month)

        # it works!
        s.map(f)
        s.apply(f)
        DataFrame(s).applymap(f)

    def test_dti_slicing(self):
        dti = DatetimeIndex(start='1/1/2005', end='12/1/2005', freq='M')
        dti2 = dti[[1, 3, 5]]

        v1 = dti2[0]
        v2 = dti2[1]
        v3 = dti2[2]

        self.assertEqual(v1, Timestamp('2/28/2005'))
        self.assertEqual(v2, Timestamp('4/30/2005'))
        self.assertEqual(v3, Timestamp('6/30/2005'))

        # don't carry freq through irregular slicing
        self.assertIsNone(dti2.freq)

    def test_woy_boundary(self):
        # make sure weeks at year boundaries are correct
        d = datetime(2013, 12, 31)
        result = Timestamp(d).week
        expected = 1  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2008, 12, 28)
        result = Timestamp(d).week
        expected = 52  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2009, 12, 31)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2010, 1, 1)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2010, 1, 3)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        self.assertEqual(result, expected)

        result = np.array([Timestamp(datetime(*args)).week
                           for args in [(2000, 1, 1), (2000, 1, 2), (
                               2005, 1, 1), (2005, 1, 2)]])
        self.assertTrue((result == [52, 52, 53, 53]).all())


class TestTsUtil(tm.TestCase):

    def test_min_valid(self):
        # Ensure that Timestamp.min is a valid Timestamp
        Timestamp(Timestamp.min)

    def test_max_valid(self):
        # Ensure that Timestamp.max is a valid Timestamp
        Timestamp(Timestamp.max)

    def test_to_datetime_bijective(self):
        # Ensure that converting to datetime and back only loses precision
        # by going from nanoseconds to microseconds.
        exp_warning = None if Timestamp.max.nanosecond == 0 else UserWarning
        with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
            self.assertEqual(
                Timestamp(Timestamp.max.to_pydatetime()).value / 1000,
                Timestamp.max.value / 1000)

        exp_warning = None if Timestamp.min.nanosecond == 0 else UserWarning
        with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
            self.assertEqual(
                Timestamp(Timestamp.min.to_pydatetime()).value / 1000,
                Timestamp.min.value / 1000)


class TestTslib(tm.TestCase):

    def test_round(self):
        stamp = Timestamp('2000-01-05 05:09:15.13')

        def _check_round(freq, expected):
            result = stamp.round(freq=freq)
            self.assertEqual(result, expected)

        for freq, expected in [('D', Timestamp('2000-01-05 00:00:00')),
                               ('H', Timestamp('2000-01-05 05:00:00')),
                               ('S', Timestamp('2000-01-05 05:09:15'))]:
            _check_round(freq, expected)

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with self.assertRaisesRegexp(ValueError, msg):
            stamp.round('foo')
