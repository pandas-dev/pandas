import nose
from distutils.version import LooseVersion
import numpy as np

from pandas import tslib
import datetime

from pandas.core.api import Timestamp, Series, Timedelta
from pandas.tslib import period_asfreq, period_ordinal, get_timezone
from pandas.tseries.index import date_range
from pandas.tseries.frequencies import get_freq
import pandas.tseries.offsets as offsets
import pandas.util.testing as tm
from pandas.util.testing import assert_series_equal

class TestTimestamp(tm.TestCase):

    def test_constructor(self):
        base_str = '2014-07-01 09:00'
        base_dt = datetime.datetime(2014, 7, 1, 9)
        base_expected = 1404205200000000000

        # confirm base representation is correct
        import calendar
        self.assertEqual(calendar.timegm(base_dt.timetuple()) * 1000000000, base_expected)

        tests = [(base_str, base_dt, base_expected),
                 ('2014-07-01 10:00', datetime.datetime(2014, 7, 1, 10),
                  base_expected + 3600 * 1000000000),
                 ('2014-07-01 09:00:00.000008000',
                  datetime.datetime(2014, 7, 1, 9, 0, 0, 8), base_expected + 8000),
                 ('2014-07-01 09:00:00.000000005',
                  Timestamp('2014-07-01 09:00:00.000000005'), base_expected + 5)]

        tm._skip_if_no_pytz()
        tm._skip_if_no_dateutil()
        import pytz
        import dateutil
        timezones = [(None, 0), ('UTC', 0), (pytz.utc, 0),
                     ('Asia/Tokyo', 9), ('US/Eastern', -4), ('dateutil/US/Pacific', -7),
                     (pytz.FixedOffset(-180), -3), (dateutil.tz.tzoffset(None, 18000), 5)]

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
                for result in [Timestamp(date_str, tz=tz), Timestamp(date, tz=tz)]:
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
        self.assertEqual(calendar.timegm(base_dt.timetuple()) * 1000000000, base_expected)

        tests = [(base_str, base_expected),
                 ('2014-07-01 12:00:00+02:00', base_expected + 3600 * 1000000000),
                 ('2014-07-01 11:00:00.000008000+02:00', base_expected + 8000),
                 ('2014-07-01 11:00:00.000000005+02:00', base_expected + 5)]

        tm._skip_if_no_pytz()
        tm._skip_if_no_dateutil()
        import pytz
        import dateutil
        timezones = [(None, 0), ('UTC', 0), (pytz.utc, 0),
                     ('Asia/Tokyo', 9), ('US/Eastern', -4),
                     ('dateutil/US/Pacific', -7),
                     (pytz.FixedOffset(-180), -3), (dateutil.tz.tzoffset(None, 18000), 5)]

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

        # This should be 2013-11-01 05:00 in UTC -> converted to Chicago tz
        result = Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')
        self.assertEqual(result.value, Timestamp('2013-11-01 05:00').value)
        expected_repr = "Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')"
        self.assertEqual(repr(result), expected_repr)
        self.assertEqual(result, eval(repr(result)))

        # This should be 2013-11-01 05:00 in UTC -> converted to Tokyo tz (+09:00)
        result = Timestamp('2013-11-01 00:00:00-0500', tz='Asia/Tokyo')
        self.assertEqual(result.value, Timestamp('2013-11-01 05:00').value)
        expected_repr = "Timestamp('2013-11-01 14:00:00+0900', tz='Asia/Tokyo')"
        self.assertEqual(repr(result), expected_repr)
        self.assertEqual(result, eval(repr(result)))

    def test_repr(self):
        tm._skip_if_no_pytz()
        tm._skip_if_no_dateutil()

        dates = ['2014-03-07', '2014-01-01 09:00', '2014-01-01 00:00:00.000000001']

        # dateutil zone change (only matters for repr)
        import dateutil
        if dateutil.__version__ >= LooseVersion('2.3'):
            timezones = ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']
        else:
            timezones = ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/America/Los_Angeles']

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

        # this can cause the tz field to be populated, but it's redundant to information in the datestring
        tm._skip_if_no_pytz()
        import pytz
        date_with_utc_offset = Timestamp('2014-03-13 00:00:00-0400', tz=None)
        self.assertIn('2014-03-13 00:00:00-0400', repr(date_with_utc_offset))
        self.assertNotIn('tzoffset', repr(date_with_utc_offset))
        self.assertIn('pytz.FixedOffset(-240)', repr(date_with_utc_offset))
        expr = repr(date_with_utc_offset).replace("'pytz.FixedOffset(-240)'",
                    'pytz.FixedOffset(-240)')
        self.assertEqual(date_with_utc_offset, eval(expr))

    def test_bounds_with_different_units(self):
        out_of_bounds_dates = (
            '1677-09-21',
            '2262-04-12',
        )

        time_units = ('D', 'h', 'm', 's', 'ms', 'us')

        for date_string in out_of_bounds_dates:
            for unit in time_units:
                self.assertRaises(
                    ValueError,
                    Timestamp,
                    np.datetime64(date_string, dtype='M8[%s]' % unit)
                )

        in_bounds_dates = (
            '1677-09-23',
            '2262-04-11',
        )

        for date_string in in_bounds_dates:
            for unit in time_units:
                Timestamp(
                    np.datetime64(date_string, dtype='M8[%s]' % unit)
                )

    def test_tz(self):
        t = '2014-02-01 09:00'
        ts = Timestamp(t)
        local = ts.tz_localize('Asia/Tokyo')
        self.assertEqual(local.hour, 9)
        self.assertEqual(local, Timestamp(t, tz='Asia/Tokyo'))
        conv = local.tz_convert('US/Eastern')
        self.assertEqual(conv,
                         Timestamp('2014-01-31 19:00', tz='US/Eastern'))
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
        self.assertRaises(ValueError, ts.tz_localize, 'US/Eastern', ambiguous='infer')

        # GH 8025
        with tm.assertRaisesRegexp(TypeError, 'Cannot localize tz-aware Timestamp, use '
                                   'tz_convert for conversions'):
            Timestamp('2011-01-01' ,tz='US/Eastern').tz_localize('Asia/Tokyo')

        with tm.assertRaisesRegexp(TypeError, 'Cannot convert tz-naive Timestamp, use '
                            'tz_localize to localize'):
            Timestamp('2011-01-01').tz_convert('Asia/Tokyo')

    def test_tz_localize_roundtrip(self):
        for tz in ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']:
            for t in ['2014-02-01 09:00', '2014-07-08 09:00', '2014-11-01 17:00',
                      '2014-11-05 00:00']:
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
            for t in ['2014-02-01 09:00', '2014-07-08 09:00', '2014-11-01 17:00',
                      '2014-11-05 00:00']:
                ts = Timestamp(t, tz='UTC')
                converted = ts.tz_convert(tz)

                reset = converted.tz_convert(None)
                self.assertEqual(reset, Timestamp(t))
                self.assertTrue(reset.tzinfo is None)
                self.assertEqual(reset, converted.tz_convert('UTC').tz_localize(None))

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
        self.assertEqual(get_timezone(Timestamp('2014-11-02 01:00Z').tzinfo), 'UTC')

    def test_now(self):
        # #9000
        ts_from_string = Timestamp('now')
        ts_from_method = Timestamp.now()
        ts_datetime = datetime.datetime.now()

        ts_from_string_tz = Timestamp('now', tz='US/Eastern')
        ts_from_method_tz = Timestamp.now(tz='US/Eastern')

        # Check that the delta between the times is less than 1s (arbitrarily small)
        delta = Timedelta(seconds=1)
        self.assertTrue((ts_from_method - ts_from_string) < delta)
        self.assertTrue((ts_from_method_tz - ts_from_string_tz) < delta)
        self.assertTrue((ts_from_string_tz.tz_localize(None) - ts_from_string) < delta)

    def test_today(self):

        ts_from_string = Timestamp('today')
        ts_from_method = Timestamp.today()
        ts_datetime = datetime.datetime.today()

        ts_from_string_tz = Timestamp('today', tz='US/Eastern')
        ts_from_method_tz = Timestamp.today(tz='US/Eastern')

        # Check that the delta between the times is less than 1s (arbitrarily small)
        delta = Timedelta(seconds=1)
        self.assertTrue((ts_from_method - ts_from_string) < delta)
        self.assertTrue((ts_datetime - ts_from_method) < delta)
        self.assertTrue((ts_datetime - ts_from_method) < delta)
        self.assertTrue((ts_from_string_tz.tz_localize(None) - ts_from_string) < delta)

class TestDatetimeParsingWrappers(tm.TestCase):
    def test_does_not_convert_mixed_integer(self):
        bad_date_strings = (
            '-50000',
            '999',
            '123.1234',
            'm',
            'T'
        )

        for bad_date_string in bad_date_strings:
            self.assertFalse(
                tslib._does_string_look_like_datetime(bad_date_string)
            )

        good_date_strings = (
            '2012-01-01',
            '01/01/2012',
            'Mon Sep 16, 2013',
            '01012012',
            '0101',
            '1-1',
        )

        for good_date_string in good_date_strings:
            self.assertTrue(
                tslib._does_string_look_like_datetime(good_date_string)
            )


class TestArrayToDatetime(tm.TestCase):
    def test_parsing_valid_dates(self):
        arr = np.array(['01-01-2013', '01-02-2013'], dtype=object)
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr),
            np.array(
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
            np.array(
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
        self.assert_numpy_array_equal(tslib.array_to_datetime(arr), arr)

        arr = np.array(['1', '2', '3', '4', '5'], dtype=object)
        self.assert_numpy_array_equal(tslib.array_to_datetime(arr), arr)

    def test_coercing_dates_outside_of_datetime64_ns_bounds(self):
        invalid_dates = [
            datetime.date(1000, 1, 1),
            datetime.datetime(1000, 1, 1),
            '1000-01-01',
            'Jan 1, 1000',
            np.datetime64('1000-01-01'),
        ]

        for invalid_date in invalid_dates:
            self.assertRaises(
                ValueError,
                tslib.array_to_datetime,
                np.array([invalid_date], dtype='object'),
                coerce=False,
                raise_=True,
            )
            self.assertTrue(
                np.array_equal(
                    tslib.array_to_datetime(
                        np.array([invalid_date], dtype='object'), coerce=True
                    ),
                    np.array([tslib.iNaT], dtype='M8[ns]')
                )
            )

        arr = np.array(['1/1/1000', '1/1/2000'], dtype=object)
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr, coerce=True),
            np.array(
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
        self.assert_numpy_array_equal(tslib.array_to_datetime(arr), arr)

        # With coercing, the invalid dates becomes iNaT
        self.assert_numpy_array_equal(
            tslib.array_to_datetime(arr, coerce=True),
            np.array(
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
            '12-31-2012 23:00:00-01:00',
        ]

        expected_output = tslib.array_to_datetime(
            np.array(['01-01-2013 00:00:00'], dtype=object)
        )

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
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'ns'), -123)

    def test_timedelta_ns_based_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(1234567898, 'ns'), 1234567898)

    def test_timedelta_us_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'us'), -123000)

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

        t = Timestamp(np.datetime64('2011-01-01 00:00:00.000000005Z'))
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

        t = Timestamp(np.datetime64('2011-01-01 00:00:00.000000010Z'))
        self.assertEqual(repr(t), "Timestamp('2011-01-01 00:00:00.000000010')")
        self.assertEqual(t.value, expected)
        self.assertEqual(t.nanosecond, 10)

    def test_nat_arithmetic(self):
        # GH 6873
        nat = tslib.NaT
        t = Timestamp('2014-01-01')
        dt = datetime.datetime(2014, 1, 1)
        delta = datetime.timedelta(3600)

        # Timestamp / datetime
        for (left, right) in [(nat, nat), (nat, t), (dt, nat)]:
            # NaT + Timestamp-like should raise TypeError
            with tm.assertRaises(TypeError):
                left + right
            with tm.assertRaises(TypeError):
                right + left

            # NaT - Timestamp-like (or inverse) returns NaT
            self.assertTrue((left - right) is tslib.NaT)
            self.assertTrue((right - left) is tslib.NaT)

        # timedelta-like
        # offsets are tested in test_offsets.py
        for (left, right) in [(nat, delta)]:
            # NaT + timedelta-like returns NaT
            self.assertTrue((left + right) is tslib.NaT)
            # timedelta-like + NaT should raise TypeError
            with tm.assertRaises(TypeError):
                right + left

            self.assertTrue((left - right) is tslib.NaT)
            with tm.assertRaises(TypeError):
                right - left


class TestTslib(tm.TestCase):

    def test_intraday_conversion_factors(self):
        self.assertEqual(period_asfreq(1, get_freq('D'), get_freq('H'), False), 24)
        self.assertEqual(period_asfreq(1, get_freq('D'), get_freq('T'), False), 1440)
        self.assertEqual(period_asfreq(1, get_freq('D'), get_freq('S'), False), 86400)
        self.assertEqual(period_asfreq(1, get_freq('D'), get_freq('L'), False), 86400000)
        self.assertEqual(period_asfreq(1, get_freq('D'), get_freq('U'), False), 86400000000)
        self.assertEqual(period_asfreq(1, get_freq('D'), get_freq('N'), False), 86400000000000)

        self.assertEqual(period_asfreq(1, get_freq('H'), get_freq('T'), False), 60)
        self.assertEqual(period_asfreq(1, get_freq('H'), get_freq('S'), False), 3600)
        self.assertEqual(period_asfreq(1, get_freq('H'), get_freq('L'), False), 3600000)
        self.assertEqual(period_asfreq(1, get_freq('H'), get_freq('U'), False), 3600000000)
        self.assertEqual(period_asfreq(1, get_freq('H'), get_freq('N'), False), 3600000000000)

        self.assertEqual(period_asfreq(1, get_freq('T'), get_freq('S'), False), 60)
        self.assertEqual(period_asfreq(1, get_freq('T'), get_freq('L'), False), 60000)
        self.assertEqual(period_asfreq(1, get_freq('T'), get_freq('U'), False), 60000000)
        self.assertEqual(period_asfreq(1, get_freq('T'), get_freq('N'), False), 60000000000)

        self.assertEqual(period_asfreq(1, get_freq('S'), get_freq('L'), False), 1000)
        self.assertEqual(period_asfreq(1, get_freq('S'), get_freq('U'), False), 1000000)
        self.assertEqual(period_asfreq(1, get_freq('S'), get_freq('N'), False), 1000000000)

        self.assertEqual(period_asfreq(1, get_freq('L'), get_freq('U'), False), 1000)
        self.assertEqual(period_asfreq(1, get_freq('L'), get_freq('N'), False), 1000000)

        self.assertEqual(period_asfreq(1, get_freq('U'), get_freq('N'), False), 1000)

    def test_period_ordinal_start_values(self):
        # information for 1.1.1970
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('Y')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('M')))
        self.assertEqual(1, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('W')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('D')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('B')))

    def test_period_ordinal_week(self):
        self.assertEqual(1, period_ordinal(1970, 1, 4, 0, 0, 0, 0, 0, get_freq('W')))
        self.assertEqual(2, period_ordinal(1970, 1, 5, 0, 0, 0, 0, 0, get_freq('W')))

        self.assertEqual(2284, period_ordinal(2013, 10, 6, 0, 0, 0, 0, 0, get_freq('W')))
        self.assertEqual(2285, period_ordinal(2013, 10, 7, 0, 0, 0, 0, 0, get_freq('W')))

    def test_period_ordinal_business_day(self):
        # Thursday
        self.assertEqual(11415, period_ordinal(2013, 10, 3, 0, 0, 0, 0, 0, get_freq('B')))
        # Friday
        self.assertEqual(11416, period_ordinal(2013, 10, 4, 0, 0, 0, 0, 0, get_freq('B')))
        # Saturday
        self.assertEqual(11417, period_ordinal(2013, 10, 5, 0, 0, 0, 0, 0, get_freq('B')))
        # Sunday
        self.assertEqual(11417, period_ordinal(2013, 10, 6, 0, 0, 0, 0, 0, get_freq('B')))
        # Monday
        self.assertEqual(11417, period_ordinal(2013, 10, 7, 0, 0, 0, 0, 0, get_freq('B')))
        # Tuesday
        self.assertEqual(11418, period_ordinal(2013, 10, 8, 0, 0, 0, 0, 0, get_freq('B')))

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
            # local tz to UTC can be differ in hourly (or higher) freqs because of DST
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

class TestTimestampOps(tm.TestCase):
    def test_timestamp_and_datetime(self):
        self.assertEqual((Timestamp(datetime.datetime(2013, 10, 13)) - datetime.datetime(2013, 10, 12)).days, 1)
        self.assertEqual((datetime.datetime(2013, 10, 12) - Timestamp(datetime.datetime(2013, 10, 13))).days, -1)

    def test_timestamp_and_series(self):
        timestamp_series = Series(date_range('2014-03-17', periods=2, freq='D', tz='US/Eastern'))
        first_timestamp = timestamp_series[0]

        delta_series = Series([np.timedelta64(0, 'D'), np.timedelta64(1, 'D')])
        assert_series_equal(timestamp_series - first_timestamp, delta_series)
        assert_series_equal(first_timestamp - timestamp_series, -delta_series)

    def test_addition_subtraction_types(self):
        # Assert on the types resulting from Timestamp +/- various date/time objects
        datetime_instance = datetime.datetime(2014, 3, 4)
        timedelta_instance = datetime.timedelta(seconds=1)
        # build a timestamp with a frequency, since then it supports addition/subtraction of integers
        timestamp_instance = date_range(datetime_instance, periods=1, freq='D')[0]

        self.assertEqual(type(timestamp_instance + 1), Timestamp)
        self.assertEqual(type(timestamp_instance - 1), Timestamp)

        # Timestamp + datetime not supported, though subtraction is supported and yields timedelta
        # more tests in tseries/base/tests/test_base.py
        self.assertEqual(type(timestamp_instance - datetime_instance), Timedelta)
        self.assertEqual(type(timestamp_instance + timedelta_instance), Timestamp)
        self.assertEqual(type(timestamp_instance - timedelta_instance), Timestamp)

        # Timestamp +/- datetime64 not supported, so not tested (could possibly assert error raised?)
        timedelta64_instance = np.timedelta64(1, 'D')
        self.assertEqual(type(timestamp_instance + timedelta64_instance), Timestamp)
        self.assertEqual(type(timestamp_instance - timedelta64_instance), Timestamp)

    def test_addition_subtraction_preserve_frequency(self):
        timestamp_instance = date_range('2014-03-05', periods=1, freq='D')[0]
        timedelta_instance = datetime.timedelta(days=1)
        original_freq = timestamp_instance.freq
        self.assertEqual((timestamp_instance + 1).freq, original_freq)
        self.assertEqual((timestamp_instance - 1).freq, original_freq)
        self.assertEqual((timestamp_instance + timedelta_instance).freq, original_freq)
        self.assertEqual((timestamp_instance - timedelta_instance).freq, original_freq)

        timedelta64_instance = np.timedelta64(1, 'D')
        self.assertEqual((timestamp_instance + timedelta64_instance).freq, original_freq)
        self.assertEqual((timestamp_instance - timedelta64_instance).freq, original_freq)

    def test_resolution(self):

        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H', 'T', 'S', 'L', 'U'],
                                  [tslib.D_RESO, tslib.D_RESO, tslib.D_RESO, tslib.D_RESO,
                                   tslib.H_RESO, tslib.T_RESO,tslib.S_RESO, tslib.MS_RESO, tslib.US_RESO]):
            for tz in [None, 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Eastern']:
                idx = date_range(start='2013-04-01', periods=30, freq=freq, tz=tz)
                result = tslib.resolution(idx.asi8, idx.tz)
                self.assertEqual(result, expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
