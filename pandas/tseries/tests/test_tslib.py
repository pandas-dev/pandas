import unittest
import nose

import numpy as np

from pandas import tslib
import datetime

from pandas.core.api import Timestamp

from pandas.tslib import period_asfreq

from pandas.tseries.frequencies import get_freq

from pandas import _np_version_under1p7


class TestTimestamp(unittest.TestCase):
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
                    tslib.Timestamp,
                    np.datetime64(date_string, dtype='M8[%s]' % unit)
                )

        in_bounds_dates = (
            '1677-09-23',
            '2262-04-11',
        )

        for date_string in in_bounds_dates:
            for unit in time_units:
                tslib.Timestamp(
                    np.datetime64(date_string, dtype='M8[%s]' % unit)
                )

    def test_barely_oob_dts(self):
        one_us = np.timedelta64(1)

        # By definition we can't go out of bounds in [ns], so we
        # convert the datetime64s to [us] so we can go out of bounds
        min_ts_us = np.datetime64(tslib.Timestamp.min).astype('M8[us]')
        max_ts_us = np.datetime64(tslib.Timestamp.max).astype('M8[us]')

        # No error for the min/max datetimes
        tslib.Timestamp(min_ts_us)
        tslib.Timestamp(max_ts_us)

        # One us less than the minimum is an error
        self.assertRaises(ValueError, tslib.Timestamp, min_ts_us - one_us)

        # One us more than the maximum is an error
        self.assertRaises(ValueError, tslib.Timestamp, max_ts_us + one_us)

class TestDatetimeParsingWrappers(unittest.TestCase):
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


class TestArrayToDatetime(unittest.TestCase):
    def test_parsing_valid_dates(self):
        arr = np.array(['01-01-2013', '01-02-2013'], dtype=object)
        self.assert_(
            np.array_equal(
                tslib.array_to_datetime(arr),
                np.array(
                    [
                        '2013-01-01T00:00:00.000000000-0000',
                        '2013-01-02T00:00:00.000000000-0000'
                    ],
                    dtype='M8[ns]'
                )
            )
        )

        arr = np.array(['Mon Sep 16 2013', 'Tue Sep 17 2013'], dtype=object)
        self.assert_(
            np.array_equal(
                tslib.array_to_datetime(arr),
                np.array(
                    [
                        '2013-09-16T00:00:00.000000000-0000',
                        '2013-09-17T00:00:00.000000000-0000'
                    ],
                    dtype='M8[ns]'
                )
            )
        )

    def test_number_looking_strings_not_into_datetime(self):
        # #4601
        # These strings don't look like datetimes so they shouldn't be
        # attempted to be converted
        arr = np.array(['-352.737091', '183.575577'], dtype=object)
        self.assert_(np.array_equal(tslib.array_to_datetime(arr), arr))

        arr = np.array(['1', '2', '3', '4', '5'], dtype=object)
        self.assert_(np.array_equal(tslib.array_to_datetime(arr), arr))

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
            self.assert_(
                np.array_equal(
                    tslib.array_to_datetime(
                        np.array([invalid_date], dtype='object'), coerce=True
                    ),
                    np.array([tslib.iNaT], dtype='M8[ns]')
                )
            )

        arr = np.array(['1/1/1000', '1/1/2000'], dtype=object)
        self.assert_(
            np.array_equal(
                tslib.array_to_datetime(arr, coerce=True),
                np.array(
                    [
                        tslib.iNaT,
                        '2000-01-01T00:00:00.000000000-0000'
                    ],
                    dtype='M8[ns]'
                )
            )
        )

    def test_coerce_of_invalid_datetimes(self):
        arr = np.array(['01-01-2013', 'not_a_date', '1'], dtype=object)

        # Without coercing, the presence of any invalid dates prevents
        # any values from being converted
        self.assert_(np.array_equal(tslib.array_to_datetime(arr), arr))

        # With coercing, the invalid dates becomes iNaT
        self.assert_(
            np.array_equal(
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
        )


class TestTimestampNsOperations(unittest.TestCase):
    def setUp(self):
        if _np_version_under1p7:
            raise nose.SkipTest('numpy >= 1.7 required')
        self.timestamp = Timestamp(datetime.datetime.utcnow())

    def assert_ns_timedelta(self, modified_timestamp, expected_value):
        value = self.timestamp.value
        modified_value = modified_timestamp.value

        self.assertEquals(modified_value - value, expected_value)

    def test_timedelta_ns_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'ns'), -123)

    def test_timedelta_ns_based_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(1234567898, 'ns'), 1234567898)

    def test_timedelta_us_arithmetic(self):
        self.assert_ns_timedelta(self.timestamp + np.timedelta64(-123, 'us'), -123000)

    def test_timedelta_ns_arithmetic(self):
        time = self.timestamp + np.timedelta64(-123, 'ms')
        self.assert_ns_timedelta(time, -123000000)

    def test_nanosecond_string_parsing(self):
        self.timestamp = Timestamp('2013-05-01 07:15:45.123456789')
        self.assertEqual(self.timestamp.value, 1367392545123456000)


class TestTslib(unittest.TestCase):

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

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
