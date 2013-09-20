import unittest

import numpy as np

from pandas import tslib
from datetime import datetime

class TestDatetimeParsingWrappers(unittest.TestCase):
    def test_verify_datetime_bounds(self):
        for year in (1, 1000, 1677, 2262, 5000):
            dt = datetime(year, 1, 1)
            self.assertRaises(
                ValueError,
                tslib.verify_datetime_bounds,
                dt
            )

        for year in (1678, 2000, 2261):
            tslib.verify_datetime_bounds(datetime(year, 1, 1))

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

    def test_dates_outside_of_datetime64_ns_bounds(self):
        # These datetimes are outside of the bounds of the
        # datetime64[ns] bounds, so they cannot be converted to
        # datetimes
        arr = np.array(['1/1/1676', '1/2/1676'], dtype=object)
        self.assert_(np.array_equal(tslib.array_to_datetime(arr), arr))

        arr = np.array(['1/1/2263', '1/2/2263'], dtype=object)
        self.assert_(np.array_equal(tslib.array_to_datetime(arr), arr))

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

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
