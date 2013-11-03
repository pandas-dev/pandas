from datetime import datetime, time, timedelta, date
import sys
import os
import unittest

import nose

import numpy as np

from pandas import Timestamp, date_range

class TestToJulianDateTimestamp(unittest.TestCase):

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

class TestDateTimeIndex(unittest.TestCase):
    def test_1700(self):
        r1 = [2345897.5, 2345898.5, 2345899.5, 2345900.5, 2345901.5]
        r2 = date_range(start=Timestamp('1710-10-01'), periods=5, freq='D').to_julian_date()
        self.assert_(isinstance(r2, np.ndarray))
        np.testing.assert_array_equal(r1, r2)

    def test_2000(self):
        r1 = [2451601.5, 2451602.5, 2451603.5, 2451604.5, 2451605.5]
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5, freq='D').to_julian_date()
        self.assert_(isinstance(r2, np.ndarray))
        np.testing.assert_array_equal(r1, r2)

    def test_hour(self):
        r1 = [2451601.5, 2451601.5416666666666666, 2451601.5833333333333333, 2451601.625, 2451601.6666666666666666]
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5, freq='H').to_julian_date()
        self.assert_(isinstance(r2, np.ndarray))
        np.testing.assert_array_equal(r1, r2)

    def test_minute(self):
        r1 = [2451601.5, 2451601.5006944444444444, 2451601.5013888888888888, 2451601.5020833333333333, 2451601.5027777777777777]
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5, freq='T').to_julian_date()
        self.assert_(isinstance(r2, np.ndarray))
        np.testing.assert_array_equal(r1, r2)

    def test_second(self):
        r1 = [2451601.5, 2451601.500011574074074, 2451601.5000231481481481, 2451601.5000347222222222, 2451601.5000462962962962]
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5, freq='S').to_julian_date()
        self.assert_(isinstance(r2, np.ndarray))
        np.testing.assert_array_equal(r1, r2)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
