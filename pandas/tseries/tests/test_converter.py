from datetime import datetime, time, timedelta, date
import sys
import os

import nose

import numpy as np
from pandas.compat import u
import pandas.util.testing as tm

try:
    import pandas.tseries.converter as converter
except ImportError:
    raise nose.SkipTest("no pandas.tseries.converter, skipping")


def test_timtetonum_accepts_unicode():
    assert(converter.time2num("00:01") == converter.time2num(u("00:01")))


class TestDateTimeConverter(tm.TestCase):

    def setUp(self):
        self.dtc = converter.DatetimeConverter()
        self.tc = converter.TimeFormatter(None)

    def test_convert_accepts_unicode(self):
        r1 = self.dtc.convert("12:22", None, None)
        r2 = self.dtc.convert(u("12:22"), None, None)
        assert(r1 == r2), "DatetimeConverter.convert should accept unicode"

    def test_conversion(self):
        rs = self.dtc.convert(['2012-1-1'], None, None)[0]
        xp = datetime(2012, 1, 1).toordinal()
        self.assertEqual(rs, xp)

        rs = self.dtc.convert('2012-1-1', None, None)
        self.assertEqual(rs, xp)

        rs = self.dtc.convert(date(2012, 1, 1), None, None)
        self.assertEqual(rs, xp)

        rs = self.dtc.convert(datetime(2012, 1, 1).toordinal(), None, None)
        self.assertEqual(rs, xp)

        rs = self.dtc.convert('2012-1-1', None, None)
        self.assertEqual(rs, xp)

    def test_time_formatter(self):
        self.tc(90000)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
