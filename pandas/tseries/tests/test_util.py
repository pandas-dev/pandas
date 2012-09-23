import nose
import unittest

import numpy as np

from pandas import Series, date_range
import pandas.util.testing as tm

from datetime import datetime, date

from pandas.tseries.tools import normalize_date
from pandas.tseries.util import pivot_annual, isleapyear

class TestPivotAnnual(unittest.TestCase):
    """
    New pandas of scikits.timeseries pivot_annual
    """
    def test_daily(self):
        rng = date_range('1/1/2000', '12/31/2004', freq='D')
        ts = Series(np.random.randn(len(rng)), index=rng)

        annual = pivot_annual(ts, 'D')

        doy = ts.index.dayofyear
        doy[(-isleapyear(ts.index.year)) & (doy >= 60)] += 1

        for i in range(1, 367):
            subset = ts[doy == i]
            subset.index = [x.year for x in subset.index]

            tm.assert_series_equal(annual[i].dropna(), subset)

        # check leap days
        leaps = ts[(ts.index.month == 2) & (ts.index.day == 29)]
        day = leaps.index.dayofyear[0]
        leaps.index = leaps.index.year
        tm.assert_series_equal(annual[day].dropna(), leaps)

    def test_weekly(self):
        pass

    def test_monthly(self):
        rng = date_range('1/1/2000', '12/31/2004', freq='M')
        ts = Series(np.random.randn(len(rng)), index=rng)

        annual = pivot_annual(ts, 'M')

        month = ts.index.month

        for i in range(1, 13):
            subset = ts[month == i]
            subset.index = [x.year for x in subset.index]
            tm.assert_series_equal(annual[i].dropna(), subset)

    def test_period_monthly(self):
        pass

    def test_period_daily(self):
        pass

    def test_period_weekly(self):
        pass


def test_normalize_date():
    value = date(2012, 9, 7)

    result = normalize_date(value)
    assert(result == datetime(2012, 9, 7))

    value = datetime(2012, 9, 7, 12)

    result = normalize_date(value)
    assert(result == datetime(2012, 9, 7))

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

