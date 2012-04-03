# pylint: disable-msg=E1101,W0612

from cStringIO import StringIO
from datetime import datetime, timedelta
import os
import operator
import unittest

import nose

from numpy import nan
import numpy as np
import numpy.ma as ma

from pandas import Index, Series, TimeSeries, DataFrame, isnull, notnull
from pandas.core.index import MultiIndex

from pandas import DatetimeIndex

import pandas.core.datetools as datetools
import pandas.core.nanops as nanops

from pandas.util import py3compat
from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm
import pandas

class TestTimeSeriesDuplicates(unittest.TestCase):

    def setUp(self):
        dates = [datetime(2000, 1, 2), datetime(2000, 1, 2),
                 datetime(2000, 1, 2), datetime(2000, 1, 3),
                 datetime(2000, 1, 3), datetime(2000, 1, 3),
                 datetime(2000, 1, 4), datetime(2000, 1, 4),
                 datetime(2000, 1, 4), datetime(2000, 1, 5)]

        self.dups = Series(np.random.randn(len(dates)), index=dates)

    def test_constructor(self):
        self.assert_(isinstance(self.dups, TimeSeries))
        self.assert_(isinstance(self.dups.index, DatetimeIndex))

    def test_is_unique_monotonic(self):
        self.assert_(not self.dups.index.is_unique)


    def test_duplicate_dates_indexing(self):

        for date in ts.index.unique():
            result = ts[date]

            mask = ts.index == date
            total = (ts.index == date).sum()
            expected = ts[mask]
            if total > 1:
                assert_series_equal(result, expected)
            else:
                assert_almost_equal(result, expected[0])



if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
