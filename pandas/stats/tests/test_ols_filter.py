from datetime import datetime
import unittest

from numpy import NaN
import numpy as np

from pandas.core.datetools import bday
from pandas.core.api import DateRange, Series, DataFrame
from pandas.stats.ols import _filter_data

class TestOLSFilter(unittest.TestCase):

    def setUp(self):
        date_index = DateRange(datetime(2009, 12, 11), periods=3, offset=bday)
        ts = Series([3, 1, 4], index=date_index)
        self.TS1 = ts

        date_index = DateRange(datetime(2009, 12, 11), periods=5, offset=bday)
        ts = Series([1, 5, 9, 2, 6], index=date_index)
        self.TS2 = ts

        date_index = DateRange(datetime(2009, 12, 11), periods=3, offset=bday)
        ts = Series([5, NaN, 3], index=date_index)
        self.TS3 = ts

        date_index = DateRange(datetime(2009, 12, 11), periods=5, offset=bday)
        ts = Series([NaN, 5, 8, 9, 7], index=date_index)
        self.TS4 = ts

        data = {'x1' : self.TS2, 'x2' : self.TS4}
        self.DF1 = DataFrame(data=data)

        data = {'x1' : self.TS2, 'x2' : self.TS4}
        self.DICT1 = data

    def testFilterWithSeriesRHS(self):
        (lhs, rhs, rhs_pre,
        index, valid) = _filter_data(self.TS1, {'x1' : self.TS2})
        self.tsAssertEqual(self.TS1, lhs)
        self.tsAssertEqual(self.TS2[:3], rhs['x1'])
        self.tsAssertEqual(self.TS2, rhs_pre['x1'])

    def testFilterWithSeriesRHS2(self):
        (lhs, rhs, rhs_pre,
        index, valid) = _filter_data(self.TS2, {'x1' : self.TS1})
        self.tsAssertEqual(self.TS2[:3], lhs)
        self.tsAssertEqual(self.TS1, rhs['x1'])
        self.tsAssertEqual(self.TS1, rhs_pre['x1'])

    def testFilterWithSeriesRHS3(self):
        (lhs, rhs, rhs_pre,
        index, valid) = _filter_data(self.TS3, {'x1' : self.TS4})
        exp_lhs = self.TS3[2:3]
        exp_rhs = self.TS4[2:3]
        exp_rhs_pre = self.TS4[1:]
        self.tsAssertEqual(exp_lhs, lhs)
        self.tsAssertEqual(exp_rhs, rhs['x1'])
        self.tsAssertEqual(exp_rhs_pre, rhs_pre['x1'])

    def testFilterWithDataFrameRHS(self):
        (lhs, rhs, rhs_pre,
        index, valid) = _filter_data(self.TS1, self.DF1)
        exp_lhs = self.TS1[1:]
        exp_rhs1 = self.TS2[1:3]
        exp_rhs2 = self.TS4[1:3]
        self.tsAssertEqual(exp_lhs, lhs)
        self.tsAssertEqual(exp_rhs1, rhs['x1'])
        self.tsAssertEqual(exp_rhs2, rhs['x2'])

    def testFilterWithDictRHS(self):
        (lhs, rhs, rhs_pre,
        index, valid) = _filter_data(self.TS1, self.DICT1)
        exp_lhs = self.TS1[1:]
        exp_rhs1 = self.TS2[1:3]
        exp_rhs2 = self.TS4[1:3]
        self.tsAssertEqual(exp_lhs, lhs)
        self.tsAssertEqual(exp_rhs1, rhs['x1'])
        self.tsAssertEqual(exp_rhs2, rhs['x2'])

    def tsAssertEqual(self, ts1, ts2):
        self.assert_(np.array_equal(ts1, ts2))

if __name__ == '__main__':
    unittest.main()
