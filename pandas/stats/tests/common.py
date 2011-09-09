# pylint: disable-msg=W0611,W0402

from datetime import datetime
import string
import unittest
import nose

import numpy as np

from pandas.core.api import DataFrame, DateRange
from pandas.util.testing import assert_almost_equal # imported in other tests
N = 100
K = 4

start = datetime(2007, 1, 1)
DATE_RANGE = DateRange(start, periods=N)

COLS = ['Col' + c for c in string.ascii_uppercase[:K]]

def makeDataFrame():
    data = DataFrame(np.random.randn(N, K),
                      columns=COLS,
                      index=DATE_RANGE)

    return data

def getBasicDatasets():
    A = makeDataFrame()
    B = makeDataFrame()
    C = makeDataFrame()

    return A, B, C

def check_for_scipy():
    try:
        import scipy
    except ImportError:
        raise nose.SkipTest('no scipy')

def check_for_statsmodels():
    try:
        import scikits.statsmodels as sm
    except Exception:
        raise nose.SkipTest('no statsmodels')


class BaseTest(unittest.TestCase):
    def setUp(self):
        check_for_scipy()
        check_for_statsmodels()


        self.A, self.B, self.C = getBasicDatasets()

        self.createData1()
        self.createData2()
        self.createData3()

    def createData1(self):
        date = datetime(2007, 1, 1)
        date2 = datetime(2007, 1, 15)
        date3 = datetime(2007, 1, 22)

        A = self.A.copy()
        B = self.B.copy()
        C = self.C.copy()

        A['ColA'][date] = np.NaN
        B['ColA'][date] = np.NaN
        C['ColA'][date] = np.NaN
        C['ColA'][date2] = np.NaN

        # truncate data to save time
        A = A[:30]
        B = B[:30]
        C = C[:30]

        self.panel_y = A
        self.panel_x = {'B' : B, 'C' : C}

        self.series_panel_y = A.filter(['ColA'])
        self.series_panel_x = {'B' : B.filter(['ColA']),
                               'C' : C.filter(['ColA'])}
        self.series_y = A['ColA']
        self.series_x = {'B' : B['ColA'],
                         'C' : C['ColA']}

    def createData2(self):
        y_data = [[1, np.NaN],
                  [2, 3],
                  [4, 5]]
        y_index = [datetime(2000, 1, 1),
                   datetime(2000, 1, 2),
                   datetime(2000, 1, 3)]
        y_cols = ['A', 'B']
        self.panel_y2 = DataFrame(np.array(y_data), index=y_index,
                                   columns=y_cols)

        x1_data = [[6, np.NaN],
                   [7, 8],
                   [9, 30],
                   [11, 12]]
        x1_index = [datetime(2000, 1, 1),
                    datetime(2000, 1, 2),
                    datetime(2000, 1, 3),
                    datetime(2000, 1, 4)]
        x1_cols = ['A', 'B']
        x1 = DataFrame(np.array(x1_data), index=x1_index,
                        columns=x1_cols)

        x2_data = [[13, 14, np.NaN],
                   [15, np.NaN, np.NaN],
                   [16, 17, 48],
                   [19, 20, 21],
                   [22, 23, 24]]
        x2_index = [datetime(2000, 1, 1),
                    datetime(2000, 1, 2),
                    datetime(2000, 1, 3),
                    datetime(2000, 1, 4),
                    datetime(2000, 1, 5)]
        x2_cols = ['C', 'A', 'B']
        x2 = DataFrame(np.array(x2_data), index=x2_index,
                        columns=x2_cols)

        self.panel_x2 = {'x1' : x1, 'x2' : x2}

    def createData3(self):
        y_data = [[1, 2],
                  [3, 4]]
        y_index = [datetime(2000, 1, 1),
                   datetime(2000, 1, 2)]
        y_cols = ['A', 'B']
        self.panel_y3 = DataFrame(np.array(y_data), index=y_index,
                                   columns=y_cols)

        x1_data = [['A', 'B'],
                   ['C', 'A']]
        x1_index = [datetime(2000, 1, 1),
                    datetime(2000, 1, 2)]
        x1_cols = ['A', 'B']
        x1 = DataFrame(np.array(x1_data), index=x1_index,
                        columns=x1_cols)

        x2_data = [['foo', 'bar'],
                   ['baz', 'foo']]
        x2_index = [datetime(2000, 1, 1),
                    datetime(2000, 1, 2)]
        x2_cols = ['A', 'B']
        x2 = DataFrame(np.array(x2_data), index=x2_index,
                        columns=x2_cols)

        self.panel_x3 = {'x1' : x1, 'x2' : x2}
