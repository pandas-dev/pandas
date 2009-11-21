from datetime import datetime
import string
import unittest

import numpy as np

from pandas.core.api import DataMatrix, DateRange


N = 100
K = 4

start = datetime(2007, 1, 1)
DATE_RANGE = DateRange(start, periods=N)

COLS = ['Col' + c for c in string.ascii_uppercase[:K]]

def makeDataMatrix():
    data = DataMatrix(np.random.randn(N, K),
                      columns=COLS,
                      index=DATE_RANGE)

    return data

def isiterable(obj):
    return getattr(obj, '__iter__', False)

def assert_almost_equal(a, b):
    if isiterable(a):
        np.testing.assert_(isiterable(b))
        np.testing.assert_equal(len(a), len(b))
        for i in xrange(len(a)):
            assert_almost_equal(a[i], b[i])
        return

    err_msg = lambda a, b: 'expected %.5f but got %.5f' % (a, b)

    if np.isnan(a):
        np.testing.assert_(np.isnan(b))
        return

    # case for zero
    if abs(a) < 1e-5:
        np.testing.assert_almost_equal(
            a, b, decimal=5, err_msg=err_msg(a, b), verbose=False)
    else:
        np.testing.assert_almost_equal(
            1, a/b, decimal=5, err_msg=err_msg(a, b), verbose=False)


def getBasicDatasets():
    A = makeDataMatrix()
    B = makeDataMatrix()
    C = makeDataMatrix()

    return A, B, C

class BaseTest(unittest.TestCase):
    def setUp(self):
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

        self.series_panel_y = A.filterItems(['ColA'])
        self.series_panel_x = {'B' : B.filterItems(['ColA']),
                               'C' : C.filterItems(['ColA'])}
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
        self.panel_y2 = DataMatrix(np.array(y_data), index=y_index, columns=y_cols)

        x1_data = [[6, np.NaN],
                   [7, 8],
                   [9, 30],
                   [11, 12]]
        x1_index = [datetime(2000, 1, 1),
                    datetime(2000, 1, 2),
                    datetime(2000, 1, 3),
                    datetime(2000, 1, 4)]
        x1_cols = ['A', 'B']
        x1 = DataMatrix(np.array(x1_data), index=x1_index, columns=x1_cols)

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
        x2 = DataMatrix(np.array(x2_data), index=x2_index, columns=x2_cols)

        self.panel_x2 = {'x1' : x1, 'x2' : x2}

    def createData3(self):
        y_data = [[1, 2],
                  [3, 4]]
        y_index = [datetime(2000, 1, 1),
                   datetime(2000, 1, 2)]
        y_cols = ['A', 'B']
        self.panel_y3 = DataMatrix(np.array(y_data), index=y_index, columns=y_cols)

        x1_data = [['A', 'B'],
                   ['C', 'A']]
        x1_index = [datetime(2000, 1, 1),
                    datetime(2000, 1, 2)]
        x1_cols = ['A', 'B']
        x1 = DataMatrix(np.array(x1_data), index=x1_index, columns=x1_cols)

        x2_data = [['3.14', '1.59'],
                   ['2.65', '3.14']]
        x2_index = [datetime(2000, 1, 1),
                    datetime(2000, 1, 2)]
        x2_cols = ['A', 'B']
        x2 = DataMatrix(np.array(x2_data), index=x2_index, columns=x2_cols)

        self.panel_x3 = {'x1' : x1, 'x2' : x2}
