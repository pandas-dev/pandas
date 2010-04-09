from datetime import datetime
import unittest

import numpy as np

from pandas.core.datetools import bday
from pandas.core.daterange import DateRange, XDateRange
from pandas.util.testing import assert_almost_equal

####
## XDateRange Tests
####


def eqXDateRange(kwargs, expected):
    assert(np.array_equal(list(XDateRange(**kwargs)), expected))

def testXDateRange1():
    eqXDateRange(dict(fromDate = datetime(2009, 3, 25),
                      nPeriods = 2),
                 [datetime(2009, 3, 25), datetime(2009, 3, 26)])

def testXDateRange2():
    eqXDateRange(dict(fromDate = datetime(2008, 1, 1),
                      toDate = datetime(2008, 1, 3)),
                 [datetime(2008, 1, 1),
                  datetime(2008, 1, 2),
                  datetime(2008, 1, 3)])

def testXDateRange3():
    eqXDateRange(dict(fromDate = datetime(2008, 1, 5),
                      toDate = datetime(2008, 1, 6)),
                 [])


class TestDateRange(unittest.TestCase):
    pass


# DateRange test

def testDateRange1():
    toDate = datetime(2009, 5, 13)
    dr = DateRange(toDate=toDate, periods=20)
    firstDate = toDate - 19 * bday

    assert len(dr) == 20
    assert dr[0] == firstDate
    assert dr[-1] == toDate

