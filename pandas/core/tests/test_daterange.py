from datetime import datetime
import unittest

import numpy as np

import pandas.core.datetools as datetools
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

START, END = datetime(2009, 1, 1), datetime(2010, 1, 1)

class TestXDateRange(unittest.TestCase):
    def test_constructor(self):
        rng = XDateRange(START, END, offset=datetools.bday)
        self.assertEquals(rng.timeRule, 'WEEKDAY')

        rng = XDateRange(START, END, timeRule='WEEKDAY')
        self.assertEquals(rng.offset, datetools.bday)

class TestDateRange(unittest.TestCase):
    def setUp(self):
        self.rng = DateRange(START, END, offset=datetools.bday)

    def test_constructor(self):
        rng = DateRange(START, END, offset=datetools.bday)

        rng = DateRange(START, periods=20, offset=datetools.bday)

        rng = DateRange(toDate=START, periods=20, offset=datetools.bday)

    def test_getCachedRange(self):
        rng = DateRange.getCachedRange(START, END, offset=datetools.bday)

        rng = DateRange.getCachedRange(START, periods=20, offset=datetools.bday)

        rng = DateRange.getCachedRange(end=START, periods=20,
                                       offset=datetools.bday)

        self.assertRaises(Exception, DateRange.getCachedRange, START, END)

        self.assertRaises(Exception, DateRange.getCachedRange, START,
                          offset=datetools.bday)

        self.assertRaises(Exception, DateRange.getCachedRange, end=END,
                          offset=datetools.bday)

        self.assertRaises(Exception, DateRange.getCachedRange, periods=20,
                          offset=datetools.bday)

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        self.assert_(comp[11])
        self.assert_(not comp[9])

    def test_repr(self):
        foo = repr(self.rng)

    def test_getitem(self):
        sliced = self.rng[10:20]

        sliced = self.rng[::5]

        fancy_indexed = self.rng[[4, 3, 2, 1, 0]]
        self.assertEquals(len(fancy_indexed), 5)
        self.assert_(not isinstance(fancy_indexed, DateRange))

    def test_shift(self):
        shifted = self.rng.shift(5)

        shifted = self.rng.shift(-5)

        shifted = self.rng.shift(0)

# DateRange test

def testDateRange1():
    toDate = datetime(2009, 5, 13)
    dr = DateRange(toDate=toDate, periods=20)
    firstDate = toDate - 19 * datetools.bday

    assert len(dr) == 20
    assert dr[0] == firstDate
    assert dr[-1] == toDate

