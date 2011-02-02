from datetime import datetime
import pickle
import unittest

import numpy as np

import pandas.core.datetools as datetools
from pandas.core.index import Index
from pandas.core.daterange import DateRange, generate_range

def eqXDateRange(kwargs, expected):
    rng = generate_range(**kwargs)
    assert(np.array_equal(list(rng), expected))

START, END = datetime(2009, 1, 1), datetime(2010, 1, 1)

class TestGeneration(unittest.TestCase):
    def test_generate(self):
        rng1 = list(generate_range(START, END, offset=datetools.bday))
        rng2 = list(generate_range(START, END, timeRule='WEEKDAY'))
        self.assert_(np.array_equal(rng1, rng2))

    def test_1(self):
        eqXDateRange(dict(start=datetime(2009, 3, 25),
                          periods=2),
                     [datetime(2009, 3, 25), datetime(2009, 3, 26)])

    def test_2(self):
        eqXDateRange(dict(start=datetime(2008, 1, 1),
                          end=datetime(2008, 1, 3)),
                     [datetime(2008, 1, 1),
                      datetime(2008, 1, 2),
                      datetime(2008, 1, 3)])

    def test_3(self):
        eqXDateRange(dict(start = datetime(2008, 1, 5),
                          end = datetime(2008, 1, 6)),
                     [])

class TestDateRange(unittest.TestCase):

    def setUp(self):
        self.rng = DateRange(START, END, offset=datetools.bday)

    def test_constructor(self):
        rng = DateRange(START, END, offset=datetools.bday)
        rng = DateRange(START, periods=20, offset=datetools.bday)
        rng = DateRange(end=START, periods=20, offset=datetools.bday)

    def test_cached_range(self):
        rng = DateRange._cached_range(START, END, offset=datetools.bday)
        rng = DateRange._cached_range(START, periods=20, offset=datetools.bday)
        rng = DateRange._cached_range(end=START, periods=20,
                                       offset=datetools.bday)

        self.assertRaises(Exception, DateRange._cached_range, START, END)

        self.assertRaises(Exception, DateRange._cached_range, START,
                          offset=datetools.bday)

        self.assertRaises(Exception, DateRange._cached_range, end=END,
                          offset=datetools.bday)

        self.assertRaises(Exception, DateRange._cached_range, periods=20,
                          offset=datetools.bday)

    def test_cached_range_bug(self):
        rng = DateRange('2010-09-01 05:00:00', periods=50,
                        offset=datetools.DateOffset(hours=6))
        self.assertEquals(len(rng), 50)
        self.assertEquals(rng[0], datetime(2010, 9, 1, 5))

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        self.assert_(comp[11])
        self.assert_(not comp[9])

    def test_repr(self):
        # only really care that it works
        repr(self.rng)

    def test_getitem(self):
        smaller = self.rng[:5]
        self.assert_(np.array_equal(smaller, self.rng.view(np.ndarray)[:5]))
        self.assertEquals(smaller.offset, self.rng.offset)

        sliced = self.rng[::5]
        self.assertEquals(sliced.offset, datetools.bday * 5)

        fancy_indexed = self.rng[[4, 3, 2, 1, 0]]
        self.assertEquals(len(fancy_indexed), 5)
        self.assert_(not isinstance(fancy_indexed, DateRange))

        # 32-bit vs. 64-bit platforms
        self.assertEquals(self.rng[4], self.rng[np.int_(4)])

    def test_shift(self):
        shifted = self.rng.shift(5)
        self.assertEquals(shifted[0], self.rng[5])
        self.assertEquals(shifted.offset, self.rng.offset)

        shifted = self.rng.shift(-5)
        self.assertEquals(shifted[5], self.rng[0])
        self.assertEquals(shifted.offset, self.rng.offset)

        shifted = self.rng.shift(0)
        self.assertEquals(shifted[0], self.rng[0])
        self.assertEquals(shifted.offset, self.rng.offset)

        rng = DateRange(START, END, offset=datetools.bmonthEnd)
        shifted = rng.shift(1, offset=datetools.bday)
        self.assertEquals(shifted[0], rng[0] + datetools.bday)

    def test_pickle_unpickle(self):
        pickled = pickle.dumps(self.rng)
        unpickled = pickle.loads(pickled)

        self.assert_(unpickled.offset is not None)

    def test_union(self):
        # overlapping
        left = self.rng[:10]
        right = self.rng[5:10]

        the_union = left.union(right)
        self.assert_(isinstance(the_union, DateRange))

        # non-overlapping, gap in middle
        left = self.rng[:5]
        right = self.rng[10:]

        the_union = left.union(right)
        self.assert_(isinstance(the_union, Index))
        self.assert_(not isinstance(the_union, DateRange))

        # non-overlapping, no gap
        left = self.rng[:5]
        right = self.rng[5:10]

        the_union = left.union(right)
        self.assert_(isinstance(the_union, DateRange))

        # order does not matter
        self.assert_(np.array_equal(right.union(left), the_union))

        # overlapping, but different offset
        rng = DateRange(START, END, offset=datetools.bmonthEnd)

        the_union = self.rng.union(rng)
        self.assert_(not isinstance(the_union, DateRange))

# DateRange test

def testDateRange1():
    end = datetime(2009, 5, 13)
    dr = DateRange(end=end, periods=20)
    firstDate = end - 19 * datetools.bday

    assert len(dr) == 20
    assert dr[0] == firstDate
    assert dr[-1] == end

