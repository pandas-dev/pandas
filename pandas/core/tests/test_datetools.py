from datetime import datetime, timedelta
import unittest

from pandas.core.datetools import (
    bday, BDay, BQuarterEnd, BMonthEnd, BYearEnd, MonthEnd,
    DateOffset, Week)

from pandas.core.daterange import XDateRange, DateRange
import pandas.core.datetools as datetools

####
## Misc function tests
####
def test_format():
    actual = datetools.format(datetime(2008, 1, 15))
    assert actual == '20080115'

def test_ole2datetime():
    actual = datetools.ole2datetime(60000)
    assert actual == datetime(2064, 4, 8)

def test_to_datetime1():
    actual = datetools.to_datetime(datetime(2008, 1, 15))
    assert actual == datetime(2008, 1, 15)

    actual = datetools.to_datetime('20080115')
    assert actual == datetime(2008, 1, 15)

    # unparseable
    s = 'Month 1, 1999'
    assert datetools.to_datetime(s) == s

def test_normalize_date():
    actual = datetools.normalize_date(datetime(2007, 10, 1, 1, 12, 5, 10))
    assert actual == datetime(2007, 10, 1)

#####
### DateOffset Tests
#####

def myAssert(actual, expected):
    assert actual == expected

class TestDateOffset(object):

    def setUp(self):
        self.d = datetime(2008, 1, 2)

    def test_repr(self):
        repr(DateOffset())
        repr(DateOffset(2))
        repr(2 * DateOffset())
        repr(2 * DateOffset(months=2))

    def test_mul(self):
        assert DateOffset(2) == 2 * DateOffset(1)
        assert DateOffset(2) == DateOffset(1) * 2

    def test_constructor(self):

        assert((self.d + DateOffset(months=2)) == datetime(2008, 3, 2))
        assert((self.d - DateOffset(months=2)) == datetime(2007, 11, 2))

        assert((self.d + DateOffset(2)) == datetime(2008, 1, 4))

        assert not DateOffset(2).isAnchored()
        assert DateOffset(1).isAnchored()

        d = datetime(2008, 1, 31)
        assert((d + DateOffset(months=1)) == datetime(2008, 2, 29))

    def test_copy(self):
        assert(DateOffset(months=2).copy() == DateOffset(months=2))

class TestBusinessDay(unittest.TestCase):

    def setUp(self):
        self.d = datetime(2008, 1, 1)

        self.offset = BDay()
        self.offset2 = BDay(2)

    def test_repr(self):
        assert repr(self.offset) == '<1 BusinessDay>'
        assert repr(self.offset2) == '<2 BusinessDays>'

    def test_with_offset(self):
        offset = self.offset + timedelta(hours=2)

        assert (self.d + offset) == datetime(2008, 1, 2, 2)

    def testEQ(self):
        myAssert(self.offset2, self.offset2)

    def test_mul(self):
        pass

    def test_hash(self):
        myAssert(hash(self.offset2), hash(self.offset2))

    def testCall(self):
        myAssert(self.offset2(self.d), datetime(2008, 1, 3))

    def testRAdd(self):
        myAssert(self.d + self.offset2, self.offset2 + self.d)

    def testSub(self):
        myAssert(self.d - self.offset2,  self.d + BDay(-2))

    def testRSub(self):
        myAssert(self.d - self.offset2, self.offset2 - self.d)

    def testMult1(self):
        myAssert(self.d + 10*self.offset, self.d + BDay(10))

    def testMult2(self):
        myAssert(self.d + (-5*BDay(-10)),
                 self.d + BDay(50))


    def testRollback1(self):
        myAssert(BDay(10).rollback(self.d), self.d)

    def testRollback2(self):
        myAssert(BDay(10).rollback(datetime(2008, 1, 5)), datetime(2008, 1, 4))

    def testRollforward1(self):
        myAssert(BDay(10).rollforward(self.d), self.d)

    def testRollforward2(self):
        myAssert(BDay(10).rollforward(datetime(2008, 1, 5)), datetime(2008, 1, 7))

    def test_onOffset(self):
        tests = [(BDay(), datetime(2008, 1, 1), True),
                 (BDay(), datetime(2008, 1, 5), False)]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_offset(self):
        tests = []

        tests.append((bday,
                      {datetime(2008, 1, 1): datetime(2008, 1, 2),
                       datetime(2008, 1, 4): datetime(2008, 1, 7),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 8)}))

        tests.append((2*bday,
                      {datetime(2008, 1, 1): datetime(2008, 1, 3),
                       datetime(2008, 1, 4): datetime(2008, 1, 8),
                       datetime(2008, 1, 5): datetime(2008, 1, 8),
                       datetime(2008, 1, 6): datetime(2008, 1, 8),
                       datetime(2008, 1, 7): datetime(2008, 1, 9)}))

        tests.append((-bday,
                      {datetime(2008, 1, 1): datetime(2007, 12, 31),
                       datetime(2008, 1, 4): datetime(2008, 1, 3),
                       datetime(2008, 1, 5): datetime(2008, 1, 4),
                       datetime(2008, 1, 6): datetime(2008, 1, 4),
                       datetime(2008, 1, 7): datetime(2008, 1, 4),
                       datetime(2008, 1, 8): datetime(2008, 1, 7)}))

        tests.append((-2*bday,
                      {datetime(2008, 1, 1): datetime(2007, 12, 28),
                       datetime(2008, 1, 4): datetime(2008, 1, 2),
                       datetime(2008, 1, 5): datetime(2008, 1, 3),
                       datetime(2008, 1, 6): datetime(2008, 1, 3),
                       datetime(2008, 1, 7): datetime(2008, 1, 3),
                       datetime(2008, 1, 8): datetime(2008, 1, 4),
                       datetime(2008, 1, 9): datetime(2008, 1, 7)}))

        tests.append((BDay(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 1, 4): datetime(2008, 1, 4),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 7)}))

        for dateOffset, cases in tests:
            for baseDate, expected in cases.iteritems():
                assertEq(dateOffset, baseDate, expected)



def assertOnOffset(offset, date, expected):
    actual = offset.onOffset(date)
    assert actual == expected

class TestWeek(unittest.TestCase):
    def test_offset(self):
        tests = []

        tests.append((datetools.week, # not business week
                      {datetime(2008, 1, 1): datetime(2008, 1, 8),
                       datetime(2008, 1, 4): datetime(2008, 1, 11),
                       datetime(2008, 1, 5): datetime(2008, 1, 12),
                       datetime(2008, 1, 6): datetime(2008, 1, 13),
                       datetime(2008, 1, 7): datetime(2008, 1, 14)}))

        tests.append((Week(dayOfWeek=0), # Mon
                      {datetime(2007, 12, 31): datetime(2008, 1, 7),
                       datetime(2008, 1, 4): datetime(2008, 1, 7),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 14)}))

        tests.append((Week(0, dayOfWeek=0), # n=0 -> roll forward. Mon
                      {datetime(2007, 12, 31): datetime(2007, 12, 31),
                       datetime(2008, 1, 4): datetime(2008, 1, 7),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 7)}))

        for dateOffset, cases in tests:
            for baseDate, expected in cases.iteritems():
                assertEq(dateOffset, baseDate, expected)

    def test_onOffset(self):

        tests = [(Week(dayOfWeek=0), datetime(2008, 1, 1), False),
                 (Week(dayOfWeek=0), datetime(2008, 1, 2), False),
                 (Week(dayOfWeek=0), datetime(2008, 1, 3), False),
                 (Week(dayOfWeek=0), datetime(2008, 1, 4), False),
                 (Week(dayOfWeek=0), datetime(2008, 1, 5), False),
                 (Week(dayOfWeek=0), datetime(2008, 1, 6), False),
                 (Week(dayOfWeek=0), datetime(2008, 1, 7), True),

                 (Week(dayOfWeek=1), datetime(2008, 1, 1), True),
                 (Week(dayOfWeek=1), datetime(2008, 1, 2), False),
                 (Week(dayOfWeek=1), datetime(2008, 1, 3), False),
                 (Week(dayOfWeek=1), datetime(2008, 1, 4), False),
                 (Week(dayOfWeek=1), datetime(2008, 1, 5), False),
                 (Week(dayOfWeek=1), datetime(2008, 1, 6), False),
                 (Week(dayOfWeek=1), datetime(2008, 1, 7), False),

                 (Week(dayOfWeek=2), datetime(2008, 1, 1), False),
                 (Week(dayOfWeek=2), datetime(2008, 1, 2), True),
                 (Week(dayOfWeek=2), datetime(2008, 1, 3), False),
                 (Week(dayOfWeek=2), datetime(2008, 1, 4), False),
                 (Week(dayOfWeek=2), datetime(2008, 1, 5), False),
                 (Week(dayOfWeek=2), datetime(2008, 1, 6), False),
                 (Week(dayOfWeek=2), datetime(2008, 1, 7), False),

                 (Week(dayOfWeek=3), datetime(2008, 1, 1), False),
                 (Week(dayOfWeek=3), datetime(2008, 1, 2), False),
                 (Week(dayOfWeek=3), datetime(2008, 1, 3), True),
                 (Week(dayOfWeek=3), datetime(2008, 1, 4), False),
                 (Week(dayOfWeek=3), datetime(2008, 1, 5), False),
                 (Week(dayOfWeek=3), datetime(2008, 1, 6), False),
                 (Week(dayOfWeek=3), datetime(2008, 1, 7), False),

                 (Week(dayOfWeek=4), datetime(2008, 1, 1), False),
                 (Week(dayOfWeek=4), datetime(2008, 1, 2), False),
                 (Week(dayOfWeek=4), datetime(2008, 1, 3), False),
                 (Week(dayOfWeek=4), datetime(2008, 1, 4), True),
                 (Week(dayOfWeek=4), datetime(2008, 1, 5), False),
                 (Week(dayOfWeek=4), datetime(2008, 1, 6), False),
                 (Week(dayOfWeek=4), datetime(2008, 1, 7), False),

                 (Week(dayOfWeek=5), datetime(2008, 1, 1), False),
                 (Week(dayOfWeek=5), datetime(2008, 1, 2), False),
                 (Week(dayOfWeek=5), datetime(2008, 1, 3), False),
                 (Week(dayOfWeek=5), datetime(2008, 1, 4), False),
                 (Week(dayOfWeek=5), datetime(2008, 1, 5), True),
                 (Week(dayOfWeek=5), datetime(2008, 1, 6), False),
                 (Week(dayOfWeek=5), datetime(2008, 1, 7), False),

                 (Week(dayOfWeek=6), datetime(2008, 1, 1), False),
                 (Week(dayOfWeek=6), datetime(2008, 1, 2), False),
                 (Week(dayOfWeek=6), datetime(2008, 1, 3), False),
                 (Week(dayOfWeek=6), datetime(2008, 1, 4), False),
                 (Week(dayOfWeek=6), datetime(2008, 1, 5), False),
                 (Week(dayOfWeek=6), datetime(2008, 1, 6), True),
                 (Week(dayOfWeek=6), datetime(2008, 1, 7), False),
             ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)
            pass

class TestBMonthEnd(unittest.TestCase):

    def test_offset(self):
        tests = []

        tests.append((BMonthEnd(),
                     {datetime(2008, 1, 1): datetime(2008, 1, 31),
                      datetime(2008, 1, 31): datetime(2008, 2, 29),
                      datetime(2006, 12, 29): datetime(2007, 1, 31),
                      datetime(2006, 12, 31): datetime(2007, 1, 31),
                      datetime(2007, 1, 1): datetime(2007, 1, 31),
                      datetime(2006, 12, 1): datetime(2006, 12, 29)}))

        tests.append((BMonthEnd(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 31),
                       datetime(2008, 1, 31): datetime(2008, 1, 31),
                       datetime(2006, 12, 29): datetime(2006, 12, 29),
                       datetime(2006, 12, 31): datetime(2007, 1, 31),
                       datetime(2007, 1, 1): datetime(2007, 1, 31)}))

        tests.append((BMonthEnd(2),
                     {datetime(2008, 1, 1): datetime(2008, 2, 29),
                      datetime(2008, 1, 31): datetime(2008, 3, 31),
                      datetime(2006, 12, 29): datetime(2007, 2, 28),
                      datetime(2006, 12, 31): datetime(2007, 2, 28),
                      datetime(2007, 1, 1): datetime(2007, 2, 28),
                      datetime(2006, 11, 1): datetime(2006, 12, 29)}))

        tests.append((BMonthEnd(-1),
                     {datetime(2007, 1, 1): datetime(2006, 12, 29),
                      datetime(2008, 6, 30): datetime(2008, 5, 30),
                      datetime(2008, 12, 31): datetime(2008, 11, 28),
                      datetime(2006, 12, 29): datetime(2006, 11, 30),
                      datetime(2006, 12, 30): datetime(2006, 12, 29),
                      datetime(2007, 1, 1): datetime(2006, 12, 29)}))

        for dateOffset, cases in tests:
            for baseDate, expected in cases.iteritems():
                assertEq(dateOffset, baseDate, expected)

    def test_onOffset(self):

        tests = [(BMonthEnd(), datetime(2007, 12, 31), True),
                 (BMonthEnd(), datetime(2008, 1, 1), False)]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

class TestBQuarterEnd(unittest.TestCase):
    def test_offset(self):
        tests = []

        tests.append((BQuarterEnd(startingMonth=1),
                      {datetime(2008, 1, 1): datetime(2008, 1, 31),
                       datetime(2008, 1, 31): datetime(2008, 4, 30),
                       datetime(2008, 2, 15): datetime(2008, 4, 30),
                       datetime(2008, 2, 29): datetime(2008, 4, 30),
                       datetime(2008, 3, 15): datetime(2008, 4, 30),
                       datetime(2008, 3, 31): datetime(2008, 4, 30),
                       datetime(2008, 4, 15): datetime(2008, 4, 30),
                       datetime(2008, 4, 30): datetime(2008, 7, 31),}))

        tests.append((BQuarterEnd(startingMonth=2),
                      {datetime(2008, 1, 1): datetime(2008, 2, 29),
                       datetime(2008, 1, 31): datetime(2008, 2, 29),
                       datetime(2008, 2, 15): datetime(2008, 2, 29),
                       datetime(2008, 2, 29): datetime(2008, 5, 30),
                       datetime(2008, 3, 15): datetime(2008, 5, 30),
                       datetime(2008, 3, 31): datetime(2008, 5, 30),
                       datetime(2008, 4, 15): datetime(2008, 5, 30),
                       datetime(2008, 4, 30): datetime(2008, 5, 30),}))

        tests.append((BQuarterEnd(startingMonth=1, n=0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 31),
                       datetime(2008, 1, 31): datetime(2008, 1, 31),
                       datetime(2008, 2, 15): datetime(2008, 4, 30),
                       datetime(2008, 2, 29): datetime(2008, 4, 30),
                       datetime(2008, 3, 15): datetime(2008, 4, 30),
                       datetime(2008, 3, 31): datetime(2008, 4, 30),
                       datetime(2008, 4, 15): datetime(2008, 4, 30),
                       datetime(2008, 4, 30): datetime(2008, 4, 30),}))

        tests.append((BQuarterEnd(startingMonth=1, n=-1),
                      {datetime(2008, 1, 1): datetime(2007, 10, 31),
                       datetime(2008, 1, 31): datetime(2007, 10, 31),
                       datetime(2008, 2, 15): datetime(2008, 1, 31),
                       datetime(2008, 2, 29): datetime(2008, 1, 31),
                       datetime(2008, 3, 15): datetime(2008, 1, 31),
                       datetime(2008, 3, 31): datetime(2008, 1, 31),
                       datetime(2008, 4, 15): datetime(2008, 1, 31),
                       datetime(2008, 4, 30): datetime(2008, 1, 31),}))

        tests.append((BQuarterEnd(startingMonth=1, n=2),
                      {datetime(2008, 1, 31): datetime(2008, 7, 31),
                       datetime(2008, 2, 15): datetime(2008, 7, 31),
                       datetime(2008, 2, 29): datetime(2008, 7, 31),
                       datetime(2008, 3, 15): datetime(2008, 7, 31),
                       datetime(2008, 3, 31): datetime(2008, 7, 31),
                       datetime(2008, 4, 15): datetime(2008, 7, 31),
                       datetime(2008, 4, 30): datetime(2008, 10, 31),}))

        for dateOffset, cases in tests:
            for baseDate, expected in cases.iteritems():
                assertEq(dateOffset, baseDate, expected)

    def test_onOffset(self):

        tests = [(BQuarterEnd(1, startingMonth=1), datetime(2008, 1, 31), True),
                 (BQuarterEnd(1, startingMonth=1), datetime(2007, 12, 31), False),
                 (BQuarterEnd(1, startingMonth=1), datetime(2008, 2, 29), False),
                 (BQuarterEnd(1, startingMonth=1), datetime(2007, 3, 30), False),
                 (BQuarterEnd(1, startingMonth=1), datetime(2007, 3, 31), False),
                 (BQuarterEnd(1, startingMonth=1), datetime(2008, 4, 30), True),
                 (BQuarterEnd(1, startingMonth=1), datetime(2008, 5, 30), False),
                 (BQuarterEnd(1, startingMonth=1), datetime(2007, 6, 29), False),
                 (BQuarterEnd(1, startingMonth=1), datetime(2007, 6, 30), False),

                 (BQuarterEnd(1, startingMonth=2), datetime(2008, 1, 31), False),
                 (BQuarterEnd(1, startingMonth=2), datetime(2007, 12, 31), False),
                 (BQuarterEnd(1, startingMonth=2), datetime(2008, 2, 29), True),
                 (BQuarterEnd(1, startingMonth=2), datetime(2007, 3, 30), False),
                 (BQuarterEnd(1, startingMonth=2), datetime(2007, 3, 31), False),
                 (BQuarterEnd(1, startingMonth=2), datetime(2008, 4, 30), False),
                 (BQuarterEnd(1, startingMonth=2), datetime(2008, 5, 30), True),
                 (BQuarterEnd(1, startingMonth=2), datetime(2007, 6, 29), False),
                 (BQuarterEnd(1, startingMonth=2), datetime(2007, 6, 30), False),

                 (BQuarterEnd(1, startingMonth=3), datetime(2008, 1, 31), False),
                 (BQuarterEnd(1, startingMonth=3), datetime(2007, 12, 31), True),
                 (BQuarterEnd(1, startingMonth=3), datetime(2008, 2, 29), False),
                 (BQuarterEnd(1, startingMonth=3), datetime(2007, 3, 30), True),
                 (BQuarterEnd(1, startingMonth=3), datetime(2007, 3, 31), False),
                 (BQuarterEnd(1, startingMonth=3), datetime(2008, 4, 30), False),
                 (BQuarterEnd(1, startingMonth=3), datetime(2008, 5, 30), False),
                 (BQuarterEnd(1, startingMonth=3), datetime(2007, 6, 29), True),
                 (BQuarterEnd(1, startingMonth=3), datetime(2007, 6, 30), False),
             ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)



class TestYearBegin(unittest.TestCase):

    def test_offset(self):
        tests = []

        tests.append((datetools.YearBegin(),
                      {datetime(2008, 1, 1): datetime(2009, 1, 1),
                       datetime(2008, 6, 30): datetime(2009, 1, 1),
                       datetime(2008, 12, 31): datetime(2009, 1, 1),
                       datetime(2005, 12, 30): datetime(2006, 1, 1),
                       datetime(2005, 12, 31): datetime(2006, 1, 1),}))

        tests.append((datetools.YearBegin(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 6, 30): datetime(2009, 1, 1),
                       datetime(2008, 12, 31): datetime(2009, 1, 1),
                       datetime(2005, 12, 30): datetime(2006, 1, 1),
                       datetime(2005, 12, 31): datetime(2006, 1, 1),}))


        tests.append((datetools.YearBegin(-1),
                      {datetime(2007, 1, 1): datetime(2006, 1, 1),
                       datetime(2008, 6, 30): datetime(2008, 1, 1),
                       datetime(2008, 12, 31): datetime(2008, 1, 1),
                       datetime(2006, 12, 29): datetime(2006, 1, 1),
                       datetime(2006, 12, 30): datetime(2006, 1, 1),
                       datetime(2007, 1, 1): datetime(2006, 1, 1),}))

        tests.append((datetools.YearBegin(-2),
                      {datetime(2007, 1, 1): datetime(2005, 1, 1),
                       datetime(2008, 6, 30): datetime(2007, 1, 1),
                       datetime(2008, 12, 31): datetime(2007, 1, 1),}))

        for dateOffset, cases in tests:
            for baseDate, expected in cases.iteritems():
                assertEq(dateOffset, baseDate, expected)



class TestBYearEnd(unittest.TestCase):

    def test_offset(self):
        tests = []

        tests.append((BYearEnd(),
                      {datetime(2008, 1, 1): datetime(2008, 12, 31),
                       datetime(2008, 6, 30): datetime(2008, 12, 31),
                       datetime(2008, 12, 31): datetime(2009, 12, 31),
                       datetime(2005, 12, 30): datetime(2006, 12, 29),
                       datetime(2005, 12, 31): datetime(2006, 12, 29),}))

        tests.append((BYearEnd(0),
                      {datetime(2008, 1, 1): datetime(2008, 12, 31),
                       datetime(2008, 6, 30): datetime(2008, 12, 31),
                       datetime(2008, 12, 31): datetime(2008, 12, 31),
                       datetime(2005, 12, 31): datetime(2006, 12, 29),}))

        tests.append((BYearEnd(-1),
                      {datetime(2007, 1, 1): datetime(2006, 12, 29),
                       datetime(2008, 6, 30): datetime(2007, 12, 31),
                       datetime(2008, 12, 31): datetime(2007, 12, 31),
                       datetime(2006, 12, 29): datetime(2005, 12, 30),
                       datetime(2006, 12, 30): datetime(2006, 12, 29),
                       datetime(2007, 1, 1): datetime(2006, 12, 29),}))

        tests.append((BYearEnd(-2),
                      {datetime(2007, 1, 1): datetime(2005, 12, 30),
                       datetime(2008, 6, 30): datetime(2006, 12, 29),
                       datetime(2008, 12, 31): datetime(2006, 12, 29),}))

        for dateOffset, cases in tests:
            for baseDate, expected in cases.iteritems():
                assertEq(dateOffset, baseDate, expected)

    def test_onOffset(self):

        tests = [
            (BYearEnd(), datetime(2007, 12, 31), True),
            (BYearEnd(), datetime(2008, 1, 1), False),
            (BYearEnd(), datetime(2006, 12, 31), False),
            (BYearEnd(), datetime(2006, 12, 29), True),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


class TestMonthEnd(unittest.TestCase):
    def test_onOffset(self):

        tests = [
            (MonthEnd(), datetime(2007, 3, 30), False),
            (MonthEnd(), datetime(2007, 3, 31), True),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

def testOnOffset():

    tests = [#(QuarterEnd(1, startingMonth=1), datetime(2008, 1, 31), True),
             #(QuarterEnd(1, startingMonth=1), datetime(2007, 12, 31), False),
             #(QuarterEnd(1, startingMonth=1), datetime(2008, 2, 29), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2007, 3, 30), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2007, 3, 31), False),
             #(QuarterEnd(1, startingMonth=1), datetime(2008, 4, 30), True),
             #(QuarterEnd(1, startingMonth=2), datetime(2008, 5, 30), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 6, 29), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 6, 30), False),

             #(QuarterEnd(1, startingMonth=2), datetime(2008, 1, 31), False),
             #(QuarterEnd(1, startingMonth=2), datetime(2007, 12, 31), False),
             #(QuarterEnd(1, startingMonth=2), datetime(2008, 2, 29), True),
             #(QuarterEnd(1, startingMonth=3), datetime(2007, 3, 30), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2007, 3, 31), False),
             #(QuarterEnd(1, startingMonth=2), datetime(2008, 4, 30), False),
             #(QuarterEnd(1, startingMonth=2), datetime(2008, 5, 30), True),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 6, 29), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 6, 30), False),

             #(QuarterEnd(1, startingMonth=3), datetime(2008, 1, 31), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2007, 12, 31), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 2, 29), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2007, 3, 30), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2007, 3, 31), True),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 4, 30), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 5, 30), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 6, 29), False),
             #(QuarterEnd(1, startingMonth=3), datetime(2008, 6, 30), True),
         ]

    for offset, date, expected in tests:
        assertOnOffset(offset, date, expected)

def assertEq(dateOffset, baseDate, expected):
    actual = dateOffset + baseDate
    assert actual == expected


def test_getOffsetName():
    try:
        datetools.getOffsetName(BDay(2))
    except Exception:
        pass
    else:
        raise Exception('failure')

    assert datetools.getOffsetName(BDay()) == 'WEEKDAY'
    assert datetools.getOffsetName(BMonthEnd()) == 'EOM'
    assert datetools.getOffsetName(Week(dayOfWeek=0)) == 'W@MON'
    assert datetools.getOffsetName(Week(dayOfWeek=1)) == 'W@TUE'
    assert datetools.getOffsetName(Week(dayOfWeek=2)) == 'W@WED'
    assert datetools.getOffsetName(Week(dayOfWeek=3)) == 'W@THU'
    assert datetools.getOffsetName(Week(dayOfWeek=4)) == 'W@FRI'


####
## XDateRange Tests
####
def eqXDateRange(kwargs, expected):
    actual = list(XDateRange(**kwargs))
    assert actual == expected

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



# DateRange test

def testDateRange1():
    toDate = datetime(2009, 5, 13)
    dr = DateRange(toDate=toDate, periods=20)
    firstDate = toDate - 19 * bday

    assert len(dr) == 20
    assert dr[0] == firstDate
    assert dr[-1] == toDate
