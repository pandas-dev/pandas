from pandas.core.datetools import BDay
from datetime import datetime
from pandas.core.datetools import *
from pandas.core.daterange import XDateRange, DateRange
import pandas.core.datetools as datetools

####
## Misc function tests
####
def testFormat():
    actual = datetools.format(datetime(2008, 1, 15))
    assert actual == '20080115'

def testOle2datetime():
    actual = datetools.ole2datetime(60000)
    assert actual == datetime(2064, 4, 8)

def testTto_datetime1():
    actual = datetools.to_datetime(datetime(2008, 1, 15))
    assert actual == datetime(2008, 1, 15)

def testTto_datetime2():
    actual = datetools.to_datetime('20080115')
    assert actual == datetime(2008, 1, 15)

def testNormalize_date():
    actual = datetools.normalize_date(datetime(2007, 10, 1, 1, 12, 5, 10))
    assert actual == datetime(2007, 10, 1)

#####
### DateOffset Tests
#####

def myAssert(actual, expected):
    assert actual == expected

def testEQ():
    myAssert(datetools.BDay(2), datetools.BDay(2))

def testHash():
    myAssert(datetools.BDay(2).__hash__(), datetools.BDay(2).__hash__())

def testCall():
    myAssert(BDay(2)(datetime(2008, 1, 1)), datetime(2008, 1, 3))

def testRAdd():
    myAssert(datetime(2008, 1, 1) + BDay(2), BDay(2) + datetime(2008, 1, 1))

def testSub():
    myAssert(datetime(2008, 1, 1) - BDay(2),  datetime(2008, 1, 1) + BDay(-2))

def testRSub():
    myAssert(datetime(2008, 1, 1) - BDay(2), BDay(2) - datetime(2008, 1, 1))

def testMult1():
    myAssert(datetime(2008, 1, 1) + 10*BDay(), datetime(2008, 1, 1) + BDay(10))

def testMult2():
    myAssert(datetime(2008, 1, 1) + (-5*BDay(-10)),
             datetime(2008, 1, 1) + BDay(50))


def testRollback1():
    myAssert(BDay(10).rollback(datetime(2008, 1, 1)), datetime(2008, 1, 1))

def testRollback2():
    myAssert(BDay(10).rollback(datetime(2008, 1, 5)), datetime(2008, 1, 4))

def testRollforward1():
    myAssert(BDay(10).rollforward(datetime(2008, 1, 1)), datetime(2008, 1, 1))

def testRollforward2():
    myAssert(BDay(10).rollforward(datetime(2008, 1, 5)), datetime(2008, 1, 7))

def assertOnOffset(offset, date, expected):
    actual = offset.onOffset(date)
    assert actual == expected

def testOnOffset():

    tests = [(BDay(), datetime(2008, 1, 1), True),
             (BDay(), datetime(2008, 1, 5), False),

             (BMonthEnd(), datetime(2007, 12, 31), True),
             (BMonthEnd(), datetime(2008, 1, 1), False),

             (BQuarterEnd(1, startingMonth=1), datetime(2008, 1, 31), True),
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

             (BYearEnd(), datetime(2007, 12, 31), True),
             (BYearEnd(), datetime(2008, 1, 1), False),
             (BYearEnd(), datetime(2006, 12, 31), False),
             (BYearEnd(), datetime(2006, 12, 29), True),

             (MonthEnd(), datetime(2007, 3, 30), False),
             (MonthEnd(), datetime(2007, 3, 31), True),

             #(QuarterEnd(1, startingMonth=1), datetime(2008, 1, 31), True),
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

             (datetools.Week(dayOfWeek=0), datetime(2008, 1, 1), False),
             (datetools.Week(dayOfWeek=0), datetime(2008, 1, 2), False),
             (datetools.Week(dayOfWeek=0), datetime(2008, 1, 3), False),
             (datetools.Week(dayOfWeek=0), datetime(2008, 1, 4), False),
             (datetools.Week(dayOfWeek=0), datetime(2008, 1, 5), False),
             (datetools.Week(dayOfWeek=0), datetime(2008, 1, 6), False),
             (datetools.Week(dayOfWeek=0), datetime(2008, 1, 7), True),

             (datetools.Week(dayOfWeek=1), datetime(2008, 1, 1), True),
             (datetools.Week(dayOfWeek=1), datetime(2008, 1, 2), False),
             (datetools.Week(dayOfWeek=1), datetime(2008, 1, 3), False),
             (datetools.Week(dayOfWeek=1), datetime(2008, 1, 4), False),
             (datetools.Week(dayOfWeek=1), datetime(2008, 1, 5), False),
             (datetools.Week(dayOfWeek=1), datetime(2008, 1, 6), False),
             (datetools.Week(dayOfWeek=1), datetime(2008, 1, 7), False),

             (datetools.Week(dayOfWeek=2), datetime(2008, 1, 1), False),
             (datetools.Week(dayOfWeek=2), datetime(2008, 1, 2), True),
             (datetools.Week(dayOfWeek=2), datetime(2008, 1, 3), False),
             (datetools.Week(dayOfWeek=2), datetime(2008, 1, 4), False),
             (datetools.Week(dayOfWeek=2), datetime(2008, 1, 5), False),
             (datetools.Week(dayOfWeek=2), datetime(2008, 1, 6), False),
             (datetools.Week(dayOfWeek=2), datetime(2008, 1, 7), False),

             (datetools.Week(dayOfWeek=3), datetime(2008, 1, 1), False),
             (datetools.Week(dayOfWeek=3), datetime(2008, 1, 2), False),
             (datetools.Week(dayOfWeek=3), datetime(2008, 1, 3), True),
             (datetools.Week(dayOfWeek=3), datetime(2008, 1, 4), False),
             (datetools.Week(dayOfWeek=3), datetime(2008, 1, 5), False),
             (datetools.Week(dayOfWeek=3), datetime(2008, 1, 6), False),
             (datetools.Week(dayOfWeek=3), datetime(2008, 1, 7), False),

             (datetools.Week(dayOfWeek=4), datetime(2008, 1, 1), False),
             (datetools.Week(dayOfWeek=4), datetime(2008, 1, 2), False),
             (datetools.Week(dayOfWeek=4), datetime(2008, 1, 3), False),
             (datetools.Week(dayOfWeek=4), datetime(2008, 1, 4), True),
             (datetools.Week(dayOfWeek=4), datetime(2008, 1, 5), False),
             (datetools.Week(dayOfWeek=4), datetime(2008, 1, 6), False),
             (datetools.Week(dayOfWeek=4), datetime(2008, 1, 7), False),

             (datetools.Week(dayOfWeek=5), datetime(2008, 1, 1), False),
             (datetools.Week(dayOfWeek=5), datetime(2008, 1, 2), False),
             (datetools.Week(dayOfWeek=5), datetime(2008, 1, 3), False),
             (datetools.Week(dayOfWeek=5), datetime(2008, 1, 4), False),
             (datetools.Week(dayOfWeek=5), datetime(2008, 1, 5), True),
             (datetools.Week(dayOfWeek=5), datetime(2008, 1, 6), False),
             (datetools.Week(dayOfWeek=5), datetime(2008, 1, 7), False),

             (datetools.Week(dayOfWeek=6), datetime(2008, 1, 1), False),
             (datetools.Week(dayOfWeek=6), datetime(2008, 1, 2), False),
             (datetools.Week(dayOfWeek=6), datetime(2008, 1, 3), False),
             (datetools.Week(dayOfWeek=6), datetime(2008, 1, 4), False),
             (datetools.Week(dayOfWeek=6), datetime(2008, 1, 5), False),
             (datetools.Week(dayOfWeek=6), datetime(2008, 1, 6), True),
             (datetools.Week(dayOfWeek=6), datetime(2008, 1, 7), False),
         ]

    for offset, date, expected in tests:
        assertOnOffset(offset, date, expected)

def assertEq(dateOffset, baseDate, expected):
    actual = dateOffset + baseDate
    assert actual == expected

def testBday():
    tests = []

    tests.append((datetools.bday,
                  {datetime(2008, 1, 1): datetime(2008, 1, 2),
                   datetime(2008, 1, 4): datetime(2008, 1, 7),
                   datetime(2008, 1, 5): datetime(2008, 1, 7),
                   datetime(2008, 1, 6): datetime(2008, 1, 7),
                   datetime(2008, 1, 7): datetime(2008, 1, 8)}))

    tests.append((2*datetools.bday,
                  {datetime(2008, 1, 1): datetime(2008, 1, 3),
                   datetime(2008, 1, 4): datetime(2008, 1, 8),
                   datetime(2008, 1, 5): datetime(2008, 1, 8),
                   datetime(2008, 1, 6): datetime(2008, 1, 8),
                   datetime(2008, 1, 7): datetime(2008, 1, 9)}))

    tests.append((-datetools.bday,
                  {datetime(2008, 1, 1): datetime(2007, 12, 31),
                   datetime(2008, 1, 4): datetime(2008, 1, 3),
                   datetime(2008, 1, 5): datetime(2008, 1, 4),
                   datetime(2008, 1, 6): datetime(2008, 1, 4),
                   datetime(2008, 1, 7): datetime(2008, 1, 4),
                   datetime(2008, 1, 8): datetime(2008, 1, 7)}))

    tests.append((-2*datetools.bday,
                  {datetime(2008, 1, 1): datetime(2007, 12, 28),
                   datetime(2008, 1, 4): datetime(2008, 1, 2),
                   datetime(2008, 1, 5): datetime(2008, 1, 3),
                   datetime(2008, 1, 6): datetime(2008, 1, 3),
                   datetime(2008, 1, 7): datetime(2008, 1, 3),
                   datetime(2008, 1, 8): datetime(2008, 1, 4),
                   datetime(2008, 1, 9): datetime(2008, 1, 7)}))

    tests.append((datetools.BDay(0),
                  {datetime(2008, 1, 1): datetime(2008, 1, 1),
                   datetime(2008, 1, 4): datetime(2008, 1, 4),
                   datetime(2008, 1, 5): datetime(2008, 1, 7),
                   datetime(2008, 1, 6): datetime(2008, 1, 7),
                   datetime(2008, 1, 7): datetime(2008, 1, 7)}))

    for dateOffset, cases in tests:
        for baseDate, expected in cases.iteritems():
            assertEq(dateOffset, baseDate, expected)

def testWeek():
    tests = []

    tests.append((datetools.week, # not business week
                  {datetime(2008, 1, 1): datetime(2008, 1, 8),
                   datetime(2008, 1, 4): datetime(2008, 1, 11),
                   datetime(2008, 1, 5): datetime(2008, 1, 12),
                   datetime(2008, 1, 6): datetime(2008, 1, 13),
                   datetime(2008, 1, 7): datetime(2008, 1, 14)}))

    tests.append((datetools.Week(dayOfWeek=0), # Mon
                  {datetime(2007, 12, 31): datetime(2008, 1, 7),
                   datetime(2008, 1, 4): datetime(2008, 1, 7),
                   datetime(2008, 1, 5): datetime(2008, 1, 7),
                   datetime(2008, 1, 6): datetime(2008, 1, 7),
                   datetime(2008, 1, 7): datetime(2008, 1, 14)}))

    tests.append((datetools.Week(0, dayOfWeek=0), # n=0 -> roll forward. Mon
                  {datetime(2007, 12, 31): datetime(2007, 12, 31),
                   datetime(2008, 1, 4): datetime(2008, 1, 7),
                   datetime(2008, 1, 5): datetime(2008, 1, 7),
                   datetime(2008, 1, 6): datetime(2008, 1, 7),
                   datetime(2008, 1, 7): datetime(2008, 1, 7)}))

    for dateOffset, cases in tests:
        for baseDate, expected in cases.iteritems():
            assertEq(dateOffset, baseDate, expected)

def testBMonthEnd():
    tests = []

    tests.append((datetools.BMonthEnd(),
                 {datetime(2008, 1, 1): datetime(2008, 1, 31),
                  datetime(2008, 1, 31): datetime(2008, 2, 29),
                  datetime(2006, 12, 29): datetime(2007, 1, 31),
                  datetime(2006, 12, 31): datetime(2007, 1, 31),
                  datetime(2007, 1, 1): datetime(2007, 1, 31),
                  datetime(2006, 12, 1): datetime(2006, 12, 29)}))

    tests.append((datetools.BMonthEnd(0),
                  {datetime(2008, 1, 1): datetime(2008, 1, 31),
                   datetime(2008, 1, 31): datetime(2008, 1, 31),
                   datetime(2006, 12, 29): datetime(2006, 12, 29),
                   datetime(2006, 12, 31): datetime(2007, 1, 31),
                   datetime(2007, 1, 1): datetime(2007, 1, 31)}))

    tests.append((datetools.BMonthEnd(2),
                 {datetime(2008, 1, 1): datetime(2008, 2, 29),
                  datetime(2008, 1, 31): datetime(2008, 3, 31),
                  datetime(2006, 12, 29): datetime(2007, 2, 28),
                  datetime(2006, 12, 31): datetime(2007, 2, 28),
                  datetime(2007, 1, 1): datetime(2007, 2, 28),
                  datetime(2006, 11, 1): datetime(2006, 12, 29)}))

    tests.append((datetools.BMonthEnd(-1),
                 {datetime(2007, 1, 1): datetime(2006, 12, 29),
                  datetime(2008, 6, 30): datetime(2008, 5, 30),
                  datetime(2008, 12, 31): datetime(2008, 11, 28),
                  datetime(2006, 12, 29): datetime(2006, 11, 30),
                  datetime(2006, 12, 30): datetime(2006, 12, 29),
                  datetime(2007, 1, 1): datetime(2006, 12, 29)}))

    for dateOffset, cases in tests:
        for baseDate, expected in cases.iteritems():
            assertEq(dateOffset, baseDate, expected)


def testBYearEnd():
    tests = []

    tests.append((datetools.BYearEnd(),
                  {datetime(2008, 1, 1): datetime(2008, 12, 31),
                   datetime(2008, 6, 30): datetime(2008, 12, 31),
                   datetime(2008, 12, 31): datetime(2009, 12, 31),
                   datetime(2005, 12, 30): datetime(2006, 12, 29),
                   datetime(2005, 12, 31): datetime(2006, 12, 29),}))

    tests.append((datetools.BYearEnd(0),
                  {datetime(2008, 1, 1): datetime(2008, 12, 31),
                   datetime(2008, 6, 30): datetime(2008, 12, 31),
                   datetime(2008, 12, 31): datetime(2008, 12, 31),
                   datetime(2005, 12, 31): datetime(2006, 12, 29),}))

    tests.append((datetools.BYearEnd(-1),
                  {datetime(2007, 1, 1): datetime(2006, 12, 29),
                   datetime(2008, 6, 30): datetime(2007, 12, 31),
                   datetime(2008, 12, 31): datetime(2007, 12, 31),
                   datetime(2006, 12, 29): datetime(2005, 12, 30),
                   datetime(2006, 12, 30): datetime(2006, 12, 29),
                   datetime(2007, 1, 1): datetime(2006, 12, 29),}))

    tests.append((datetools.BYearEnd(-2),
                  {datetime(2007, 1, 1): datetime(2005, 12, 30),
                   datetime(2008, 6, 30): datetime(2006, 12, 29),
                   datetime(2008, 12, 31): datetime(2006, 12, 29),}))

    for dateOffset, cases in tests:
        for baseDate, expected in cases.iteritems():
            assertEq(dateOffset, baseDate, expected)

def testYearBegin():
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


def testBQuarterEnd():
    tests = []

    tests.append((datetools.BQuarterEnd(),
                  {datetime(2008, 1, 1): datetime(2008, 3, 31),
                   datetime(2008, 1, 31): datetime(2008, 3, 31),
                   datetime(2008, 2, 15): datetime(2008, 3, 31),
                   datetime(2008, 2, 29): datetime(2008, 3, 31),
                   datetime(2008, 3, 15): datetime(2008, 3, 31),
                   datetime(2008, 3, 31): datetime(2008, 6, 30),
                   datetime(2008, 4, 15): datetime(2008, 6, 30),
                   datetime(2008, 4, 30): datetime(2008, 6, 30),}))

    tests.append((datetools.BQuarterEnd(n = 0),
                  {datetime(2008, 1, 1): datetime(2008, 3, 31),
                   datetime(2008, 1, 31): datetime(2008, 3, 31),
                   datetime(2008, 2, 15): datetime(2008, 3, 31),
                   datetime(2008, 2, 29): datetime(2008, 3, 31),
                   datetime(2008, 3, 15): datetime(2008, 3, 31),
                   datetime(2008, 3, 31): datetime(2008, 3, 31),
                   datetime(2008, 4, 15): datetime(2008, 6, 30),
                   datetime(2008, 4, 30): datetime(2008, 6, 30),}))

    tests.append((datetools.BQuarterEnd(n = -1),
                  {datetime(2008, 1, 1): datetime(2007, 12, 31),
                   datetime(2008, 1, 31): datetime(2007, 12, 31),
                   datetime(2008, 2, 15): datetime(2007, 12, 31),
                   datetime(2008, 2, 29): datetime(2007, 12, 31),
                   datetime(2008, 3, 15): datetime(2007, 12, 31),
                   datetime(2008, 3, 31): datetime(2007, 12, 31),
                   datetime(2008, 4, 15): datetime(2008, 3, 31),
                  datetime(2008, 4, 30): datetime(2008, 3, 31),}))

    tests.append((datetools.BQuarterEnd(n = 2),
                  {datetime(2008, 1, 1): datetime(2008, 6, 30),
                   datetime(2008, 1, 31): datetime(2008, 6, 30),
                   datetime(2008, 2, 15): datetime(2008, 6, 30),
                   datetime(2008, 2, 29): datetime(2008, 6, 30),
                   datetime(2008, 3, 15): datetime(2008, 6, 30),
                   datetime(2008, 3, 31): datetime(2008, 9, 30),
                   datetime(2008, 4, 15): datetime(2008, 9, 30),
                   datetime(2008, 4, 30): datetime(2008, 9, 30),}))

    for dateOffset, cases in tests:
        for baseDate, expected in cases.iteritems():
            assertEq(dateOffset, baseDate, expected)

def testBQuarterEndOffsets():
    tests = []

    tests.append((datetools.BQuarterEnd(startingMonth=1),
                  {datetime(2008, 1, 1): datetime(2008, 1, 31),
                   datetime(2008, 1, 31): datetime(2008, 4, 30),
                   datetime(2008, 2, 15): datetime(2008, 4, 30),
                   datetime(2008, 2, 29): datetime(2008, 4, 30),
                   datetime(2008, 3, 15): datetime(2008, 4, 30),
                   datetime(2008, 3, 31): datetime(2008, 4, 30),
                   datetime(2008, 4, 15): datetime(2008, 4, 30),
                   datetime(2008, 4, 30): datetime(2008, 7, 31),}))

    tests.append((datetools.BQuarterEnd(startingMonth=2),
                  {datetime(2008, 1, 1): datetime(2008, 2, 29),
                   datetime(2008, 1, 31): datetime(2008, 2, 29),
                   datetime(2008, 2, 15): datetime(2008, 2, 29),
                   datetime(2008, 2, 29): datetime(2008, 5, 30),
                   datetime(2008, 3, 15): datetime(2008, 5, 30),
                   datetime(2008, 3, 31): datetime(2008, 5, 30),
                   datetime(2008, 4, 15): datetime(2008, 5, 30),
                   datetime(2008, 4, 30): datetime(2008, 5, 30),}))

    tests.append((datetools.BQuarterEnd(startingMonth=1, n=0),
                  {datetime(2008, 1, 1): datetime(2008, 1, 31),
                   datetime(2008, 1, 31): datetime(2008, 1, 31),
                   datetime(2008, 2, 15): datetime(2008, 4, 30),
                   datetime(2008, 2, 29): datetime(2008, 4, 30),
                   datetime(2008, 3, 15): datetime(2008, 4, 30),
                   datetime(2008, 3, 31): datetime(2008, 4, 30),
                   datetime(2008, 4, 15): datetime(2008, 4, 30),
                   datetime(2008, 4, 30): datetime(2008, 4, 30),}))

    tests.append((datetools.BQuarterEnd(startingMonth=1, n=-1),
                  {datetime(2008, 1, 1): datetime(2007, 10, 31),
                   datetime(2008, 1, 31): datetime(2007, 10, 31),
                   datetime(2008, 2, 15): datetime(2008, 1, 31),
                   datetime(2008, 2, 29): datetime(2008, 1, 31),
                   datetime(2008, 3, 15): datetime(2008, 1, 31),
                   datetime(2008, 3, 31): datetime(2008, 1, 31),
                   datetime(2008, 4, 15): datetime(2008, 1, 31),
                   datetime(2008, 4, 30): datetime(2008, 1, 31),}))

    tests.append((datetools.BQuarterEnd(startingMonth=1, n=2),
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

def assertEqual(a, b):
    actual = dateOffset + baseDate
    assert actual == expected

def testDateRange1():
    toDate = datetime(2009, 5, 13)
    dr = DateRange(toDate=toDate, periods=20)
    firstDate = toDate - 19 * datetools.bday

    assert len(dr) == 20
    assert dr[0] == firstDate
    assert dr[-1] == toDate
