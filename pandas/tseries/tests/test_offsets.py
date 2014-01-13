from datetime import date, datetime, timedelta
from dateutil.relativedelta import relativedelta
from pandas.compat import range
from pandas import compat
import nose
from nose.tools import assert_raises

import numpy as np

from pandas.core.datetools import (
    bday, BDay, cday, CDay, BQuarterEnd, BMonthEnd, BYearEnd, MonthEnd,
    MonthBegin, BYearBegin, QuarterBegin, BQuarterBegin, BMonthBegin,
    DateOffset, Week, YearBegin, YearEnd, Hour, Minute, Second, Day, Micro,
    Milli, Nano,
    WeekOfMonth, format, ole2datetime, QuarterEnd, to_datetime, normalize_date,
    get_offset, get_offset_name, get_standard_freq)

from pandas.tseries.frequencies import _offset_map
from pandas.tseries.index import _to_m8, DatetimeIndex, _daterange_cache
from pandas.tseries.tools import parse_time_string
import pandas.tseries.offsets as offsets

from pandas.tslib import monthrange, OutOfBoundsDatetime
from pandas.lib import Timestamp
from pandas.util.testing import assertRaisesRegexp
import pandas.util.testing as tm
from pandas.tseries.offsets import BusinessMonthEnd, CacheableOffset, \
    LastWeekOfMonth, FY5253, FY5253Quarter, WeekDay

from pandas import _np_version_under1p7

_multiprocess_can_split_ = True


def test_monthrange():
    import calendar
    for y in range(2000, 2013):
        for m in range(1, 13):
            assert monthrange(y, m) == calendar.monthrange(y, m)


def _skip_if_no_cday():
    if cday is None:
        raise nose.SkipTest("CustomBusinessDay not available.")


####
## Misc function tests
####


def test_format():
    actual = format(datetime(2008, 1, 15))
    assert actual == '20080115'


def test_ole2datetime():
    actual = ole2datetime(60000)
    assert actual == datetime(2064, 4, 8)

    assert_raises(ValueError, ole2datetime, 60)


def test_to_datetime1():
    actual = to_datetime(datetime(2008, 1, 15))
    assert actual == datetime(2008, 1, 15)

    actual = to_datetime('20080115')
    assert actual == datetime(2008, 1, 15)

    # unparseable
    s = 'Month 1, 1999'
    assert to_datetime(s) == s


def test_normalize_date():
    actual = normalize_date(datetime(2007, 10, 1, 1, 12, 5, 10))
    assert actual == datetime(2007, 10, 1)


def test_to_m8():
    valb = datetime(2007, 10, 1)
    valu = _to_m8(valb)
    tm.assert_isinstance(valu, np.datetime64)
    # assert valu == np.datetime64(datetime(2007,10,1))

# def test_datetime64_box():
#    valu = np.datetime64(datetime(2007,10,1))
#    valb = _dt_box(valu)
#    assert type(valb) == datetime
#    assert valb == datetime(2007,10,1)

#####
### DateOffset Tests
#####

class TestBase(tm.TestCase):
    _offset = None

    def test_apply_out_of_range(self):
        if self._offset is None:
            raise nose.SkipTest("_offset not defined")

        # try to create an out-of-bounds result timestamp; if we can't create the offset
        # skip
        try:
            offset = self._offset(10000)

            result = Timestamp('20080101') + offset
            self.assert_(isinstance(result, datetime))
        except (OutOfBoundsDatetime):
            raise
        except (ValueError, KeyError):
            raise nose.SkipTest("cannot create out_of_range offset")

class TestDateOffset(TestBase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.d = Timestamp(datetime(2008, 1, 2))
        _offset_map.clear()

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

    def test_eq(self):
        offset1 = DateOffset(days=1)
        offset2 = DateOffset(days=365)

        self.assert_(offset1 != offset2)
        self.assert_(not (offset1 == offset2))


class TestBusinessDay(TestBase):
    _multiprocess_can_split_ = True
    _offset = BDay

    def setUp(self):
        self.d = datetime(2008, 1, 1)

        self.offset = BDay()
        self.offset2 = BDay(2)

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = BDay()
        offset2 = BDay()
        offset2.normalize = True
        self.assertEqual(offset, offset2)

    def test_repr(self):
        self.assertEqual(repr(self.offset), '<BusinessDay>')
        assert repr(self.offset2) == '<2 * BusinessDays>'

        expected = '<BusinessDay: offset=datetime.timedelta(1)>'
        assert repr(self.offset + timedelta(1)) == expected

    def test_with_offset(self):
        offset = self.offset + timedelta(hours=2)

        assert (self.d + offset) == datetime(2008, 1, 2, 2)

    def testEQ(self):
        self.assertEqual(self.offset2, self.offset2)

    def test_mul(self):
        pass

    def test_hash(self):
        self.assertEqual(hash(self.offset2), hash(self.offset2))

    def testCall(self):
        self.assertEqual(self.offset2(self.d), datetime(2008, 1, 3))

    def testRAdd(self):
        self.assertEqual(self.d + self.offset2, self.offset2 + self.d)

    def testSub(self):
        off = self.offset2
        self.assertRaises(Exception, off.__sub__, self.d)
        self.assertEqual(2 * off - off, off)

        self.assertEqual(self.d - self.offset2, self.d + BDay(-2))

    def testRSub(self):
        self.assertEqual(self.d - self.offset2, (-self.offset2).apply(self.d))

    def testMult1(self):
        self.assertEqual(self.d + 10 * self.offset, self.d + BDay(10))

    def testMult2(self):
        self.assertEqual(self.d + (-5 * BDay(-10)),
                         self.d + BDay(50))

    def testRollback1(self):
        self.assertEqual(BDay(10).rollback(self.d), self.d)

    def testRollback2(self):
        self.assertEqual(
            BDay(10).rollback(datetime(2008, 1, 5)), datetime(2008, 1, 4))

    def testRollforward1(self):
        self.assertEqual(BDay(10).rollforward(self.d), self.d)

    def testRollforward2(self):
        self.assertEqual(
            BDay(10).rollforward(datetime(2008, 1, 5)), datetime(2008, 1, 7))

    def test_roll_date_object(self):
        offset = BDay()

        dt = date(2012, 9, 15)

        result = offset.rollback(dt)
        self.assertEqual(result, datetime(2012, 9, 14))

        result = offset.rollforward(dt)
        self.assertEqual(result, datetime(2012, 9, 17))

        offset = offsets.Day()
        result = offset.rollback(dt)
        self.assertEqual(result, datetime(2012, 9, 15))

        result = offset.rollforward(dt)
        self.assertEqual(result, datetime(2012, 9, 15))

    def test_onOffset(self):
        tests = [(BDay(), datetime(2008, 1, 1), True),
                 (BDay(), datetime(2008, 1, 5), False)]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_apply(self):
        tests = []

        tests.append((bday,
                      {datetime(2008, 1, 1): datetime(2008, 1, 2),
                       datetime(2008, 1, 4): datetime(2008, 1, 7),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 8)}))

        tests.append((2 * bday,
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

        tests.append((-2 * bday,
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

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_apply_large_n(self):
        dt = datetime(2012, 10, 23)

        result = dt + BDay(10)
        self.assertEqual(result, datetime(2012, 11, 6))

        result = dt + BDay(100) - BDay(100)
        self.assertEqual(result, dt)

        off = BDay() * 6
        rs = datetime(2012, 1, 1) - off
        xp = datetime(2011, 12, 23)
        self.assertEqual(rs, xp)

        st = datetime(2011, 12, 18)
        rs = st + off
        xp = datetime(2011, 12, 26)
        self.assertEqual(rs, xp)
        
        off = BDay() * 10
        rs = datetime(2014, 1, 5) + off # see #5890
        xp = datetime(2014, 1, 17)
        self.assertEqual(rs, xp)

    def test_apply_corner(self):
        self.assertRaises(TypeError, BDay().apply, BMonthEnd())

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = BDay()
        offset2 = BDay()
        self.assertFalse(offset1 != offset2)


class TestCustomBusinessDay(TestBase):
    _multiprocess_can_split_ = True
    _offset = CDay

    def setUp(self):
        self.d = datetime(2008, 1, 1)

        _skip_if_no_cday()
        self.offset = CDay()
        self.offset2 = CDay(2)

    def test_different_normalize_equals(self):
        # equivalent in this special case
        offset = CDay()
        offset2 = CDay()
        offset2.normalize = True
        self.assertEqual(offset, offset2)

    def test_repr(self):
        assert repr(self.offset) == '<CustomBusinessDay>'
        assert repr(self.offset2) == '<2 * CustomBusinessDays>'

        expected = '<BusinessDay: offset=datetime.timedelta(1)>'
        assert repr(self.offset + timedelta(1)) == expected

    def test_with_offset(self):
        offset = self.offset + timedelta(hours=2)

        assert (self.d + offset) == datetime(2008, 1, 2, 2)

    def testEQ(self):
        self.assertEqual(self.offset2, self.offset2)

    def test_mul(self):
        pass

    def test_hash(self):
        self.assertEqual(hash(self.offset2), hash(self.offset2))

    def testCall(self):
        self.assertEqual(self.offset2(self.d), datetime(2008, 1, 3))

    def testRAdd(self):
        self.assertEqual(self.d + self.offset2, self.offset2 + self.d)

    def testSub(self):
        off = self.offset2
        self.assertRaises(Exception, off.__sub__, self.d)
        self.assertEqual(2 * off - off, off)

        self.assertEqual(self.d - self.offset2, self.d + CDay(-2))

    def testRSub(self):
        self.assertEqual(self.d - self.offset2, (-self.offset2).apply(self.d))

    def testMult1(self):
        self.assertEqual(self.d + 10 * self.offset, self.d + CDay(10))

    def testMult2(self):
        self.assertEqual(self.d + (-5 * CDay(-10)),
                         self.d + CDay(50))

    def testRollback1(self):
        self.assertEqual(CDay(10).rollback(self.d), self.d)

    def testRollback2(self):
        self.assertEqual(
            CDay(10).rollback(datetime(2008, 1, 5)), datetime(2008, 1, 4))

    def testRollforward1(self):
        self.assertEqual(CDay(10).rollforward(self.d), self.d)

    def testRollforward2(self):
        self.assertEqual(
            CDay(10).rollforward(datetime(2008, 1, 5)), datetime(2008, 1, 7))

    def test_roll_date_object(self):
        offset = CDay()

        dt = date(2012, 9, 15)

        result = offset.rollback(dt)
        self.assertEqual(result, datetime(2012, 9, 14))

        result = offset.rollforward(dt)
        self.assertEqual(result, datetime(2012, 9, 17))

        offset = offsets.Day()
        result = offset.rollback(dt)
        self.assertEqual(result, datetime(2012, 9, 15))

        result = offset.rollforward(dt)
        self.assertEqual(result, datetime(2012, 9, 15))

    def test_onOffset(self):
        tests = [(CDay(), datetime(2008, 1, 1), True),
                 (CDay(), datetime(2008, 1, 5), False)]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_apply(self):
        from pandas.core.datetools import cday
        tests = []

        tests.append((cday,
                      {datetime(2008, 1, 1): datetime(2008, 1, 2),
                       datetime(2008, 1, 4): datetime(2008, 1, 7),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 8)}))

        tests.append((2 * cday,
                      {datetime(2008, 1, 1): datetime(2008, 1, 3),
                       datetime(2008, 1, 4): datetime(2008, 1, 8),
                       datetime(2008, 1, 5): datetime(2008, 1, 8),
                       datetime(2008, 1, 6): datetime(2008, 1, 8),
                       datetime(2008, 1, 7): datetime(2008, 1, 9)}))

        tests.append((-cday,
                      {datetime(2008, 1, 1): datetime(2007, 12, 31),
                       datetime(2008, 1, 4): datetime(2008, 1, 3),
                       datetime(2008, 1, 5): datetime(2008, 1, 4),
                       datetime(2008, 1, 6): datetime(2008, 1, 4),
                       datetime(2008, 1, 7): datetime(2008, 1, 4),
                       datetime(2008, 1, 8): datetime(2008, 1, 7)}))

        tests.append((-2 * cday,
                      {datetime(2008, 1, 1): datetime(2007, 12, 28),
                       datetime(2008, 1, 4): datetime(2008, 1, 2),
                       datetime(2008, 1, 5): datetime(2008, 1, 3),
                       datetime(2008, 1, 6): datetime(2008, 1, 3),
                       datetime(2008, 1, 7): datetime(2008, 1, 3),
                       datetime(2008, 1, 8): datetime(2008, 1, 4),
                       datetime(2008, 1, 9): datetime(2008, 1, 7)}))

        tests.append((CDay(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 1, 4): datetime(2008, 1, 4),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 7)}))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_apply_large_n(self):
        dt = datetime(2012, 10, 23)

        result = dt + CDay(10)
        self.assertEqual(result, datetime(2012, 11, 6))

        result = dt + CDay(100) - CDay(100)
        self.assertEqual(result, dt)

        off = CDay() * 6
        rs = datetime(2012, 1, 1) - off
        xp = datetime(2011, 12, 23)
        self.assertEqual(rs, xp)

        st = datetime(2011, 12, 18)
        rs = st + off
        xp = datetime(2011, 12, 26)
        self.assertEqual(rs, xp)

    def test_apply_corner(self):
        self.assertRaises(Exception, CDay().apply, BMonthEnd())

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = CDay()
        offset2 = CDay()
        self.assertFalse(offset1 != offset2)

    def test_holidays(self):
        # Define a TradingDay offset
        holidays = ['2012-05-01', datetime(2013, 5, 1),
                    np.datetime64('2014-05-01')]
        tday = CDay(holidays=holidays)
        for year in range(2012, 2015):
            dt = datetime(year, 4, 30)
            xp = datetime(year, 5, 2)
            rs = dt + tday
            self.assertEqual(rs, xp)

    def test_weekmask(self):
        weekmask_saudi = 'Sat Sun Mon Tue Wed'  # Thu-Fri Weekend
        weekmask_uae = '1111001'                # Fri-Sat Weekend
        weekmask_egypt = [1,1,1,1,0,0,1]        # Fri-Sat Weekend
        bday_saudi = CDay(weekmask=weekmask_saudi)
        bday_uae = CDay(weekmask=weekmask_uae)
        bday_egypt = CDay(weekmask=weekmask_egypt)
        dt = datetime(2013, 5, 1)
        xp_saudi = datetime(2013, 5, 4)
        xp_uae = datetime(2013, 5, 2)
        xp_egypt = datetime(2013, 5, 2)
        self.assertEqual(xp_saudi, dt + bday_saudi)
        self.assertEqual(xp_uae, dt + bday_uae)
        self.assertEqual(xp_egypt, dt + bday_egypt)
        xp2 = datetime(2013, 5, 5)
        self.assertEqual(xp2, dt + 2 * bday_saudi)
        self.assertEqual(xp2, dt + 2 * bday_uae)
        self.assertEqual(xp2, dt + 2 * bday_egypt)

    def test_weekmask_and_holidays(self):
        weekmask_egypt = 'Sun Mon Tue Wed Thu'  # Fri-Sat Weekend
        holidays = ['2012-05-01', datetime(2013, 5, 1),
                    np.datetime64('2014-05-01')]
        bday_egypt = CDay(holidays=holidays, weekmask=weekmask_egypt)
        dt = datetime(2013, 4, 30)
        xp_egypt = datetime(2013, 5, 5)
        self.assertEqual(xp_egypt, dt + 2 * bday_egypt)


def assertOnOffset(offset, date, expected):
    actual = offset.onOffset(date)
    assert actual == expected, ("\nExpected: %s\nActual: %s\nFor Offset: %s)"
                                "\nAt Date: %s" %
                                (expected, actual, offset, date))


class TestWeek(TestBase):
    _offset = Week

    def test_repr(self):
        self.assertEqual(repr(Week(weekday=0)), "<Week: weekday=0>")
        self.assertEqual(repr(Week(n=-1, weekday=0)), "<-1 * Week: weekday=0>")
        self.assertEqual(repr(Week(n=-2, weekday=0)), "<-2 * Weeks: weekday=0>")

    def test_corner(self):
        self.assertRaises(ValueError, Week, weekday=7)
        assertRaisesRegexp(ValueError, "Day must be", Week, weekday=-1)

    def test_isAnchored(self):
        self.assert_(Week(weekday=0).isAnchored())
        self.assert_(not Week().isAnchored())
        self.assert_(not Week(2, weekday=2).isAnchored())
        self.assert_(not Week(2).isAnchored())

    def test_offset(self):
        tests = []

        tests.append((Week(),  # not business week
                      {datetime(2008, 1, 1): datetime(2008, 1, 8),
                       datetime(2008, 1, 4): datetime(2008, 1, 11),
                       datetime(2008, 1, 5): datetime(2008, 1, 12),
                       datetime(2008, 1, 6): datetime(2008, 1, 13),
                       datetime(2008, 1, 7): datetime(2008, 1, 14)}))

        tests.append((Week(weekday=0),  # Mon
                      {datetime(2007, 12, 31): datetime(2008, 1, 7),
                       datetime(2008, 1, 4): datetime(2008, 1, 7),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 14)}))

        tests.append((Week(0, weekday=0),  # n=0 -> roll forward. Mon
                      {datetime(2007, 12, 31): datetime(2007, 12, 31),
                       datetime(2008, 1, 4): datetime(2008, 1, 7),
                       datetime(2008, 1, 5): datetime(2008, 1, 7),
                       datetime(2008, 1, 6): datetime(2008, 1, 7),
                       datetime(2008, 1, 7): datetime(2008, 1, 7)}))

        tests.append((Week(-2, weekday=1),  # n=0 -> roll forward. Mon
                      {datetime(2010, 4, 6): datetime(2010, 3, 23),
                       datetime(2010, 4, 8): datetime(2010, 3, 30),
                       datetime(2010, 4, 5): datetime(2010, 3, 23)}))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_onOffset(self):
        for weekday in range(7):
            offset = Week(weekday=weekday)

            for day in range(1, 8):
                date = datetime(2008, 1, day)

                if day % 7 == weekday:
                    expected = True
                else:
                    expected = False
            assertOnOffset(offset, date, expected)

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = Week()
        offset2 = Week()
        self.assertFalse(offset1 != offset2)


class TestWeekOfMonth(TestBase):
    _offset = WeekOfMonth

    def test_constructor(self):
        assertRaisesRegexp(ValueError, "^N cannot be 0", WeekOfMonth, n=0, week=1, weekday=1)
        assertRaisesRegexp(ValueError, "^Week", WeekOfMonth, n=1, week=4, weekday=0)
        assertRaisesRegexp(ValueError, "^Week", WeekOfMonth, n=1, week=-1, weekday=0)
        assertRaisesRegexp(ValueError, "^Day", WeekOfMonth, n=1, week=0, weekday=-1)
        assertRaisesRegexp(ValueError, "^Day", WeekOfMonth, n=1, week=0, weekday=7)

    def test_repr(self):
        self.assertEqual(repr(WeekOfMonth(weekday=1,week=2)), "<WeekOfMonth: week=2, weekday=1>")

    def test_offset(self):
        date1 = datetime(2011, 1, 4)  # 1st Tuesday of Month
        date2 = datetime(2011, 1, 11)  # 2nd Tuesday of Month
        date3 = datetime(2011, 1, 18)  # 3rd Tuesday of Month
        date4 = datetime(2011, 1, 25)  # 4th Tuesday of Month

        # see for loop for structure
        test_cases = [
            (-2, 2, 1, date1, datetime(2010, 11, 16)),
            (-2, 2, 1, date2, datetime(2010, 11, 16)),
            (-2, 2, 1, date3, datetime(2010, 11, 16)),
            (-2, 2, 1, date4, datetime(2010, 12, 21)),

            (-1, 2, 1, date1, datetime(2010, 12, 21)),
            (-1, 2, 1, date2, datetime(2010, 12, 21)),
            (-1, 2, 1, date3, datetime(2010, 12, 21)),
            (-1, 2, 1, date4, datetime(2011, 1, 18)),

            (1, 0, 0, date1, datetime(2011, 2, 7)),
            (1, 0, 0, date2, datetime(2011, 2, 7)),
            (1, 0, 0, date3, datetime(2011, 2, 7)),
            (1, 0, 0, date4, datetime(2011, 2, 7)),
            (1, 0, 1, date1, datetime(2011, 2, 1)),
            (1, 0, 1, date2, datetime(2011, 2, 1)),
            (1, 0, 1, date3, datetime(2011, 2, 1)),
            (1, 0, 1, date4, datetime(2011, 2, 1)),
            (1, 0, 2, date1, datetime(2011, 1, 5)),
            (1, 0, 2, date2, datetime(2011, 2, 2)),
            (1, 0, 2, date3, datetime(2011, 2, 2)),
            (1, 0, 2, date4, datetime(2011, 2, 2)),

            (1, 2, 1, date1, datetime(2011, 1, 18)),
            (1, 2, 1, date2, datetime(2011, 1, 18)),
            (1, 2, 1, date3, datetime(2011, 2, 15)),
            (1, 2, 1, date4, datetime(2011, 2, 15)),

            (2, 2, 1, date1, datetime(2011, 2, 15)),
            (2, 2, 1, date2, datetime(2011, 2, 15)),
            (2, 2, 1, date3, datetime(2011, 3, 15)),
            (2, 2, 1, date4, datetime(2011, 3, 15)),
        ]

        for n, week, weekday, date, expected in test_cases:
            offset = WeekOfMonth(n, week=week, weekday=weekday)
            assertEq(offset, date, expected)

        # try subtracting
        result = datetime(2011, 2, 1) - WeekOfMonth(week=1, weekday=2)
        self.assertEqual(result, datetime(2011, 1, 12))
        result = datetime(2011, 2, 3) - WeekOfMonth(week=0, weekday=2)
        self.assertEqual(result, datetime(2011, 2, 2))

    def test_onOffset(self):
        test_cases = [
            (0, 0, datetime(2011, 2, 7), True),
            (0, 0, datetime(2011, 2, 6), False),
            (0, 0, datetime(2011, 2, 14), False),
            (1, 0, datetime(2011, 2, 14), True),
            (0, 1, datetime(2011, 2, 1), True),
            (0, 1, datetime(2011, 2, 8), False),
        ]

        for week, weekday, date, expected in test_cases:
            offset = WeekOfMonth(week=week, weekday=weekday)
            self.assert_(offset.onOffset(date) == expected)

class TestLastWeekOfMonth(TestBase):
    _offset = LastWeekOfMonth

    def test_constructor(self):
        assertRaisesRegexp(ValueError, "^N cannot be 0", \
                           LastWeekOfMonth, n=0, weekday=1)

        assertRaisesRegexp(ValueError, "^Day", LastWeekOfMonth, n=1, weekday=-1)
        assertRaisesRegexp(ValueError, "^Day", LastWeekOfMonth, n=1, weekday=7)

    def test_offset(self):
        #### Saturday
        last_sat = datetime(2013,8,31)
        next_sat = datetime(2013,9,28)
        offset_sat = LastWeekOfMonth(n=1, weekday=5)

        one_day_before = (last_sat + timedelta(days=-1))
        self.assert_(one_day_before + offset_sat == last_sat)

        one_day_after = (last_sat + timedelta(days=+1))
        self.assert_(one_day_after + offset_sat == next_sat)

        #Test On that day
        self.assert_(last_sat + offset_sat == next_sat)

        #### Thursday

        offset_thur = LastWeekOfMonth(n=1, weekday=3)
        last_thurs = datetime(2013,1,31)
        next_thurs = datetime(2013,2,28)

        one_day_before = last_thurs + timedelta(days=-1)
        self.assert_(one_day_before + offset_thur == last_thurs)

        one_day_after = last_thurs + timedelta(days=+1)
        self.assert_(one_day_after + offset_thur == next_thurs)

        # Test on that day
        self.assert_(last_thurs + offset_thur == next_thurs)

        three_before = last_thurs + timedelta(days=-3)
        self.assert_(three_before + offset_thur == last_thurs)

        two_after = last_thurs + timedelta(days=+2)
        self.assert_(two_after + offset_thur == next_thurs)

        offset_sunday = LastWeekOfMonth(n=1, weekday=WeekDay.SUN)
        self.assert_(datetime(2013,7,31) + offset_sunday == datetime(2013,8,25))

    def test_onOffset(self):
        test_cases = [
            (WeekDay.SUN, datetime(2013, 1, 27), True),
            (WeekDay.SAT, datetime(2013, 3, 30), True),
            (WeekDay.MON, datetime(2013, 2, 18), False), #Not the last Mon
            (WeekDay.SUN, datetime(2013, 2, 25), False), #Not a SUN
            (WeekDay.MON, datetime(2013, 2, 25), True),
            (WeekDay.SAT, datetime(2013, 11, 30), True),

            (WeekDay.SAT, datetime(2006, 8, 26), True),
            (WeekDay.SAT, datetime(2007, 8, 25), True),
            (WeekDay.SAT, datetime(2008, 8, 30), True),
            (WeekDay.SAT, datetime(2009, 8, 29), True),
            (WeekDay.SAT, datetime(2010, 8, 28), True),
            (WeekDay.SAT, datetime(2011, 8, 27), True),
            (WeekDay.SAT, datetime(2019, 8, 31), True),
        ]

        for weekday, date, expected in test_cases:
            offset = LastWeekOfMonth(weekday=weekday)
            self.assert_(offset.onOffset(date) == expected, date)


class TestBMonthBegin(TestBase):
    _offset = BMonthBegin

    def test_offset(self):
        tests = []

        tests.append((BMonthBegin(),
                     {datetime(2008, 1, 1): datetime(2008, 2, 1),
                      datetime(2008, 1, 31): datetime(2008, 2, 1),
                      datetime(2006, 12, 29): datetime(2007, 1, 1),
                      datetime(2006, 12, 31): datetime(2007, 1, 1),
                      datetime(2006, 9, 1): datetime(2006, 10, 2),
                      datetime(2007, 1, 1): datetime(2007, 2, 1),
                      datetime(2006, 12, 1): datetime(2007, 1, 1)}))

        tests.append((BMonthBegin(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2006, 10, 2): datetime(2006, 10, 2),
                       datetime(2008, 1, 31): datetime(2008, 2, 1),
                       datetime(2006, 12, 29): datetime(2007, 1, 1),
                       datetime(2006, 12, 31): datetime(2007, 1, 1),
                       datetime(2006, 9, 15): datetime(2006, 10, 2)}))

        tests.append((BMonthBegin(2),
                     {datetime(2008, 1, 1): datetime(2008, 3, 3),
                      datetime(2008, 1, 15): datetime(2008, 3, 3),
                      datetime(2006, 12, 29): datetime(2007, 2, 1),
                      datetime(2006, 12, 31): datetime(2007, 2, 1),
                      datetime(2007, 1, 1): datetime(2007, 3, 1),
                      datetime(2006, 11, 1): datetime(2007, 1, 1)}))

        tests.append((BMonthBegin(-1),
                     {datetime(2007, 1, 1): datetime(2006, 12, 1),
                      datetime(2008, 6, 30): datetime(2008, 6, 2),
                      datetime(2008, 6, 1): datetime(2008, 5, 1),
                      datetime(2008, 3, 10): datetime(2008, 3, 3),
                      datetime(2008, 12, 31): datetime(2008, 12, 1),
                      datetime(2006, 12, 29): datetime(2006, 12, 1),
                      datetime(2006, 12, 30): datetime(2006, 12, 1),
                      datetime(2007, 1, 1): datetime(2006, 12, 1)}))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_onOffset(self):

        tests = [(BMonthBegin(), datetime(2007, 12, 31), False),
                 (BMonthBegin(), datetime(2008, 1, 1), True),
                 (BMonthBegin(), datetime(2001, 4, 2), True),
                 (BMonthBegin(), datetime(2008, 3, 3), True)]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = BMonthBegin()
        offset2 = BMonthBegin()
        self.assertFalse(offset1 != offset2)


class TestBMonthEnd(TestBase):
    _offset = BMonthEnd

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

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_normalize(self):
        dt = datetime(2007, 1, 1, 3)

        result = dt + BMonthEnd()
        expected = dt.replace(hour=0) + BMonthEnd()
        self.assertEqual(result, expected)

    def test_onOffset(self):

        tests = [(BMonthEnd(), datetime(2007, 12, 31), True),
                 (BMonthEnd(), datetime(2008, 1, 1), False)]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_offsets_compare_equal(self):
        # root cause of #456
        offset1 = BMonthEnd()
        offset2 = BMonthEnd()
        self.assertFalse(offset1 != offset2)


class TestMonthBegin(TestBase):
    _offset = MonthBegin

    def test_offset(self):
        tests = []

        # NOTE: I'm not entirely happy with the logic here for Begin -ss
        # see thread 'offset conventions' on the ML
        tests.append((MonthBegin(),
                     {datetime(2008, 1, 31): datetime(2008, 2, 1),
                      datetime(2008, 2, 1): datetime(2008, 3, 1),
                      datetime(2006, 12, 31): datetime(2007, 1, 1),
                      datetime(2006, 12, 1): datetime(2007, 1, 1),
                      datetime(2007, 1, 31): datetime(2007, 2, 1)}))

        tests.append((MonthBegin(0),
                      {datetime(2008, 1, 31): datetime(2008, 2, 1),
                       datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2006, 12, 3): datetime(2007, 1, 1),
                       datetime(2007, 1, 31): datetime(2007, 2, 1)}))

        tests.append((MonthBegin(2),
                     {datetime(2008, 2, 29): datetime(2008, 4, 1),
                      datetime(2008, 1, 31): datetime(2008, 3, 1),
                      datetime(2006, 12, 31): datetime(2007, 2, 1),
                      datetime(2007, 12, 28): datetime(2008, 2, 1),
                      datetime(2007, 1, 1): datetime(2007, 3, 1),
                      datetime(2006, 11, 1): datetime(2007, 1, 1)}))

        tests.append((MonthBegin(-1),
                     {datetime(2007, 1, 1): datetime(2006, 12, 1),
                      datetime(2008, 5, 31): datetime(2008, 5, 1),
                      datetime(2008, 12, 31): datetime(2008, 12, 1),
                      datetime(2006, 12, 29): datetime(2006, 12, 1),
                      datetime(2006, 1, 2): datetime(2006, 1, 1)}))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)


class TestMonthEnd(TestBase):
    _offset = MonthEnd

    def test_offset(self):
        tests = []

        tests.append((MonthEnd(),
                     {datetime(2008, 1, 1): datetime(2008, 1, 31),
                      datetime(2008, 1, 31): datetime(2008, 2, 29),
                      datetime(2006, 12, 29): datetime(2006, 12, 31),
                      datetime(2006, 12, 31): datetime(2007, 1, 31),
                      datetime(2007, 1, 1): datetime(2007, 1, 31),
                      datetime(2006, 12, 1): datetime(2006, 12, 31)}))

        tests.append((MonthEnd(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 31),
                       datetime(2008, 1, 31): datetime(2008, 1, 31),
                       datetime(2006, 12, 29): datetime(2006, 12, 31),
                       datetime(2006, 12, 31): datetime(2006, 12, 31),
                       datetime(2007, 1, 1): datetime(2007, 1, 31)}))

        tests.append((MonthEnd(2),
                     {datetime(2008, 1, 1): datetime(2008, 2, 29),
                      datetime(2008, 1, 31): datetime(2008, 3, 31),
                      datetime(2006, 12, 29): datetime(2007, 1, 31),
                      datetime(2006, 12, 31): datetime(2007, 2, 28),
                      datetime(2007, 1, 1): datetime(2007, 2, 28),
                      datetime(2006, 11, 1): datetime(2006, 12, 31)}))

        tests.append((MonthEnd(-1),
                     {datetime(2007, 1, 1): datetime(2006, 12, 31),
                      datetime(2008, 6, 30): datetime(2008, 5, 31),
                      datetime(2008, 12, 31): datetime(2008, 11, 30),
                      datetime(2006, 12, 29): datetime(2006, 11, 30),
                      datetime(2006, 12, 30): datetime(2006, 11, 30),
                      datetime(2007, 1, 1): datetime(2006, 12, 31)}))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    # def test_day_of_month(self):
    #     dt = datetime(2007, 1, 1)

    #     offset = MonthEnd(day=20)

    #     result = dt + offset
    #     self.assertEqual(result, datetime(2007, 1, 20))

    #     result = result + offset
    #     self.assertEqual(result, datetime(2007, 2, 20))

    def test_normalize(self):
        dt = datetime(2007, 1, 1, 3)

        result = dt + MonthEnd()
        expected = dt.replace(hour=0) + MonthEnd()
        self.assertEqual(result, expected)

    def test_onOffset(self):

        tests = [(MonthEnd(), datetime(2007, 12, 31), True),
                 (MonthEnd(), datetime(2008, 1, 1), False)]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


class TestBQuarterBegin(TestBase):
    _offset = BQuarterBegin

    def test_repr(self):
        self.assertEqual(repr(BQuarterBegin()),"<BusinessQuarterBegin: startingMonth=3>")
        self.assertEqual(repr(BQuarterBegin(startingMonth=3)), "<BusinessQuarterBegin: startingMonth=3>")
        self.assertEqual(repr(BQuarterBegin(startingMonth=1)), "<BusinessQuarterBegin: startingMonth=1>")

    def test_isAnchored(self):
        self.assert_(BQuarterBegin(startingMonth=1).isAnchored())
        self.assert_(BQuarterBegin().isAnchored())
        self.assert_(not BQuarterBegin(2, startingMonth=1).isAnchored())

    def test_offset(self):
        tests = []

        tests.append((BQuarterBegin(startingMonth=1),
                      {datetime(2008, 1, 1): datetime(2008, 4, 1),
                       datetime(2008, 1, 31): datetime(2008, 4, 1),
                       datetime(2008, 2, 15): datetime(2008, 4, 1),
                       datetime(2008, 2, 29): datetime(2008, 4, 1),
                       datetime(2008, 3, 15): datetime(2008, 4, 1),
                       datetime(2008, 3, 31): datetime(2008, 4, 1),
                       datetime(2008, 4, 15): datetime(2008, 7, 1),
                       datetime(2007, 3, 15): datetime(2007, 4, 2),
                       datetime(2007, 2, 28): datetime(2007, 4, 2),
                       datetime(2007, 1, 1): datetime(2007, 4, 2),
                       datetime(2007, 4, 15): datetime(2007, 7, 2),
                       datetime(2007, 7, 1): datetime(2007, 7, 2),
                       datetime(2007, 4, 1): datetime(2007, 4, 2),
                       datetime(2007, 4, 2): datetime(2007, 7, 2),
                       datetime(2008, 4, 30): datetime(2008, 7, 1), }))

        tests.append((BQuarterBegin(startingMonth=2),
                      {datetime(2008, 1, 1): datetime(2008, 2, 1),
                       datetime(2008, 1, 31): datetime(2008, 2, 1),
                       datetime(2008, 1, 15): datetime(2008, 2, 1),
                       datetime(2008, 2, 29): datetime(2008, 5, 1),
                       datetime(2008, 3, 15): datetime(2008, 5, 1),
                       datetime(2008, 3, 31): datetime(2008, 5, 1),
                       datetime(2008, 4, 15): datetime(2008, 5, 1),
                       datetime(2008, 8, 15): datetime(2008, 11, 3),
                       datetime(2008, 9, 15): datetime(2008, 11, 3),
                       datetime(2008, 11, 1): datetime(2008, 11, 3),
                       datetime(2008, 4, 30): datetime(2008, 5, 1), }))

        tests.append((BQuarterBegin(startingMonth=1, n=0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2007, 12, 31): datetime(2008, 1, 1),
                       datetime(2008, 2, 15): datetime(2008, 4, 1),
                       datetime(2008, 2, 29): datetime(2008, 4, 1),
                       datetime(2008, 1, 15): datetime(2008, 4, 1),
                       datetime(2008, 2, 27): datetime(2008, 4, 1),
                       datetime(2008, 3, 15): datetime(2008, 4, 1),
                       datetime(2007, 4, 1): datetime(2007, 4, 2),
                       datetime(2007, 4, 2): datetime(2007, 4, 2),
                       datetime(2007, 7, 1): datetime(2007, 7, 2),
                       datetime(2007, 4, 15): datetime(2007, 7, 2),
                       datetime(2007, 7, 2): datetime(2007, 7, 2), }))

        tests.append((BQuarterBegin(startingMonth=1, n=-1),
                      {datetime(2008, 1, 1): datetime(2007, 10, 1),
                       datetime(2008, 1, 31): datetime(2008, 1, 1),
                       datetime(2008, 2, 15): datetime(2008, 1, 1),
                       datetime(2008, 2, 29): datetime(2008, 1, 1),
                       datetime(2008, 3, 15): datetime(2008, 1, 1),
                       datetime(2008, 3, 31): datetime(2008, 1, 1),
                       datetime(2008, 4, 15): datetime(2008, 4, 1),
                       datetime(2007, 7, 3): datetime(2007, 7, 2),
                       datetime(2007, 4, 3): datetime(2007, 4, 2),
                       datetime(2007, 7, 2): datetime(2007, 4, 2),
                       datetime(2008, 4, 1): datetime(2008, 1, 1), }))

        tests.append((BQuarterBegin(startingMonth=1, n=2),
                      {datetime(2008, 1, 1): datetime(2008, 7, 1),
                       datetime(2008, 1, 15): datetime(2008, 7, 1),
                       datetime(2008, 2, 29): datetime(2008, 7, 1),
                       datetime(2008, 3, 15): datetime(2008, 7, 1),
                       datetime(2007, 3, 31): datetime(2007, 7, 2),
                       datetime(2007, 4, 15): datetime(2007, 10, 1),
                       datetime(2008, 4, 30): datetime(2008, 10, 1), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

        # corner
        offset = BQuarterBegin(n=-1, startingMonth=1)
        self.assertEqual(datetime(2007, 4, 3) + offset, datetime(2007, 4, 2))


class TestBQuarterEnd(TestBase):
    _offset = BQuarterEnd

    def test_repr(self):
        self.assertEqual(repr(BQuarterEnd()),"<BusinessQuarterEnd: startingMonth=3>")
        self.assertEqual(repr(BQuarterEnd(startingMonth=3)), "<BusinessQuarterEnd: startingMonth=3>")
        self.assertEqual(repr(BQuarterEnd(startingMonth=1)), "<BusinessQuarterEnd: startingMonth=1>")

    def test_isAnchored(self):
        self.assert_(BQuarterEnd(startingMonth=1).isAnchored())
        self.assert_(BQuarterEnd().isAnchored())
        self.assert_(not BQuarterEnd(2, startingMonth=1).isAnchored())

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
                       datetime(2008, 4, 30): datetime(2008, 7, 31), }))

        tests.append((BQuarterEnd(startingMonth=2),
                      {datetime(2008, 1, 1): datetime(2008, 2, 29),
                       datetime(2008, 1, 31): datetime(2008, 2, 29),
                       datetime(2008, 2, 15): datetime(2008, 2, 29),
                       datetime(2008, 2, 29): datetime(2008, 5, 30),
                       datetime(2008, 3, 15): datetime(2008, 5, 30),
                       datetime(2008, 3, 31): datetime(2008, 5, 30),
                       datetime(2008, 4, 15): datetime(2008, 5, 30),
                       datetime(2008, 4, 30): datetime(2008, 5, 30), }))

        tests.append((BQuarterEnd(startingMonth=1, n=0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 31),
                       datetime(2008, 1, 31): datetime(2008, 1, 31),
                       datetime(2008, 2, 15): datetime(2008, 4, 30),
                       datetime(2008, 2, 29): datetime(2008, 4, 30),
                       datetime(2008, 3, 15): datetime(2008, 4, 30),
                       datetime(2008, 3, 31): datetime(2008, 4, 30),
                       datetime(2008, 4, 15): datetime(2008, 4, 30),
                       datetime(2008, 4, 30): datetime(2008, 4, 30), }))

        tests.append((BQuarterEnd(startingMonth=1, n=-1),
                      {datetime(2008, 1, 1): datetime(2007, 10, 31),
                       datetime(2008, 1, 31): datetime(2007, 10, 31),
                       datetime(2008, 2, 15): datetime(2008, 1, 31),
                       datetime(2008, 2, 29): datetime(2008, 1, 31),
                       datetime(2008, 3, 15): datetime(2008, 1, 31),
                       datetime(2008, 3, 31): datetime(2008, 1, 31),
                       datetime(2008, 4, 15): datetime(2008, 1, 31),
                       datetime(2008, 4, 30): datetime(2008, 1, 31), }))

        tests.append((BQuarterEnd(startingMonth=1, n=2),
                      {datetime(2008, 1, 31): datetime(2008, 7, 31),
                       datetime(2008, 2, 15): datetime(2008, 7, 31),
                       datetime(2008, 2, 29): datetime(2008, 7, 31),
                       datetime(2008, 3, 15): datetime(2008, 7, 31),
                       datetime(2008, 3, 31): datetime(2008, 7, 31),
                       datetime(2008, 4, 15): datetime(2008, 7, 31),
                       datetime(2008, 4, 30): datetime(2008, 10, 31), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

        # corner
        offset = BQuarterEnd(n=-1, startingMonth=1)
        self.assertEqual(datetime(2010, 1, 31) + offset, datetime(2010, 1, 29))

    def test_onOffset(self):

        tests = [
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
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

def makeFY5253LastOfMonthQuarter(*args, **kwds):
    return FY5253Quarter(*args, variation="last", **kwds)

def makeFY5253NearestEndMonthQuarter(*args, **kwds):
    return FY5253Quarter(*args, variation="nearest", **kwds)

def makeFY5253NearestEndMonth(*args, **kwds):
    return FY5253(*args, variation="nearest", **kwds)

def makeFY5253LastOfMonth(*args, **kwds):
    return FY5253(*args, variation="last", **kwds)

class TestFY5253LastOfMonth(TestBase):

    def test_onOffset(self):

        offset_lom_sat_aug = makeFY5253LastOfMonth(1, startingMonth=8, weekday=WeekDay.SAT)
        offset_lom_sat_sep = makeFY5253LastOfMonth(1, startingMonth=9, weekday=WeekDay.SAT)

        tests = [
            #From Wikipedia (see: http://en.wikipedia.org/wiki/4%E2%80%934%E2%80%935_calendar#Last_Saturday_of_the_month_at_fiscal_year_end)
            (offset_lom_sat_aug, datetime(2006, 8, 26), True),
            (offset_lom_sat_aug, datetime(2007, 8, 25), True),
            (offset_lom_sat_aug, datetime(2008, 8, 30), True),
            (offset_lom_sat_aug, datetime(2009, 8, 29), True),
            (offset_lom_sat_aug, datetime(2010, 8, 28), True),
            (offset_lom_sat_aug, datetime(2011, 8, 27), True),
            (offset_lom_sat_aug, datetime(2012, 8, 25), True),
            (offset_lom_sat_aug, datetime(2013, 8, 31), True),
            (offset_lom_sat_aug, datetime(2014, 8, 30), True),
            (offset_lom_sat_aug, datetime(2015, 8, 29), True),
            (offset_lom_sat_aug, datetime(2016, 8, 27), True),
            (offset_lom_sat_aug, datetime(2017, 8, 26), True),
            (offset_lom_sat_aug, datetime(2018, 8, 25), True),
            (offset_lom_sat_aug, datetime(2019, 8, 31), True),

            (offset_lom_sat_aug, datetime(2006, 8, 27), False),
            (offset_lom_sat_aug, datetime(2007, 8, 28), False),
            (offset_lom_sat_aug, datetime(2008, 8, 31), False),
            (offset_lom_sat_aug, datetime(2009, 8, 30), False),
            (offset_lom_sat_aug, datetime(2010, 8, 29), False),
            (offset_lom_sat_aug, datetime(2011, 8, 28), False),

            (offset_lom_sat_aug, datetime(2006, 8, 25), False),
            (offset_lom_sat_aug, datetime(2007, 8, 24), False),
            (offset_lom_sat_aug, datetime(2008, 8, 29), False),
            (offset_lom_sat_aug, datetime(2009, 8, 28), False),
            (offset_lom_sat_aug, datetime(2010, 8, 27), False),
            (offset_lom_sat_aug, datetime(2011, 8, 26), False),
            (offset_lom_sat_aug, datetime(2019, 8, 30), False),

            #From GMCR (see for example: http://yahoo.brand.edgar-online.com/Default.aspx?companyid=3184&formtypeID=7)
            (offset_lom_sat_sep, datetime(2010, 9, 25), True),
            (offset_lom_sat_sep, datetime(2011, 9, 24), True),
            (offset_lom_sat_sep, datetime(2012, 9, 29), True),

        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_apply(self):
        offset_lom_aug_sat = makeFY5253LastOfMonth(startingMonth=8, weekday=WeekDay.SAT)
        offset_lom_aug_sat_1 = makeFY5253LastOfMonth(n=1, startingMonth=8, weekday=WeekDay.SAT)

        date_seq_lom_aug_sat = [datetime(2006, 8, 26), datetime(2007, 8, 25),
                                datetime(2008, 8, 30), datetime(2009, 8, 29),
                                datetime(2010, 8, 28), datetime(2011, 8, 27),
                                datetime(2012, 8, 25), datetime(2013, 8, 31),
                                datetime(2014, 8, 30), datetime(2015, 8, 29),
                                datetime(2016, 8, 27)]

        tests = [
                 (offset_lom_aug_sat, date_seq_lom_aug_sat),
                 (offset_lom_aug_sat_1, date_seq_lom_aug_sat),
                 (offset_lom_aug_sat, [datetime(2006, 8, 25)] + date_seq_lom_aug_sat),
                 (offset_lom_aug_sat_1, [datetime(2006, 8, 27)] + date_seq_lom_aug_sat[1:]),
                 (makeFY5253LastOfMonth(n=-1, startingMonth=8, weekday=WeekDay.SAT), list(reversed(date_seq_lom_aug_sat))),
                ]
        for test in tests:
            offset, data = test
            current = data[0]
            for datum in data[1:]:
                current = current + offset
                self.assertEqual(current, datum)

class TestFY5253NearestEndMonth(TestBase):

    def test_get_target_month_end(self):
        self.assertEqual(makeFY5253NearestEndMonth(startingMonth=8, weekday=WeekDay.SAT).get_target_month_end(datetime(2013,1,1)), datetime(2013,8,31))
        self.assertEqual(makeFY5253NearestEndMonth(startingMonth=12, weekday=WeekDay.SAT).get_target_month_end(datetime(2013,1,1)), datetime(2013,12,31))
        self.assertEqual(makeFY5253NearestEndMonth(startingMonth=2, weekday=WeekDay.SAT).get_target_month_end(datetime(2013,1,1)), datetime(2013,2,28))

    def test_get_year_end(self):
        self.assertEqual(makeFY5253NearestEndMonth(startingMonth=8, weekday=WeekDay.SAT).get_year_end(datetime(2013,1,1)), datetime(2013,8,31))
        self.assertEqual(makeFY5253NearestEndMonth(startingMonth=8, weekday=WeekDay.SUN).get_year_end(datetime(2013,1,1)), datetime(2013,9,1))
        self.assertEqual(makeFY5253NearestEndMonth(startingMonth=8, weekday=WeekDay.FRI).get_year_end(datetime(2013,1,1)), datetime(2013,8,30))

        offset_n = FY5253(weekday=WeekDay.TUE, startingMonth=12,
                      variation="nearest")
        self.assertEqual(offset_n.get_year_end(datetime(2012,1,1)), datetime(2013,1,1))
        self.assertEqual(offset_n.get_year_end(datetime(2012,1,10)), datetime(2013,1,1))

        self.assertEqual(offset_n.get_year_end(datetime(2013,1,1)), datetime(2013,12,31))
        self.assertEqual(offset_n.get_year_end(datetime(2013,1,2)), datetime(2013,12,31))
        self.assertEqual(offset_n.get_year_end(datetime(2013,1,3)), datetime(2013,12,31))
        self.assertEqual(offset_n.get_year_end(datetime(2013,1,10)), datetime(2013,12,31))

        JNJ = FY5253(n=1, startingMonth=12, weekday=6, variation="nearest")
        self.assertEqual(JNJ.get_year_end(datetime(2006, 1, 1)), datetime(2006, 12, 31))

    def test_onOffset(self):
        offset_lom_aug_sat = makeFY5253NearestEndMonth(1, startingMonth=8, weekday=WeekDay.SAT)
        offset_lom_aug_thu = makeFY5253NearestEndMonth(1, startingMonth=8, weekday=WeekDay.THU)
        offset_n = FY5253(weekday=WeekDay.TUE, startingMonth=12,
                      variation="nearest")

        tests = [
#             From Wikipedia (see: http://en.wikipedia.org/wiki/4%E2%80%934%E2%80%935_calendar#Saturday_nearest_the_end_of_month)
#             2006-09-02   2006 September 2
#             2007-09-01   2007 September 1
#             2008-08-30   2008 August 30    (leap year)
#             2009-08-29   2009 August 29
#             2010-08-28   2010 August 28
#             2011-09-03   2011 September 3
#             2012-09-01   2012 September 1  (leap year)
#             2013-08-31   2013 August 31
#             2014-08-30   2014 August 30
#             2015-08-29   2015 August 29
#             2016-09-03   2016 September 3  (leap year)
#             2017-09-02   2017 September 2
#             2018-09-01   2018 September 1
#             2019-08-31   2019 August 31
            (offset_lom_aug_sat, datetime(2006, 9, 2), True),
            (offset_lom_aug_sat, datetime(2007, 9, 1), True),
            (offset_lom_aug_sat, datetime(2008, 8, 30), True),
            (offset_lom_aug_sat, datetime(2009, 8, 29), True),
            (offset_lom_aug_sat, datetime(2010, 8, 28), True),
            (offset_lom_aug_sat, datetime(2011, 9, 3), True),

            (offset_lom_aug_sat, datetime(2016, 9, 3), True),
            (offset_lom_aug_sat, datetime(2017, 9, 2), True),
            (offset_lom_aug_sat, datetime(2018, 9, 1), True),
            (offset_lom_aug_sat, datetime(2019, 8, 31), True),

            (offset_lom_aug_sat, datetime(2006, 8, 27), False),
            (offset_lom_aug_sat, datetime(2007, 8, 28), False),
            (offset_lom_aug_sat, datetime(2008, 8, 31), False),
            (offset_lom_aug_sat, datetime(2009, 8, 30), False),
            (offset_lom_aug_sat, datetime(2010, 8, 29), False),
            (offset_lom_aug_sat, datetime(2011, 8, 28), False),

            (offset_lom_aug_sat, datetime(2006, 8, 25), False),
            (offset_lom_aug_sat, datetime(2007, 8, 24), False),
            (offset_lom_aug_sat, datetime(2008, 8, 29), False),
            (offset_lom_aug_sat, datetime(2009, 8, 28), False),
            (offset_lom_aug_sat, datetime(2010, 8, 27), False),
            (offset_lom_aug_sat, datetime(2011, 8, 26), False),
            (offset_lom_aug_sat, datetime(2019, 8, 30), False),

            #From Micron, see: http://google.brand.edgar-online.com/?sym=MU&formtypeID=7
            (offset_lom_aug_thu, datetime(2012, 8, 30), True),
            (offset_lom_aug_thu, datetime(2011, 9, 1), True),

            (offset_n, datetime(2012, 12, 31), False),
            (offset_n, datetime(2013, 1, 1), True),
            (offset_n, datetime(2013, 1, 2), False),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_apply(self):
        date_seq_nem_8_sat = [datetime(2006, 9, 2), datetime(2007, 9, 1),
                              datetime(2008, 8, 30), datetime(2009, 8, 29),
                              datetime(2010, 8, 28), datetime(2011, 9, 3)]

        JNJ = [datetime(2005, 1, 2), datetime(2006, 1, 1),
               datetime(2006, 12, 31), datetime(2007, 12, 30),
               datetime(2008, 12, 28), datetime(2010, 1, 3),
               datetime(2011, 1, 2), datetime(2012, 1, 1),
               datetime(2012, 12, 30)]

        DEC_SAT = FY5253(n=-1, startingMonth=12, weekday=5, variation="nearest")

        tests = [
                (makeFY5253NearestEndMonth(startingMonth=8, weekday=WeekDay.SAT), date_seq_nem_8_sat),
                (makeFY5253NearestEndMonth(n=1, startingMonth=8, weekday=WeekDay.SAT), date_seq_nem_8_sat),
                (makeFY5253NearestEndMonth(startingMonth=8, weekday=WeekDay.SAT), [datetime(2006, 9, 1)] + date_seq_nem_8_sat),
                (makeFY5253NearestEndMonth(n=1, startingMonth=8, weekday=WeekDay.SAT), [datetime(2006, 9, 3)] + date_seq_nem_8_sat[1:]),
                (makeFY5253NearestEndMonth(n=-1, startingMonth=8, weekday=WeekDay.SAT), list(reversed(date_seq_nem_8_sat))),
                (makeFY5253NearestEndMonth(n=1, startingMonth=12, weekday=WeekDay.SUN), JNJ),
                (makeFY5253NearestEndMonth(n=-1, startingMonth=12, weekday=WeekDay.SUN), list(reversed(JNJ))),
                (makeFY5253NearestEndMonth(n=1, startingMonth=12, weekday=WeekDay.SUN), [datetime(2005,1,2), datetime(2006, 1, 1)]),
                (makeFY5253NearestEndMonth(n=1, startingMonth=12, weekday=WeekDay.SUN), [datetime(2006,1,2), datetime(2006, 12, 31)]),
                (DEC_SAT, [datetime(2013,1,15), datetime(2012,12,29)])
                ]
        for test in tests:
            offset, data = test
            current = data[0]
            for datum in data[1:]:
                current = current + offset
                self.assertEqual(current, datum)

class TestFY5253LastOfMonthQuarter(TestBase):

    def test_isAnchored(self):
        self.assert_(makeFY5253LastOfMonthQuarter(startingMonth=1, weekday=WeekDay.SAT, qtr_with_extra_week=4).isAnchored())
        self.assert_(makeFY5253LastOfMonthQuarter(weekday=WeekDay.SAT, startingMonth=3, qtr_with_extra_week=4).isAnchored())
        self.assert_(not makeFY5253LastOfMonthQuarter(2, startingMonth=1, weekday=WeekDay.SAT, qtr_with_extra_week=4).isAnchored())

    def test_equality(self):
        self.assertEqual(makeFY5253LastOfMonthQuarter(startingMonth=1, weekday=WeekDay.SAT, qtr_with_extra_week=4), makeFY5253LastOfMonthQuarter(startingMonth=1, weekday=WeekDay.SAT, qtr_with_extra_week=4))
        self.assertNotEqual(makeFY5253LastOfMonthQuarter(startingMonth=1, weekday=WeekDay.SAT, qtr_with_extra_week=4), makeFY5253LastOfMonthQuarter(startingMonth=1, weekday=WeekDay.SUN, qtr_with_extra_week=4))
        self.assertNotEqual(makeFY5253LastOfMonthQuarter(startingMonth=1, weekday=WeekDay.SAT, qtr_with_extra_week=4), makeFY5253LastOfMonthQuarter(startingMonth=2, weekday=WeekDay.SAT, qtr_with_extra_week=4))

    def test_offset(self):
        offset = makeFY5253LastOfMonthQuarter(1, startingMonth=9, weekday=WeekDay.SAT, qtr_with_extra_week=4)
        offset2 = makeFY5253LastOfMonthQuarter(2, startingMonth=9, weekday=WeekDay.SAT, qtr_with_extra_week=4)
        offset4 = makeFY5253LastOfMonthQuarter(4, startingMonth=9, weekday=WeekDay.SAT, qtr_with_extra_week=4)

        offset_neg1 = makeFY5253LastOfMonthQuarter(-1, startingMonth=9, weekday=WeekDay.SAT, qtr_with_extra_week=4)
        offset_neg2 = makeFY5253LastOfMonthQuarter(-2, startingMonth=9, weekday=WeekDay.SAT, qtr_with_extra_week=4)

        GMCR = [datetime(2010, 3, 27),
                datetime(2010, 6, 26),
                datetime(2010, 9, 25),
                datetime(2010, 12, 25),
                datetime(2011, 3, 26),
                datetime(2011, 6, 25),
                datetime(2011, 9, 24),
                datetime(2011, 12, 24),
                datetime(2012, 3, 24),
                datetime(2012, 6, 23),
                datetime(2012, 9, 29),
                datetime(2012, 12, 29),
                datetime(2013, 3, 30),
                datetime(2013, 6, 29)]


        assertEq(offset, base=GMCR[0], expected=GMCR[1])
        assertEq(offset, base=GMCR[0] + relativedelta(days=-1), expected=GMCR[0])
        assertEq(offset, base=GMCR[1], expected=GMCR[2])

        assertEq(offset2, base=GMCR[0], expected=GMCR[2])
        assertEq(offset4, base=GMCR[0], expected=GMCR[4])

        assertEq(offset_neg1, base=GMCR[-1], expected=GMCR[-2])
        assertEq(offset_neg1, base=GMCR[-1] + relativedelta(days=+1), expected=GMCR[-1])
        assertEq(offset_neg2, base=GMCR[-1], expected=GMCR[-3])

        date = GMCR[0] + relativedelta(days=-1)
        for expected in GMCR:
            assertEq(offset, date, expected)
            date = date + offset

        date = GMCR[-1] + relativedelta(days=+1)
        for expected in reversed(GMCR):
            assertEq(offset_neg1, date, expected)
            date = date + offset_neg1


    def test_onOffset(self):
        lomq_aug_sat_4 = makeFY5253LastOfMonthQuarter(1, startingMonth=8, weekday=WeekDay.SAT, qtr_with_extra_week=4)
        lomq_sep_sat_4 = makeFY5253LastOfMonthQuarter(1, startingMonth=9, weekday=WeekDay.SAT, qtr_with_extra_week=4)

        tests = [
            #From Wikipedia
            (lomq_aug_sat_4, datetime(2006, 8, 26), True),
            (lomq_aug_sat_4, datetime(2007, 8, 25), True),
            (lomq_aug_sat_4, datetime(2008, 8, 30), True),
            (lomq_aug_sat_4, datetime(2009, 8, 29), True),
            (lomq_aug_sat_4, datetime(2010, 8, 28), True),
            (lomq_aug_sat_4, datetime(2011, 8, 27), True),
            (lomq_aug_sat_4, datetime(2019, 8, 31), True),

            (lomq_aug_sat_4, datetime(2006, 8, 27), False),
            (lomq_aug_sat_4, datetime(2007, 8, 28), False),
            (lomq_aug_sat_4, datetime(2008, 8, 31), False),
            (lomq_aug_sat_4, datetime(2009, 8, 30), False),
            (lomq_aug_sat_4, datetime(2010, 8, 29), False),
            (lomq_aug_sat_4, datetime(2011, 8, 28), False),

            (lomq_aug_sat_4, datetime(2006, 8, 25), False),
            (lomq_aug_sat_4, datetime(2007, 8, 24), False),
            (lomq_aug_sat_4, datetime(2008, 8, 29), False),
            (lomq_aug_sat_4, datetime(2009, 8, 28), False),
            (lomq_aug_sat_4, datetime(2010, 8, 27), False),
            (lomq_aug_sat_4, datetime(2011, 8, 26), False),
            (lomq_aug_sat_4, datetime(2019, 8, 30), False),

            #From GMCR
            (lomq_sep_sat_4, datetime(2010, 9, 25), True),
            (lomq_sep_sat_4, datetime(2011, 9, 24), True),
            (lomq_sep_sat_4, datetime(2012, 9, 29), True),

            (lomq_sep_sat_4, datetime(2013, 6, 29), True),
            (lomq_sep_sat_4, datetime(2012, 6, 23), True),
            (lomq_sep_sat_4, datetime(2012, 6, 30), False),

            (lomq_sep_sat_4, datetime(2013, 3, 30), True),
            (lomq_sep_sat_4, datetime(2012, 3, 24), True),

            (lomq_sep_sat_4, datetime(2012, 12, 29), True),
            (lomq_sep_sat_4, datetime(2011, 12, 24), True),

            #INTC (extra week in Q1)
            #See: http://www.intc.com/releasedetail.cfm?ReleaseID=542844
            (makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1), datetime(2011, 4, 2), True),

            #see: http://google.brand.edgar-online.com/?sym=INTC&formtypeID=7
            (makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1), datetime(2012, 12, 29), True),
            (makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1), datetime(2011, 12, 31), True),
            (makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1), datetime(2010, 12, 25), True),

        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_year_has_extra_week(self):
        #End of long Q1
        self.assertTrue(makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1).year_has_extra_week(datetime(2011, 4, 2)))

        #Start of long Q1
        self.assertTrue(makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1).year_has_extra_week(datetime(2010, 12, 26)))

        #End of year before year with long Q1
        self.assertFalse(makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1).year_has_extra_week(datetime(2010, 12, 25)))

        for year in [x for x in range(1994, 2011+1) if x not in [2011, 2005, 2000, 1994]]:
            self.assertFalse(makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1).year_has_extra_week(datetime(year, 4, 2)))

        #Other long years
        self.assertTrue(makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1).year_has_extra_week(datetime(2005, 4, 2)))
        self.assertTrue(makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1).year_has_extra_week(datetime(2000, 4, 2)))
        self.assertTrue(makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1).year_has_extra_week(datetime(1994, 4, 2)))

    def test_get_weeks(self):
        sat_dec_1 = makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=1)
        sat_dec_4 = makeFY5253LastOfMonthQuarter(1, startingMonth=12, weekday=WeekDay.SAT, qtr_with_extra_week=4)

        self.assertEqual(sat_dec_1.get_weeks(datetime(2011, 4, 2)), [14, 13, 13, 13])
        self.assertEqual(sat_dec_4.get_weeks(datetime(2011, 4, 2)), [13, 13, 13, 14])
        self.assertEqual(sat_dec_1.get_weeks(datetime(2010, 12, 25)), [13, 13, 13, 13])

class TestFY5253NearestEndMonthQuarter(TestBase):

    def test_onOffset(self):

        offset_nem_sat_aug_4 = makeFY5253NearestEndMonthQuarter(1, startingMonth=8, weekday=WeekDay.SAT, qtr_with_extra_week=4)
        offset_nem_thu_aug_4 = makeFY5253NearestEndMonthQuarter(1, startingMonth=8, weekday=WeekDay.THU, qtr_with_extra_week=4)
        offset_n = FY5253(weekday=WeekDay.TUE, startingMonth=12,
                      variation="nearest", qtr_with_extra_week=4)

        tests = [
            #From Wikipedia
            (offset_nem_sat_aug_4, datetime(2006, 9, 2), True),
            (offset_nem_sat_aug_4, datetime(2007, 9, 1), True),
            (offset_nem_sat_aug_4, datetime(2008, 8, 30), True),
            (offset_nem_sat_aug_4, datetime(2009, 8, 29), True),
            (offset_nem_sat_aug_4, datetime(2010, 8, 28), True),
            (offset_nem_sat_aug_4, datetime(2011, 9, 3), True),

            (offset_nem_sat_aug_4, datetime(2016, 9, 3), True),
            (offset_nem_sat_aug_4, datetime(2017, 9, 2), True),
            (offset_nem_sat_aug_4, datetime(2018, 9, 1), True),
            (offset_nem_sat_aug_4, datetime(2019, 8, 31), True),

            (offset_nem_sat_aug_4, datetime(2006, 8, 27), False),
            (offset_nem_sat_aug_4, datetime(2007, 8, 28), False),
            (offset_nem_sat_aug_4, datetime(2008, 8, 31), False),
            (offset_nem_sat_aug_4, datetime(2009, 8, 30), False),
            (offset_nem_sat_aug_4, datetime(2010, 8, 29), False),
            (offset_nem_sat_aug_4, datetime(2011, 8, 28), False),

            (offset_nem_sat_aug_4, datetime(2006, 8, 25), False),
            (offset_nem_sat_aug_4, datetime(2007, 8, 24), False),
            (offset_nem_sat_aug_4, datetime(2008, 8, 29), False),
            (offset_nem_sat_aug_4, datetime(2009, 8, 28), False),
            (offset_nem_sat_aug_4, datetime(2010, 8, 27), False),
            (offset_nem_sat_aug_4, datetime(2011, 8, 26), False),
            (offset_nem_sat_aug_4, datetime(2019, 8, 30), False),

            #From Micron, see: http://google.brand.edgar-online.com/?sym=MU&formtypeID=7
            (offset_nem_thu_aug_4, datetime(2012, 8, 30), True),
            (offset_nem_thu_aug_4, datetime(2011, 9, 1), True),

            #See: http://google.brand.edgar-online.com/?sym=MU&formtypeID=13
            (offset_nem_thu_aug_4, datetime(2013, 5, 30), True),
            (offset_nem_thu_aug_4, datetime(2013, 2, 28), True),
            (offset_nem_thu_aug_4, datetime(2012, 11, 29), True),
            (offset_nem_thu_aug_4, datetime(2012, 5, 31), True),
            (offset_nem_thu_aug_4, datetime(2007, 3, 1), True),
            (offset_nem_thu_aug_4, datetime(1994, 3, 3), True),

            (offset_n, datetime(2012, 12, 31), False),
            (offset_n, datetime(2013, 1, 1), True),
            (offset_n, datetime(2013, 1, 2), False)
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)

    def test_offset(self):
        offset = makeFY5253NearestEndMonthQuarter(1, startingMonth=8, weekday=WeekDay.THU, qtr_with_extra_week=4)

        MU = [datetime(2012, 5, 31), datetime(2012, 8, 30), datetime(2012, 11, 29), datetime(2013, 2, 28), datetime(2013, 5, 30)]

        date = MU[0] + relativedelta(days=-1)
        for expected in MU:
            assertEq(offset, date, expected)
            date = date + offset

        assertEq(offset, datetime(2012, 5, 31), datetime(2012, 8, 30))
        assertEq(offset, datetime(2012, 5, 30), datetime(2012, 5, 31))

        offset2 = FY5253Quarter(weekday=5, startingMonth=12,
                     variation="last", qtr_with_extra_week=4)

        assertEq(offset2, datetime(2013,1,15), datetime(2013, 3, 30))

class TestQuarterBegin(TestBase):

    def test_repr(self):
        self.assertEqual(repr(QuarterBegin()), "<QuarterBegin: startingMonth=3>")
        self.assertEqual(repr(QuarterBegin(startingMonth=3)), "<QuarterBegin: startingMonth=3>")
        self.assertEqual(repr(QuarterBegin(startingMonth=1)),"<QuarterBegin: startingMonth=1>")

    def test_isAnchored(self):
        self.assert_(QuarterBegin(startingMonth=1).isAnchored())
        self.assert_(QuarterBegin().isAnchored())
        self.assert_(not QuarterBegin(2, startingMonth=1).isAnchored())

    def test_offset(self):
        tests = []

        tests.append((QuarterBegin(startingMonth=1),
                      {datetime(2007, 12, 1): datetime(2008, 1, 1),
                       datetime(2008, 1, 1): datetime(2008, 4, 1),
                       datetime(2008, 2, 15): datetime(2008, 4, 1),
                       datetime(2008, 2, 29): datetime(2008, 4, 1),
                       datetime(2008, 3, 15): datetime(2008, 4, 1),
                       datetime(2008, 3, 31): datetime(2008, 4, 1),
                       datetime(2008, 4, 15): datetime(2008, 7, 1),
                       datetime(2008, 4, 1): datetime(2008, 7, 1), }))

        tests.append((QuarterBegin(startingMonth=2),
                      {datetime(2008, 1, 1): datetime(2008, 2, 1),
                       datetime(2008, 1, 31): datetime(2008, 2, 1),
                       datetime(2008, 1, 15): datetime(2008, 2, 1),
                       datetime(2008, 2, 29): datetime(2008, 5, 1),
                       datetime(2008, 3, 15): datetime(2008, 5, 1),
                       datetime(2008, 3, 31): datetime(2008, 5, 1),
                       datetime(2008, 4, 15): datetime(2008, 5, 1),
                       datetime(2008, 4, 30): datetime(2008, 5, 1), }))

        tests.append((QuarterBegin(startingMonth=1, n=0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 12, 1): datetime(2009, 1, 1),
                       datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 2, 15): datetime(2008, 4, 1),
                       datetime(2008, 2, 29): datetime(2008, 4, 1),
                       datetime(2008, 3, 15): datetime(2008, 4, 1),
                       datetime(2008, 3, 31): datetime(2008, 4, 1),
                       datetime(2008, 4, 15): datetime(2008, 4, 1),
                       datetime(2008, 4, 30): datetime(2008, 4, 1), }))

        tests.append((QuarterBegin(startingMonth=1, n=-1),
                      {datetime(2008, 1, 1): datetime(2007, 10, 1),
                       datetime(2008, 1, 31): datetime(2008, 1, 1),
                       datetime(2008, 2, 15): datetime(2008, 1, 1),
                       datetime(2008, 2, 29): datetime(2008, 1, 1),
                       datetime(2008, 3, 15): datetime(2008, 1, 1),
                       datetime(2008, 3, 31): datetime(2008, 1, 1),
                       datetime(2008, 4, 15): datetime(2008, 4, 1),
                       datetime(2008, 4, 30): datetime(2008, 4, 1),
                       datetime(2008, 7, 1): datetime(2008, 4, 1)}))

        tests.append((QuarterBegin(startingMonth=1, n=2),
                      {datetime(2008, 1, 1): datetime(2008, 7, 1),
                       datetime(2008, 2, 15): datetime(2008, 7, 1),
                       datetime(2008, 2, 29): datetime(2008, 7, 1),
                       datetime(2008, 3, 15): datetime(2008, 7, 1),
                       datetime(2008, 3, 31): datetime(2008, 7, 1),
                       datetime(2008, 4, 15): datetime(2008, 10, 1),
                       datetime(2008, 4, 1): datetime(2008, 10, 1), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

        # corner
        offset = QuarterBegin(n=-1, startingMonth=1)
        self.assertEqual(datetime(2010, 2, 1) + offset, datetime(2010, 1, 1))


class TestQuarterEnd(TestBase):
    _offset = QuarterEnd

    def test_repr(self):
        self.assertEqual(repr(QuarterEnd()), "<QuarterEnd: startingMonth=3>")
        self.assertEqual(repr(QuarterEnd(startingMonth=3)), "<QuarterEnd: startingMonth=3>")
        self.assertEqual(repr(QuarterEnd(startingMonth=1)), "<QuarterEnd: startingMonth=1>")

    def test_isAnchored(self):
        self.assert_(QuarterEnd(startingMonth=1).isAnchored())
        self.assert_(QuarterEnd().isAnchored())
        self.assert_(not QuarterEnd(2, startingMonth=1).isAnchored())

    def test_offset(self):
        tests = []

        tests.append((QuarterEnd(startingMonth=1),
                      {datetime(2008, 1, 1): datetime(2008, 1, 31),
                       datetime(2008, 1, 31): datetime(2008, 4, 30),
                       datetime(2008, 2, 15): datetime(2008, 4, 30),
                       datetime(2008, 2, 29): datetime(2008, 4, 30),
                       datetime(2008, 3, 15): datetime(2008, 4, 30),
                       datetime(2008, 3, 31): datetime(2008, 4, 30),
                       datetime(2008, 4, 15): datetime(2008, 4, 30),
                       datetime(2008, 4, 30): datetime(2008, 7, 31), }))

        tests.append((QuarterEnd(startingMonth=2),
                      {datetime(2008, 1, 1): datetime(2008, 2, 29),
                       datetime(2008, 1, 31): datetime(2008, 2, 29),
                       datetime(2008, 2, 15): datetime(2008, 2, 29),
                       datetime(2008, 2, 29): datetime(2008, 5, 31),
                       datetime(2008, 3, 15): datetime(2008, 5, 31),
                       datetime(2008, 3, 31): datetime(2008, 5, 31),
                       datetime(2008, 4, 15): datetime(2008, 5, 31),
                       datetime(2008, 4, 30): datetime(2008, 5, 31), }))

        tests.append((QuarterEnd(startingMonth=1, n=0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 31),
                       datetime(2008, 1, 31): datetime(2008, 1, 31),
                       datetime(2008, 2, 15): datetime(2008, 4, 30),
                       datetime(2008, 2, 29): datetime(2008, 4, 30),
                       datetime(2008, 3, 15): datetime(2008, 4, 30),
                       datetime(2008, 3, 31): datetime(2008, 4, 30),
                       datetime(2008, 4, 15): datetime(2008, 4, 30),
                       datetime(2008, 4, 30): datetime(2008, 4, 30), }))

        tests.append((QuarterEnd(startingMonth=1, n=-1),
                      {datetime(2008, 1, 1): datetime(2007, 10, 31),
                       datetime(2008, 1, 31): datetime(2007, 10, 31),
                       datetime(2008, 2, 15): datetime(2008, 1, 31),
                       datetime(2008, 2, 29): datetime(2008, 1, 31),
                       datetime(2008, 3, 15): datetime(2008, 1, 31),
                       datetime(2008, 3, 31): datetime(2008, 1, 31),
                       datetime(2008, 4, 15): datetime(2008, 1, 31),
                       datetime(2008, 4, 30): datetime(2008, 1, 31),
                       datetime(2008, 7, 1): datetime(2008, 4, 30)}))

        tests.append((QuarterEnd(startingMonth=1, n=2),
                      {datetime(2008, 1, 31): datetime(2008, 7, 31),
                       datetime(2008, 2, 15): datetime(2008, 7, 31),
                       datetime(2008, 2, 29): datetime(2008, 7, 31),
                       datetime(2008, 3, 15): datetime(2008, 7, 31),
                       datetime(2008, 3, 31): datetime(2008, 7, 31),
                       datetime(2008, 4, 15): datetime(2008, 7, 31),
                       datetime(2008, 4, 30): datetime(2008, 10, 31), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

        # corner
        offset = QuarterEnd(n=-1, startingMonth=1)
        self.assertEqual(datetime(2010, 2, 1) + offset, datetime(2010, 1, 31))

    def test_onOffset(self):

        tests = [(QuarterEnd(1, startingMonth=1), datetime(2008, 1, 31), True),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2007, 12, 31), False),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2008, 2, 29), False),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2007, 3, 30), False),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2007, 3, 31), False),
                 (QuarterEnd(1, startingMonth=1), datetime(2008, 4, 30), True),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2008, 5, 30), False),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2008, 5, 31), False),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2007, 6, 29), False),
                 (QuarterEnd(
                     1, startingMonth=1), datetime(2007, 6, 30), False),

                 (QuarterEnd(
                     1, startingMonth=2), datetime(2008, 1, 31), False),
                 (QuarterEnd(
                     1, startingMonth=2), datetime(2007, 12, 31), False),
                 (QuarterEnd(1, startingMonth=2), datetime(2008, 2, 29), True),
                 (QuarterEnd(
                     1, startingMonth=2), datetime(2007, 3, 30), False),
                 (QuarterEnd(
                     1, startingMonth=2), datetime(2007, 3, 31), False),
                 (QuarterEnd(
                     1, startingMonth=2), datetime(2008, 4, 30), False),
                 (QuarterEnd(
                     1, startingMonth=2), datetime(2008, 5, 30), False),
                 (QuarterEnd(1, startingMonth=2), datetime(2008, 5, 31), True),
                 (QuarterEnd(
                     1, startingMonth=2), datetime(2007, 6, 29), False),
                 (QuarterEnd(
                     1, startingMonth=2), datetime(2007, 6, 30), False),

                 (QuarterEnd(
                     1, startingMonth=3), datetime(2008, 1, 31), False),
                 (QuarterEnd(
                     1, startingMonth=3), datetime(2007, 12, 31), True),
                 (QuarterEnd(
                     1, startingMonth=3), datetime(2008, 2, 29), False),
                 (QuarterEnd(
                     1, startingMonth=3), datetime(2007, 3, 30), False),
                 (QuarterEnd(1, startingMonth=3), datetime(2007, 3, 31), True),
                 (QuarterEnd(
                     1, startingMonth=3), datetime(2008, 4, 30), False),
                 (QuarterEnd(
                     1, startingMonth=3), datetime(2008, 5, 30), False),
                 (QuarterEnd(
                     1, startingMonth=3), datetime(2008, 5, 31), False),
                 (QuarterEnd(
                     1, startingMonth=3), datetime(2007, 6, 29), False),
                 (QuarterEnd(1, startingMonth=3), datetime(2007, 6, 30), True),
                 ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


class TestBYearBegin(TestBase):
    _offset = BYearBegin

    def test_misspecified(self):
        self.assertRaises(ValueError, BYearBegin, month=13)
        self.assertRaises(ValueError, BYearEnd, month=13)

    def test_offset(self):
        tests = []

        tests.append((BYearBegin(),
                      {datetime(2008, 1, 1): datetime(2009, 1, 1),
                       datetime(2008, 6, 30): datetime(2009, 1, 1),
                       datetime(2008, 12, 31): datetime(2009, 1, 1),
                       datetime(2011, 1, 1): datetime(2011, 1, 3),
                       datetime(2011, 1, 3): datetime(2012, 1, 2),
                       datetime(2005, 12, 30): datetime(2006, 1, 2),
                       datetime(2005, 12, 31): datetime(2006, 1, 2)
                       }
                      ))

        tests.append((BYearBegin(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 6, 30): datetime(2009, 1, 1),
                       datetime(2008, 12, 31): datetime(2009, 1, 1),
                       datetime(2005, 12, 30): datetime(2006, 1, 2),
                       datetime(2005, 12, 31): datetime(2006, 1, 2), }))

        tests.append((BYearBegin(-1),
                      {datetime(2007, 1, 1): datetime(2006, 1, 2),
                       datetime(2009, 1, 4): datetime(2009, 1, 1),
                       datetime(2009, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 6, 30): datetime(2008, 1, 1),
                       datetime(2008, 12, 31): datetime(2008, 1, 1),
                       datetime(2006, 12, 29): datetime(2006, 1, 2),
                       datetime(2006, 12, 30): datetime(2006, 1, 2),
                       datetime(2006, 1, 1): datetime(2005, 1, 3), }))

        tests.append((BYearBegin(-2),
                      {datetime(2007, 1, 1): datetime(2005, 1, 3),
                       datetime(2007, 6, 30): datetime(2006, 1, 2),
                       datetime(2008, 12, 31): datetime(2007, 1, 1), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)


class TestYearBegin(TestBase):
    _offset = YearBegin

    def test_misspecified(self):
        self.assertRaises(ValueError, YearBegin, month=13)

    def test_offset(self):
        tests = []

        tests.append((YearBegin(),
                      {datetime(2008, 1, 1): datetime(2009, 1, 1),
                       datetime(2008, 6, 30): datetime(2009, 1, 1),
                       datetime(2008, 12, 31): datetime(2009, 1, 1),
                       datetime(2005, 12, 30): datetime(2006, 1, 1),
                       datetime(2005, 12, 31): datetime(2006, 1, 1), }))

        tests.append((YearBegin(0),
                      {datetime(2008, 1, 1): datetime(2008, 1, 1),
                       datetime(2008, 6, 30): datetime(2009, 1, 1),
                       datetime(2008, 12, 31): datetime(2009, 1, 1),
                       datetime(2005, 12, 30): datetime(2006, 1, 1),
                       datetime(2005, 12, 31): datetime(2006, 1, 1), }))

        tests.append((YearBegin(-1),
                      {datetime(2007, 1, 1): datetime(2006, 1, 1),
                       datetime(2007, 1, 15): datetime(2007, 1, 1),
                       datetime(2008, 6, 30): datetime(2008, 1, 1),
                       datetime(2008, 12, 31): datetime(2008, 1, 1),
                       datetime(2006, 12, 29): datetime(2006, 1, 1),
                       datetime(2006, 12, 30): datetime(2006, 1, 1),
                       datetime(2007, 1, 1): datetime(2006, 1, 1), }))

        tests.append((YearBegin(-2),
                      {datetime(2007, 1, 1): datetime(2005, 1, 1),
                       datetime(2008, 6, 30): datetime(2007, 1, 1),
                       datetime(2008, 12, 31): datetime(2007, 1, 1), }))

        tests.append((YearBegin(month=4),
                      {datetime(2007, 4, 1): datetime(2008, 4, 1),
                       datetime(2007, 4, 15): datetime(2008, 4, 1),
                       datetime(2007, 3, 1): datetime(2007, 4, 1),
                       datetime(2007, 12, 15): datetime(2008, 4, 1),
                       datetime(2012, 1, 31): datetime(2012, 4, 1), }))

        tests.append((YearBegin(0, month=4),
                      {datetime(2007, 4, 1): datetime(2007, 4, 1),
                       datetime(2007, 3, 1): datetime(2007, 4, 1),
                       datetime(2007, 12, 15): datetime(2008, 4, 1),
                       datetime(2012, 1, 31): datetime(2012, 4, 1), }))

        tests.append((YearBegin(-1, month=4),
                      {datetime(2007, 4, 1): datetime(2006, 4, 1),
                       datetime(2007, 3, 1): datetime(2006, 4, 1),
                       datetime(2007, 12, 15): datetime(2007, 4, 1),
                       datetime(2012, 1, 31): datetime(2011, 4, 1), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_onOffset(self):

        tests = [
            (YearBegin(), datetime(2007, 1, 3), False),
            (YearBegin(), datetime(2008, 1, 1), True),
            (YearBegin(), datetime(2006, 12, 31), False),
            (YearBegin(), datetime(2006, 1, 2), False),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


class TestBYearEndLagged(TestBase):

    def test_bad_month_fail(self):
        self.assertRaises(Exception, BYearEnd, month=13)
        self.assertRaises(Exception, BYearEnd, month=0)

    def test_offset(self):
        tests = []

        tests.append((BYearEnd(month=6),
                      {datetime(2008, 1, 1): datetime(2008, 6, 30),
                      datetime(2007, 6, 30): datetime(2008, 6, 30)},
                      ))

        tests.append((BYearEnd(n=-1, month=6),
                      {datetime(2008, 1, 1): datetime(2007, 6, 29),
                      datetime(2007, 6, 30): datetime(2007, 6, 29)},
                      ))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                self.assertEqual(base + offset, expected)

    def test_roll(self):
        offset = BYearEnd(month=6)
        date = datetime(2009, 11, 30)

        self.assertEqual(offset.rollforward(date), datetime(2010, 6, 30))
        self.assertEqual(offset.rollback(date), datetime(2009, 6, 30))

    def test_onOffset(self):

        tests = [
            (BYearEnd(month=2), datetime(2007, 2, 28), True),
            (BYearEnd(month=6), datetime(2007, 6, 30), False),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


class TestBYearEnd(TestBase):
    _offset = BYearEnd

    def test_offset(self):
        tests = []

        tests.append((BYearEnd(),
                      {datetime(2008, 1, 1): datetime(2008, 12, 31),
                       datetime(2008, 6, 30): datetime(2008, 12, 31),
                       datetime(2008, 12, 31): datetime(2009, 12, 31),
                       datetime(2005, 12, 30): datetime(2006, 12, 29),
                       datetime(2005, 12, 31): datetime(2006, 12, 29), }))

        tests.append((BYearEnd(0),
                      {datetime(2008, 1, 1): datetime(2008, 12, 31),
                       datetime(2008, 6, 30): datetime(2008, 12, 31),
                       datetime(2008, 12, 31): datetime(2008, 12, 31),
                       datetime(2005, 12, 31): datetime(2006, 12, 29), }))

        tests.append((BYearEnd(-1),
                      {datetime(2007, 1, 1): datetime(2006, 12, 29),
                       datetime(2008, 6, 30): datetime(2007, 12, 31),
                       datetime(2008, 12, 31): datetime(2007, 12, 31),
                       datetime(2006, 12, 29): datetime(2005, 12, 30),
                       datetime(2006, 12, 30): datetime(2006, 12, 29),
                       datetime(2007, 1, 1): datetime(2006, 12, 29), }))

        tests.append((BYearEnd(-2),
                      {datetime(2007, 1, 1): datetime(2005, 12, 30),
                       datetime(2008, 6, 30): datetime(2006, 12, 29),
                       datetime(2008, 12, 31): datetime(2006, 12, 29), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_onOffset(self):

        tests = [
            (BYearEnd(), datetime(2007, 12, 31), True),
            (BYearEnd(), datetime(2008, 1, 1), False),
            (BYearEnd(), datetime(2006, 12, 31), False),
            (BYearEnd(), datetime(2006, 12, 29), True),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


class TestYearEnd(TestBase):
    _offset = YearEnd

    def test_misspecified(self):
        self.assertRaises(ValueError, YearEnd, month=13)

    def test_offset(self):
        tests = []

        tests.append((YearEnd(),
                      {datetime(2008, 1, 1): datetime(2008, 12, 31),
                       datetime(2008, 6, 30): datetime(2008, 12, 31),
                       datetime(2008, 12, 31): datetime(2009, 12, 31),
                       datetime(2005, 12, 30): datetime(2005, 12, 31),
                       datetime(2005, 12, 31): datetime(2006, 12, 31), }))

        tests.append((YearEnd(0),
                      {datetime(2008, 1, 1): datetime(2008, 12, 31),
                       datetime(2008, 6, 30): datetime(2008, 12, 31),
                       datetime(2008, 12, 31): datetime(2008, 12, 31),
                       datetime(2005, 12, 30): datetime(2005, 12, 31), }))

        tests.append((YearEnd(-1),
                      {datetime(2007, 1, 1): datetime(2006, 12, 31),
                       datetime(2008, 6, 30): datetime(2007, 12, 31),
                       datetime(2008, 12, 31): datetime(2007, 12, 31),
                       datetime(2006, 12, 29): datetime(2005, 12, 31),
                       datetime(2006, 12, 30): datetime(2005, 12, 31),
                       datetime(2007, 1, 1): datetime(2006, 12, 31), }))

        tests.append((YearEnd(-2),
                      {datetime(2007, 1, 1): datetime(2005, 12, 31),
                       datetime(2008, 6, 30): datetime(2006, 12, 31),
                       datetime(2008, 12, 31): datetime(2006, 12, 31), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_onOffset(self):

        tests = [
            (YearEnd(), datetime(2007, 12, 31), True),
            (YearEnd(), datetime(2008, 1, 1), False),
            (YearEnd(), datetime(2006, 12, 31), True),
            (YearEnd(), datetime(2006, 12, 29), False),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


class TestYearEndDiffMonth(TestBase):

    def test_offset(self):
        tests = []

        tests.append((YearEnd(month=3),
                      {datetime(2008, 1, 1): datetime(2008, 3, 31),
                       datetime(2008, 2, 15): datetime(2008, 3, 31),
                       datetime(2008, 3, 31): datetime(2009, 3, 31),
                       datetime(2008, 3, 30): datetime(2008, 3, 31),
                       datetime(2005, 3, 31): datetime(2006, 3, 31),
                       datetime(2006, 7, 30): datetime(2007, 3, 31)}))

        tests.append((YearEnd(0, month=3),
                      {datetime(2008, 1, 1): datetime(2008, 3, 31),
                       datetime(2008, 2, 28): datetime(2008, 3, 31),
                       datetime(2008, 3, 31): datetime(2008, 3, 31),
                       datetime(2005, 3, 30): datetime(2005, 3, 31), }))

        tests.append((YearEnd(-1, month=3),
                      {datetime(2007, 1, 1): datetime(2006, 3, 31),
                       datetime(2008, 2, 28): datetime(2007, 3, 31),
                       datetime(2008, 3, 31): datetime(2007, 3, 31),
                       datetime(2006, 3, 29): datetime(2005, 3, 31),
                       datetime(2006, 3, 30): datetime(2005, 3, 31),
                       datetime(2007, 3, 1): datetime(2006, 3, 31), }))

        tests.append((YearEnd(-2, month=3),
                      {datetime(2007, 1, 1): datetime(2005, 3, 31),
                       datetime(2008, 6, 30): datetime(2007, 3, 31),
                       datetime(2008, 3, 31): datetime(2006, 3, 31), }))

        for offset, cases in tests:
            for base, expected in compat.iteritems(cases):
                assertEq(offset, base, expected)

    def test_onOffset(self):

        tests = [
            (YearEnd(month=3), datetime(2007, 3, 31), True),
            (YearEnd(month=3), datetime(2008, 1, 1), False),
            (YearEnd(month=3), datetime(2006, 3, 31), True),
            (YearEnd(month=3), datetime(2006, 3, 29), False),
        ]

        for offset, date, expected in tests:
            assertOnOffset(offset, date, expected)


def assertEq(offset, base, expected):
    actual = offset + base
    actual_swapped = base + offset
    try:
        assert actual == expected
        assert actual_swapped == expected
    except AssertionError:
        raise AssertionError("\nExpected: %s\nActual: %s\nFor Offset: %s)"
                             "\nAt Date: %s" %
                            (expected, actual, offset, base))


def test_Hour():
    assertEq(Hour(), datetime(2010, 1, 1), datetime(2010, 1, 1, 1))
    assertEq(Hour(-1), datetime(2010, 1, 1, 1), datetime(2010, 1, 1))
    assertEq(2 * Hour(), datetime(2010, 1, 1), datetime(2010, 1, 1, 2))
    assertEq(-1 * Hour(), datetime(2010, 1, 1, 1), datetime(2010, 1, 1))

    assert (Hour(3) + Hour(2)) == Hour(5)
    assert (Hour(3) - Hour(2)) == Hour()

    assert(Hour(4) != Hour(1))

    assert not Hour().isAnchored()


def test_Minute():
    assertEq(Minute(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 1))
    assertEq(Minute(-1), datetime(2010, 1, 1, 0, 1), datetime(2010, 1, 1))
    assertEq(2 * Minute(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 2))
    assertEq(-1 * Minute(), datetime(2010, 1, 1, 0, 1), datetime(2010, 1, 1))

    assert (Minute(3) + Minute(2)) == Minute(5)
    assert (Minute(3) - Minute(2)) == Minute()
    assert(Minute(5) != Minute())

    assert not Minute().isAnchored()


def test_Second():
    assertEq(Second(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 0, 1))
    assertEq(Second(-1), datetime(2010, 1, 1, 0, 0, 1), datetime(2010, 1, 1))
    assertEq(2 * Second(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 0, 2))
    assertEq(
        -1 * Second(), datetime(2010, 1, 1, 0, 0, 1), datetime(2010, 1, 1))

    assert (Second(3) + Second(2)) == Second(5)
    assert (Second(3) - Second(2)) == Second()

    assert not Second().isAnchored()


def test_Millisecond():
    assertEq(Milli(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 0, 0, 1000))
    assertEq(Milli(-1), datetime(2010, 1, 1, 0, 0, 0, 1000), datetime(2010, 1, 1))
    assertEq(Milli(2), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 0, 0, 2000))
    assertEq(2 * Milli(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 0, 0, 2000))
    assertEq(-1 * Milli(), datetime(2010, 1, 1, 0, 0, 0, 1000), datetime(2010, 1, 1))

    assert (Milli(3) + Milli(2)) == Milli(5)
    assert (Milli(3) - Milli(2)) == Milli()


def test_MillisecondTimestampArithmetic():
    assertEq(Milli(), Timestamp('2010-01-01'), Timestamp('2010-01-01 00:00:00.001'))
    assertEq(Milli(-1), Timestamp('2010-01-01 00:00:00.001'), Timestamp('2010-01-01'))


def test_Microsecond():
    assertEq(Micro(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 0, 0, 1))
    assertEq(Micro(-1), datetime(2010, 1, 1, 0, 0, 0, 1), datetime(2010, 1, 1))
    assertEq(2 * Micro(), datetime(2010, 1, 1), datetime(2010, 1, 1, 0, 0, 0, 2))
    assertEq(-1 * Micro(), datetime(2010, 1, 1, 0, 0, 0, 1), datetime(2010, 1, 1))

    assert (Micro(3) + Micro(2)) == Micro(5)
    assert (Micro(3) - Micro(2)) == Micro()


def test_NanosecondGeneric():
    if _np_version_under1p7:
        raise nose.SkipTest('numpy >= 1.7 required')
    timestamp = Timestamp(datetime(2010, 1, 1))
    assert timestamp.nanosecond == 0

    result = timestamp + Nano(10)
    assert result.nanosecond == 10

    reverse_result = Nano(10) + timestamp
    assert reverse_result.nanosecond == 10


def test_Nanosecond():
    if _np_version_under1p7:
        raise nose.SkipTest('numpy >= 1.7 required')

    timestamp = Timestamp(datetime(2010, 1, 1))
    assertEq(Nano(), timestamp, timestamp + np.timedelta64(1, 'ns'))
    assertEq(Nano(-1), timestamp + np.timedelta64(1, 'ns'), timestamp)
    assertEq(2 * Nano(), timestamp, timestamp + np.timedelta64(2, 'ns'))
    assertEq(-1 * Nano(), timestamp + np.timedelta64(1, 'ns'), timestamp)

    assert (Nano(3) + Nano(2)) == Nano(5)
    assert (Nano(3) - Nano(2)) == Nano()


def test_tick_offset():
    assert not Day().isAnchored()
    assert not Milli().isAnchored()
    assert not Micro().isAnchored()
    assert not Nano().isAnchored()


def test_compare_ticks():
    offsets = [Hour, Minute, Second, Milli, Micro]

    for kls in offsets:
        three = kls(3)
        four = kls(4)

        for _ in range(10):
            assert(three < kls(4))
            assert(kls(3) < four)
            assert(four > kls(3))
            assert(kls(4) > three)
            assert(kls(3) == kls(3))
            assert(kls(3) != kls(4))


class TestOffsetNames(tm.TestCase):
    def test_get_offset_name(self):
        assertRaisesRegexp(ValueError, 'Bad rule.*BusinessDays', get_offset_name, BDay(2))

        assert get_offset_name(BDay()) == 'B'
        assert get_offset_name(BMonthEnd()) == 'BM'
        assert get_offset_name(Week(weekday=0)) == 'W-MON'
        assert get_offset_name(Week(weekday=1)) == 'W-TUE'
        assert get_offset_name(Week(weekday=2)) == 'W-WED'
        assert get_offset_name(Week(weekday=3)) == 'W-THU'
        assert get_offset_name(Week(weekday=4)) == 'W-FRI'

        self.assertEqual(get_offset_name(LastWeekOfMonth(weekday=WeekDay.SUN)), "LWOM-SUN")
        self.assertEqual(get_offset_name(makeFY5253LastOfMonthQuarter(weekday=1, startingMonth=3, qtr_with_extra_week=4)),"REQ-L-MAR-TUE-4")
        self.assertEqual(get_offset_name(makeFY5253NearestEndMonthQuarter(weekday=1, startingMonth=3, qtr_with_extra_week=3)), "REQ-N-MAR-TUE-3")

def test_get_offset():
    assertRaisesRegexp(ValueError, "rule.*GIBBERISH", get_offset, 'gibberish')
    assertRaisesRegexp(ValueError, "rule.*QS-JAN-B", get_offset, 'QS-JAN-B')
    pairs = [
             ('B', BDay()), ('b', BDay()), ('bm', BMonthEnd()),
             ('Bm', BMonthEnd()), ('W-MON', Week(weekday=0)),
             ('W-TUE', Week(weekday=1)), ('W-WED', Week(weekday=2)),
             ('W-THU', Week(weekday=3)), ('W-FRI', Week(weekday=4)),
             ('w@Sat', Week(weekday=5)),
             ("RE-N-DEC-MON", makeFY5253NearestEndMonth(weekday=0, startingMonth=12)),
             ("RE-L-DEC-TUE", makeFY5253LastOfMonth(weekday=1, startingMonth=12)),
             ("REQ-L-MAR-TUE-4", makeFY5253LastOfMonthQuarter(weekday=1, startingMonth=3, qtr_with_extra_week=4)),
             ("REQ-L-DEC-MON-3", makeFY5253LastOfMonthQuarter(weekday=0, startingMonth=12, qtr_with_extra_week=3)),
             ("REQ-N-DEC-MON-3", makeFY5253NearestEndMonthQuarter(weekday=0, startingMonth=12, qtr_with_extra_week=3)),
             ]

    for name, expected in pairs:
        offset = get_offset(name)
        assert offset == expected, ("Expected %r to yield %r (actual: %r)" %
                                    (name, expected, offset))


def test_parse_time_string():
    (date, parsed, reso) = parse_time_string('4Q1984')
    (date_lower, parsed_lower, reso_lower) = parse_time_string('4q1984')
    assert date == date_lower
    assert parsed == parsed_lower
    assert reso == reso_lower


def test_get_standard_freq():
    fstr = get_standard_freq('W')
    assert fstr == get_standard_freq('w')
    assert fstr == get_standard_freq('1w')
    assert fstr == get_standard_freq(('W', 1))
    assert fstr == get_standard_freq('WeEk')

    fstr = get_standard_freq('5Q')
    assert fstr == get_standard_freq('5q')
    assert fstr == get_standard_freq('5QuarTer')
    assert fstr == get_standard_freq(('q', 5))


def test_quarterly_dont_normalize():
    date = datetime(2012, 3, 31, 5, 30)

    offsets = (QuarterBegin, QuarterEnd, BQuarterEnd, BQuarterBegin)

    for klass in offsets:
        result = date + klass()
        assert(result.time() == date.time())


class TestOffsetAliases(tm.TestCase):

    def setUp(self):
        _offset_map.clear()

    def test_alias_equality(self):
        for k, v in compat.iteritems(_offset_map):
            if v is None:
                continue
            self.assertEqual(k, v.copy())

    def test_rule_code(self):
        lst = ['M', 'MS', 'BM', 'BMS', 'D', 'B', 'H', 'T', 'S', 'L', 'U']
        for k in lst:
            self.assertEqual(k, get_offset(k).rule_code)
            # should be cached - this is kind of an internals test...
            assert k in _offset_map
            self.assertEqual(k, (get_offset(k) * 3).rule_code)

        suffix_lst = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
        base = 'W'
        for v in suffix_lst:
            alias = '-'.join([base, v])
            self.assertEqual(alias, get_offset(alias).rule_code)
            self.assertEqual(alias, (get_offset(alias) * 5).rule_code)

        suffix_lst = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG',
                      'SEP', 'OCT', 'NOV', 'DEC']
        base_lst = ['A', 'AS', 'BA', 'BAS', 'Q', 'QS', 'BQ', 'BQS']
        for base in base_lst:
            for v in suffix_lst:
                alias = '-'.join([base, v])
                self.assertEqual(alias, get_offset(alias).rule_code)
                self.assertEqual(alias, (get_offset(alias) * 5).rule_code)


def test_apply_ticks():
    result = offsets.Hour(3).apply(offsets.Hour(4))
    exp = offsets.Hour(7)
    assert(result == exp)


def test_delta_to_tick():
    delta = timedelta(3)

    tick = offsets._delta_to_tick(delta)
    assert(tick == offsets.Day(3))


def test_dateoffset_misc():
    oset = offsets.DateOffset(months=2, days=4)
    # it works
    result = oset.freqstr

    assert(not offsets.DateOffset(months=2) == 2)


def test_freq_offsets():
    off = BDay(1, offset=timedelta(0, 1800))
    assert(off.freqstr == 'B+30Min')

    off = BDay(1, offset=timedelta(0, -1800))
    assert(off.freqstr == 'B-30Min')


def get_all_subclasses(cls):
    ret = set()
    this_subclasses = cls.__subclasses__()
    ret = ret | set(this_subclasses)
    for this_subclass in this_subclasses:
        ret | get_all_subclasses(this_subclass)
    return ret

class TestCaching(tm.TestCase):
    no_simple_ctr = [WeekOfMonth, FY5253,
                     FY5253Quarter,
                     LastWeekOfMonth]

    def test_should_cache_month_end(self):
        self.assertTrue(MonthEnd()._should_cache())

    def test_should_cache_bmonth_end(self):
        self.assertTrue(BusinessMonthEnd()._should_cache())

    def test_should_cache_week_month(self):
        self.assertTrue(WeekOfMonth(weekday=1, week=2)._should_cache())

    def test_all_cacheableoffsets(self):
        for subclass in get_all_subclasses(CacheableOffset):
            if subclass.__name__[0] == "_" \
                or subclass in TestCaching.no_simple_ctr:
                continue
            self.run_X_index_creation(subclass)

    def setUp(self):
        _daterange_cache.clear()
        _offset_map.clear()

    def run_X_index_creation(self, cls):
        inst1 = cls()
        if not inst1.isAnchored():
            self.assertFalse(inst1._should_cache(), cls)
            return

        self.assertTrue(inst1._should_cache(), cls)

        DatetimeIndex(start=datetime(2013,1,31), end=datetime(2013,3,31), freq=inst1, normalize=True)
        self.assertTrue(cls() in _daterange_cache, cls)

    def test_month_end_index_creation(self):
        DatetimeIndex(start=datetime(2013,1,31), end=datetime(2013,3,31), freq=MonthEnd(), normalize=True)
        self.assertTrue(MonthEnd() in _daterange_cache)

    def test_bmonth_end_index_creation(self):
        DatetimeIndex(start=datetime(2013,1,31), end=datetime(2013,3,29), freq=BusinessMonthEnd(), normalize=True)
        self.assertTrue(BusinessMonthEnd() in _daterange_cache)

    def test_week_of_month_index_creation(self):
        inst1 = WeekOfMonth(weekday=1, week=2)
        DatetimeIndex(start=datetime(2013,1,31), end=datetime(2013,3,29), freq=inst1, normalize=True)
        inst2 = WeekOfMonth(weekday=1, week=2)
        self.assertTrue(inst2 in _daterange_cache)


class TestReprNames(tm.TestCase):
    def test_str_for_named_is_name(self):
        # look at all the amazing combinations!
        month_prefixes = ['A', 'AS', 'BA', 'BAS', 'Q', 'BQ', 'BQS', 'QS']
        names = [prefix + '-' + month for prefix in month_prefixes
                    for month in ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                                  'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']]
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
        names += ['W-' + day for day in days]
        names += ['WOM-' + week + day for week in ('1', '2', '3', '4')
                                       for day in days]
        #singletons
        names += ['S', 'T', 'U', 'BM', 'BMS', 'BQ', 'QS'] # No 'Q'
        _offset_map.clear()
        for name in names:
            offset = get_offset(name)
            self.assertEqual(repr(offset), name)
            self.assertEqual(str(offset), name)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
