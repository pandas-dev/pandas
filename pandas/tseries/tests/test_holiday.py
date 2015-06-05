
from datetime import datetime
import pandas.util.testing as tm
from pandas import DatetimeIndex
from pandas.tseries.holiday import (
    USFederalHolidayCalendar, USMemorialDay, USThanksgivingDay,
    nearest_workday, next_monday_or_tuesday, next_monday,
    previous_friday, sunday_to_monday, Holiday, DateOffset,
    MO, Timestamp, AbstractHolidayCalendar, get_calendar,
    HolidayCalendarFactory, next_workday, previous_workday,
    before_nearest_workday, EasterMonday, GoodFriday,
    after_nearest_workday, weekend_to_monday)
from pytz import utc
import nose

class TestCalendar(tm.TestCase):

    def setUp(self):
        self.holiday_list = [
                       datetime(2012, 1, 2),
                       datetime(2012, 1, 16),
                       datetime(2012, 2, 20),
                       datetime(2012, 5, 28),
                       datetime(2012, 7, 4),
                       datetime(2012, 9, 3),
                       datetime(2012, 10, 8),
                       datetime(2012, 11, 12),
                       datetime(2012, 11, 22),
                       datetime(2012, 12, 25)]

        self.start_date = datetime(2012, 1, 1)
        self.end_date = datetime(2012, 12, 31)

    def test_calendar(self):

        calendar = USFederalHolidayCalendar()
        holidays = calendar.holidays(self.start_date,
                                     self.end_date)

        holidays_1 = calendar.holidays(
                        self.start_date.strftime('%Y-%m-%d'),
                        self.end_date.strftime('%Y-%m-%d'))
        holidays_2 = calendar.holidays(
                        Timestamp(self.start_date),
                        Timestamp(self.end_date))

        self.assertEqual(list(holidays.to_pydatetime()),
                         self.holiday_list)
        self.assertEqual(list(holidays_1.to_pydatetime()),
                         self.holiday_list)
        self.assertEqual(list(holidays_2.to_pydatetime()),
                         self.holiday_list)

    def test_calendar_caching(self):
        # Test for issue #9552

        class TestCalendar(AbstractHolidayCalendar):
            def __init__(self, name=None, rules=None):
                super(TestCalendar, self).__init__(
                    name=name,
                    rules=rules
                )

        jan1 = TestCalendar(rules=[Holiday('jan1', year=2015, month=1, day=1)])
        jan2 = TestCalendar(rules=[Holiday('jan2', year=2015, month=1, day=2)])

        tm.assert_index_equal(
            jan1.holidays(),
            DatetimeIndex(['01-Jan-2015'])
        )
        tm.assert_index_equal(
            jan2.holidays(),
            DatetimeIndex(['02-Jan-2015'])
        )


class TestHoliday(tm.TestCase):

    def setUp(self):
        self.start_date = datetime(2011, 1, 1)
        self.end_date   = datetime(2020, 12, 31)

    def check_results(self, holiday, start, end, expected):
        self.assertEqual(list(holiday.dates(start, end)), expected)
        # Verify that timezone info is preserved.
        self.assertEqual(
            list(
                holiday.dates(
                    utc.localize(Timestamp(start)),
                    utc.localize(Timestamp(end)),
                )
            ),
            [utc.localize(dt) for dt in expected],
        )

    def test_usmemorialday(self):
        self.check_results(
            holiday=USMemorialDay,
            start=self.start_date,
            end=self.end_date,
            expected=[
                datetime(2011, 5, 30),
                datetime(2012, 5, 28),
                datetime(2013, 5, 27),
                datetime(2014, 5, 26),
                datetime(2015, 5, 25),
                datetime(2016, 5, 30),
                datetime(2017, 5, 29),
                datetime(2018, 5, 28),
                datetime(2019, 5, 27),
                datetime(2020, 5, 25),
            ],
        )

    def test_non_observed_holiday(self):

        self.check_results(
            Holiday('July 4th Eve', month=7, day=3),
            start="2001-01-01",
            end="2003-03-03",
            expected=[
                Timestamp('2001-07-03 00:00:00'),
                Timestamp('2002-07-03 00:00:00')
            ]
        )

        self.check_results(
            Holiday('July 4th Eve', month=7, day=3, days_of_week=(0, 1, 2, 3)),
            start="2001-01-01",
            end="2008-03-03",
            expected=[
                Timestamp('2001-07-03 00:00:00'),
                Timestamp('2002-07-03 00:00:00'),
                Timestamp('2003-07-03 00:00:00'),
                Timestamp('2006-07-03 00:00:00'),
                Timestamp('2007-07-03 00:00:00'),
            ]
        )

    def test_easter(self):

        self.check_results(
            EasterMonday,
            start=self.start_date,
            end=self.end_date,
            expected=[
                Timestamp('2011-04-25 00:00:00'),
                Timestamp('2012-04-09 00:00:00'),
                Timestamp('2013-04-01 00:00:00'),
                Timestamp('2014-04-21 00:00:00'),
                Timestamp('2015-04-06 00:00:00'),
                Timestamp('2016-03-28 00:00:00'),
                Timestamp('2017-04-17 00:00:00'),
                Timestamp('2018-04-02 00:00:00'),
                Timestamp('2019-04-22 00:00:00'),
                Timestamp('2020-04-13 00:00:00'),
            ],
        )
        self.check_results(
            GoodFriday,
            start=self.start_date,
            end=self.end_date,
            expected=[
                Timestamp('2011-04-22 00:00:00'),
                Timestamp('2012-04-06 00:00:00'),
                Timestamp('2013-03-29 00:00:00'),
                Timestamp('2014-04-18 00:00:00'),
                Timestamp('2015-04-03 00:00:00'),
                Timestamp('2016-03-25 00:00:00'),
                Timestamp('2017-04-14 00:00:00'),
                Timestamp('2018-03-30 00:00:00'),
                Timestamp('2019-04-19 00:00:00'),
                Timestamp('2020-04-10 00:00:00'),
            ],
        )

    def test_usthanksgivingday(self):

        self.check_results(
            USThanksgivingDay,
            start=self.start_date,
            end=self.end_date,
            expected=[
                datetime(2011, 11, 24),
                datetime(2012, 11, 22),
                datetime(2013, 11, 28),
                datetime(2014, 11, 27),
                datetime(2015, 11, 26),
                datetime(2016, 11, 24),
                datetime(2017, 11, 23),
                datetime(2018, 11, 22),
                datetime(2019, 11, 28),
                datetime(2020, 11, 26),
            ],
        )

    def test_argument_types(self):
        holidays = USThanksgivingDay.dates(self.start_date,
                                           self.end_date)

        holidays_1 = USThanksgivingDay.dates(
                        self.start_date.strftime('%Y-%m-%d'),
                        self.end_date.strftime('%Y-%m-%d'))

        holidays_2 = USThanksgivingDay.dates(
                        Timestamp(self.start_date),
                        Timestamp(self.end_date))

        self.assertEqual(holidays, holidays_1)
        self.assertEqual(holidays, holidays_2)

    def test_special_holidays(self):
        base_date = [datetime(2012, 5, 28)]
        holiday_1 = Holiday('One-Time', year=2012, month=5, day=28)
        holiday_2 = Holiday('Range', month=5, day=28,
                            start_date=datetime(2012, 1, 1),
                            end_date=datetime(2012, 12, 31),
                            offset=DateOffset(weekday=MO(1)))

        self.assertEqual(base_date,
                         holiday_1.dates(self.start_date, self.end_date))
        self.assertEqual(base_date,
                         holiday_2.dates(self.start_date, self.end_date))

    def test_get_calendar(self):
        class TestCalendar(AbstractHolidayCalendar):
            rules = []

        calendar = get_calendar('TestCalendar')
        self.assertEqual(TestCalendar, calendar.__class__)

    def test_factory(self):
        class_1 = HolidayCalendarFactory('MemorialDay', AbstractHolidayCalendar,
                                         USMemorialDay)
        class_2 = HolidayCalendarFactory('Thansksgiving', AbstractHolidayCalendar,
                                         USThanksgivingDay)
        class_3 = HolidayCalendarFactory('Combined', class_1, class_2)

        self.assertEqual(len(class_1.rules), 1)
        self.assertEqual(len(class_2.rules), 1)
        self.assertEqual(len(class_3.rules), 2)


class TestObservanceRules(tm.TestCase):

    def setUp(self):
        self.we =   datetime(2014, 4, 9)
        self.th =   datetime(2014, 4, 10)
        self.fr =   datetime(2014, 4, 11)
        self.sa =   datetime(2014, 4, 12)
        self.su =   datetime(2014, 4, 13)
        self.mo =   datetime(2014, 4, 14)
        self.tu =   datetime(2014, 4, 15)

    def test_next_monday(self):
        self.assertEqual(next_monday(self.sa), self.mo)
        self.assertEqual(next_monday(self.su), self.mo)

    def test_next_monday_or_tuesday(self):
        self.assertEqual(next_monday_or_tuesday(self.sa), self.mo)
        self.assertEqual(next_monday_or_tuesday(self.su), self.tu)
        self.assertEqual(next_monday_or_tuesday(self.mo), self.tu)

    def test_previous_friday(self):
        self.assertEqual(previous_friday(self.sa), self.fr)
        self.assertEqual(previous_friday(self.su), self.fr)

    def test_sunday_to_monday(self):
        self.assertEqual(sunday_to_monday(self.su), self.mo)

    def test_nearest_workday(self):
        self.assertEqual(nearest_workday(self.sa), self.fr)
        self.assertEqual(nearest_workday(self.su), self.mo)
        self.assertEqual(nearest_workday(self.mo), self.mo)

    def test_weekend_to_monday(self):
        self.assertEqual(weekend_to_monday(self.sa), self.mo)
        self.assertEqual(weekend_to_monday(self.su), self.mo)
        self.assertEqual(weekend_to_monday(self.mo), self.mo)

    def test_next_workday(self):
        self.assertEqual(next_workday(self.sa), self.mo)
        self.assertEqual(next_workday(self.su), self.mo)
        self.assertEqual(next_workday(self.mo), self.tu)

    def test_previous_workday(self):
        self.assertEqual(previous_workday(self.sa), self.fr)
        self.assertEqual(previous_workday(self.su), self.fr)
        self.assertEqual(previous_workday(self.tu), self.mo)

    def test_before_nearest_workday(self):
        self.assertEqual(before_nearest_workday(self.sa), self.th)
        self.assertEqual(before_nearest_workday(self.su), self.fr)
        self.assertEqual(before_nearest_workday(self.tu), self.mo)
    
    def test_after_nearest_workday(self):
        self.assertEqual(after_nearest_workday(self.sa), self.mo)
        self.assertEqual(after_nearest_workday(self.su), self.tu)
        self.assertEqual(after_nearest_workday(self.fr), self.mo)
    

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)

