
from datetime import datetime
import pandas.util.testing as tm
from pandas.tseries.holiday import (
    USFederalHolidayCalendar, USMemorialDay, USThanksgivingDay,
    nearest_workday, next_monday_or_tuesday, next_monday,
    previous_friday, sunday_to_monday, Holiday, DateOffset,
    MO, Timestamp, AbstractHolidayCalendar, get_calendar,
    HolidayCalendarFactory)

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

class TestHoliday(tm.TestCase):

    def setUp(self):
        self.start_date = datetime(2011, 1, 1)
        self.end_date   = datetime(2020, 12, 31)

    def test_usmemorialday(self):
        holidays = USMemorialDay.dates(self.start_date,
                                       self.end_date)
        holidayList = [
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
                       ]
        self.assertEqual(list(holidays), holidayList)

    def test_usthanksgivingday(self):
        holidays = USThanksgivingDay.dates(self.start_date,
                                           self.end_date)
        holidayList = [
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
                       ]

        self.assertEqual(list(holidays), holidayList)

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

