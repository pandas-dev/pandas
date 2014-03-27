
from datetime import datetime
import pandas.util.testing as tm
from pandas.tseries.holiday import (
    USFederalHolidayCalendar, USMemorialDay, USThanksgivingDay)

class TestCalendar(tm.TestCase):
    
    def test_calendar(self):

        calendar = USFederalHolidayCalendar()
        holidays = calendar.holidays(datetime(2012, 1, 1), datetime(2012, 12, 31))
        
        holidayList = [
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
        
        self.assertEqual(list(holidays.to_pydatetime()), holidayList)
        
class TestHoliday(tm.TestCase):
    
    def setUp(self):
        self.start_date = datetime(2011, 1, 1)
        self.end_date   = datetime(2020, 12, 31)
    
    def test_usmemorialday(self):
        holidays = USMemorialDay.dates(self.start_date, self.end_date)
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
        self.assertEqual(holidays, holidayList)
        
    def test_usthanksgivingday(self):
        holidays = USThanksgivingDay.dates(self.start_date, self.end_date)
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
        self.assertEqual(holidays, holidayList)