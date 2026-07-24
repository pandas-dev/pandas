from dateparser.calendars.jalali_parser import jalali_parser

from . import CalendarBase

class JalaliCalendar(CalendarBase):
    parser: type[jalali_parser]
