from dateparser.calendars import CalendarBase
from dateparser.calendars.hijri_parser import hijri_parser

class HijriCalendar(CalendarBase):
    parser: type[hijri_parser]
