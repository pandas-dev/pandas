import datetime

from ..core import IslamicCalendar

class Turkey(IslamicCalendar):
    def get_delta_islamic_holidays(self, year: int) -> datetime.timedelta: ...
