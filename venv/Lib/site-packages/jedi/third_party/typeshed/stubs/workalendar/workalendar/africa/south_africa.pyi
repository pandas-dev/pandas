import datetime

from ..core import WesternCalendar

class SouthAfrica(WesternCalendar):
    def get_easter_monday_or_family_day(self, year: int) -> tuple[datetime.date, str]: ...
