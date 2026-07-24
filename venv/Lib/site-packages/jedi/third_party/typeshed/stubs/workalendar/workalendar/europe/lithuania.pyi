import datetime

from ..core import WesternCalendar

class Lithuania(WesternCalendar):
    def get_mothers_day(self, year: int) -> tuple[datetime.date, str]: ...
    def get_fathers_day(self, year: int) -> tuple[datetime.date, str]: ...
