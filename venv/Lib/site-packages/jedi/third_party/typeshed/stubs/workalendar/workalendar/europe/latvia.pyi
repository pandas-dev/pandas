import datetime

from ..core import WesternCalendar

class Latvia(WesternCalendar):
    def get_independence_days(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_republic_days(self, year: int) -> list[tuple[datetime.date, str]]: ...
