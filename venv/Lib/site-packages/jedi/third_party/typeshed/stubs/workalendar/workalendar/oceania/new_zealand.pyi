import datetime

from ..core import WesternCalendar

class NewZealand(WesternCalendar):
    def get_queens_birthday(self, year: int) -> tuple[datetime.date, str]: ...
    def get_labour_day(self, year: int) -> tuple[datetime.date, str]: ...
