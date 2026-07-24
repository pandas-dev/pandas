import datetime

from ..core import WesternCalendar

class Iceland(WesternCalendar):
    def get_first_day_of_summer(self, year: int) -> datetime.date: ...
    def get_commerce_day(self, year: int) -> datetime.date: ...
