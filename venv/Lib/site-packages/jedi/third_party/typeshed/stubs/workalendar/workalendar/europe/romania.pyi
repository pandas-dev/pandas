import datetime

from ..core import OrthodoxCalendar

class Romania(OrthodoxCalendar):
    def get_childrens_day(self, year: int) -> list[tuple[datetime.date, str]]: ...
    def get_liberation_day(self, year: int) -> list[tuple[datetime.date, str]]: ...
