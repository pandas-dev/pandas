import datetime

from ..core import OrthodoxCalendar

class Belarus(OrthodoxCalendar):
    def get_radonitsa(self, year: int) -> datetime.date: ...
