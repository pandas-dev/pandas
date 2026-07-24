import datetime

from ..core import WesternCalendar

class Denmark(WesternCalendar):
    def get_store_bededag(self, year: int) -> datetime.date: ...
