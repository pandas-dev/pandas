import datetime
from collections.abc import Iterable

from ..core import ChineseNewYearCalendar

class Taiwan(ChineseNewYearCalendar):
    def is_working_day(
        self,
        day: datetime.date | datetime.datetime,
        extra_working_days: Iterable[datetime.date | datetime.datetime] | None = None,
        extra_holidays: Iterable[datetime.date | datetime.datetime] | None = None,
    ) -> bool: ...
