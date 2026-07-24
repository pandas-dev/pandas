import datetime
from typing import Final

from ..core import ChineseNewYearCalendar

holidays: Final[dict[int, dict[str, list[tuple[int, int]]]]]
workdays: Final[dict[int, dict[str, list[tuple[int, int]]]]]

class China(ChineseNewYearCalendar):
    extra_working_days: list[datetime.date]
    def __init__(self, *args, **kwargs) -> None: ...
