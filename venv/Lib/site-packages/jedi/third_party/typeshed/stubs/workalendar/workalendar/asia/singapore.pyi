import datetime
from typing import ClassVar

from ..core import ChineseNewYearCalendar, IslamicMixin, WesternMixin

class Singapore(WesternMixin, IslamicMixin, ChineseNewYearCalendar):
    DEEPAVALI: ClassVar[dict[int, datetime.date]]
