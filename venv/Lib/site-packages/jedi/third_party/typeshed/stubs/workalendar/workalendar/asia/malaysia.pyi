import datetime
from typing import ClassVar

from ..core import ChineseNewYearCalendar, IslamicMixin

class Malaysia(IslamicMixin, ChineseNewYearCalendar):
    MSIA_DEEPAVALI: ClassVar[dict[int, datetime.date]]
    MSIA_THAIPUSAM: ClassVar[dict[int, datetime.date]]
