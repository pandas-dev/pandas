from typing import ClassVar

from ..core import IslamoWesternCalendar

class Nigeria(IslamoWesternCalendar):
    shift_sunday_holidays: ClassVar[bool]
