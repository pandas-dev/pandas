from typing import ClassVar

from ..core import WesternCalendar

class Cyprus(WesternCalendar):
    include_christmas_day: ClassVar[bool]
