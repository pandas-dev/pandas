from typing import ClassVar

from ..core import OrthodoxCalendar

class Ukraine(OrthodoxCalendar):
    shift_sunday_holidays: ClassVar[bool]
