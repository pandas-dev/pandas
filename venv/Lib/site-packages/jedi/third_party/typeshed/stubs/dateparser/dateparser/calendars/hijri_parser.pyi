from typing import ClassVar

from dateparser.calendars import non_gregorian_parser

class hijri:
    @classmethod
    def to_gregorian(cls, year: int | None = None, month: int | None = None, day: int | None = None) -> tuple[int, int, int]: ...
    @classmethod
    def from_gregorian(
        cls, year: int | None = None, month: int | None = None, day: int | None = None
    ) -> tuple[int, int, int]: ...
    @classmethod
    def month_length(cls, year: int | None, month: int | None) -> int: ...

class HijriDate:
    year: int
    month: int
    day: int
    def __init__(self, year: int, month: int, day: int) -> None: ...
    def weekday(self) -> int | None: ...

class hijri_parser(non_gregorian_parser):
    calendar_converter: ClassVar[type[hijri]]
    default_year: ClassVar[int]
    default_month: ClassVar[int]
    default_day: ClassVar[int]
    non_gregorian_date_cls: ClassVar[type[HijriDate]]
    def handle_two_digit_year(self, year: int) -> int: ...
