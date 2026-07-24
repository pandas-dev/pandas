from typing import ClassVar

from dateparser.calendars import non_gregorian_parser

class PersianDate:
    year: int
    month: int
    day: int
    def __init__(self, year: int, month: int, day: int) -> None: ...
    def weekday(self) -> int | None: ...

class jalali_parser(non_gregorian_parser):
    # `calendar_converter` is `convertdate.persian` module
    default_year: ClassVar[int]
    default_month: ClassVar[int]
    default_day: ClassVar[int]
    non_gregorian_date_cls: ClassVar[type[PersianDate]]
    def handle_two_digit_year(self, year: int) -> int: ...
