import datetime
from typing import Any
from typing_extensions import TypeAlias

from ..core import Calendar

_HebrewDate: TypeAlias = Any | datetime.date  # from `pyluach.dates` package

class Israel(Calendar):
    def get_hebrew_independence_day(self, jewish_year: int) -> list[tuple[_HebrewDate, str]]: ...
