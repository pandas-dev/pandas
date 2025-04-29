from datetime import date, timedelta
from typing import SupportsFloat, TypeVar, overload
from typing_extensions import Self, TypeAlias

# See #9817 for why we reexport this here
from ._common import weekday as weekday

_DateT = TypeVar("_DateT", bound=date)
# Work around attribute and type having the same name.
_Weekday: TypeAlias = weekday

MO: weekday
TU: weekday
WE: weekday
TH: weekday
FR: weekday
SA: weekday
SU: weekday

class relativedelta:
    years: int
    months: int
    days: int
    leapdays: int
    hours: int
    minutes: int
    seconds: int
    microseconds: int
    year: int | None
    month: int | None
    weekday: _Weekday | None
    day: int | None
    hour: int | None
    minute: int | None
    second: int | None
    microsecond: int | None
    def __init__(
        self,
        dt1: date | None = None,
        dt2: date | None = None,
        years: int = 0,
        months: int = 0,
        days: int = 0,
        leapdays: int = 0,
        weeks: int = 0,
        hours: int = 0,
        minutes: int = 0,
        seconds: int = 0,
        microseconds: int = 0,
        year: int | None = None,
        month: int | None = None,
        day: int | None = None,
        weekday: int | _Weekday | None = None,
        yearday: int | None = None,
        nlyearday: int | None = None,
        hour: int | None = None,
        minute: int | None = None,
        second: int | None = None,
        microsecond: int | None = None,
    ) -> None: ...
    @property
    def weeks(self) -> int: ...
    @weeks.setter
    def weeks(self, value: int) -> None: ...
    def normalized(self) -> Self: ...
    @overload
    def __add__(self, other: timedelta | relativedelta) -> Self: ...
    @overload
    def __add__(self, other: _DateT) -> _DateT: ...
    @overload
    def __radd__(self, other: timedelta | relativedelta) -> Self: ...
    @overload
    def __radd__(self, other: _DateT) -> _DateT: ...
    @overload
    def __rsub__(self, other: timedelta | relativedelta) -> Self: ...
    @overload
    def __rsub__(self, other: _DateT) -> _DateT: ...
    def __sub__(self, other: relativedelta) -> Self: ...
    def __neg__(self) -> Self: ...
    def __bool__(self) -> bool: ...
    def __nonzero__(self) -> bool: ...
    def __mul__(self, other: SupportsFloat) -> Self: ...
    def __rmul__(self, other: SupportsFloat) -> Self: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __div__(self, other: SupportsFloat) -> Self: ...
    def __truediv__(self, other: SupportsFloat) -> Self: ...
    def __abs__(self) -> Self: ...
    def __hash__(self) -> int: ...
