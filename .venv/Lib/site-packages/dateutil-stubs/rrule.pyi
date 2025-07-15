import datetime
from _typeshed import Incomplete
from collections.abc import Generator, Iterable, Iterator, Sequence
from typing import Final, Literal
from typing_extensions import Self, TypeAlias

from ._common import weekday as weekdaybase

__all__ = [
    "rrule",
    "rruleset",
    "rrulestr",
    "YEARLY",
    "MONTHLY",
    "WEEKLY",
    "DAILY",
    "HOURLY",
    "MINUTELY",
    "SECONDLY",
    "MO",
    "TU",
    "WE",
    "TH",
    "FR",
    "SA",
    "SU",
]

M366MASK: Final[tuple[int, ...]]
MDAY366MASK: Final[tuple[int, ...]]
MDAY365MASK: Final[tuple[int, ...]]
NMDAY366MASK: Final[tuple[int, ...]]
NMDAY365MASK: Final[list[int]]
M366RANGE: Final[tuple[int, ...]]
M365RANGE: Final[tuple[int, ...]]
WDAYMASK: Final[list[int]]
M365MASK: Final[tuple[int, ...]]
FREQNAMES: Final[list[str]]
YEARLY: Final = 0
MONTHLY: Final = 1
WEEKLY: Final = 2
DAILY: Final = 3
HOURLY: Final = 4
MINUTELY: Final = 5
SECONDLY: Final = 6

class weekday(weekdaybase): ...

weekdays: tuple[weekday, weekday, weekday, weekday, weekday, weekday, weekday]
MO: weekday
TU: weekday
WE: weekday
TH: weekday
FR: weekday
SA: weekday
SU: weekday

class rrulebase:
    def __init__(self, cache: bool | None = False) -> None: ...
    def __iter__(self) -> Iterator[datetime.datetime]: ...
    def __getitem__(self, item): ...
    def __contains__(self, item) -> bool: ...
    def count(self) -> int | None: ...
    def before(self, dt, inc: bool = False): ...
    def after(self, dt, inc: bool = False): ...
    def xafter(self, dt, count=None, inc: bool = False) -> Generator[Incomplete]: ...
    def between(self, after, before, inc: bool = False, count: int = 1) -> list[Incomplete]: ...

class rrule(rrulebase):
    def __init__(
        self,
        freq: Literal[0, 1, 2, 3, 4, 5, 6],
        dtstart: datetime.date | None = None,
        interval: int = 1,
        wkst: weekday | int | None = None,
        count: int | None = None,
        until: datetime.date | int | None = None,
        bysetpos: int | Iterable[int] | None = None,
        bymonth: int | Iterable[int] | None = None,
        bymonthday: int | Iterable[int] | None = None,
        byyearday: int | Iterable[int] | None = None,
        byeaster: int | Iterable[int] | None = None,
        byweekno: int | Iterable[int] | None = None,
        byweekday: int | weekday | Iterable[int] | Iterable[weekday] | None = None,
        byhour: int | Iterable[int] | None = None,
        byminute: int | Iterable[int] | None = None,
        bysecond: int | Iterable[int] | None = None,
        cache: bool | None = False,
    ) -> None: ...
    def replace(
        self,
        *,
        freq: Literal[0, 1, 2, 3, 4, 5, 6] = ...,
        dtstart: datetime.date | None = ...,
        interval: int = ...,
        wkst: weekday | int | None = ...,
        count: int | None = ...,
        until: datetime.date | int | None = ...,
        bysetpos: int | Iterable[int] | None = None,
        bymonth: int | Iterable[int] | None = None,
        bymonthday: int | Iterable[int] | None = None,
        byyearday: int | Iterable[int] | None = None,
        byeaster: int | Iterable[int] | None = None,
        byweekno: int | Iterable[int] | None = None,
        byweekday: int | weekday | Iterable[int] | Iterable[weekday] | None = None,
        byhour: int | Iterable[int] | None = None,
        byminute: int | Iterable[int] | None = None,
        bysecond: int | Iterable[int] | None = None,
        cache: bool | None = ...,
    ) -> Self: ...

_RRule: TypeAlias = rrule

class _iterinfo:
    rrule: _RRule
    def __init__(self, rrule: _RRule) -> None: ...
    yearlen: int | None
    nextyearlen: int | None
    yearordinal: int | None
    yearweekday: int | None
    mmask: Sequence[int] | None
    mdaymask: Sequence[int] | None
    nmdaymask: Sequence[int] | None
    wdaymask: Sequence[int] | None
    mrange: Sequence[int] | None
    wnomask: Sequence[int] | None
    nwdaymask: Sequence[int] | None
    eastermask: Sequence[int] | None
    lastyear: int | None
    lastmonth: int | None
    def rebuild(self, year: int, month: int) -> None: ...
    def ydayset(self, year: int, month: int, day: int) -> tuple[Iterable[int | None], int, int]: ...
    def mdayset(self, year: int, month: int, day: int) -> tuple[Iterable[int | None], int, int]: ...
    def wdayset(self, year: int, month: int, day: int) -> tuple[Iterable[int | None], int, int]: ...
    def ddayset(self, year: int, month: int, day: int) -> tuple[Iterable[int | None], int, int]: ...
    def htimeset(self, hour: int, minute: int, second: int) -> list[datetime.time]: ...
    def mtimeset(self, hour: int, minute: int, second: int) -> list[datetime.time]: ...
    def stimeset(self, hour: int, minute: int, second: int) -> tuple[datetime.time, ...]: ...

class rruleset(rrulebase):
    class _genitem:
        dt: Incomplete
        genlist: list[Incomplete]
        gen: Incomplete
        def __init__(self, genlist, gen) -> None: ...
        def __next__(self) -> None: ...
        next = __next__
        def __lt__(self, other) -> bool: ...
        def __gt__(self, other) -> bool: ...
        def __eq__(self, other) -> bool: ...
        def __ne__(self, other) -> bool: ...

    def __init__(self, cache: bool | None = False) -> None: ...
    def rrule(self, rrule: _RRule) -> None: ...
    def rdate(self, rdate) -> None: ...
    def exrule(self, exrule) -> None: ...
    def exdate(self, exdate) -> None: ...

class _rrulestr:
    def __call__(
        self,
        s: str,
        *,
        dtstart: datetime.date | None = None,
        cache: bool | None = False,
        unfold: bool = False,
        forceset: bool = False,
        compatible: bool = False,
        ignoretz: bool = False,
        tzids=None,
        tzinfos=None,
    ) -> rrule | rruleset: ...

rrulestr: _rrulestr
