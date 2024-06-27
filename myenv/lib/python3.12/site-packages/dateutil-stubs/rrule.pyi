import datetime
from _typeshed import Incomplete
from collections.abc import Iterable, Iterator, Sequence
from typing_extensions import TypeAlias

from ._common import weekday as weekdaybase

YEARLY: int
MONTHLY: int
WEEKLY: int
DAILY: int
HOURLY: int
MINUTELY: int
SECONDLY: int

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
    def __init__(self, cache: bool = False) -> None: ...
    def __iter__(self) -> Iterator[datetime.datetime]: ...
    def __getitem__(self, item): ...
    def __contains__(self, item): ...
    def count(self): ...
    def before(self, dt, inc: bool = False): ...
    def after(self, dt, inc: bool = False): ...
    def xafter(self, dt, count: Incomplete | None = None, inc: bool = False): ...
    def between(self, after, before, inc: bool = False, count: int = 1): ...

class rrule(rrulebase):
    def __init__(
        self,
        freq,
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
        cache: bool = False,
    ) -> None: ...
    def replace(self, **kwargs): ...

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
    def rebuild(self, year, month): ...
    def ydayset(self, year, month, day): ...
    def mdayset(self, year, month, day): ...
    def wdayset(self, year, month, day): ...
    def ddayset(self, year, month, day): ...
    def htimeset(self, hour, minute, second): ...
    def mtimeset(self, hour, minute, second): ...
    def stimeset(self, hour, minute, second): ...

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

    def __init__(self, cache: bool = False) -> None: ...
    def rrule(self, rrule: _RRule): ...
    def rdate(self, rdate): ...
    def exrule(self, exrule): ...
    def exdate(self, exdate): ...

class _rrulestr:
    def __call__(self, s, **kwargs) -> rrule | rruleset: ...

rrulestr: _rrulestr
