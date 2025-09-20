import re
from _typeshed import Incomplete, SupportsRead
from collections.abc import Callable, Mapping
from datetime import _TzInfo, datetime
from io import StringIO
from typing import IO, Any
from typing_extensions import Self, TypeAlias

_FileOrStr: TypeAlias = bytes | str | IO[str] | IO[Any]
_TzData: TypeAlias = _TzInfo | int | str | None
_TzInfos: TypeAlias = Mapping[str, _TzData] | Callable[[str, int], _TzData]

__all__ = ["parse", "parserinfo", "ParserError"]

class _timelex:
    _split_decimal: re.Pattern[str]
    instream: StringIO | SupportsRead[str]
    charstack: list[str]
    tokenstack: list[str]
    eof: bool
    def __init__(self, instream: str | bytes | bytearray | SupportsRead[str]) -> None: ...
    def get_token(self) -> str | None: ...
    def __iter__(self) -> Self: ...
    def __next__(self) -> str: ...
    def next(self) -> str: ...
    @classmethod
    def split(cls, s: str) -> list[str]: ...
    @classmethod
    def isword(cls, nextchar: str) -> bool: ...
    @classmethod
    def isnum(cls, nextchar: str) -> bool: ...
    @classmethod
    def isspace(cls, nextchar: str) -> bool: ...

class _resultbase:
    def __init__(self) -> None: ...
    def _repr(self, classname: str) -> str: ...
    def __len__(self) -> int: ...

class parserinfo:
    JUMP: list[str]
    WEEKDAYS: list[tuple[str, str]]
    MONTHS: list[tuple[str, str] | tuple[str, str, str]]
    HMS: list[tuple[str, str, str]]
    AMPM: list[tuple[str, str]]
    UTCZONE: list[str]
    PERTAIN: list[str]
    TZOFFSET: dict[str, int]
    def __init__(self, dayfirst: bool = False, yearfirst: bool = False) -> None: ...
    def jump(self, name: str) -> bool: ...
    def weekday(self, name: str) -> int | None: ...
    def month(self, name: str) -> int | None: ...
    def hms(self, name: str) -> int | None: ...
    def ampm(self, name: str) -> int | None: ...
    def pertain(self, name: str) -> bool: ...
    def utczone(self, name: str) -> bool: ...
    def tzoffset(self, name: str) -> int | None: ...
    def convertyear(self, year: int, century_specified: bool = False) -> int: ...
    def validate(self, res: datetime) -> bool: ...

class _ymd(list[Incomplete]):
    century_specified: bool
    dstridx: int | None
    mstridx: int | None
    ystridx: int | None
    def __init__(self, *args, **kwargs) -> None: ...
    @property
    def has_year(self) -> bool: ...
    @property
    def has_month(self) -> bool: ...
    @property
    def has_day(self) -> bool: ...
    def could_be_day(self, value): ...
    def append(self, val, label=None): ...
    def _resolve_from_stridxs(self, strids): ...
    def resolve_ymd(self, yearfirst: bool | None, dayfirst: bool | None): ...

class parser:
    info: parserinfo
    def __init__(self, info: parserinfo | None = None) -> None: ...
    def parse(
        self,
        timestr: _FileOrStr,
        default: datetime | None = None,
        ignoretz: bool = False,
        tzinfos: _TzInfos | None = None,
        *,
        dayfirst: bool | None = ...,
        yearfirst: bool | None = ...,
        fuzzy: bool = ...,
        fuzzy_with_tokens: bool = ...,
    ) -> datetime: ...

DEFAULTPARSER: parser

def parse(
    timestr: _FileOrStr,
    parserinfo: parserinfo | None = None,
    *,
    dayfirst: bool | None = ...,
    yearfirst: bool | None = ...,
    ignoretz: bool = ...,
    fuzzy: bool = ...,
    fuzzy_with_tokens: bool = ...,
    default: datetime | None = ...,
    tzinfos: _TzInfos | None = ...,
) -> datetime: ...

class _tzparser:
    class _result(_resultbase):
        __slots__ = ["stdabbr", "stdoffset", "dstabbr", "dstoffset", "start", "end"]
        stdabbr: str | None
        stdoffset: int | None
        dstabbr: str | None
        dstoffset: int | None
        start: _attr
        end: _attr

        class _attr(_resultbase):
            __slots__ = ["month", "week", "weekday", "yday", "jyday", "day", "time"]
            month: int | None
            week: int | None
            weekday: int | None
            yday: int | None
            jyday: int | None
            day: int | None
            time: int | None

        def __init__(self): ...

    def parse(self, tzstr: str | re.Pattern[str]) -> _result | None: ...

DEFAULTTZPARSER: _tzparser

class ParserError(ValueError): ...
class UnknownTimezoneWarning(RuntimeWarning): ...
