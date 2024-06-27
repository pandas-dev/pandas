from collections.abc import Callable, Mapping
from datetime import datetime, tzinfo
from typing import IO, Any
from typing_extensions import TypeAlias

from .isoparser import isoparse as isoparse, isoparser as isoparser

_FileOrStr: TypeAlias = bytes | str | IO[str] | IO[Any]
_TzData: TypeAlias = tzinfo | int | str | None
_TzInfo: TypeAlias = Mapping[str, _TzData] | Callable[[str, int], _TzData]

class parserinfo:
    JUMP: list[str]
    WEEKDAYS: list[tuple[str, ...]]
    MONTHS: list[tuple[str, ...]]
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
    def convertyear(self, year: int) -> int: ...
    def validate(self, res: datetime) -> bool: ...

class parser:
    def __init__(self, info: parserinfo | None = None) -> None: ...
    def parse(
        self,
        timestr: _FileOrStr,
        default: datetime | None = None,
        ignoretz: bool = False,
        tzinfos: _TzInfo | None = None,
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
    tzinfos: _TzInfo | None = ...,
) -> datetime: ...

class _tzparser: ...

DEFAULTTZPARSER: _tzparser

class ParserError(ValueError): ...
