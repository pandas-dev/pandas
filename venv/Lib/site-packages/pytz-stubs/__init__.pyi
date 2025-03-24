import datetime
from _typeshed import Unused
from collections.abc import Mapping
from typing import ClassVar

from .exceptions import (
    AmbiguousTimeError as AmbiguousTimeError,
    InvalidTimeError as InvalidTimeError,
    NonExistentTimeError as NonExistentTimeError,
    UnknownTimeZoneError as UnknownTimeZoneError,
)
from .tzinfo import BaseTzInfo as BaseTzInfo, DstTzInfo, StaticTzInfo

# Actually named UTC and then masked with a singleton with the same name
class _UTCclass(BaseTzInfo):
    def localize(self, dt: datetime.datetime, is_dst: bool | None = False) -> datetime.datetime: ...
    def normalize(self, dt: datetime.datetime, is_dst: bool | None = False) -> datetime.datetime: ...
    def tzname(self, dt: datetime.datetime | None) -> str: ...
    def utcoffset(self, dt: datetime.datetime | None) -> datetime.timedelta: ...
    def dst(self, dt: datetime.datetime | None) -> datetime.timedelta: ...

utc: _UTCclass
UTC: _UTCclass

def timezone(zone: str) -> _UTCclass | StaticTzInfo | DstTzInfo: ...

class _FixedOffset(datetime.tzinfo):
    zone: ClassVar[None]
    def __init__(self, minutes: int) -> None: ...
    def utcoffset(self, dt: Unused) -> datetime.timedelta | None: ...
    def dst(self, dt: Unused) -> datetime.timedelta: ...
    def tzname(self, dt: Unused) -> None: ...
    def localize(self, dt: datetime.datetime, is_dst: bool | None = False) -> datetime.datetime: ...
    def normalize(self, dt: datetime.datetime, is_dst: bool | None = False) -> datetime.datetime: ...

def FixedOffset(offset: int, _tzinfos: dict[int, _FixedOffset] = {}) -> _UTCclass | _FixedOffset: ...

all_timezones: list[str]
all_timezones_set: set[str]
common_timezones: list[str]
common_timezones_set: set[str]
country_timezones: Mapping[str, list[str]]
country_names: Mapping[str, str]
ZERO: datetime.timedelta
HOUR: datetime.timedelta
VERSION: str
