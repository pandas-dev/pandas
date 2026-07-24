from collections.abc import Generator
from datetime import date, datetime, timedelta, tzinfo

total_seconds = timedelta.total_seconds

def dt_to_timestamp(dt: datetime) -> int: ...
def isoparse(iso_str: str) -> datetime: ...
def parse_timedelta(text: str) -> timedelta: ...

parse_td = parse_timedelta

def decimal_relative_time(
    d: datetime, other: datetime | None = None, ndigits: int = 0, cardinalize: bool = True
) -> tuple[float, str]: ...
def relative_time(d: datetime, other: datetime | None = None, ndigits: int = 0) -> str: ...
def strpdate(string: str, format: str) -> date: ...
def daterange(start: date, stop: date, step: int = 1, inclusive: bool = False) -> Generator[date]: ...

ZERO: timedelta
HOUR: timedelta

class ConstantTZInfo(tzinfo):
    name: str
    offset: timedelta
    def __init__(self, name: str = "ConstantTZ", offset: timedelta = ...) -> None: ...
    @property
    def utcoffset_hours(self) -> str: ...
    def utcoffset(self, dt: datetime | None) -> timedelta: ...
    def tzname(self, dt: datetime | None) -> str: ...
    def dst(self, dt: datetime | None) -> timedelta: ...

UTC: ConstantTZInfo
EPOCH_AWARE: datetime

class LocalTZInfo(tzinfo):
    def is_dst(self, dt: datetime) -> bool: ...
    def utcoffset(self, dt: datetime) -> timedelta: ...  # type: ignore[override] # Doesn't support None
    def dst(self, dt: datetime) -> timedelta: ...  # type: ignore[override] # Doesn't support None
    def tzname(self, dt: datetime) -> str: ...  # type: ignore[override] # Doesn't support None

LocalTZ: LocalTZInfo
DSTSTART_2007: datetime
DSTEND_2007: datetime
DSTSTART_1987_2006: datetime
DSTEND_1987_2006: datetime
DSTSTART_1967_1986: datetime
DSTEND_1967_1986: datetime

class USTimeZone(tzinfo):
    stdoffset: timedelta
    reprname: str
    stdname: str
    dstname: str
    def __init__(self, hours: int, reprname: str, stdname: str, dstname: str) -> None: ...
    def tzname(self, dt: datetime | None) -> str: ...
    def utcoffset(self, dt: datetime | None) -> timedelta: ...
    def dst(self, dt: datetime | None) -> timedelta: ...

Eastern: USTimeZone
Central: USTimeZone
Mountain: USTimeZone
Pacific: USTimeZone
