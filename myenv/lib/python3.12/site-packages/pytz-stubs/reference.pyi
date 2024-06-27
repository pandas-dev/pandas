import datetime

from pytz import UTC as UTC

class FixedOffset(datetime.tzinfo):
    def __init__(self, offset: float, name: str) -> None: ...
    def utcoffset(self, dt: datetime.datetime | None) -> datetime.timedelta: ...
    def tzname(self, dt: datetime.datetime | None) -> str: ...
    def dst(self, dt: datetime.datetime | None) -> datetime.timedelta: ...

STDOFFSET: datetime.timedelta
DSTOFFSET: datetime.timedelta

class LocalTimezone(datetime.tzinfo):
    def utcoffset(self, dt: datetime.datetime) -> datetime.timedelta: ...  # type: ignore[override]
    def dst(self, dt: datetime.datetime) -> datetime.timedelta: ...  # type: ignore[override]
    def tzname(self, dt: datetime.datetime) -> str: ...  # type: ignore[override]

Local: LocalTimezone
DSTSTART: datetime.datetime
DSTEND: datetime.datetime

def first_sunday_on_or_after(dt: datetime.datetime) -> datetime.datetime: ...

class USTimeZone(datetime.tzinfo):
    stdoffset: datetime.timedelta
    reprname: str
    stdname: str
    dstname: str
    def __init__(self, hours: float, reprname: str, stdname: str, dstname: str) -> None: ...
    def tzname(self, dt: datetime.datetime | None) -> str: ...
    def utcoffset(self, dt: datetime.datetime | None) -> datetime.timedelta: ...
    def dst(self, dt: datetime.datetime | None) -> datetime.timedelta: ...

Eastern: USTimeZone
Central: USTimeZone
Mountain: USTimeZone
Pacific: USTimeZone
