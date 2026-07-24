import datetime
import sys
import winreg
from _typeshed import Incomplete
from typing import Any, SupportsIndex

if sys.platform == "win32":
    handle: winreg.HKEYType
    tzparent: winreg.HKEYType
    parentsize: int
    localkey: winreg.HKEYType
    WEEKS: datetime.timedelta
    def list_timezones() -> list[str]: ...

    class win32tz(datetime.tzinfo):
        data: win32tz_data
        def __init__(self, name: str | None) -> None: ...
        def utcoffset(self, dt: datetime.datetime) -> datetime.timedelta: ...  # type: ignore[override]
        def dst(self, dt: datetime.datetime) -> datetime.timedelta: ...  # type: ignore[override]
        def tzname(self, dt: datetime.datetime) -> str | None: ...  # type: ignore[override]

    def pickNthWeekday(
        year: SupportsIndex, month: SupportsIndex, dayofweek: int, hour: SupportsIndex, minute: SupportsIndex, whichweek: int
    ) -> datetime.datetime | None: ...

    class win32tz_data:
        display: Incomplete
        dstname: Incomplete
        stdname: Incomplete
        stdoffset: Incomplete
        dstoffset: Incomplete
        stdmonth: Incomplete
        stddayofweek: Incomplete
        stdweeknumber: Incomplete
        stdhour: Incomplete
        stdminute: Incomplete
        dstmonth: Incomplete
        dstdayofweek: Incomplete
        dstweeknumber: Incomplete
        dsthour: Incomplete
        dstminute: Incomplete
        def __init__(self, path: str | None) -> None: ...

    def valuesToDict(key: winreg._KeyType) -> dict[str, Any]: ...
