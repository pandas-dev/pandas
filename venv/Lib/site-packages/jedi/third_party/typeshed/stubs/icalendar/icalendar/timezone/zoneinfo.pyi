import datetime
from typing import Final, Literal
from zoneinfo import ZoneInfo

from dateutil.rrule import rrule

from ..cal import Timezone
from ..prop import vRecur
from .provider import TZProvider

__all__ = ["ZONEINFO"]

class ZONEINFO(TZProvider):
    @property
    def name(self) -> Literal["zoneinfo"]: ...
    utc: Final[ZoneInfo]
    def localize(self, dt: datetime.datetime, tz: ZoneInfo) -> datetime.datetime: ...  # type: ignore[override]
    def localize_utc(self, dt: datetime.datetime) -> datetime.datetime: ...
    def timezone(self, name: str) -> datetime.tzinfo | None: ...
    def knows_timezone_id(self, id: str) -> bool: ...
    def fix_rrule_until(self, rrule: rrule, ical_rrule: vRecur) -> None: ...
    def create_timezone(self, tz: Timezone) -> datetime.tzinfo: ...  # type: ignore[override]
    def uses_pytz(self) -> Literal[False]: ...
    def uses_zoneinfo(self) -> Literal[True]: ...
