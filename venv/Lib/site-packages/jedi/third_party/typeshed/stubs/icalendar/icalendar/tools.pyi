import datetime
from typing import Final
from typing_extensions import TypeGuard, TypeIs, deprecated

from pytz.tzinfo import BaseTzInfo

from .prop import vText

__all__ = ["UIDGenerator", "is_date", "is_datetime", "to_datetime", "is_pytz", "is_pytz_dt", "normalize_pytz"]

class UIDGenerator:
    chars: Final[list[str]]
    @staticmethod
    @deprecated("Use the Python standard library's :func:`uuid.uuid4` instead.")
    def rnd_string(length: int = 16) -> str: ...
    @staticmethod
    @deprecated("Use the Python standard library's :func:`uuid.uuid5` instead.")
    def uid(host_name: str = "example.com", unique: str = "") -> vText: ...

def is_date(dt: datetime.date) -> bool: ...  # and not datetime.date
def is_datetime(dt: datetime.date) -> TypeIs[datetime.datetime]: ...
def to_datetime(dt: datetime.date) -> datetime.datetime: ...
def is_pytz(tz: datetime.tzinfo) -> TypeIs[BaseTzInfo]: ...
def is_pytz_dt(dt: datetime.date) -> TypeGuard[datetime.datetime]: ...  # and dt.tzinfo is BaseTZInfo
def normalize_pytz(dt: datetime.date) -> datetime.datetime: ...
