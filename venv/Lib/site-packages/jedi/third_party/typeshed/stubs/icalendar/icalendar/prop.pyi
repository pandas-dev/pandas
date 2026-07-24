import datetime
from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, SupportsKeysAndGetItem, Unused
from collections.abc import Iterable, Iterator
from enum import Enum
from re import Pattern
from typing import Any, ClassVar, Final, Literal, Protocol, SupportsIndex, overload, type_check_only
from typing_extensions import Self, TypeAlias

from .caselessdict import CaselessDict
from .parser import Parameters
from .parser_tools import ICAL_TYPE
from .timezone import tzid_from_dt as tzid_from_dt, tzid_from_tzinfo as tzid_from_tzinfo

__all__ = [
    "DURATION_REGEX",
    "TimeBase",
    "TypesFactory",
    "WEEKDAY_RULE",
    "tzid_from_dt",
    "vBinary",
    "vBoolean",
    "vCalAddress",
    "vCategory",
    "vDDDLists",
    "vDDDTypes",
    "vDate",
    "vDatetime",
    "vDuration",
    "vFloat",
    "vFrequency",
    "vGeo",
    "vInline",
    "vInt",
    "vMonth",
    "vPeriod",
    "vRecur",
    "vText",
    "vTime",
    "vUTCOffset",
    "vUri",
    "vWeekday",
    "tzid_from_tzinfo",
    "vSkip",
]

_PropType: TypeAlias = type[Any]  # any of the v* classes in this file
_PeriodTuple: TypeAlias = tuple[datetime.datetime, datetime.datetime | datetime.timedelta]
_AnyTimeType: TypeAlias = datetime.datetime | datetime.date | datetime.timedelta | datetime.time | _PeriodTuple

@type_check_only
class _vType(Protocol):
    def to_ical(self) -> bytes | str: ...

DURATION_REGEX: Final[Pattern[str]]
WEEKDAY_RULE: Final[Pattern[str]]

class vBinary:
    obj: str
    params: Parameters
    def __init__(self, obj: str | bytes) -> None: ...
    def to_ical(self) -> bytes: ...
    @staticmethod
    def from_ical(ical: ICAL_TYPE) -> bytes: ...
    def __eq__(self, other: object) -> bool: ...

class vBoolean(int):
    BOOL_MAP: Final[CaselessDict[bool]]
    params: Parameters
    def __new__(cls, x: ConvertibleToInt = ..., /, *, params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> Literal[b"TRUE", b"FALSE"]: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> bool: ...

class vText(str):
    __slots__ = ("encoding", "params")
    encoding: str
    params: Parameters
    def __new__(cls, value: ICAL_TYPE, encoding: str = "utf-8", params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...
    ALTREP: property
    LANGUAGE: property
    RELTYPE: property

class vCalAddress(str):
    __slots__ = ("params",)
    params: Parameters
    def __new__(cls, value: ICAL_TYPE, encoding: str = "utf-8", params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...
    @property
    def email(self) -> str: ...
    @property
    def name(self) -> str: ...
    @name.setter
    def name(self, value: str) -> None: ...
    @name.deleter
    def name(self) -> None: ...
    CN: property
    CUTYPE: property
    DELEGATED_FROM: property
    DELEGATED_TO: property
    DIR: property
    LANGUAGE: property
    PARTSTAT: property
    ROLE: property
    RSVP: property
    SENT_BY: property

class vFloat(float):
    params: Parameters
    def __new__(cls, x: ConvertibleToFloat = ..., /, *, params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...

class vInt(int):
    params: Parameters
    def __new__(cls, x: ConvertibleToInt = ..., /, *, params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...

class vDDDLists:
    params: Parameters
    dts: list[vDDDTypes]
    def __init__(self, dt_list: Iterable[_AnyTimeType] | _AnyTimeType) -> None: ...
    def to_ical(self) -> bytes: ...
    @staticmethod
    def from_ical(ical: str, timezone: str | datetime.timezone | None = None) -> list[Incomplete]: ...
    def __eq__(self, other: object) -> bool: ...

class vCategory:
    cats: list[vText]
    params: Parameters
    def __init__(self, c_list: Iterable[ICAL_TYPE] | ICAL_TYPE, params: SupportsKeysAndGetItem[str, str] = {}) -> None: ...
    def __iter__(self) -> Iterator[str]: ...
    def to_ical(self) -> bytes: ...
    @staticmethod
    def from_ical(ical: ICAL_TYPE) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    RANGE: property
    RELATED: property
    TZID: property

class TimeBase:
    params: Parameters
    ignore_for_equality: set[str]
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    RANGE: property
    RELATED: property
    TZID: property

class vDDDTypes(TimeBase):
    params: Parameters
    dt: _AnyTimeType
    def __init__(self, dt: _AnyTimeType) -> None: ...
    def to_ical(self) -> bytes: ...
    @overload
    @classmethod
    def from_ical(cls, ical: Self, timezone: Unused | None = None) -> _AnyTimeType: ...
    # Return type is one of vDuration, vPeriod, vDatetime, vDate, or vTime,
    # depending on the ical string.
    @overload
    @classmethod
    def from_ical(cls, ical: str, timezone: datetime.timezone | str | None = None) -> Any: ...

class vDate(TimeBase):
    dt: datetime.date
    params: Parameters
    def __init__(self, dt: datetime.date) -> None: ...
    def to_ical(self) -> bytes: ...
    @staticmethod
    def from_ical(ical: ICAL_TYPE) -> datetime.date: ...

class vDatetime(TimeBase):
    dt: datetime.datetime
    params: Parameters
    def __init__(self, dt: datetime.datetime, params: SupportsKeysAndGetItem[str, str] = {}) -> None: ...
    def to_ical(self) -> bytes: ...
    @staticmethod
    def from_ical(ical: ICAL_TYPE, timezone: datetime.timezone | str | None = None) -> datetime.datetime: ...

class vDuration(TimeBase):
    td: datetime.timedelta
    params: Parameters
    def __init__(self, td: datetime.timedelta, params: SupportsKeysAndGetItem[str, str] = {}) -> None: ...
    def to_ical(self) -> bytes: ...
    @staticmethod
    def from_ical(ical: str) -> datetime.timedelta: ...
    @property
    def dt(self) -> datetime.timedelta: ...

class vPeriod(TimeBase):
    params: Parameters
    start: datetime.datetime
    end: datetime.datetime
    by_duration: bool
    duration: datetime.timedelta
    def __init__(self, per: _PeriodTuple) -> None: ...
    def overlaps(self, other: vPeriod) -> bool: ...
    def to_ical(self) -> bytes: ...
    # Return type is a tuple of vDuration, vPeriod, vDatetime, vDate, or vTime,
    # depending on the ical string. If the ical string is formed according to
    # the iCalendar specification, this should always return a
    # (datetime, datetime) or a (datetime, timedelta) tuple, but this is not
    # enforced.
    @staticmethod
    def from_ical(ical: str, timezone: datetime.timezone | str | None = None) -> tuple[Any, Any]: ...
    @property
    def dt(self) -> _PeriodTuple: ...
    FBTYPE: property

class vWeekday(str):
    __slots__ = ("params", "relative", "weekday")
    week_days: Final[CaselessDict[int]]
    weekday: Literal["SU", "MO", "TU", "WE", "TH", "FR", "SA"] | None
    relative: int | None
    params: Parameters
    def __new__(cls, value: ICAL_TYPE, encoding: str = "utf-8", params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...

class vFrequency(str):
    __slots__ = ("params",)
    frequencies: Final[CaselessDict[str]]
    params: Parameters
    def __new__(cls, value: ICAL_TYPE, encoding: str = "utf-8", params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...

class vMonth(int):
    params: Parameters
    def __new__(cls, month: vMonth | str | int, params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: vMonth | str | int) -> Self: ...
    @property
    def leap(self) -> bool: ...
    @leap.setter
    def leap(self, value: bool) -> None: ...

class vSkip(vText, Enum):
    OMIT = "OMIT"
    FORWARD = "FORWARD"
    BACKWARD = "BACKWARD"

    def __reduce_ex__(self, _p: Unused) -> tuple[Self, tuple[str]]: ...

# The type of the values depend on the key. Each key maps to a v* class, and
# the allowed types are the types that the corresponding v* class can parse.
class vRecur(CaselessDict[Iterable[Any] | Any]):
    params: Parameters
    frequencies: Final[list[str]]
    canonical_order: ClassVar[tuple[str, ...]]
    types: Final[CaselessDict[_PropType]]
    def __init__(
        self, *args, params: SupportsKeysAndGetItem[str, str] = {}, **kwargs: list[Any] | tuple[Any, ...] | Any
    ) -> None: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def parse_type(cls, key: str, values: str) -> list[Any]: ...  # Returns a list of v* objects
    @classmethod
    def from_ical(cls, ical: vRecur | str) -> Self: ...

class vTime(TimeBase):
    dt: datetime.time | datetime.datetime
    params: Parameters
    @overload
    def __init__(self, dt: datetime.time | datetime.datetime, /) -> None: ...
    # args are passed to the datetime.time() constructor
    @overload
    def __init__(
        self,
        hour: SupportsIndex = ...,
        minute: SupportsIndex = ...,
        second: SupportsIndex = ...,
        microsecond: SupportsIndex = ...,
        tzinfo: datetime.tzinfo | None = ...,
        /,
    ) -> None: ...
    def to_ical(self) -> str: ...
    @staticmethod
    def from_ical(ical: ICAL_TYPE) -> datetime.time: ...

class vUri(str):
    __slots__ = ("params",)
    params: Parameters
    def __new__(cls, value: ICAL_TYPE, encoding: str = "utf-8", params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...

class vGeo:
    latitude: float
    longitude: float
    params: Parameters
    def __init__(self, geo: tuple[float | str, float | str], params: SupportsKeysAndGetItem[str, str] = {}) -> None: ...
    def to_ical(self) -> str: ...
    @staticmethod
    def from_ical(ical: str) -> tuple[float, float]: ...
    def __eq__(self, other: _vType) -> bool: ...  # type: ignore[override]

class vUTCOffset:
    ignore_exceptions: bool
    td: datetime.timedelta
    params: Parameters
    def __init__(self, td: datetime.timedelta, params: SupportsKeysAndGetItem[str, str] = {}) -> None: ...
    def to_ical(self) -> str: ...
    @classmethod
    def from_ical(cls, ical: Self | ICAL_TYPE) -> datetime.timedelta: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

class vInline(str):
    __slots__ = ("params",)
    params: Parameters
    def __new__(cls, value: ICAL_TYPE, encoding: str = "utf-8", params: SupportsKeysAndGetItem[str, str] = {}) -> Self: ...
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, ical: ICAL_TYPE) -> Self: ...

class TypesFactory(CaselessDict[_PropType]):
    all_types: tuple[_PropType, ...]
    types_map: CaselessDict[str]
    def for_property(self, name: str) -> _PropType: ...
    # value is str | bytes, depending on what the v* class supports
    def to_ical(self, name: str, value: Any) -> bytes: ...
    # value and return type depend on what the v* class supports
    def from_ical(self, name: str, value: Any) -> Any: ...
