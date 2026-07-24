import datetime
from _typeshed import Incomplete
from collections.abc import Iterable
from datetime import timedelta, tzinfo
from typing import Final, Literal, overload

from .base import Component
from .behavior import Behavior

DATENAMES: Final[tuple[str, ...]]
RULENAMES: Final[tuple[str, ...]]
DATESANDRULES: Final[tuple[str, ...]]
PRODID: Final[str]
WEEKDAYS: Final[tuple[str, ...]]
FREQUENCIES: Final[tuple[str, ...]]
ZERO_DELTA: Final[timedelta]
twoHours: Final[timedelta]

@overload
def toUnicode(s: None) -> None: ...
@overload
def toUnicode(s: str | bytes) -> str: ...
def registerTzid(tzid: str | bytes, tzinfo) -> None: ...
def getTzid(tzid: str, smart: bool = True): ...

utc: tzinfo  # dateutil.tz.tz.tzutc (subclass of tzinfo)

class TimezoneComponent(Component):
    isNative: bool
    behavior: Incomplete
    tzinfo: Incomplete
    name: str
    useBegin: bool
    def __init__(self, tzinfo=None, *args, **kwds) -> None: ...
    @classmethod
    def registerTzinfo(cls, tzinfo) -> str | None: ...
    def gettzinfo(self): ...
    tzid: Incomplete
    daylight: Incomplete
    standard: Incomplete
    def settzinfo(self, tzinfo, start: int = 2000, end: int = 2030) -> None: ...
    normal_attributes: Incomplete
    @staticmethod
    def pickTzid(tzinfo, allowUTC: bool = False) -> str | None: ...
    def prettyPrint(self, level: int, tabwidth: int) -> None: ...  # type: ignore[override]

class RecurringComponent(Component):
    isNative: bool
    def __init__(self, *args, **kwds) -> None: ...
    def getrruleset(self, addRDate: bool = False): ...
    def setrruleset(self, rruleset): ...
    rruleset: Incomplete
    def __setattr__(self, name, value) -> None: ...

class TextBehavior(Behavior):
    base64string: str
    @classmethod
    def decode(cls, line) -> None: ...
    @classmethod
    def encode(cls, line) -> None: ...

class VCalendarComponentBehavior(Behavior):
    defaultBehavior: Incomplete
    isComponent: bool

class RecurringBehavior(VCalendarComponentBehavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...
    @staticmethod
    def generateImplicitParameters(obj) -> None: ...

class DateTimeBehavior(Behavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @classmethod
    def transformFromNative(cls, obj): ...

class UTCDateTimeBehavior(DateTimeBehavior):
    forceUTC: bool

class DateOrDateTimeBehavior(Behavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...

class MultiDateBehavior(Behavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...

class MultiTextBehavior(Behavior):
    listSeparator: str
    @classmethod
    def decode(cls, line) -> None: ...
    @classmethod
    def encode(cls, line) -> None: ...

class SemicolonMultiTextBehavior(MultiTextBehavior):
    listSeparator: str

class VCalendar2_0(VCalendarComponentBehavior):
    name: str
    description: str
    versionString: str
    sortFirst: Incomplete
    @classmethod
    def generateImplicitParameters(cls, obj) -> None: ...
    @classmethod
    def serialize(cls, obj, buf, lineLength, validate: bool = True): ...

class VTimezone(VCalendarComponentBehavior):
    name: str
    hasNative: bool
    description: str
    sortFirst: Incomplete
    @classmethod
    def validate(cls, obj, raiseException: bool, *args) -> bool: ...  # type: ignore[override]
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...

class TZID(Behavior): ...

class DaylightOrStandard(VCalendarComponentBehavior):
    hasNative: bool

class VEvent(RecurringBehavior):
    name: str
    sortFirst: Incomplete
    description: str
    @classmethod
    def validate(cls, obj, raiseException: bool, *args) -> bool: ...  # type: ignore[override]

class VTodo(RecurringBehavior):
    name: str
    description: str
    @classmethod
    def validate(cls, obj, raiseException: bool, *args) -> bool: ...  # type: ignore[override]

class VJournal(RecurringBehavior):
    name: str

class VFreeBusy(VCalendarComponentBehavior):
    name: str
    description: str
    sortFirst: Incomplete

class VAlarm(VCalendarComponentBehavior):
    name: str
    description: str
    @staticmethod
    def generateImplicitParameters(obj) -> None: ...
    @classmethod
    def validate(cls, obj, raiseException: bool, *args) -> bool: ...  # type: ignore[override]

class VAvailability(VCalendarComponentBehavior):
    name: str
    description: str
    sortFirst: Incomplete
    @classmethod
    def validate(cls, obj, raiseException: bool, *args) -> bool: ...  # type: ignore[override]

class Available(RecurringBehavior):
    name: str
    sortFirst: Incomplete
    description: str
    @classmethod
    def validate(cls, obj, raiseException: bool, *args) -> bool: ...  # type: ignore[override]

class Duration(Behavior):
    name: str
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...

class Trigger(Behavior):
    name: str
    description: str
    hasNative: bool
    forceUTC: bool
    @staticmethod
    def transformToNative(obj): ...
    @staticmethod
    def transformFromNative(obj): ...

class PeriodBehavior(Behavior):
    hasNative: bool
    @staticmethod
    def transformToNative(obj): ...
    @classmethod
    def transformFromNative(cls, obj): ...

class FreeBusy(PeriodBehavior):
    name: str
    forceUTC: bool

class RRule(Behavior): ...

utcDateTimeList: list[str]
dateTimeOrDateList: list[str]
textList: list[str]

def numToDigits(num: float | None, places: int) -> str: ...
def timedeltaToString(delta: datetime.timedelta) -> str: ...
def timeToString(dateOrDateTime: datetime.date | datetime.datetime) -> str: ...
def dateToString(date: datetime.date) -> str: ...
def dateTimeToString(dateTime: datetime.datetime, convertToUTC: bool = False) -> str: ...
def deltaToOffset(delta: datetime.timedelta) -> str: ...
def periodToString(period, convertToUTC: bool = False): ...
def isDuration(s: str) -> bool: ...
def stringToDate(s: str) -> datetime.date: ...
def stringToDateTime(s: str, tzinfo: datetime._TzInfo | None = None, strict: bool = False) -> datetime.datetime: ...

escapableCharList: str

def stringToTextValues(
    s: str, listSeparator: str = ",", charList: Iterable[str] | None = None, strict: bool = False
) -> list[str]: ...
def stringToDurations(s: str, strict: bool = False) -> list[datetime.timedelta]: ...
def parseDtstart(contentline, allowSignatureMismatch: bool = False): ...
def stringToPeriod(s: str, tzinfo=None): ...
def getTransition(transitionTo: Literal["daylight", "standard"], year: int, tzinfo: datetime._TzInfo): ...
def tzinfo_eq(
    tzinfo1: datetime._TzInfo | None, tzinfo2: datetime._TzInfo | None, startYear: int = 2000, endYear: int = 2020
) -> bool: ...
