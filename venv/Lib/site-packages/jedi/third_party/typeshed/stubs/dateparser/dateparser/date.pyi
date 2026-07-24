import re
from _typeshed import Incomplete
from collections import OrderedDict
from collections.abc import Callable, Iterable, Iterator, Set as AbstractSet
from datetime import date, datetime, tzinfo
from typing import Any, ClassVar, Final, Literal, NamedTuple, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from dateparser.conf import Settings
from dateparser.languages.loader import LocaleDataLoader
from dateparser.languages.locale import Locale

_DateT = TypeVar("_DateT", bound=date)

_DetectLanguagesFunction: TypeAlias = Callable[[str, float], list[str]]
_Period: TypeAlias = Literal["time", "day", "week", "month", "year"]
# Work around attribute and type having the same name.
_Weekday: TypeAlias = Incomplete  # Actually it's dateutil._common.weekday class

@type_check_only
class _DateData(NamedTuple):
    date_obj: datetime | None
    locale: str | None
    period: _Period | None

APOSTROPHE_LOOK_ALIKE_CHARS: Final[list[str]]
RE_NBSP: Final[re.Pattern[str]]
RE_SPACES: Final[re.Pattern[str]]
RE_TRIM_SPACES: Final[re.Pattern[str]]
RE_TRIM_COLONS: Final[re.Pattern[str]]
RE_SANITIZE_SKIP: Final[re.Pattern[str]]
RE_SANITIZE_RUSSIAN: Final[re.Pattern[str]]
RE_SANITIZE_PERIOD: Final[re.Pattern[str]]
RE_SANITIZE_ON: Final[re.Pattern[str]]
RE_SANITIZE_APOSTROPHE: Final[re.Pattern[str]]
RE_SEARCH_TIMESTAMP: Final[re.Pattern[str]]
RE_SANITIZE_CROATIAN: Final[re.Pattern[str]]
RE_SEARCH_NEGATIVE_TIMESTAMP: Final[re.Pattern[str]]

def sanitize_spaces(date_string: str) -> str: ...
def date_range(
    begin: _DateT,
    end: _DateT,
    *,
    dt1: date | None = None,
    dt2: date | None = None,
    years: int = 0,
    months: int = 0,
    days: int = 0,
    leapdays: int = 0,
    weeks: int = 0,
    hours: int = 0,
    minutes: int = 0,
    seconds: int = 0,
    microseconds: int = 0,
    weekday: int | _Weekday | None = None,
    yearday: int | None = None,
    nlyearday: int | None = None,
    microsecond: int | None = None,
) -> Iterator[_DateT]: ...
def get_intersecting_periods(
    low: _DateT, high: _DateT, period: Literal["year", "month", "week", "day", "hour", "minute", "second", "microsecond"] = "day"
) -> Iterator[_DateT]: ...
def sanitize_date(date_string: str) -> str: ...
def get_date_from_timestamp(date_string: str, settings: Settings, negative: bool | None = False) -> datetime | None: ...
def parse_with_formats(date_string: str, date_formats: Iterable[str], settings: Settings) -> DateData: ...

class _DateLocaleParser:
    locale: Locale
    date_string: str
    date_formats: list[str] | tuple[str, ...] | AbstractSet[str] | None
    def __init__(
        self,
        locale: Locale,
        date_string: str,
        date_formats: list[str] | tuple[str, ...] | AbstractSet[str] | None,
        settings: Settings | None = None,
    ) -> None: ...
    @classmethod
    def parse(
        cls,
        locale: Locale,
        date_string: str,
        date_formats: list[str] | tuple[str, ...] | AbstractSet[str] | None = None,
        settings: Settings | None = None,
    ) -> DateData: ...
    def _parse(self) -> DateData | None: ...
    def _try_timestamp(self) -> DateData: ...
    def _try_freshness_parser(self) -> DateData | None: ...
    def _try_absolute_parser(self) -> DateData | None: ...
    def _try_nospaces_parser(self) -> DateData | None: ...
    def _try_parser(self, parse_method: Callable[[str, Settings, tzinfo | None], tuple[datetime, str]]) -> DateData | None: ...
    def _try_given_formats(self) -> DateData | None: ...
    def _get_translated_date(self) -> str: ...
    def _get_translated_date_with_formatting(self) -> str: ...
    def _is_valid_date_data(self, date_data: DateData) -> bool: ...

class DateData:
    date_obj: datetime | None
    locale: str | None
    period: _Period | None
    def __init__(self, *, date_obj: datetime | None = None, period: _Period | None = None, locale: str | None = None) -> None: ...
    @overload
    def __getitem__(self, k: Literal["date_obj"]) -> datetime | None: ...
    @overload
    def __getitem__(self, k: Literal["locale"]) -> str | None: ...
    @overload
    def __getitem__(self, k: Literal["period"]) -> _Period | None: ...
    @overload
    def __setitem__(self, k: Literal["date_obj"], v: datetime) -> None: ...
    @overload
    def __setitem__(self, k: Literal["locale"], v: str) -> None: ...
    @overload
    def __setitem__(self, k: Literal["period"], v: _Period) -> None: ...

class DateDataParser:
    _settings: Settings
    locale_loader: ClassVar[LocaleDataLoader | None]
    try_previous_locales: bool
    use_given_order: bool
    languages: list[str] | None
    locales: list[str] | tuple[str, ...] | AbstractSet[str] | None
    region: str
    detect_languages_function: _DetectLanguagesFunction | None
    previous_locales: OrderedDict[Locale, None]
    def __init__(
        self,
        languages: list[str] | tuple[str, ...] | AbstractSet[str] | None = None,
        locales: list[str] | tuple[str, ...] | AbstractSet[str] | None = None,
        region: str | None = None,
        try_previous_locales: bool = False,
        use_given_order: bool = False,
        settings: Settings | dict[str, Any] | None = None,
        detect_languages_function: _DetectLanguagesFunction | None = None,
    ) -> None: ...
    def get_date_data(
        self, date_string: str, date_formats: list[str] | tuple[str, ...] | AbstractSet[str] | None = None
    ) -> DateData: ...
    def get_date_tuple(
        self, date_string: str, date_formats: list[str] | tuple[str, ...] | AbstractSet[str] | None = None
    ) -> _DateData: ...
    def _get_applicable_locales(self, date_string: str) -> Iterator[Locale]: ...
    def _is_applicable_locale(self, locale: Locale, date_string: str) -> bool: ...
    @classmethod
    def _get_locale_loader(cls: type[DateDataParser]) -> LocaleDataLoader: ...
