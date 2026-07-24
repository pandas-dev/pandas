import datetime
from typing import Any, Final, Literal, TypedDict, type_check_only

from dateparser.conf import Settings

from .date import DateDataParser, _DetectLanguagesFunction

__version__: Final[str]

_default_parser: DateDataParser

@type_check_only
class _Settings(TypedDict, total=False):  # noqa: Y049
    DATE_ORDER: str
    PREFER_LOCALE_DATE_ORDER: bool
    TIMEZONE: str
    TO_TIMEZONE: str
    RETURN_AS_TIMEZONE_AWARE: bool
    PREFER_MONTH_OF_YEAR: Literal["current", "first", "last"]
    PREFER_DAY_OF_MONTH: Literal["current", "first", "last"]
    PREFER_DATES_FROM: Literal["current_period", "future", "past"]
    RELATIVE_BASE: datetime.datetime
    STRICT_PARSING: bool
    REQUIRE_PARTS: list[Literal["day", "month", "year"]]
    SKIP_TOKENS: list[str]
    NORMALIZE: bool
    RETURN_TIME_AS_PERIOD: bool
    RETURN_TIME_SPAN: bool
    DEFAULT_START_OF_WEEK: Literal["monday", "sunday"]
    DEFAULT_DAYS_IN_MONTH: int
    PARSERS: list[Literal["timestamp", "relative-time", "custom-formats", "absolute-time", "no-spaces-time"]]
    DEFAULT_LANGUAGES: list[str]
    LANGUAGE_DETECTION_CONFIDENCE_THRESHOLD: float
    CACHE_SIZE_LIMIT: int

def parse(
    date_string: str,
    date_formats: list[str] | tuple[str, ...] | set[str] | None = None,
    languages: list[str] | tuple[str, ...] | set[str] | None = None,
    locales: list[str] | tuple[str, ...] | set[str] | None = None,
    region: str | None = None,
    settings: Settings | dict[str, Any] | None = None,
    detect_languages_function: _DetectLanguagesFunction | None = None,
) -> datetime.datetime | None: ...
