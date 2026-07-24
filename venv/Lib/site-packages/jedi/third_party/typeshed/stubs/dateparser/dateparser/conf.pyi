import datetime
from collections.abc import Callable
from typing import Any, Literal, TypeVar
from typing_extensions import ParamSpec, Self

_P = ParamSpec("_P")
_R = TypeVar("_R")

class Settings:
    # Next attributes are optional and may be missing.
    # Please keep in sync with _Settings TypedDict
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

    def __new__(cls, *args, **kwargs) -> Self: ...
    def __init__(self, settings: dict[str, Any] | None = None) -> None: ...
    @classmethod
    def get_key(cls, settings: dict[str, Any] | None = None) -> str: ...
    def replace(self, mod_settings: dict[str, Any] | None = None, **kwds) -> Self: ...

settings: Settings

def apply_settings(f: Callable[_P, _R]) -> Callable[_P, _R]: ...

class SettingValidationError(ValueError): ...

def check_settings(settings: Settings) -> None: ...
