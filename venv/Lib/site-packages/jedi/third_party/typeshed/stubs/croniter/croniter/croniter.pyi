import datetime
from _typeshed import Unused
from collections.abc import Generator, Iterable
from re import Match, Pattern
from typing import Any, Final, Generic, Literal, Protocol, TypeVar, overload, type_check_only
from typing_extensions import Never, Self, TypeAlias

_R_co = TypeVar("_R_co", float, datetime.datetime, default=float, covariant=True)
_R2_co = TypeVar("_R2_co", float, datetime.datetime, covariant=True)
_Expressions: TypeAlias = list[str]  # fixed-length list of 5 or 6 strings
ExpandedExpression: TypeAlias = list[int | Literal["*", "l"]]

@type_check_only
class _AllIter(Protocol[_R_co]):
    @overload
    def __call__(
        self, ret_type: type[_R2_co], start_time: float | datetime.datetime | None = None, update_current: bool | None = None
    ) -> Generator[_R2_co]: ...
    @overload
    def __call__(
        self, ret_type: None = None, start_time: float | datetime.datetime | None = None, update_current: bool | None = None
    ) -> Generator[_R_co]: ...

def is_32bit() -> bool: ...

OVERFLOW32B_MODE: Final[bool]

UTC_DT: Final[datetime.timezone]
EPOCH: Final[datetime.datetime]
M_ALPHAS: Final[dict[str, int]]
DOW_ALPHAS: Final[dict[str, int]]

MINUTE_FIELD: Final = 0
HOUR_FIELD: Final = 1
DAY_FIELD: Final = 2
MONTH_FIELD: Final = 3
DOW_FIELD: Final = 4
SECOND_FIELD: Final = 5
YEAR_FIELD: Final = 6

UNIX_FIELDS: Final[tuple[int, int, int, int, int]]
SECOND_FIELDS: Final[tuple[int, int, int, int, int, int]]
YEAR_FIELDS: Final[tuple[int, int, int, int, int, int, int]]

step_search_re: Final[Pattern[str]]
only_int_re: Final[Pattern[str]]

DAYS: Final[tuple[int, int, int, int, int, int, int, int, int, int, int, int]]
WEEKDAYS: Final[str]
MONTHS: Final[str]
star_or_int_re: Final[Pattern[str]]
special_dow_re: Final[Pattern[str]]
nearest_weekday_re: Final[Pattern[str]]
re_star: Final[Pattern[str]]
hash_expression_re: Final[Pattern[str]]

CRON_FIELDS: Final[dict[str | int, tuple[int, ...]]]
UNIX_CRON_LEN: Final = 5
SECOND_CRON_LEN: Final = 6
YEAR_CRON_LEN: Final = 7
VALID_LEN_EXPRESSION: Final[set[int]]
MARKER: object

def datetime_to_timestamp(d: datetime.datetime) -> float: ...

class CroniterError(ValueError): ...
class CroniterBadTypeRangeError(TypeError): ...
class CroniterBadCronError(CroniterError): ...
class CroniterUnsupportedSyntaxError(CroniterBadCronError): ...
class CroniterBadDateError(CroniterError): ...
class CroniterNotAlphaError(CroniterError): ...

class croniter(Generic[_R_co]):
    MONTHS_IN_YEAR: Final = 12
    RANGES: Final[
        tuple[
            tuple[int, int], tuple[int, int], tuple[int, int], tuple[int, int], tuple[int, int], tuple[int, int], tuple[int, int]
        ]
    ]
    ALPHACONV: Final[
        tuple[
            dict[Never, Never],
            dict[Never, Never],
            dict[str, str],
            dict[str, int],
            dict[str, int],
            dict[Never, Never],
            dict[Never, Never],
        ]
    ]
    LOWMAP: Final[
        tuple[
            dict[Never, Never],
            dict[Never, Never],
            dict[int, int],
            dict[int, int],
            dict[int, int],
            dict[Never, Never],
            dict[Never, Never],
        ]
    ]
    LEN_MEANS_ALL: Final[tuple[int, int, int, int, int, int, int]]

    second_at_beginning: bool
    tzinfo: datetime.tzinfo | None

    start_time: float
    dst_start_time: float
    cur: float

    expanded: list[list[str]]
    nth_weekday_of_month: dict[str, set[int]]
    expressions: _Expressions
    nearest_weekday: set[int]
    fields: tuple[int, ...]

    @overload
    def __new__(
        cls,
        expr_format: str,
        start_time: float | datetime.datetime | None = None,
        ret_type: type[float] = ...,
        day_or: bool = True,
        max_years_between_matches: int | None = None,
        is_prev: bool = False,
        hash_id: str | bytes | None = None,
        implement_cron_bug: bool = False,
        second_at_beginning: bool = False,
        expand_from_start_time: bool = False,
    ) -> croniter[float]: ...
    @overload
    def __new__(
        cls,
        expr_format: str,
        start_time: float | datetime.datetime | None,
        ret_type: type[datetime.datetime],
        day_or: bool = True,
        max_years_between_matches: int | None = None,
        is_prev: bool = False,
        hash_id: str | bytes | None = None,
        implement_cron_bug: bool = False,
        second_at_beginning: bool = False,
        expand_from_start_time: bool = False,
    ) -> croniter[datetime.datetime]: ...
    @overload
    def __new__(
        cls,
        expr_format: str,
        *,
        ret_type: type[datetime.datetime],
        day_or: bool = True,
        max_years_between_matches: int | None = None,
        is_prev: bool = False,
        hash_id: str | bytes | None = None,
        implement_cron_bug: bool = False,
        second_at_beginning: bool = False,
        expand_from_start_time: bool = False,
    ) -> croniter[datetime.datetime]: ...
    def __init__(
        self,
        expr_format: str,
        start_time: float | datetime.datetime | None = None,
        ret_type: type[_R_co] = ...,
        day_or: bool = True,
        max_years_between_matches: int | None = None,
        is_prev: bool = False,
        hash_id: str | bytes | None = None,
        implement_cron_bug: bool = False,
        second_at_beginning: bool = False,
        expand_from_start_time: bool = False,
    ) -> None: ...
    @overload
    def get_next(
        self, ret_type: type[_R2_co], start_time: float | datetime.datetime | None = None, update_current: bool = True
    ) -> _R2_co: ...
    @overload
    def get_next(
        self, ret_type: None = None, start_time: float | datetime.datetime | None = None, update_current: bool = True
    ) -> _R_co: ...
    @overload
    def get_prev(
        self, ret_type: type[_R2_co], start_time: float | datetime.datetime | None = None, update_current: bool = True
    ) -> _R2_co: ...
    @overload
    def get_prev(
        self, ret_type: None = None, start_time: float | datetime.datetime | None = None, update_current: bool = True
    ) -> _R_co: ...
    @overload
    def get_current(self, ret_type: type[_R2_co]) -> _R2_co: ...
    @overload
    def get_current(self, ret_type: None = None) -> _R_co: ...
    def set_current(self, start_time: float | datetime.datetime | None, force: bool = True) -> float: ...
    @staticmethod
    def datetime_to_timestamp(d: datetime.datetime) -> float: ...
    def timestamp_to_datetime(self, timestamp: float, tzinfo: datetime.tzinfo | None = ...) -> datetime.datetime: ...
    @overload
    def all_next(
        self, ret_type: type[_R2_co], start_time: float | datetime.datetime | None = None, update_current: bool | None = None
    ) -> Generator[_R2_co]: ...
    @overload
    def all_next(
        self, ret_type: None = None, start_time: float | datetime.datetime | None = None, update_current: bool | None = None
    ) -> Generator[_R_co]: ...
    @overload
    def all_prev(
        self, ret_type: type[_R2_co], start_time: float | datetime.datetime | None = None, update_current: bool | None = None
    ) -> Generator[_R2_co]: ...
    @overload
    def all_prev(
        self, ret_type: None = None, start_time: float | datetime.datetime | None = None, update_current: bool | None = None
    ) -> Generator[_R_co]: ...
    def iter(self, *args: Unused, **kwargs: Unused) -> _AllIter[_R_co]: ...
    def __iter__(self) -> Self: ...
    @overload
    def next(
        self,
        ret_type: type[_R2_co],
        start_time: float | datetime.datetime | None = None,
        is_prev: bool | None = None,
        update_current: bool | None = None,
    ) -> _R2_co: ...
    @overload
    def next(
        self,
        ret_type: None = None,
        start_time: float | datetime.datetime | None = None,
        is_prev: bool | None = None,
        update_current: bool | None = None,
    ) -> _R_co: ...
    __next__ = next
    @classmethod
    def value_alias(
        cls,
        val: int,
        field_index: Literal[0, 1, 2, 3, 4, 5, 6],
        len_expressions: int | list[Any] | dict[Any, Any] | tuple[Any, ...] | set[Any] = 5,
    ) -> int: ...
    DAYS_IN_MONTH: Final[dict[int, int]]
    @classmethod
    def expand(
        cls,
        expr_format: str,
        hash_id: bytes | None = None,
        second_at_beginning: bool = False,
        from_timestamp: float | None = None,
        strict: bool = False,
        strict_year: int | Iterable[int] | None = None,
    ) -> tuple[list[ExpandedExpression], dict[str, set[int]]]: ...
    @classmethod
    def is_valid(
        cls,
        expression: str,
        hash_id: bytes | None = None,
        encoding: str = "UTF-8",
        second_at_beginning: bool = False,
        strict: bool = False,
        strict_year: int | Iterable[int] | None = None,
    ) -> bool: ...
    @classmethod
    def match(
        cls,
        cron_expression: str,
        testdate: float | datetime.datetime | None,
        day_or: bool = True,
        second_at_beginning: bool = False,
        precision_in_seconds: int | None = None,
    ) -> bool: ...
    @classmethod
    def match_range(
        cls,
        cron_expression: str,
        from_datetime: datetime.datetime,
        to_datetime: datetime.datetime,
        day_or: bool = True,
        second_at_beginning: bool = False,
        precision_in_seconds: int | None = None,
    ) -> bool: ...

@overload
def croniter_range(
    start: float | datetime.datetime,
    stop: float | datetime.datetime,
    expr_format: str,
    ret_type: type[_R2_co],
    day_or: bool = True,
    exclude_ends: bool = False,
    _croniter: type[croniter] | None = None,
    second_at_beginning: bool = False,
    expand_from_start_time: bool = False,
) -> Generator[_R2_co]: ...
@overload
def croniter_range(
    start: float,
    stop: float | datetime.datetime,
    expr_format: str,
    ret_type: None = None,
    day_or: bool = True,
    exclude_ends: bool = False,
    _croniter: type[croniter] | None = None,
    second_at_beginning: bool = False,
    expand_from_start_time: bool = False,
) -> Generator[float]: ...
@overload
def croniter_range(
    start: datetime.datetime,
    stop: float | datetime.datetime,
    expr_format: str,
    ret_type: None = None,
    day_or: bool = True,
    exclude_ends: bool = False,
    _croniter: type[croniter] | None = None,
    second_at_beginning: bool = False,
    expand_from_start_time: bool = False,
) -> Generator[datetime.datetime]: ...

class HashExpander:
    cron: croniter
    def __init__(self, cronit: croniter) -> None: ...
    @overload
    def do(
        self,
        idx: int,
        hash_type: Literal["r"],
        hash_id: None = None,
        range_end: int | None = None,
        range_begin: int | None = None,
    ) -> int: ...
    @overload
    def do(
        self, idx: int, hash_type: str, hash_id: bytes, range_end: int | None = None, range_begin: int | None = None
    ) -> int: ...
    @overload
    def do(
        self, idx: int, hash_type: str = "h", *, hash_id: bytes, range_end: int | None = None, range_begin: int | None = None
    ) -> int: ...
    def match(self, efl: Unused, idx: Unused, expr: str, hash_id: bytes | None = None, **kw: Unused) -> Match[str] | None: ...
    def expand(
        self,
        efl: object,
        idx: int,
        expr: str,
        hash_id: bytes | None = None,
        match: Match[str] | None | Literal[""] = "",
        **kw: object,
    ) -> str: ...

EXPANDERS: dict[str, type[HashExpander]]
