from _typeshed import ConvertibleToInt, Incomplete, Unused
from re import Pattern
from typing import ClassVar, Final, Literal, SupportsIndex, overload
from typing_extensions import TypeGuard

from openpyxl.descriptors import Strict, String
from openpyxl.descriptors.base import Integer
from openpyxl.descriptors.serialisable import Serialisable

BUILTIN_FORMATS: Final[dict[int, str]]
BUILTIN_FORMATS_MAX_SIZE: Final = 164
BUILTIN_FORMATS_REVERSE: Final[dict[str, int]]
FORMAT_GENERAL: Final = "General"
FORMAT_TEXT: Final = "@"
FORMAT_NUMBER: Final = "0"
FORMAT_NUMBER_00: Final = "0.00"
FORMAT_NUMBER_COMMA_SEPARATED1: Final = "#,##0.00"
FORMAT_NUMBER_COMMA_SEPARATED2: Final = "#,##0.00_-"
FORMAT_PERCENTAGE: Final = "0%"
FORMAT_PERCENTAGE_00: Final = "0.00%"
FORMAT_DATE_YYYYMMDD2: Final = "yyyy-mm-dd"
FORMAT_DATE_YYMMDD: Final = "yy-mm-dd"
FORMAT_DATE_DDMMYY: Final = "dd/mm/yy"
FORMAT_DATE_DMYSLASH: Final = "d/m/y"
FORMAT_DATE_DMYMINUS: Final = "d-m-y"
FORMAT_DATE_DMMINUS: Final = "d-m"
FORMAT_DATE_MYMINUS: Final = "m-y"
FORMAT_DATE_XLSX14: Final = "mm-dd-yy"
FORMAT_DATE_XLSX15: Final = "d-mmm-yy"
FORMAT_DATE_XLSX16: Final = "d-mmm"
FORMAT_DATE_XLSX17: Final = "mmm-yy"
FORMAT_DATE_XLSX22: Final = "m/d/yy h:mm"
FORMAT_DATE_DATETIME: Final = "yyyy-mm-dd h:mm:ss"
FORMAT_DATE_TIME1: Final = "h:mm AM/PM"
FORMAT_DATE_TIME2: Final = "h:mm:ss AM/PM"
FORMAT_DATE_TIME3: Final = "h:mm"
FORMAT_DATE_TIME4: Final = "h:mm:ss"
FORMAT_DATE_TIME5: Final = "mm:ss"
FORMAT_DATE_TIME6: Final = "h:mm:ss"
FORMAT_DATE_TIME7: Final = "i:s.S"
FORMAT_DATE_TIME8: Final = "h:mm:ss@"
FORMAT_DATE_TIMEDELTA: Final = "[hh]:mm:ss"
FORMAT_DATE_YYMMDDSLASH: Final = "yy/mm/dd@"
FORMAT_CURRENCY_USD_SIMPLE: Final = '"$"#,##0.00_-'
FORMAT_CURRENCY_USD: Final = "$#,##0_-"
FORMAT_CURRENCY_EUR_SIMPLE: Final = "[$EUR ]#,##0.00_-"

COLORS: Final[str]
LITERAL_GROUP: Final = r'".*?"'
LOCALE_GROUP: Final = r"\[(?!hh?\]|mm?\]|ss?\])[^\]]*\]"
STRIP_RE: Final[Pattern[str]]
TIMEDELTA_RE: Final[Pattern[str]]

def is_date_format(fmt: str | None) -> TypeGuard[str]: ...
def is_timedelta_format(fmt: str | None) -> TypeGuard[str]: ...
@overload
def is_datetime(fmt: None) -> None: ...
@overload
def is_datetime(fmt: str) -> Literal["datetime", "date", "time"] | None: ...
def is_builtin(fmt: str | None) -> TypeGuard[str]: ...
def builtin_format_code(index: int) -> str | None: ...
@overload
def builtin_format_id(fmt: None) -> None: ...
@overload
def builtin_format_id(fmt: str) -> int | None: ...

class NumberFormatDescriptor(String[Incomplete]):
    def __set__(self, instance: Serialisable | Strict, value) -> None: ...

class NumberFormat(Serialisable):
    numFmtId: Integer[Literal[False]]
    formatCode: String[Literal[False]]
    def __init__(self, numFmtId: ConvertibleToInt, formatCode: str) -> None: ...

class NumberFormatList(Serialisable):
    # Overwritten by property below
    # count: Integer
    numFmt: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, count: Unused = None, numFmt=()) -> None: ...
    @property
    def count(self) -> int: ...
    def __getitem__(self, idx: SupportsIndex) -> NumberFormat: ...
