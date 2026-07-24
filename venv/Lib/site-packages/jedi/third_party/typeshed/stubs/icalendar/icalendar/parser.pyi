from _collections_abc import dict_keys
from _typeshed import Incomplete
from collections.abc import Iterable
from re import Pattern
from typing import AnyStr, ClassVar, Final, overload
from typing_extensions import Self

from .caselessdict import CaselessDict
from .parser_tools import ICAL_TYPE
from .prop import _vType

__all__ = [
    "Contentline",
    "Contentlines",
    "FOLD",
    "NAME",
    "NEWLINE",
    "Parameters",
    "QUNSAFE_CHAR",
    "QUOTABLE",
    "UNSAFE_CHAR",
    "dquote",
    "escape_char",
    "escape_string",
    "foldline",
    "param_value",
    "q_join",
    "q_split",
    "rfc_6868_escape",
    "rfc_6868_unescape",
    "uFOLD",
    "unescape_char",
    "unescape_list_or_string",
    "unescape_string",
    "validate_param_value",
    "validate_token",
]

def escape_char(text: str) -> str: ...
def unescape_char(text: AnyStr) -> AnyStr: ...
def foldline(line: str, limit: int = 75, fold_sep: str = "\r\n ") -> str: ...
def param_value(value: str | list[str] | tuple[str, ...] | Incomplete, always_quote: bool = False) -> str: ...

NAME: Final[Pattern[str]]
UNSAFE_CHAR: Final[Pattern[str]]
QUNSAFE_CHAR: Final[Pattern[str]]
FOLD: Final[Pattern[bytes]]
uFOLD: Final[Pattern[str]]
NEWLINE: Final[Pattern[str]]

def validate_token(name: str) -> None: ...
def validate_param_value(value: str, quoted: bool = True) -> None: ...

QUOTABLE: Final[Pattern[str]]

def dquote(val: str, always_quote: bool = False) -> str: ...
def q_split(st: str, sep: str = ",", maxsplit: int = -1) -> list[str]: ...
def q_join(lst: Iterable[str], sep: str = ",", always_quote: bool = False) -> str: ...

class Parameters(CaselessDict[str]):
    always_quoted: ClassVar[tuple[str, ...]]
    quote_also: ClassVar[dict[str, str]]
    def params(self) -> dict_keys[str, str]: ...
    def to_ical(self, sorted: bool = True) -> bytes: ...
    @classmethod
    def from_ical(cls, st: str, strict: bool = False) -> Self: ...

def escape_string(val: str) -> str: ...
def unescape_string(val: str) -> str: ...

RFC_6868_UNESCAPE_REGEX: Final[Pattern[str]]

def rfc_6868_unescape(param_value: str) -> str: ...

RFC_6868_ESCAPE_REGEX: Final[Pattern[str]]

def rfc_6868_escape(param_value: str) -> str: ...
@overload
def unescape_list_or_string(val: list[str]) -> list[str]: ...
@overload
def unescape_list_or_string(val: str) -> str: ...

class Contentline(str):
    __slots__ = ("strict",)
    strict: bool
    def __new__(cls, value: str | bytes, strict: bool = False, encoding: str = "utf-8") -> Self: ...
    @classmethod
    def from_parts(cls, name: ICAL_TYPE, params: Parameters, values: _vType | ICAL_TYPE, sorted: bool = True) -> Self: ...
    def parts(self) -> tuple[str, Parameters, str]: ...
    @classmethod
    def from_ical(cls, ical: str | bytes, strict: bool = False) -> Self: ...
    def to_ical(self) -> bytes: ...

class Contentlines(list[Contentline]):
    def to_ical(self) -> bytes: ...
    @classmethod
    def from_ical(cls, st: str | bytes) -> Self: ...
