from abc import ABCMeta
from collections.abc import ItemsView, Mapping, ValuesView
from io import StringIO as StringIO
from re import Pattern
from typing import Any, Final, SupportsIndex, TypeVar
from typing_extensions import TypeGuard
from urllib.parse import parse_qs, quote, unquote, urlencode as urlencode, urlparse as urlparse

_KT = TypeVar("_KT")
_VT_co = TypeVar("_VT_co", covariant=True)

url_quote = quote
url_unquote = unquote
url_parse_qs = parse_qs

PY2: Final = False
PY3: Final = True
RE_NUM: Final[Pattern[str]]
ON_LINUX: Final[bool]
ON_OSX: Final[bool]
ON_WINDOWS: Final[bool]

class AbstractBase(metaclass=ABCMeta): ...

SOCKET_ERROR = OSError
SOL_TCP: Final[int]
basestring: Final[tuple[type[str]]]
str_or_bytes: Final[tuple[type[str], type[bytes]]]
xrange = range
unicode_type = str

def time_now() -> float: ...
def dictkeys(dct: Mapping[_KT, Any]) -> list[_KT]: ...
def dictvalues(dct: Mapping[Any, _VT_co]) -> list[_VT_co]: ...
def dict_iteritems(dct: Mapping[_KT, _VT_co]) -> ItemsView[_KT, _VT_co]: ...
def dict_itervalues(dct: Mapping[Any, _VT_co]) -> ValuesView[_VT_co]: ...
def byte(*args: SupportsIndex) -> bytes: ...

class long(int): ...

def canonical_str(value: object) -> str: ...
def is_integer(value: object) -> TypeGuard[int]: ...
def as_bytes(value: str | bytes) -> bytes: ...
def to_digit(value: str) -> int: ...
def get_linux_version(release_str: str) -> tuple[int, int, int]: ...

HAVE_SIGNAL: Final[bool]
EINTR_IS_EXPOSED: Final = False
LINUX_VERSION: tuple[int, int, int] | None
