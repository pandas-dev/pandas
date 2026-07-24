from typing import Any, Final, TypeVar, overload
from typing_extensions import TypeAlias

_T = TypeVar("_T")

__all__ = ["DEFAULT_ENCODING", "SEQUENCE_TYPES", "ICAL_TYPE", "data_encode", "from_unicode", "to_unicode"]

SEQUENCE_TYPES: Final[tuple[type[Any], ...]]
DEFAULT_ENCODING: str
ICAL_TYPE: TypeAlias = str | bytes

def from_unicode(value: ICAL_TYPE, encoding: str = "utf-8") -> bytes: ...
def to_unicode(value: ICAL_TYPE, encoding: str = "utf-8-sig") -> str: ...
@overload
def data_encode(data: ICAL_TYPE, encoding: str = "utf-8") -> bytes: ...
@overload
def data_encode(data: dict[Any, Any], encoding: str = "utf-8") -> dict[Any, Any]: ...
@overload
def data_encode(data: list[Any] | tuple[Any, ...], encoding: str = "utf-8") -> list[Any]: ...
@overload
def data_encode(data: _T, encoding: str = "utf-8") -> _T: ...
