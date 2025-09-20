"""Collection of functions for building custom `json_default` functions.

In general functions come in pairs of `use_x_default` and `x_default`, where the former is used
to determine if you should call the latter.

Most `use_x_default` functions also act as a [`TypeGuard`](https://mypy.readthedocs.io/en/stable/type_narrowing.html#user-defined-type-guards).
"""

### IMPORTS
### ============================================================================
## Future
from __future__ import annotations

## Standard Library
import base64
import dataclasses
import datetime
import enum
import sys
from types import TracebackType
from typing import Any
import traceback
import uuid

if sys.version_info >= (3, 10):
    from typing import TypeGuard
else:
    from typing_extensions import TypeGuard

## Installed

## Application


### FUNCTIONS
### ============================================================================
def unknown_default(obj: Any) -> str:
    """Backup default function for any object type.

    Will attempt to use `str` or `repr`. If both functions error will return
    the string `"__could_not_encode__"`.

    Args:
        obj: object to handle
    """
    try:
        return str(obj)
    except Exception:  # pylint: disable=broad-exception-caught
        pass
    try:
        return repr(obj)
    except Exception:  # pylint: disable=broad-exception-caught
        pass
    return "__could_not_encode__"


## Types
## -----------------------------------------------------------------------------
def use_type_default(obj: Any) -> TypeGuard[type]:
    """Default check function for `type` objects (aka classes)."""
    return isinstance(obj, type)


def type_default(obj: type) -> str:
    """Default function for `type` objects.

    Args:
        obj: object to handle
    """
    return obj.__name__


## Dataclasses
## -----------------------------------------------------------------------------
def use_dataclass_default(obj: Any) -> bool:
    """Default check function for dataclass instances"""
    return dataclasses.is_dataclass(obj) and not isinstance(obj, type)


def dataclass_default(obj) -> dict[str, Any]:
    """Default function for dataclass instances

    Args:
        obj: object to handle
    """
    return dataclasses.asdict(obj)


## Dates and Times
## -----------------------------------------------------------------------------
def use_time_default(obj: Any) -> TypeGuard[datetime.time]:
    """Default check function for `datetime.time` instances"""
    return isinstance(obj, datetime.time)


def time_default(obj: datetime.time) -> str:
    """Default function for `datetime.time` instances

    Args:
        obj: object to handle
    """
    return obj.isoformat()


def use_date_default(obj: Any) -> TypeGuard[datetime.date]:
    """Default check function for `datetime.date` instances"""
    return isinstance(obj, datetime.date)


def date_default(obj: datetime.date) -> str:
    """Default function for `datetime.date` instances

    Args:
        obj: object to handle
    """
    return obj.isoformat()


def use_datetime_default(obj: Any) -> TypeGuard[datetime.datetime]:
    """Default check function for `datetime.datetime` instances"""
    return isinstance(obj, datetime.datetime)


def datetime_default(obj: datetime.datetime) -> str:
    """Default function for `datetime.datetime` instances

    Args:
        obj: object to handle
    """
    return obj.isoformat()


def use_datetime_any(obj: Any) -> TypeGuard[datetime.time | datetime.date | datetime.datetime]:
    """Default check function for `datetime` related instances"""
    return isinstance(obj, (datetime.time, datetime.date, datetime.datetime))


def datetime_any(obj: datetime.time | datetime.date | datetime.date) -> str:
    """Default function for `datetime` related instances

    Args:
        obj: object to handle
    """
    return obj.isoformat()


## Exception and Tracebacks
## -----------------------------------------------------------------------------
def use_exception_default(obj: Any) -> TypeGuard[BaseException]:
    """Default check function for exception instances.

    Exception classes are not treated specially and should be handled by the
    `[use_]type_default` functions.
    """
    return isinstance(obj, BaseException)


def exception_default(obj: BaseException) -> str:
    """Default function for exception instances

    Args:
        obj: object to handle
    """
    return f"{obj.__class__.__name__}: {obj}"


def use_traceback_default(obj: Any) -> TypeGuard[TracebackType]:
    """Default check function for tracebacks"""
    return isinstance(obj, TracebackType)


def traceback_default(obj: TracebackType) -> str:
    """Default function for tracebacks

    Args:
        obj: object to handle
    """
    return "".join(traceback.format_tb(obj)).strip()


## Enums
## -----------------------------------------------------------------------------
def use_enum_default(obj: Any) -> TypeGuard[enum.Enum | enum.EnumMeta]:
    """Default check function for enums.

    Supports both enum classes and enum values.
    """
    return isinstance(obj, (enum.Enum, enum.EnumMeta))


def enum_default(obj: enum.Enum | enum.EnumMeta) -> Any | list[Any]:
    """Default function for enums.

    Supports both enum classes and enum values.

    Args:
        obj: object to handle
    """
    if isinstance(obj, enum.Enum):
        return obj.value
    return [e.value for e in obj]  # type: ignore[var-annotated]


## UUIDs
## -----------------------------------------------------------------------------
def use_uuid_default(obj: Any) -> TypeGuard[uuid.UUID]:
    """Default check function for `uuid.UUID` instances"""
    return isinstance(obj, uuid.UUID)


def uuid_default(obj: uuid.UUID) -> str:
    """Default function for `uuid.UUID` instances

    Formats the UUID using "hyphen" format.

    Args:
        obj: object to handle
    """
    return str(obj)


## Bytes
## -----------------------------------------------------------------------------
def use_bytes_default(obj: Any) -> TypeGuard[bytes | bytearray]:
    """Default check function for bytes"""
    return isinstance(obj, (bytes, bytearray))


def bytes_default(obj: bytes | bytearray, url_safe: bool = True) -> str:
    """Default function for bytes

    Args:
        obj: object to handle
        url_safe: use URL safe base 64 character set.

    Returns:
        The byte data as a base 64 string.
    """
    if url_safe:
        return base64.urlsafe_b64encode(obj).decode("utf8")
    return base64.b64encode(obj).decode("utf8")
