from collections.abc import Callable
from typing import TypeVar

from .parser import Parameters

class IcalendarProperty:
    params: Parameters

_T = TypeVar("_T")

def string_parameter(
    name: str,
    doc: str,
    default: Callable[..., str | None] = ...,
    convert: Callable[[str], _T] | None = None,
    convert_to: Callable[[_T], str] | None = None,
) -> property: ...

ALTREP: property
CN: property
CUTYPE: property

def quoted_list_parameter(name: str, doc: str) -> property: ...

DELEGATED_FROM: property
DELEGATED_TO: property
DIR: property
FBTYPE: property
LANGUAGE: property
MEMBER: property
PARTSTAT: property
RANGE: property
RELATED: property
ROLE: property

def boolean_parameter(name: str, default: bool, doc: str) -> property: ...

RSVP: property
SENT_BY: property
TZID: property
RELTYPE: property

__all__ = [
    "string_parameter",
    "quoted_list_parameter",
    "ALTREP",
    "CN",
    "CUTYPE",
    "DELEGATED_FROM",
    "DELEGATED_TO",
    "DIR",
    "FBTYPE",
    "LANGUAGE",
    "MEMBER",
    "PARTSTAT",
    "RANGE",
    "RELATED",
    "ROLE",
    "RSVP",
    "SENT_BY",
    "TZID",
]
