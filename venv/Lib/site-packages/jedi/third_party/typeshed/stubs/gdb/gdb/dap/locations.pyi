from _typeshed import Unused
from typing import TypedDict, type_check_only

from .sources import Source

@type_check_only
class _Line(TypedDict):
    line: int

@type_check_only
class _BreakpointLocationsResult(TypedDict):
    breakpoints: list[_Line]

def breakpoint_locations(
    *, source: Source, line: int, endLine: int | None = None, **extra: Unused
) -> _BreakpointLocationsResult: ...
