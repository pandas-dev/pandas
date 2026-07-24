from _typeshed import Incomplete, Unused
from collections.abc import Sequence
from contextlib import AbstractContextManager
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

import gdb

from .sources import Source

@type_check_only
class _SourceBreakpoint(TypedDict):
    source: str
    line: int
    condition: NotRequired[str | None]
    hitCondition: NotRequired[str | None]
    logMessage: NotRequired[str | None]

@type_check_only
class _ExceptionFilterOptions(TypedDict):
    filderId: str
    condition: NotRequired[str | None]

@type_check_only
class _BreakpointDescriptor(TypedDict):
    id: int
    verified: bool
    reason: NotRequired[str]  # only present when verified is False. Possibly only literal "pending" or "failed"
    message: NotRequired[str]  # only present when reason == "failed"
    source: NotRequired[Source]
    line: NotRequired[int]
    instructionReference: NotRequired[str]

@type_check_only
class _SetBreakpointResult(TypedDict):
    breakpoints: list[_BreakpointDescriptor]

# frozenset entries are tuples from _SourceBreakpoint.items() or _ExceptionFilterOptions.items()
breakpoint_map: dict[str, dict[frozenset[Incomplete], gdb.Breakpoint]]

def suppress_new_breakpoint_event() -> AbstractContextManager[None]: ...
def set_breakpoint(*, source: Source, breakpoints: Sequence[_SourceBreakpoint] = (), **args: Unused) -> _SetBreakpointResult: ...
def set_fn_breakpoint(*, breakpoints: Sequence[_SourceBreakpoint], **args: Unused) -> _SetBreakpointResult: ...
def set_insn_breakpoints(
    *, breakpoints: Sequence[_SourceBreakpoint], offset: int | None = None, **args: Unused
) -> _SetBreakpointResult: ...
def set_exception_breakpoints(
    *, filters: Sequence[str], filterOptions: Sequence[_ExceptionFilterOptions] = (), **args: Unused
) -> _SetBreakpointResult: ...
