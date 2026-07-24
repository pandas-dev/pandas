from _typeshed import Incomplete
from collections.abc import Callable, Iterable
from typing import Literal, SupportsIndex
from typing_extensions import TypeAlias

import pythoncom

com_error = pythoncom.com_error

# Type of pythoncom._univgw.WriteFromOutTuple
# The two tuples must be of equal length
_WriteFromOutTupleType: TypeAlias = Callable[
    [tuple[Incomplete, ...] | None, tuple[Incomplete, ...] | None, int], Incomplete | None
]

def RegisterInterfaces(
    typelibGUID, lcid, major, minor, interface_names: Iterable[str] | None = None
) -> list[tuple[Incomplete, Incomplete, Incomplete]]: ...

class Arg:
    name: Incomplete
    vt: Incomplete
    inOut: Incomplete
    default: Incomplete
    clsid: Incomplete
    size: Incomplete
    offset: int
    def __init__(self, arg_info, name=None) -> None: ...

class Method:
    dispid: Incomplete
    invkind: Incomplete
    name: Incomplete
    args: list[Arg]
    cbArgs: Incomplete
    def __init__(self, method_info, isEventSink: bool | Literal[0, 1] = 0) -> None: ...

class Definition:
    def __init__(self, iid, is_dispatch, method_defs) -> None: ...
    def iid(self): ...
    def vtbl_argsizes(self) -> list[Incomplete]: ...
    def vtbl_argcounts(self) -> list[int]: ...
    def dispatch(
        self, ob, index: SupportsIndex, argPtr, ReadFromInTuple=..., WriteFromOutTuple: _WriteFromOutTupleType = ...
    ): ...
