from _typeshed import Unused
from collections.abc import Iterable
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

import gdb

from ..FrameDecorator import FrameDecorator, SymValueWrapper
from .varref import BaseReference, _ReferenceDescriptor

frame_to_scope: dict[int, _ScopeReference]

@type_check_only
class _ScopeReferenceDescriptor(_ReferenceDescriptor):
    presentationHint: str
    expensive: bool
    namedVariables: int
    line: NotRequired[int]

@type_check_only
class _ScopesResult(TypedDict):
    scopes: list[_ScopeReferenceDescriptor]

def clear_scopes(event: Unused) -> None: ...
def set_finish_value(val: gdb.Value) -> None: ...
def symbol_value(sym: SymValueWrapper, frame: FrameDecorator) -> tuple[str, gdb.Value]: ...

class _ScopeReference(BaseReference):
    hint: str
    frame: FrameDecorator
    inf_frame: gdb.Frame
    function: str | None
    line: int | None
    var_list: tuple[SymValueWrapper, ...]
    def __init__(self, name: str, hint: str, frame: FrameDecorator, var_list: Iterable[SymValueWrapper]) -> None: ...
    def to_object(self) -> _ScopeReferenceDescriptor: ...
    def has_children(self) -> bool: ...
    def child_count(self) -> int: ...
    # note: parameter named changed from 'index' to 'idx'
    def fetch_one_child(self, idx: int) -> tuple[str, gdb.Value]: ...

class _FinishScopeReference(_ScopeReference): ...

class _RegisterReference(_ScopeReference):
    def __init__(self, name: str, frame: FrameDecorator) -> None: ...

def scopes(*, frameId: int, **extra: Unused) -> _ScopesResult: ...
