import abc
from _typeshed import Unused
from collections import defaultdict
from collections.abc import Generator
from contextlib import AbstractContextManager
from typing import TypedDict, type_check_only
from typing_extensions import NotRequired

import gdb

@type_check_only
class _ValueFormat(TypedDict, total=False):
    hex: bool

@type_check_only
class _ReferenceDescriptor(TypedDict):
    # Result of BaseReference.to_object()
    variableReference: int
    name: NotRequired[str]

@type_check_only
class _VariableReferenceDescriptor(_ReferenceDescriptor):
    # Result of VariableReference.to_object()
    indexedVariables: NotRequired[int]
    namedVariables: NotRequired[int]
    memoryReference: NotRequired[str]
    type: NotRequired[str]
    # Below key name set by VariableReference.result_name
    # Could be modelled with extra_items=str if PEP 728 is accepted.
    value: NotRequired[str]

all_variables: list[BaseReference]

def clear_vars(event: Unused) -> None: ...
def apply_format(value_format: _ValueFormat | None) -> AbstractContextManager[None]: ...

class BaseReference(abc.ABC):
    ref: int
    name: str
    children: list[VariableReference | None] | None
    by_name: dict[str, VariableReference]
    name_counts: defaultdict[str, int]
    def __init__(self, name: str) -> None: ...
    def to_object(self) -> _ReferenceDescriptor: ...
    @abc.abstractmethod
    def has_children(self) -> bool: ...
    def reset_children(self): ...
    @abc.abstractmethod
    def fetch_one_child(self, index: int) -> tuple[str, gdb.Value]: ...
    @abc.abstractmethod
    def child_count(self) -> int: ...
    def fetch_children(self, start: int, count: int) -> Generator[VariableReference]: ...
    def find_child_by_name(self, name: str) -> VariableReference: ...

class VariableReference(BaseReference):
    result_name: str
    value: gdb.Value
    child_cache: list[tuple[int | str, gdb.Value]] | None
    count: int | None
    printer: gdb._PrettyPrinter
    def __init__(self, name: str, value: gdb.Value, result_name: str = "value") -> None: ...
    def assign(self, value: gdb.Value) -> None: ...
    def has_children(self) -> bool: ...
    def cache_children(self) -> list[tuple[int | str, gdb.Value]]: ...
    def child_count(self) -> int: ...
    def to_object(self) -> _VariableReferenceDescriptor: ...
    # note: parameter named changed from 'index' to 'idx'
    def fetch_one_child(self, idx: int) -> tuple[str, gdb.Value]: ...

def find_variable(ref: int) -> BaseReference: ...
