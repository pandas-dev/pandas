from _typeshed import Incomplete, Unused
from collections.abc import Generator, Iterable, Sized
from typing import Any, Protocol, TypeVar, type_check_only
from typing_extensions import Self

from openpyxl.descriptors import Strict
from openpyxl.descriptors.serialisable import Serialisable, _SerialisableTreeElement
from openpyxl.xml._functions_overloads import _HasGet
from openpyxl.xml.functions import Element

from .base import Alias, Descriptor

_T = TypeVar("_T")
_ContainerT = TypeVar("_ContainerT")

@type_check_only
class _SupportsFromTree(Protocol):
    @classmethod
    def from_tree(cls, node: _SerialisableTreeElement) -> Any: ...

@type_check_only
class _SupportsToTree(Protocol):
    def to_tree(self) -> Element: ...

# `_ContainerT` is the internal container type (which defaults to `list`), or
# `IndexedList` if unique is `True`.
class Sequence(Descriptor[_ContainerT]):
    expected_type: type[Any]  # expected type of the sequence elements
    seq_types: tuple[type, ...]  # allowed settable sequence types, defaults to `list`, `tuple`
    idx_base: int
    unique: bool
    container: type  # internal container type, defaults to `list`
    # seq must be an instance of any of the declared `seq_types`.
    def __set__(self, instance: Serialisable | Strict, seq: Any) -> None: ...
    def to_tree(self, tagname: str | None, obj: Iterable[object], namespace: str | None = None) -> Generator[Element]: ...

# `_T` is the type of the elements in the sequence.
class UniqueSequence(Sequence[set[_T]]):
    seq_types: tuple[type, ...]  # defaults to `list`, `tuple`, `set`
    container: type[set[_T]]

# See `Sequence` for the meaning of `_ContainerT`.
class ValueSequence(Sequence[_ContainerT]):
    attribute: str
    def to_tree(
        self, tagname: str, obj: Iterable[object], namespace: str | None = None  # type: ignore[override]
    ) -> Generator[Element]: ...
    def from_tree(self, node: _HasGet[_T]) -> _T: ...

@type_check_only
class _NestedSequenceToTreeObj(Sized, Iterable[_SupportsToTree], Protocol): ...

# See `Sequence` for the meaning of `_ContainerT`.
class NestedSequence(Sequence[_ContainerT]):
    count: bool
    expected_type: type[_SupportsFromTree]
    def to_tree(  # type: ignore[override]
        self, tagname: str, obj: _NestedSequenceToTreeObj, namespace: str | None = None
    ) -> Element: ...
    # returned list generic type should be same as the return type of expected_type.from_tree(node)
    # Which can really be anything given the wildly different, and sometimes generic, from_tree return types
    def from_tree(self, node: Iterable[_SerialisableTreeElement]) -> list[Any]: ...

# `_T` is the type of the elements in the sequence.
class MultiSequence(Sequence[list[_T]]):
    def __set__(self, instance: Serialisable | Strict, seq: tuple[_T, ...] | list[_T]) -> None: ...
    def to_tree(
        self, tagname: Unused, obj: Iterable[_SupportsToTree], namespace: str | None = None  # type: ignore[override]
    ) -> Generator[Element]: ...

class MultiSequencePart(Alias):
    expected_type: type[Incomplete]
    store: Incomplete
    def __init__(self, expected_type, store) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value) -> None: ...
    def __get__(self, instance: Unused, cls: Unused) -> Self: ...
