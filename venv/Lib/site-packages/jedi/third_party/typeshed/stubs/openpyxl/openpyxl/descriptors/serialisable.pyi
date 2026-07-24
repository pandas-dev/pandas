from _typeshed import ConvertibleToInt, Incomplete, SupportsIter
from collections.abc import Iterator
from typing import Any, ClassVar, Final, Protocol, type_check_only
from typing_extensions import Self

from openpyxl.descriptors import MetaSerialisable
from openpyxl.xml.functions import Element

from ..xml._functions_overloads import _HasAttrib, _HasGet, _HasTagAndGet, _HasText, _SupportsFindChartLines

# For any override directly re-using Serialisable.from_tree
@type_check_only
class _ChildSerialisableTreeElement(_HasAttrib, _HasText, SupportsIter[Incomplete], Protocol): ...

@type_check_only
class _SerialisableTreeElement(_HasGet[object], _SupportsFindChartLines, _ChildSerialisableTreeElement, Protocol): ...

KEYWORDS: Final[frozenset[str]]
seq_types: Final[tuple[type[list[Any]], type[tuple[Any, ...]]]]

class Serialisable(metaclass=MetaSerialisable):
    # These dunders are always set at runtime by MetaSerialisable so they can't be None
    __attrs__: ClassVar[tuple[str, ...]]
    __nested__: ClassVar[tuple[str, ...]]
    __elements__: ClassVar[tuple[str, ...]]
    __namespaced__: ClassVar[tuple[tuple[str, str], ...]]
    idx_base: int
    # Needs overrides in many sub-classes. But a lot of subclasses are instantiated without overriding it, so can't be abstract
    # Subclasses "overrides" this property with a ClassVar, and Serialisable is too widely used,
    # so it can't be typed as NoReturn either without introducing many false-positives.
    @property
    def tagname(self) -> str: ...
    namespace: ClassVar[str | None]
    # Note: To respect the Liskov substitution principle, the protocol for node includes all child class requirements.
    # Same with the return type to avoid override issues.
    # See comment in xml/functions.pyi as to why use a protocol instead of Element
    # Child classes should be more precise than _SerialisableTreeElement !
    # Use _ChildSerialisableTreeElement instead for child classes that reuse Serialisable.from_tree directly.
    @classmethod
    def from_tree(cls, node: _SerialisableTreeElement) -> Self | None: ...
    # Note: To respect the Liskov substitution principle, idx is a type union of all child class requirements.
    # Use Unused instead for child classes that reuse Serialisable.to_tree directly.
    def to_tree(
        self,
        tagname: str | None = None,
        idx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt | None = None,
        namespace: str | None = None,
    ) -> Element: ...
    def __iter__(self) -> Iterator[tuple[str, str]]: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __add__(self, other): ...
    def __copy__(self): ...
