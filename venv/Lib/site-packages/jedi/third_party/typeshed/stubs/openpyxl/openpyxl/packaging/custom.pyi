from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete
from collections.abc import Iterator
from datetime import datetime
from typing import Any, Final, Generic, Literal, TypeVar
from typing_extensions import Self, TypeAlias

from openpyxl.descriptors import Sequence, Strict
from openpyxl.descriptors.base import Bool, DateTime, Float, Integer, String, _ConvertibleToBool
from openpyxl.descriptors.nested import NestedText
from openpyxl.descriptors.serialisable import _ChildSerialisableTreeElement
from openpyxl.xml.functions import Element

_T = TypeVar("_T")

# Does not reimplement anything, so runtime also has incompatible supertypes
class NestedBoolText(Bool[Incomplete], NestedText[Incomplete, Incomplete]): ...  # type: ignore[misc]

class _TypedProperty(Strict, Generic[_T]):
    name: String[Literal[False]]
    # Since this is internal, just list all possible values
    value: (
        Integer[Literal[False]]
        | Float[Literal[False]]
        | String[Literal[True]]
        | DateTime[Literal[False]]
        | Bool[Literal[False]]
        | String[Literal[False]]
    )
    def __init__(self, name: str, value: _T) -> None: ...
    def __eq__(self, other: _TypedProperty[Any]) -> bool: ...  # type: ignore[override]

class IntProperty(_TypedProperty[ConvertibleToInt]):
    value: Integer[Literal[False]]

class FloatProperty(_TypedProperty[ConvertibleToFloat]):
    value: Float[Literal[False]]

class StringProperty(_TypedProperty[str | None]):
    value: String[Literal[True]]

class DateTimeProperty(_TypedProperty[datetime]):
    value: DateTime[Literal[False]]

class BoolProperty(_TypedProperty[_ConvertibleToBool]):
    value: Bool[Literal[False]]

class LinkProperty(_TypedProperty[str]):
    value: String[Literal[False]]

_MappingPropertyType: TypeAlias = StringProperty | IntProperty | FloatProperty | DateTimeProperty | BoolProperty | LinkProperty
CLASS_MAPPING: Final[dict[type[_MappingPropertyType], str]]
XML_MAPPING: Final[dict[str, type[_MappingPropertyType]]]

class CustomPropertyList(Strict, Generic[_T]):
    props: Sequence[list[_TypedProperty[_T]]]
    def __init__(self) -> None: ...
    @classmethod
    def from_tree(cls, tree: _ChildSerialisableTreeElement) -> Self: ...
    def append(self, prop) -> None: ...
    def to_tree(self) -> Element: ...
    def __len__(self) -> int: ...
    @property
    def names(self) -> list[str]: ...
    def __getitem__(self, name): ...
    def __delitem__(self, name) -> None: ...
    def __iter__(self) -> Iterator[_TypedProperty[_T]]: ...
