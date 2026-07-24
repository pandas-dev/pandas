from typing import TypeVar
from typing_extensions import Self

from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.xml.functions import Element

_T = TypeVar("_T", bound=Serialisable)

# Abstract base class.
class ElementList(list[_T]):
    @property
    def tagname(self) -> str: ...  # abstract
    @property
    def expected_type(self) -> type[_T]: ...  # abstract
    @classmethod
    def from_tree(cls, tree: Element) -> Self: ...
    def to_tree(self) -> Element: ...
    def append(self, value: _T) -> None: ...
