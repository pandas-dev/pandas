from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal

from openpyxl.descriptors.base import Typed
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.xml.functions import Element

class Record(Serialisable):
    tagname: ClassVar[str]
    m: Incomplete
    n: Incomplete
    b: Incomplete
    e: Incomplete
    s: Incomplete
    d: Incomplete
    x: Incomplete
    def __init__(self, _fields=(), m=None, n=None, b=None, e=None, s=None, d=None, x=None) -> None: ...

class RecordList(Serialisable):
    mime_type: str
    rel_type: str
    tagname: ClassVar[str]
    r: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, count: Unused = None, r=(), extLst: ExtensionList | None = None) -> None: ...
    @property
    def count(self) -> int: ...
    def to_tree(self) -> Element: ...  # type: ignore[override]
    @property
    def path(self) -> str: ...
