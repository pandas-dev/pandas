from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal, SupportsIndex

from openpyxl.descriptors.base import Alias, Typed
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.styles import Alignment, Border, Fill, Font, Protection
from openpyxl.styles.numbers import NumberFormat

class DifferentialStyle(Serialisable):
    tagname: ClassVar[str]
    __elements__: ClassVar[tuple[str, ...]]
    font: Typed[Font, Literal[True]]
    numFmt: Typed[NumberFormat, Literal[True]]
    fill: Typed[Fill, Literal[True]]
    alignment: Typed[Alignment, Literal[True]]
    border: Typed[Border, Literal[True]]
    protection: Typed[Protection, Literal[True]]
    extLst: ExtensionList | None
    def __init__(
        self,
        font: Font | None = None,
        numFmt: NumberFormat | None = None,
        fill: Fill | None = None,
        alignment: Alignment | None = None,
        border: Border | None = None,
        protection: Protection | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class DifferentialStyleList(Serialisable):
    tagname: ClassVar[str]
    dxf: Incomplete
    styles: Alias
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, dxf=(), count: Unused = None) -> None: ...
    def append(self, dxf: DifferentialStyle) -> None: ...
    def add(self, dxf: DifferentialStyle) -> int: ...
    def __bool__(self) -> bool: ...
    def __getitem__(self, idx: SupportsIndex) -> DifferentialStyle: ...
    @property
    def count(self) -> int: ...
