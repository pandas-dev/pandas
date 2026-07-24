from _typeshed import ConvertibleToInt, Unused
from typing import ClassVar, Literal, overload

from openpyxl.chart.label import DataLabel as _DataLabel
from openpyxl.chart.marker import Marker
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText
from openpyxl.descriptors.base import Alias, Typed
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedInteger, NestedText
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

class PivotSource(Serialisable):
    tagname: ClassVar[str]
    name: NestedText[str, Literal[False]]
    fmtId: NestedInteger[Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self, name: object, fmtId: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt, extLst: Unused = None
    ) -> None: ...
    @overload
    def __init__(
        self, name: object = None, *, fmtId: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt, extLst: Unused = None
    ) -> None: ...

class PivotFormat(Serialisable):
    tagname: ClassVar[str]
    idx: NestedInteger[Literal[False]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    txPr: Typed[RichText, Literal[True]]
    TextBody: Alias
    marker: Typed[Marker, Literal[True]]
    dLbl: Typed[_DataLabel, Literal[True]]
    DataLabel: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        idx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt = 0,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        marker: Marker | None = None,
        dLbl: _DataLabel | None = None,
        extLst: Unused = None,
    ) -> None: ...
