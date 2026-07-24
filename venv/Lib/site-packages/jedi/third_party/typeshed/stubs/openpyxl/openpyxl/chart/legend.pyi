from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.layout import Layout
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedInteger, NestedSet
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

_LegendLegendPos: TypeAlias = Literal["b", "tr", "l", "r", "t"]

class LegendEntry(Serialisable):
    tagname: ClassVar[str]
    idx: NestedInteger[Literal[False]]
    delete: NestedBool[Literal[False]]
    txPr: Typed[RichText, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        idx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt = 0,
        delete: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = False,
        txPr: RichText | None = None,
        extLst: Unused = None,
    ) -> None: ...

class Legend(Serialisable):
    tagname: ClassVar[str]
    legendPos: NestedSet[_LegendLegendPos]
    position: Alias
    legendEntry: Incomplete
    layout: Typed[Layout, Literal[True]]
    overlay: NestedBool[Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    txPr: Typed[RichText, Literal[True]]
    textProperties: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        legendPos: _HasTagAndGet[_LegendLegendPos] | _LegendLegendPos = "r",
        legendEntry=(),
        layout: Layout | None = None,
        overlay: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        extLst: Unused = None,
    ) -> None: ...
