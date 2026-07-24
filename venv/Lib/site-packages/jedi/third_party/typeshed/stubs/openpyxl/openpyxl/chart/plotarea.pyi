from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import Self

from openpyxl.chart.layout import Layout
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool
from openpyxl.descriptors.serialisable import Serialisable, _ChildSerialisableTreeElement
from openpyxl.xml.functions import Element

from ..xml._functions_overloads import _HasTagAndGet

class DataTable(Serialisable):
    tagname: ClassVar[str]
    showHorzBorder: NestedBool[Literal[True]]
    showVertBorder: NestedBool[Literal[True]]
    showOutline: NestedBool[Literal[True]]
    showKeys: NestedBool[Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    txPr: Typed[RichText, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        showHorzBorder: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        showVertBorder: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        showOutline: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        showKeys: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        extLst: Unused = None,
    ) -> None: ...

class PlotArea(Serialisable):
    tagname: ClassVar[str]
    layout: Typed[Layout, Literal[True]]
    dTable: Typed[DataTable, Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    areaChart: Incomplete
    area3DChart: Incomplete
    lineChart: Incomplete
    line3DChart: Incomplete
    stockChart: Incomplete
    radarChart: Incomplete
    scatterChart: Incomplete
    pieChart: Incomplete
    pie3DChart: Incomplete
    doughnutChart: Incomplete
    barChart: Incomplete
    bar3DChart: Incomplete
    ofPieChart: Incomplete
    surfaceChart: Incomplete
    surface3DChart: Incomplete
    bubbleChart: Incomplete
    valAx: Incomplete
    catAx: Incomplete
    dateAx: Incomplete
    serAx: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        layout: Layout | None = None,
        dTable: DataTable | None = None,
        spPr: GraphicalProperties | None = None,
        _charts=(),
        _axes=(),
        extLst: Unused = None,
    ) -> None: ...
    def to_tree(self, tagname: str | None = None, idx: Unused = None, namespace: Unused = None) -> Element: ...
    @classmethod
    def from_tree(cls, node: _ChildSerialisableTreeElement) -> Self: ...
