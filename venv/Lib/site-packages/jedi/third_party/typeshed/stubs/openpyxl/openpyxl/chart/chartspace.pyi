from _typeshed import ConvertibleToFloat, Incomplete, Unused
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.chart._3d import Surface, View3D
from openpyxl.chart.legend import Legend
from openpyxl.chart.pivot import PivotSource
from openpyxl.chart.plotarea import PlotArea
from openpyxl.chart.print_settings import PrintSettings
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText
from openpyxl.chart.title import Title
from openpyxl.descriptors.base import Alias, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedMinMax, NestedNoneSet, NestedString, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.colors import ColorMapping
from openpyxl.xml.functions import Element

from ..xml._functions_overloads import _HasTagAndGet

_ChartContainerDispBlanksAs: TypeAlias = Literal["span", "gap", "zero"]

class ChartContainer(Serialisable):
    tagname: ClassVar[str]
    title: Typed[Title, Literal[True]]
    autoTitleDeleted: NestedBool[Literal[True]]
    pivotFmts: Incomplete

    # Same as _3DBase
    # https://github.com/python/mypy/issues/6700
    view3D: Typed[View3D, Literal[True]]
    floor: Typed[Surface, Literal[True]]
    sideWall: Typed[Surface, Literal[True]]
    backWall: Typed[Surface, Literal[True]]

    plotArea: Typed[PlotArea, Literal[False]]
    legend: Typed[Legend, Literal[True]]
    plotVisOnly: NestedBool[Literal[False]]
    dispBlanksAs: NestedNoneSet[_ChartContainerDispBlanksAs]
    showDLblsOverMax: NestedBool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        title: Title | None = None,
        autoTitleDeleted: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        pivotFmts=(),
        view3D: View3D | None = None,
        floor: Surface | None = None,
        sideWall: Surface | None = None,
        backWall: Surface | None = None,
        plotArea: PlotArea | None = None,
        legend: Legend | None = None,
        plotVisOnly: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = True,
        dispBlanksAs: _NestedNoneSetParam[_ChartContainerDispBlanksAs] = "gap",
        showDLblsOverMax: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        extLst: Unused = None,
    ) -> None: ...

class Protection(Serialisable):
    tagname: ClassVar[str]
    chartObject: NestedBool[Literal[True]]
    data: NestedBool[Literal[True]]
    formatting: NestedBool[Literal[True]]
    selection: NestedBool[Literal[True]]
    userInterface: NestedBool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        chartObject: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        data: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        formatting: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        selection: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        userInterface: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
    ) -> None: ...

class ExternalData(Serialisable):
    tagname: ClassVar[str]
    autoUpdate: NestedBool[Literal[True]]
    id: String[Literal[False]]
    @overload
    def __init__(
        self, autoUpdate: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None, *, id: str
    ) -> None: ...
    @overload
    def __init__(self, autoUpdate: Incomplete | None, id: str) -> None: ...

class ChartSpace(Serialisable):
    tagname: ClassVar[str]
    date1904: NestedBool[Literal[True]]
    lang: NestedString[Literal[True]]
    roundedCorners: NestedBool[Literal[True]]
    style: NestedMinMax[float, Literal[True]]
    clrMapOvr: Typed[ColorMapping, Literal[True]]
    pivotSource: Typed[PivotSource, Literal[True]]
    protection: Typed[Protection, Literal[True]]
    chart: Typed[ChartContainer, Literal[False]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphical_properties: Alias
    txPr: Typed[RichText, Literal[True]]
    textProperties: Alias
    externalData: Typed[ExternalData, Literal[True]]
    printSettings: Typed[PrintSettings, Literal[True]]
    userShapes: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        date1904: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        lang: object = None,
        roundedCorners: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        style: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        clrMapOvr: ColorMapping | None = None,
        pivotSource: PivotSource | None = None,
        protection: Protection | None = None,
        *,
        chart: ChartContainer,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        externalData: ExternalData | None = None,
        printSettings: PrintSettings | None = None,
        userShapes=None,
        extLst: Unused = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        date1904: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None,
        lang: object,
        roundedCorners: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None,
        style: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None,
        clrMapOvr: ColorMapping | None,
        pivotSource: PivotSource | None,
        protection: Protection | None,
        chart: ChartContainer,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        externalData: ExternalData | None = None,
        printSettings: PrintSettings | None = None,
        userShapes=None,
        extLst: Unused = None,
    ) -> None: ...
    def to_tree(self, tagname: Unused = None, idx: Unused = None, namespace: Unused = None) -> Element: ...
