from _typeshed import ConvertibleToFloat, ConvertibleToInt, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.data_source import NumFmt
from openpyxl.chart.layout import Layout
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText, Text
from openpyxl.descriptors.base import Alias, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedFloat, NestedInteger, NestedSet
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.xml._functions_overloads import _HasTagAndGet

_TrendlineTrendlineType: TypeAlias = Literal["exp", "linear", "log", "movingAvg", "poly", "power"]

class TrendlineLabel(Serialisable):
    tagname: ClassVar[str]
    layout: Typed[Layout, Literal[True]]
    tx: Typed[Text, Literal[True]]
    numFmt: Typed[NumFmt, Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    txPr: Typed[RichText, Literal[True]]
    textProperties: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        layout: Layout | None = None,
        tx: Text | None = None,
        numFmt: NumFmt | None = None,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        extLst: Unused = None,
    ) -> None: ...

class Trendline(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[True]]
    spPr: Typed[ExtensionList, Literal[True]]
    graphicalProperties: Alias
    trendlineType: NestedSet[_TrendlineTrendlineType]
    order: NestedInteger[Literal[True]]
    period: NestedInteger[Literal[True]]
    forward: NestedFloat[Literal[True]]
    backward: NestedFloat[Literal[True]]
    intercept: NestedFloat[Literal[True]]
    dispRSqr: NestedBool[Literal[True]]
    dispEq: NestedBool[Literal[True]]
    trendlineLbl: Typed[ExtensionList, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        name: str | None = None,
        spPr: ExtensionList | None = None,
        trendlineType: _HasTagAndGet[_TrendlineTrendlineType] | _TrendlineTrendlineType = "linear",
        order: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        period: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        forward: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        backward: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        intercept: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        dispRSqr: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        dispEq: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        trendlineLbl: ExtensionList | None = None,
        extLst: Unused = None,
    ) -> None: ...
