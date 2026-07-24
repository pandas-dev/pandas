from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal, overload
from typing_extensions import Self, TypeAlias

from openpyxl.chart.layout import Layout
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText, Text
from openpyxl.chart.title import Title, TitleDescriptor
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import (
    NestedBool,
    NestedFloat,
    NestedInteger,
    NestedMinMax,
    NestedNoneSet,
    NestedSet,
    _NestedNoneSetParam,
)
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet, _SupportsFindAndIterAndAttribAndText

_ScalingOrientation: TypeAlias = Literal["maxMin", "minMax"]
_BaseAxisAxPos: TypeAlias = Literal["b", "l", "r", "t"]
_BaseAxisTickMark: TypeAlias = Literal["cross", "in", "out"]
_BaseAxisTickLblPos: TypeAlias = Literal["high", "low", "nextTo"]
_BaseAxisCrosses: TypeAlias = Literal["autoZero", "max", "min"]
_DisplayUnitsLabelListBuiltInUnit: TypeAlias = Literal[
    "hundreds",
    "thousands",
    "tenThousands",
    "hundredThousands",
    "millions",
    "tenMillions",
    "hundredMillions",
    "billions",
    "trillions",
]
_NumericAxisCrossBetween: TypeAlias = Literal["between", "midCat"]
_TextAxisLblAlgn: TypeAlias = Literal["ctr", "l", "r"]
_DateAxisTimeUnit: TypeAlias = Literal["days", "months", "years"]

class ChartLines(Serialisable):
    tagname: ClassVar[str]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    def __init__(self, spPr: GraphicalProperties | None = None) -> None: ...

class Scaling(Serialisable):
    tagname: ClassVar[str]
    logBase: NestedFloat[Literal[True]]
    orientation: NestedSet[_ScalingOrientation]
    max: NestedFloat[Literal[True]]
    min: NestedFloat[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        logBase: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        orientation: _HasTagAndGet[_ScalingOrientation] | _ScalingOrientation = "minMax",
        max: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        min: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        extLst: Unused = None,
    ) -> None: ...

class _BaseAxis(Serialisable):
    axId: NestedInteger[Literal[False]]
    scaling: Typed[Scaling, Literal[False]]
    delete: NestedBool[Literal[True]]
    axPos: NestedSet[_BaseAxisAxPos]
    majorGridlines: Typed[ChartLines, Literal[True]]
    minorGridlines: Typed[ChartLines, Literal[True]]
    title: TitleDescriptor
    numFmt: Incomplete
    number_format: Alias
    majorTickMark: NestedNoneSet[_BaseAxisTickMark]
    minorTickMark: NestedNoneSet[_BaseAxisTickMark]
    tickLblPos: NestedNoneSet[_BaseAxisTickLblPos]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    txPr: Typed[RichText, Literal[True]]
    textProperties: Alias
    crossAx: NestedInteger[Literal[False]]
    crosses: NestedNoneSet[_BaseAxisCrosses]
    crossesAt: NestedFloat[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        axId: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt,
        scaling: Scaling | None,
        delete: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None,
        axPos: _HasTagAndGet[_BaseAxisAxPos] | _BaseAxisAxPos,
        majorGridlines: ChartLines | None,
        minorGridlines: ChartLines | None,
        title: str | Title | None,
        numFmt: Incomplete | None,
        majorTickMark: _NestedNoneSetParam[_BaseAxisTickMark],
        minorTickMark: _NestedNoneSetParam[_BaseAxisTickMark],
        tickLblPos: _NestedNoneSetParam[_BaseAxisTickLblPos],
        spPr: GraphicalProperties | None,
        txPr: RichText | None,
        crossAx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt,
        crosses: _NestedNoneSetParam[_BaseAxisCrosses] = None,
        crossesAt: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        axId: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt,
        scaling: Scaling | None = None,
        delete: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        axPos: _HasTagAndGet[_BaseAxisAxPos] | _BaseAxisAxPos = "l",
        majorGridlines: ChartLines | None = None,
        minorGridlines: ChartLines | None = None,
        title: str | Title | None = None,
        numFmt=None,
        majorTickMark=None,
        minorTickMark=None,
        tickLblPos=None,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
        *,
        crossAx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt,
        crosses=None,
        crossesAt: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
    ) -> None: ...

class DisplayUnitsLabel(Serialisable):
    tagname: ClassVar[str]
    layout: Typed[Layout, Literal[True]]
    tx: Typed[Text, Literal[True]]
    text: Alias
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    txPr: Typed[RichText, Literal[True]]
    textPropertes: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        layout: Layout | None = None,
        tx: Text | None = None,
        spPr: GraphicalProperties | None = None,
        txPr: RichText | None = None,
    ) -> None: ...

class DisplayUnitsLabelList(Serialisable):
    tagname: ClassVar[str]
    custUnit: NestedFloat[Literal[True]]
    builtInUnit: NestedNoneSet[_DisplayUnitsLabelListBuiltInUnit]
    dispUnitsLbl: Typed[DisplayUnitsLabel, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        custUnit: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        builtInUnit: (
            _HasTagAndGet[_DisplayUnitsLabelListBuiltInUnit] | _DisplayUnitsLabelListBuiltInUnit | Literal["none"] | None
        ) = None,
        dispUnitsLbl: DisplayUnitsLabel | None = None,
        extLst: Unused = None,
    ) -> None: ...

class NumericAxis(_BaseAxis):
    tagname: ClassVar[str]
    # Same as parent
    # axId = _BaseAxis.axId
    # scaling = _BaseAxis.scaling
    # delete = _BaseAxis.delete
    # axPos = _BaseAxis.axPos
    # majorGridlines = _BaseAxis.majorGridlines
    # minorGridlines = _BaseAxis.minorGridlines
    # title = _BaseAxis.title
    # numFmt = _BaseAxis.numFmt
    # majorTickMark = _BaseAxis.majorTickMark
    # minorTickMark = _BaseAxis.minorTickMark
    # tickLblPos = _BaseAxis.tickLblPos
    # spPr = _BaseAxis.spPr
    # txPr = _BaseAxis.txPr
    # crossAx = _BaseAxis.crossAx
    # crosses = _BaseAxis.crosses
    # crossesAt = _BaseAxis.crossesAt
    crossBetween: NestedNoneSet[_NumericAxisCrossBetween]
    majorUnit: NestedFloat[Literal[True]]
    minorUnit: NestedFloat[Literal[True]]
    dispUnits: Typed[DisplayUnitsLabelList, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        crossBetween: _NestedNoneSetParam[_NumericAxisCrossBetween] = None,
        majorUnit: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        minorUnit: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        dispUnits: DisplayUnitsLabelList | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...
    @classmethod
    def from_tree(cls, node: _SupportsFindAndIterAndAttribAndText) -> Self: ...

class TextAxis(_BaseAxis):
    tagname: ClassVar[str]
    # Same as parent
    # axId = _BaseAxis.axId
    # scaling = _BaseAxis.scaling
    # delete = _BaseAxis.delete
    # axPos = _BaseAxis.axPos
    # majorGridlines = _BaseAxis.majorGridlines
    # minorGridlines = _BaseAxis.minorGridlines
    # title = _BaseAxis.title
    # numFmt = _BaseAxis.numFmt
    # majorTickMark = _BaseAxis.majorTickMark
    # minorTickMark = _BaseAxis.minorTickMark
    # tickLblPos = _BaseAxis.tickLblPos
    # spPr = _BaseAxis.spPr
    # txPr = _BaseAxis.txPr
    # crossAx = _BaseAxis.crossAx
    # crosses = _BaseAxis.crosses
    # crossesAt = _BaseAxis.crossesAt
    auto: NestedBool[Literal[True]]
    lblAlgn: NestedNoneSet[_TextAxisLblAlgn]
    lblOffset: NestedMinMax[float, Literal[False]]
    tickLblSkip: NestedInteger[Literal[True]]
    tickMarkSkip: NestedInteger[Literal[True]]
    noMultiLvlLbl: NestedBool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        auto: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        lblAlgn: _NestedNoneSetParam[_TextAxisLblAlgn] = None,
        lblOffset: _HasTagAndGet[ConvertibleToFloat] | ConvertibleToFloat = 100,
        tickLblSkip: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        tickMarkSkip: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        noMultiLvlLbl: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...

class DateAxis(TextAxis):
    tagname: ClassVar[str]
    # Same as parent and grandparent
    # axId = _BaseAxis.axId
    # scaling = _BaseAxis.scaling
    # delete = _BaseAxis.delete
    # axPos = _BaseAxis.axPos
    # majorGridlines = _BaseAxis.majorGridlines
    # minorGridlines = _BaseAxis.minorGridlines
    # title = _BaseAxis.title
    # numFmt = _BaseAxis.numFmt
    # majorTickMark = _BaseAxis.majorTickMark
    # minorTickMark = _BaseAxis.minorTickMark
    # tickLblPos = _BaseAxis.tickLblPos
    # spPr = _BaseAxis.spPr
    # txPr = _BaseAxis.txPr
    # crossAx = _BaseAxis.crossAx
    # crosses = _BaseAxis.crosses
    # crossesAt = _BaseAxis.crossesAt
    auto: NestedBool[Literal[True]]
    lblOffset: NestedInteger[Literal[True]]  # type: ignore[assignment]
    baseTimeUnit: NestedNoneSet[_DateAxisTimeUnit]
    majorUnit: NestedFloat[Literal[True]]
    majorTimeUnit: NestedNoneSet[_DateAxisTimeUnit]
    minorUnit: NestedFloat[Literal[True]]
    minorTimeUnit: NestedNoneSet[_DateAxisTimeUnit]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        auto: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        lblOffset: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        baseTimeUnit: _NestedNoneSetParam[_DateAxisTimeUnit] = None,
        majorUnit: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        majorTimeUnit: _NestedNoneSetParam[_DateAxisTimeUnit] = None,
        minorUnit: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        minorTimeUnit: _NestedNoneSetParam[_DateAxisTimeUnit] = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...

class SeriesAxis(_BaseAxis):
    tagname: ClassVar[str]
    # Same as parent
    # axId = _BaseAxis.axId
    # scaling = _BaseAxis.scaling
    # delete = _BaseAxis.delete
    # axPos = _BaseAxis.axPos
    # majorGridlines = _BaseAxis.majorGridlines
    # minorGridlines = _BaseAxis.minorGridlines
    # title = _BaseAxis.title
    # numFmt = _BaseAxis.numFmt
    # majorTickMark = _BaseAxis.majorTickMark
    # minorTickMark = _BaseAxis.minorTickMark
    # tickLblPos = _BaseAxis.tickLblPos
    # spPr = _BaseAxis.spPr
    # txPr = _BaseAxis.txPr
    # crossAx = _BaseAxis.crossAx
    # crosses = _BaseAxis.crosses
    # crossesAt = _BaseAxis.crossesAt
    tickLblSkip: NestedInteger[Literal[True]]
    tickMarkSkip: NestedInteger[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        tickLblSkip: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        tickMarkSkip: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...
