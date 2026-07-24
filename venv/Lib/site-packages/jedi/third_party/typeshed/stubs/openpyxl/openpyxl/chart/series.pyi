from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.data_source import AxDataSource, NumDataSource, StrRef
from openpyxl.chart.error_bar import ErrorBars
from openpyxl.chart.label import DataLabelList
from openpyxl.chart.marker import Marker
from openpyxl.chart.picture import PictureOptions
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.trendline import Trendline
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedInteger, NestedNoneSet, NestedText, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.xml.functions import Element

from ..xml._functions_overloads import _HasTagAndGet

_SeriesShape: TypeAlias = Literal["cone", "coneToMax", "box", "cylinder", "pyramid", "pyramidToMax"]

attribute_mapping: Incomplete

class SeriesLabel(Serialisable):
    tagname: ClassVar[str]
    strRef: Typed[StrRef, Literal[True]]
    v: NestedText[str, Literal[True]]
    value: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, strRef: StrRef | None = None, v: object = None) -> None: ...

class Series(Serialisable):
    tagname: ClassVar[str]
    idx: NestedInteger[Literal[False]]
    order: NestedInteger[Literal[False]]
    tx: Typed[SeriesLabel, Literal[True]]
    title: Alias
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Incomplete
    pictureOptions: Typed[PictureOptions, Literal[True]]
    dPt: Incomplete
    data_points: Alias
    dLbls: Typed[DataLabelList, Literal[True]]
    labels: Alias
    trendline: Typed[Trendline, Literal[True]]
    errBars: Typed[ErrorBars, Literal[True]]
    cat: Typed[AxDataSource, Literal[True]]
    identifiers: Alias
    val: Typed[NumDataSource, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    invertIfNegative: NestedBool[Literal[True]]
    shape: NestedNoneSet[_SeriesShape]
    xVal: Typed[AxDataSource, Literal[True]]
    yVal: Typed[NumDataSource, Literal[True]]
    bubbleSize: Typed[NumDataSource, Literal[True]]
    zVal: Alias
    bubble3D: NestedBool[Literal[True]]
    marker: Typed[Marker, Literal[True]]
    smooth: NestedBool[Literal[True]]
    explosion: NestedInteger[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        idx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt = 0,
        order: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt = 0,
        tx: SeriesLabel | None = None,
        spPr: GraphicalProperties | None = None,
        pictureOptions: PictureOptions | None = None,
        dPt=(),
        dLbls: DataLabelList | None = None,
        trendline: Trendline | None = None,
        errBars: ErrorBars | None = None,
        cat: AxDataSource | None = None,
        val: NumDataSource | None = None,
        invertIfNegative: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        shape: _NestedNoneSetParam[_SeriesShape] = None,
        xVal: AxDataSource | None = None,
        yVal: NumDataSource | None = None,
        bubbleSize: NumDataSource | None = None,
        bubble3D: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        marker: Marker | None = None,
        smooth: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        explosion: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        extLst: Unused = None,
    ) -> None: ...
    def to_tree(  # type: ignore[override]
        self, tagname: str | None = None, idx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt | None = None
    ) -> Element: ...

class XYSeries(Series):
    # Same as parent
    # idx = Series.idx
    # order = Series.order
    # tx = Series.tx
    # spPr = Series.spPr
    # dPt = Series.dPt
    # dLbls = Series.dLbls
    # trendline = Series.trendline
    # errBars = Series.errBars
    # xVal = Series.xVal
    # yVal = Series.yVal
    # invertIfNegative = Series.invertIfNegative
    # bubbleSize = Series.bubbleSize
    # bubble3D = Series.bubble3D
    # marker = Series.marker
    # smooth = Series.smooth
    ...
