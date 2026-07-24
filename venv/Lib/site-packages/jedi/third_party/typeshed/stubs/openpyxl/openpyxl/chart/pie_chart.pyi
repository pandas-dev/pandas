from _typeshed import ConvertibleToFloat, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.axis import ChartLines
from openpyxl.chart.label import DataLabelList
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedFloat, NestedMinMax, NestedNoneSet, NestedSet, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet
from ._chart import ChartBase

_ProjectedPieChartOfPieType: TypeAlias = Literal["pie", "bar"]
_ProjectedPieChartSplitType: TypeAlias = Literal["auto", "cust", "percent", "pos", "val"]

class _PieChartBase(ChartBase):
    varyColors: NestedBool[Literal[True]]
    ser: Incomplete
    dLbls: Typed[DataLabelList, Literal[True]]
    dataLabels: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        varyColors: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = True,
        ser=(),
        dLbls: DataLabelList | None = None,
    ) -> None: ...

class PieChart(_PieChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # varyColors = _PieChartBase.varyColors
    # ser = _PieChartBase.ser
    # dLbls = _PieChartBase.dLbls
    firstSliceAng: NestedMinMax[float, Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self, firstSliceAng: _HasTagAndGet[ConvertibleToFloat] | ConvertibleToFloat = 0, extLst: Unused = None, **kw
    ) -> None: ...

class PieChart3D(_PieChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # varyColors = _PieChartBase.varyColors
    # ser = _PieChartBase.ser
    # dLbls = _PieChartBase.dLbls
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]

class DoughnutChart(_PieChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # varyColors = _PieChartBase.varyColors
    # ser = _PieChartBase.ser
    # dLbls = _PieChartBase.dLbls
    firstSliceAng: NestedMinMax[float, Literal[False]]
    holeSize: NestedMinMax[float, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        firstSliceAng: _HasTagAndGet[ConvertibleToFloat] | ConvertibleToFloat = 0,
        holeSize: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = 10,
        extLst: Unused = None,
        **kw,
    ) -> None: ...

class CustomSplit(Serialisable):
    tagname: ClassVar[str]
    secondPiePt: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, secondPiePt=()) -> None: ...

class ProjectedPieChart(_PieChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # varyColors = _PieChartBase.varyColors
    # ser = _PieChartBase.ser
    # dLbls = _PieChartBase.dLbls
    ofPieType: NestedSet[_ProjectedPieChartOfPieType]
    type: Alias
    gapWidth: Incomplete
    splitType: NestedNoneSet[_ProjectedPieChartSplitType]
    splitPos: NestedFloat[Literal[True]]
    custSplit: Typed[CustomSplit, Literal[True]]
    secondPieSize: NestedMinMax[float, Literal[True]]
    serLines: Typed[ChartLines, Literal[True]]
    join_lines: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        ofPieType: _HasTagAndGet[_ProjectedPieChartOfPieType] | _ProjectedPieChartOfPieType = "pie",
        gapWidth=None,
        splitType: _NestedNoneSetParam[_ProjectedPieChartSplitType] = "auto",
        splitPos: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        custSplit: CustomSplit | None = None,
        secondPieSize: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = 75,
        serLines: ChartLines | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...
