from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.axis import ChartLines, NumericAxis, SeriesAxis, TextAxis
from openpyxl.chart.label import DataLabelList
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedNoneSet, NestedSet, _NestedNoneSetParam

from ..xml._functions_overloads import _HasTagAndGet
from ._3d import _3DBase
from ._chart import ChartBase

_BarChartBaseBarDir: TypeAlias = Literal["bar", "col"]
_BarChartBaseGrouping: TypeAlias = Literal["percentStacked", "clustered", "standard", "stacked"]
_BarChart3DShape: TypeAlias = Literal["cone", "coneToMax", "box", "cylinder", "pyramid", "pyramidToMax"]

class _BarChartBase(ChartBase):
    barDir: NestedSet[_BarChartBaseBarDir]
    type: Alias
    grouping: NestedSet[_BarChartBaseGrouping]
    varyColors: NestedBool[Literal[True]]
    ser: Incomplete
    dLbls: Typed[DataLabelList, Literal[True]]
    dataLabels: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        barDir: _HasTagAndGet[_BarChartBaseBarDir] | _BarChartBaseBarDir = "col",
        grouping: _HasTagAndGet[_BarChartBaseGrouping] | _BarChartBaseGrouping = "clustered",
        varyColors: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        ser=(),
        dLbls: DataLabelList | None = None,
        **kw,
    ) -> None: ...

class BarChart(_BarChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # barDir = _BarChartBase.barDir
    # grouping = _BarChartBase.grouping
    # varyColors = _BarChartBase.varyColors
    # ser = _BarChartBase.ser
    # dLbls = _BarChartBase.dLbls
    gapWidth: Incomplete
    overlap: Incomplete
    serLines: Typed[ChartLines, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[TextAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    legend: Incomplete
    def __init__(
        self, gapWidth: int = 150, overlap=None, serLines: ChartLines | None = None, extLst: Unused = None, **kw
    ) -> None: ...

class BarChart3D(_BarChartBase, _3DBase):
    tagname: ClassVar[str]
    # Same as parents
    # barDir = _BarChartBase.barDir
    # grouping = _BarChartBase.grouping
    # varyColors = _BarChartBase.varyColors
    # ser = _BarChartBase.ser
    # dLbls = _BarChartBase.dLbls
    # view3D = _3DBase.view3D
    # floor = _3DBase.floor
    # sideWall = _3DBase.sideWall
    # backWall = _3DBase.backWall
    gapWidth: Incomplete
    gapDepth: Incomplete
    shape: NestedNoneSet[_BarChart3DShape]
    serLines: Typed[ChartLines, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[TextAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    z_axis: Typed[SeriesAxis, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        gapWidth: int = 150,
        gapDepth: int = 150,
        shape: _NestedNoneSetParam[_BarChart3DShape] = None,
        serLines: ChartLines | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...
