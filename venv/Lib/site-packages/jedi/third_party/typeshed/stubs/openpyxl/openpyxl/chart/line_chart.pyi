from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.axis import ChartLines, NumericAxis, _BaseAxis
from openpyxl.chart.label import DataLabelList
from openpyxl.chart.updown_bars import UpDownBars
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedSet

from ..xml._functions_overloads import _HasTagAndGet
from ._chart import ChartBase

_LineChartBaseGrouping: TypeAlias = Literal["percentStacked", "standard", "stacked"]

class _LineChartBase(ChartBase):
    grouping: NestedSet[_LineChartBaseGrouping]
    varyColors: NestedBool[Literal[True]]
    ser: Incomplete
    dLbls: Typed[DataLabelList, Literal[True]]
    dataLabels: Alias
    dropLines: Typed[ChartLines, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        grouping: _HasTagAndGet[_LineChartBaseGrouping] | _LineChartBaseGrouping = "standard",
        varyColors: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        ser=(),
        dLbls: DataLabelList | None = None,
        dropLines: ChartLines | None = None,
        **kw,
    ) -> None: ...

class LineChart(_LineChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # grouping = _LineChartBase.grouping
    # varyColors = _LineChartBase.varyColors
    # ser = _LineChartBase.ser
    # dLbls = _LineChartBase.dLbls
    # dropLines = _LineChartBase.dropLines
    hiLowLines: Typed[ChartLines, Literal[True]]
    upDownBars: Typed[UpDownBars, Literal[True]]
    marker: NestedBool[Literal[True]]
    smooth: NestedBool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[_BaseAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        hiLowLines: ChartLines | None = None,
        upDownBars: UpDownBars | None = None,
        marker: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        smooth: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...

class LineChart3D(_LineChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # grouping = _LineChartBase.grouping
    # varyColors = _LineChartBase.varyColors
    # ser = _LineChartBase.ser
    # dLbls = _LineChartBase.dLbls
    # dropLines = _LineChartBase.dropLines
    gapDepth: Incomplete
    hiLowLines: Typed[ChartLines, Literal[True]]
    upDownBars: Typed[UpDownBars, Literal[True]]
    marker: NestedBool[Literal[True]]
    smooth: NestedBool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[ExtensionList, Literal[False]]
    y_axis: Typed[ExtensionList, Literal[False]]
    z_axis: Typed[ExtensionList, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        gapDepth=None,
        hiLowLines: ChartLines | None = None,
        upDownBars: UpDownBars | None = None,
        marker: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        smooth: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        **kw,
    ) -> None: ...
