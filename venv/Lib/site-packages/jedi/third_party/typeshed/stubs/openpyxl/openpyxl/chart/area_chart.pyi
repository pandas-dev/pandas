from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.axis import ChartLines, NumericAxis, SeriesAxis, TextAxis
from openpyxl.chart.label import DataLabelList
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedSet

from ..xml._functions_overloads import _HasTagAndGet
from ._chart import ChartBase

_AreaChartBaseGrouping: TypeAlias = Literal["percentStacked", "standard", "stacked"]

class _AreaChartBase(ChartBase):
    grouping: NestedSet[_AreaChartBaseGrouping]
    varyColors: NestedBool[Literal[True]]
    ser: Incomplete
    dLbls: Typed[DataLabelList, Literal[True]]
    dataLabels: Alias
    dropLines: Typed[ChartLines, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        grouping: _HasTagAndGet[_AreaChartBaseGrouping] | _AreaChartBaseGrouping = "standard",
        varyColors: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        ser=(),
        dLbls: DataLabelList | None = None,
        dropLines: ChartLines | None = None,
    ) -> None: ...

class AreaChart(_AreaChartBase):
    tagname: ClassVar[str]
    # Same as parent
    # grouping = _AreaChartBase.grouping
    # varyColors = _AreaChartBase.varyColors
    # ser = _AreaChartBase.ser
    # dLbls = _AreaChartBase.dLbls
    # dropLines = _AreaChartBase.dropLines
    x_axis: Typed[TextAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, axId: Unused = None, extLst: Unused = None, **kw) -> None: ...

class AreaChart3D(AreaChart):
    tagname: ClassVar[str]
    # Same as parent and grandparent
    # grouping = _AreaChartBase.grouping
    # varyColors = _AreaChartBase.varyColors
    # ser = _AreaChartBase.ser
    # dLbls = _AreaChartBase.dLbls
    # dropLines = _AreaChartBase.dropLines
    gapDepth: Incomplete
    x_axis: Typed[TextAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    z_axis: Typed[SeriesAxis, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, gapDepth=None, **kw) -> None: ...
