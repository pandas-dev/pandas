from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.axis import NumericAxis, TextAxis
from openpyxl.chart.label import DataLabelList
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedSet

from ..xml._functions_overloads import _HasTagAndGet
from ._chart import ChartBase

_RadarChartRadarStyle: TypeAlias = Literal["standard", "marker", "filled"]

class RadarChart(ChartBase):
    tagname: ClassVar[str]
    radarStyle: NestedSet[_RadarChartRadarStyle]
    type: Alias
    varyColors: NestedBool[Literal[True]]
    ser: Incomplete
    dLbls: Typed[DataLabelList, Literal[True]]
    dataLabels: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[TextAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        radarStyle: _HasTagAndGet[_RadarChartRadarStyle] | _RadarChartRadarStyle = "standard",
        varyColors: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        ser=(),
        dLbls: DataLabelList | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...
