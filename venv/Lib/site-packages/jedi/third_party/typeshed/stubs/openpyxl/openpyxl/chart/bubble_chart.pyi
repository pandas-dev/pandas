from _typeshed import ConvertibleToFloat, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.axis import NumericAxis
from openpyxl.chart.label import DataLabelList
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedMinMax, NestedNoneSet, _NestedNoneSetParam

from ..xml._functions_overloads import _HasTagAndGet
from ._chart import ChartBase

_BubbleChartSizeRepresents: TypeAlias = Literal["area", "w"]

class BubbleChart(ChartBase):
    tagname: ClassVar[str]
    varyColors: NestedBool[Literal[True]]
    ser: Incomplete
    dLbls: Typed[DataLabelList, Literal[True]]
    dataLabels: Alias
    bubble3D: NestedBool[Literal[True]]
    bubbleScale: NestedMinMax[float, Literal[True]]
    showNegBubbles: NestedBool[Literal[True]]
    sizeRepresents: NestedNoneSet[_BubbleChartSizeRepresents]
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[NumericAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        varyColors: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        ser=(),
        dLbls: DataLabelList | None = None,
        bubble3D: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        bubbleScale: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        showNegBubbles: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        sizeRepresents: _NestedNoneSetParam[_BubbleChartSizeRepresents] = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...
