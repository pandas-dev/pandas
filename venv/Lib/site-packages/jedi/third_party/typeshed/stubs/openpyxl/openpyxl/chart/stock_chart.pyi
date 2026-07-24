from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal

from openpyxl.chart.axis import ChartLines, NumericAxis, TextAxis
from openpyxl.chart.label import DataLabelList
from openpyxl.chart.updown_bars import UpDownBars
from openpyxl.descriptors.base import Alias, Typed
from openpyxl.descriptors.excel import ExtensionList

from ._chart import ChartBase

class StockChart(ChartBase):
    tagname: ClassVar[str]
    ser: Incomplete
    dLbls: Typed[DataLabelList, Literal[True]]
    dataLabels: Alias
    dropLines: Typed[ChartLines, Literal[True]]
    hiLowLines: Typed[ChartLines, Literal[True]]
    upDownBars: Typed[UpDownBars, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[TextAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        ser=(),
        dLbls: DataLabelList | None = None,
        dropLines: ChartLines | None = None,
        hiLowLines: ChartLines | None = None,
        upDownBars: UpDownBars | None = None,
        extLst: Unused = None,
        **kw,
    ) -> None: ...
