from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal

from openpyxl.chart.axis import ChartLines
from openpyxl.descriptors.base import Typed
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable

class UpDownBars(Serialisable):
    tagname: ClassVar[str]
    gapWidth: Incomplete
    upBars: Typed[ChartLines, Literal[True]]
    downBars: Typed[ChartLines, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self, gapWidth: int = 150, upBars: ChartLines | None = None, downBars: ChartLines | None = None, extLst: Unused = None
    ) -> None: ...
