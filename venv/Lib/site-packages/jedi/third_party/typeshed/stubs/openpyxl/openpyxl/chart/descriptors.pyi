from typing import Literal

from openpyxl.chart.data_source import NumFmt
from openpyxl.descriptors import Strict, Typed
from openpyxl.descriptors.nested import NestedMinMax
from openpyxl.descriptors.serialisable import Serialisable

class NestedGapAmount(NestedMinMax[float, bool]):
    allow_none: bool
    min: float
    max: float

class NestedOverlap(NestedMinMax[float, bool]):
    allow_none: bool
    min: float
    max: float

class NumberFormatDescriptor(Typed[NumFmt, Literal[True]]):
    expected_type: type[NumFmt]
    allow_none: Literal[True]
    def __set__(self, instance: Serialisable | Strict, value) -> None: ...
