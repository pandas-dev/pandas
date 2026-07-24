from _typeshed import Incomplete, Unused
from typing import ClassVar

from openpyxl.cell import _CellOrMergedCell
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.worksheet.worksheet import Worksheet

from .cell_range import CellRange

class MergeCell(CellRange):
    tagname: ClassVar[str]
    # Same as CellRange.coord
    # https://github.com/python/mypy/issues/6700
    @property
    def ref(self) -> str: ...
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, ref=None) -> None: ...
    def __copy__(self): ...

class MergeCells(Serialisable):
    tagname: ClassVar[str]
    # Overwritten by property below
    # count: Integer
    mergeCell: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, count: Unused = None, mergeCell=()) -> None: ...
    @property
    def count(self) -> int: ...

class MergedCellRange(CellRange):
    ws: Worksheet
    start_cell: _CellOrMergedCell
    def __init__(self, worksheet: Worksheet, coord) -> None: ...
    def format(self) -> None: ...
    def __contains__(self, coord: str) -> bool: ...
    def __copy__(self): ...
