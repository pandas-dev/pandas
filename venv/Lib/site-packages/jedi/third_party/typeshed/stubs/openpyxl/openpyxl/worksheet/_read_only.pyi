from _typeshed import SupportsGetItem
from collections.abc import Generator

from openpyxl import _VisibilityType
from openpyxl.cell import _CellGetValue, _CellOrMergedCell
from openpyxl.utils.cell import _RangeBoundariesTuple
from openpyxl.workbook.workbook import Workbook
from openpyxl.worksheet.worksheet import Worksheet

def read_dimension(source) -> _RangeBoundariesTuple | None: ...

class ReadOnlyWorksheet:
    cell = Worksheet.cell
    iter_rows = Worksheet.iter_rows
    # Same as Worksheet.values
    # https://github.com/python/mypy/issues/6700
    @property
    def values(self) -> Generator[tuple[_CellGetValue, ...]]: ...
    # Same as Worksheet.rows
    # https://github.com/python/mypy/issues/6700
    @property
    def rows(self) -> Generator[tuple[_CellOrMergedCell, ...]]: ...
    __getitem__ = Worksheet.__getitem__
    __iter__ = Worksheet.__iter__
    parent: Workbook
    title: str
    sheet_state: _VisibilityType
    def __init__(
        self, parent_workbook: Workbook, title: str, worksheet_path, shared_strings: SupportsGetItem[int, str]
    ) -> None: ...
    def calculate_dimension(self, force: bool = False): ...
    def reset_dimensions(self) -> None: ...
    @property
    def min_row(self) -> int: ...
    @property
    def max_row(self) -> int | None: ...
    @property
    def min_column(self) -> int: ...
    @property
    def max_column(self) -> int | None: ...
