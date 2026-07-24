from _typeshed import ConvertibleToInt, Unused
from re import Pattern
from typing import Final, Literal, overload
from typing_extensions import Self

from openpyxl.descriptors import Integer, Strict, String
from openpyxl.descriptors.base import Typed
from openpyxl.utils.cell import SHEETRANGE_RE as SHEETRANGE_RE

from .cell_range import MultiCellRange

COL_RANGE: Final[str]
COL_RANGE_RE: Final[Pattern[str]]
ROW_RANGE: Final[str]
ROW_RANGE_RE: Final[Pattern[str]]
TITLES_REGEX: Final[Pattern[str]]
PRINT_AREA_RE: Final[Pattern[str]]

class ColRange(Strict):
    min_col: String[Literal[False]]
    max_col: String[Literal[False]]
    @overload
    def __init__(self, range_string: None = None, *, min_col: str, max_col: str) -> None: ...
    @overload
    def __init__(self, range_string, min_col: Unused = None, max_col: Unused = None) -> None: ...
    def __eq__(self, other: object) -> bool: ...

class RowRange(Strict):
    min_row: Integer[Literal[False]]
    max_row: Integer[Literal[False]]
    @overload
    def __init__(self, range_string: None, min_row: ConvertibleToInt, max_row: ConvertibleToInt) -> None: ...
    @overload
    def __init__(self, range_string, min_row: Unused = None, max_row: Unused = None) -> None: ...
    def __eq__(self, other: object) -> bool: ...

class PrintTitles(Strict):
    cols: Typed[ColRange, Literal[True]]
    rows: Typed[RowRange, Literal[True]]
    title: String[Literal[False]]
    def __init__(self, cols: ColRange | None = None, rows: RowRange | None = None, title: str = "") -> None: ...
    @classmethod
    def from_string(cls, value: str) -> Self: ...
    def __eq__(self, other: object) -> bool: ...

class PrintArea(MultiCellRange):
    title: str
    @classmethod
    def from_string(cls, value) -> Self: ...
    def __init__(self, ranges=(), title: Unused = "") -> None: ...
    def __eq__(self, other: str | MultiCellRange) -> bool: ...  # type: ignore[override]
