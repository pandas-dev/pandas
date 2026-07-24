from _typeshed import ConvertibleToInt, Unused
from collections.abc import Generator
from typing import Literal, overload

from openpyxl.descriptors import Strict
from openpyxl.descriptors.base import MinMax, String
from openpyxl.workbook.child import _WorkbookChild
from openpyxl.worksheet._read_only import ReadOnlyWorksheet

class DummyWorksheet:
    title: str
    def __init__(self, title: str) -> None: ...

class Reference(Strict):
    min_row: MinMax[int, Literal[False]]
    max_row: MinMax[int, Literal[False]]
    min_col: MinMax[int, Literal[False]]
    max_col: MinMax[int, Literal[False]]
    range_string: String[Literal[True]]
    worksheet: _WorkbookChild | ReadOnlyWorksheet | DummyWorksheet
    @overload
    def __init__(
        self,
        *,
        worksheet: _WorkbookChild | ReadOnlyWorksheet | DummyWorksheet | None = None,
        min_col: Unused = None,
        min_row: Unused = None,
        max_col: Unused = None,
        max_row: Unused = None,
        range_string: str,
    ) -> None: ...
    @overload
    def __init__(
        self,
        worksheet: _WorkbookChild | ReadOnlyWorksheet,
        min_col: ConvertibleToInt,
        min_row: ConvertibleToInt,
        max_col: ConvertibleToInt | None = None,
        max_row: ConvertibleToInt | None = None,
        range_string: str | None = None,
    ) -> None: ...
    def __len__(self) -> int: ...
    def __eq__(self, other: object) -> bool: ...
    @property
    def rows(self) -> Generator[Reference]: ...
    @property
    def cols(self) -> Generator[Reference]: ...
    def pop(self): ...
    @property
    def sheetname(self) -> str: ...
