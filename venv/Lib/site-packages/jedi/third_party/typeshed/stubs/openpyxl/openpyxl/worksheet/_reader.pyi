from _typeshed import Incomplete, SupportsGetItem, Unused
from collections.abc import Container, Generator
from datetime import datetime
from typing import Final
from xml.etree.ElementTree import _FileRead

from openpyxl.cell import _AnyCellValue, _CellOrMergedCell
from openpyxl.cell.rich_text import CellRichText
from openpyxl.descriptors.serialisable import _ChildSerialisableTreeElement, _SerialisableTreeElement
from openpyxl.formula.translate import Translator
from openpyxl.utils.cell import _RangeBoundariesTuple

from ..xml._functions_overloads import _HasAttrib, _SupportsIterAndAttrib
from .hyperlink import HyperlinkList
from .pagebreak import ColBreak, RowBreak
from .protection import SheetProtection
from .table import TablePartList
from .worksheet import Worksheet

CELL_TAG: Final[str]
VALUE_TAG: Final[str]
FORMULA_TAG: Final[str]
MERGE_TAG: Final[str]
INLINE_STRING: Final[str]
COL_TAG: Final[str]
ROW_TAG: Final[str]
CF_TAG: Final[str]
LEGACY_TAG: Final[str]
PROT_TAG: Final[str]
EXT_TAG: Final[str]
HYPERLINK_TAG: Final[str]
TABLE_TAG: Final[str]
PRINT_TAG: Final[str]
MARGINS_TAG: Final[str]
PAGE_TAG: Final[str]
HEADER_TAG: Final[str]
FILTER_TAG: Final[str]
VALIDATION_TAG: Final[str]
PROPERTIES_TAG: Final[str]
VIEWS_TAG: Final[str]
FORMAT_TAG: Final[str]
ROW_BREAK_TAG: Final[str]
COL_BREAK_TAG: Final[str]
SCENARIOS_TAG: Final[str]
DATA_TAG: Final[str]
DIMENSION_TAG: Final[str]
CUSTOM_VIEWS_TAG: Final[str]

def parse_richtext_string(element: _ChildSerialisableTreeElement) -> CellRichText | str: ...

class WorkSheetParser:
    min_row: Incomplete | None
    min_col: Incomplete | None
    epoch: datetime
    source: _FileRead
    shared_strings: SupportsGetItem[int, str]
    data_only: bool
    shared_formulae: dict[Incomplete, Translator]
    row_counter: int
    col_counter: int
    tables: TablePartList
    date_formats: Container[int]
    timedelta_formats: Container[int]
    row_dimensions: dict[Incomplete, Incomplete]
    column_dimensions: dict[Incomplete, Incomplete]
    number_formats: list[Incomplete]
    keep_vba: bool
    hyperlinks: HyperlinkList
    formatting: list[Incomplete]
    legacy_drawing: Incomplete | None
    merged_cells: Incomplete | None
    row_breaks: RowBreak
    col_breaks: ColBreak
    rich_text: bool
    protection: SheetProtection  # initialized after call to parse_sheet_protection()

    def __init__(
        self,
        src: _FileRead,
        shared_strings: SupportsGetItem[int, str],
        data_only: bool = False,
        epoch: datetime = ...,
        date_formats: Container[int] = ...,
        timedelta_formats: Container[int] = ...,
        rich_text: bool = False,
    ) -> None: ...
    def parse(self) -> Generator[Incomplete]: ...
    def parse_dimensions(self) -> _RangeBoundariesTuple | None: ...
    def parse_cell(self, element) -> dict[str, _AnyCellValue]: ...
    def parse_formula(self, element): ...
    def parse_column_dimensions(self, col: _HasAttrib) -> None: ...
    def parse_row(self, row: _SupportsIterAndAttrib) -> tuple[int, list[dict[str, _AnyCellValue]]]: ...
    def parse_formatting(self, element: _ChildSerialisableTreeElement) -> None: ...
    def parse_sheet_protection(self, element: _SerialisableTreeElement) -> None: ...
    def parse_extensions(self, element: _ChildSerialisableTreeElement) -> None: ...
    def parse_legacy(self, element: _ChildSerialisableTreeElement) -> None: ...
    def parse_row_breaks(self, element: _ChildSerialisableTreeElement) -> None: ...
    def parse_col_breaks(self, element: _ChildSerialisableTreeElement) -> None: ...
    def parse_custom_views(self, element: Unused) -> None: ...

class WorksheetReader:
    ws: Worksheet
    parser: WorkSheetParser
    tables: list[Incomplete]
    def __init__(
        self, ws: Worksheet, xml_source: _FileRead, shared_strings: SupportsGetItem[int, str], data_only: bool, rich_text: bool
    ) -> None: ...
    def bind_cells(self) -> None: ...
    def bind_formatting(self) -> None: ...
    def bind_tables(self) -> None: ...
    def bind_merged_cells(self) -> None: ...
    def bind_hyperlinks(self) -> None: ...
    def normalize_merged_cell_link(self, coord: str) -> _CellOrMergedCell | None: ...
    def bind_col_dimensions(self) -> None: ...
    def bind_row_dimensions(self) -> None: ...
    def bind_properties(self) -> None: ...
    def bind_all(self) -> None: ...
