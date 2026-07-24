from _typeshed import ConvertibleToInt, Incomplete
from collections.abc import Generator, Iterable, Iterator
from types import GeneratorType
from typing import Any, Final, Literal, NoReturn, overload
from typing_extensions import deprecated

from openpyxl import _Decodable, _VisibilityType
from openpyxl.cell import _AnyCellValue, _CellGetValue, _CellOrMergedCell, _CellSetValue
from openpyxl.cell.cell import Cell
from openpyxl.chart._chart import ChartBase
from openpyxl.drawing.image import Image
from openpyxl.formatting.formatting import ConditionalFormattingList
from openpyxl.workbook.child import _WorkbookChild
from openpyxl.workbook.defined_name import DefinedNameDict
from openpyxl.workbook.workbook import Workbook
from openpyxl.worksheet.cell_range import CellRange, MultiCellRange
from openpyxl.worksheet.datavalidation import DataValidation, DataValidationList
from openpyxl.worksheet.dimensions import ColumnDimension, DimensionHolder, RowDimension, SheetFormatProperties
from openpyxl.worksheet.filters import AutoFilter
from openpyxl.worksheet.page import PageMargins, PrintOptions, PrintPageSetup
from openpyxl.worksheet.pagebreak import ColBreak, RowBreak
from openpyxl.worksheet.properties import WorksheetProperties
from openpyxl.worksheet.protection import SheetProtection
from openpyxl.worksheet.scenario import ScenarioList
from openpyxl.worksheet.table import Table, TableList
from openpyxl.worksheet.views import SheetView, SheetViewList

class Worksheet(_WorkbookChild):
    mime_type: str
    BREAK_NONE: Final = 0
    BREAK_ROW: Final = 1
    BREAK_COLUMN: Final = 2

    SHEETSTATE_VISIBLE: Final = "visible"
    SHEETSTATE_HIDDEN: Final = "hidden"
    SHEETSTATE_VERYHIDDEN: Final = "veryHidden"

    PAPERSIZE_LETTER: Final = "1"
    PAPERSIZE_LETTER_SMALL: Final = "2"
    PAPERSIZE_TABLOID: Final = "3"
    PAPERSIZE_LEDGER: Final = "4"
    PAPERSIZE_LEGAL: Final = "5"
    PAPERSIZE_STATEMENT: Final = "6"
    PAPERSIZE_EXECUTIVE: Final = "7"
    PAPERSIZE_A3: Final = "8"
    PAPERSIZE_A4: Final = "9"
    PAPERSIZE_A4_SMALL: Final = "10"
    PAPERSIZE_A5: Final = "11"

    ORIENTATION_PORTRAIT: Final = "portrait"
    ORIENTATION_LANDSCAPE: Final = "landscape"

    _cells: dict[tuple[int, int], _CellOrMergedCell]  # private but very useful to understand typing
    row_dimensions: DimensionHolder[int, RowDimension]
    column_dimensions: DimensionHolder[str, ColumnDimension]
    row_breaks: RowBreak
    col_breaks: ColBreak
    merged_cells: MultiCellRange
    data_validations: DataValidationList
    sheet_state: _VisibilityType
    page_setup: PrintPageSetup
    print_options: PrintOptions
    page_margins: PageMargins
    views: SheetViewList
    protection: SheetProtection
    defined_names: DefinedNameDict
    auto_filter: AutoFilter
    conditional_formatting: ConditionalFormattingList
    legacy_drawing: Incomplete | None
    sheet_properties: WorksheetProperties
    sheet_format: SheetFormatProperties
    scenarios: ScenarioList

    def __init__(self, parent: Workbook | None, title: str | _Decodable | None = None) -> None: ...
    @property
    def sheet_view(self) -> SheetView: ...
    @property
    def selected_cell(self) -> str | None: ...
    @property
    def active_cell(self) -> str | None: ...
    @property
    def array_formulae(self) -> dict[str, str]: ...
    @property
    def show_gridlines(self) -> bool | None: ...
    @property
    def freeze_panes(self) -> str | None: ...
    @freeze_panes.setter
    def freeze_panes(self, topLeftCell: str | Cell | None = None) -> None: ...
    # A MergedCell value should be kept to None
    @overload
    def cell(self, row: int, column: int, value: None = None) -> _CellOrMergedCell: ...
    @overload
    def cell(self, row: int, column: int, value: _CellSetValue = None) -> Cell: ...
    # An int is necessarily a row selection
    @overload
    def __getitem__(self, key: int) -> tuple[_CellOrMergedCell, ...]: ...
    # A slice is necessarily a row or rows, even if targetting a single cell
    @overload
    def __getitem__(self, key: slice) -> tuple[Any, ...]: ...  # tuple[AnyOf[_CellOrMergedCell, tuple[_CellOrMergedCell, ...]]]
    # A str could be an individual cell, row, column or full range
    @overload
    def __getitem__(
        self, key: str
    ) -> Any: ...  # AnyOf[_CellOrMergedCell, tuple[_CellOrMergedCell, ...], tuple[tuple[_CellOrMergedCell, ...], ...]]
    def __setitem__(self, key: str, value: _CellSetValue) -> None: ...
    def __iter__(self) -> Iterator[tuple[_CellOrMergedCell, ...]]: ...
    def __delitem__(self, key: str) -> None: ...
    @property
    def min_row(self) -> int: ...
    @property
    def max_row(self) -> int: ...
    @property
    def min_column(self) -> int: ...
    @property
    def max_column(self) -> int: ...
    def calculate_dimension(self) -> str: ...
    @property
    def dimensions(self) -> str: ...
    @overload
    def iter_rows(
        self, min_row: int | None, max_row: int | None, min_col: int | None, max_col: int | None, values_only: Literal[True]
    ) -> Generator[tuple[_CellGetValue, ...]]: ...
    @overload
    def iter_rows(
        self,
        min_row: int | None = None,
        max_row: int | None = None,
        min_col: int | None = None,
        max_col: int | None = None,
        *,
        values_only: Literal[True],
    ) -> Generator[tuple[_CellGetValue, ...]]: ...
    @overload
    def iter_rows(
        self,
        min_row: int | None = None,
        max_row: int | None = None,
        min_col: int | None = None,
        max_col: int | None = None,
        values_only: Literal[False] = False,
    ) -> Generator[tuple[_CellOrMergedCell, ...]]: ...
    @overload
    def iter_rows(
        self, min_row: int | None, max_row: int | None, min_col: int | None, max_col: int | None, values_only: bool
    ) -> Generator[tuple[_CellOrMergedCell, ...]] | Generator[tuple[_CellGetValue, ...]]: ...
    @overload
    def iter_rows(
        self,
        min_row: int | None = None,
        max_row: int | None = None,
        min_col: int | None = None,
        max_col: int | None = None,
        *,
        values_only: bool,
    ) -> Generator[tuple[_CellOrMergedCell, ...]] | Generator[tuple[_CellGetValue, ...]]: ...
    @property
    def rows(self) -> Generator[tuple[_CellOrMergedCell, ...]]: ...
    @property
    def values(self) -> Generator[tuple[_CellGetValue, ...]]: ...
    @overload
    def iter_cols(
        self, min_col: int | None, max_col: int | None, min_row: int | None, max_row: int | None, values_only: Literal[True]
    ) -> Generator[tuple[_CellGetValue, ...]]: ...
    @overload
    def iter_cols(
        self,
        min_col: int | None = None,
        max_col: int | None = None,
        min_row: int | None = None,
        max_row: int | None = None,
        *,
        values_only: Literal[True],
    ) -> Generator[tuple[_CellGetValue, ...]]: ...
    @overload
    def iter_cols(
        self,
        min_col: int | None = None,
        max_col: int | None = None,
        min_row: int | None = None,
        max_row: int | None = None,
        values_only: Literal[False] = False,
    ) -> Generator[tuple[_CellOrMergedCell, ...]]: ...
    @overload
    def iter_cols(
        self, min_col: int | None, max_col: int | None, min_row: int | None, max_row: int | None, values_only: bool
    ) -> Generator[tuple[_CellOrMergedCell, ...]] | Generator[tuple[_CellGetValue, ...]]: ...
    @overload
    def iter_cols(
        self,
        min_col: int | None = None,
        max_col: int | None = None,
        min_row: int | None = None,
        max_row: int | None = None,
        *,
        values_only: bool,
    ) -> Generator[tuple[_CellOrMergedCell, ...]] | Generator[tuple[_CellGetValue, ...]]: ...
    @property
    def columns(self) -> Generator[tuple[_CellOrMergedCell, ...]]: ...
    @property
    def column_groups(self) -> list[str]: ...
    def set_printer_settings(
        self, paper_size: int | None, orientation: Literal["default", "portrait", "landscape"] | None
    ) -> None: ...
    def add_data_validation(self, data_validation: DataValidation) -> None: ...
    def add_chart(self, chart: ChartBase, anchor: str | None = None) -> None: ...
    def add_image(self, img: Image, anchor: str | None = None) -> None: ...
    def add_table(self, table: Table) -> None: ...
    @property
    def tables(self) -> TableList: ...
    def add_pivot(self, pivot) -> None: ...
    # Same overload as CellRange.__init__
    @overload
    def merge_cells(
        self, range_string: str, start_row: None = None, start_column: None = None, end_row: None = None, end_column: None = None
    ) -> None: ...
    @overload
    def merge_cells(
        self,
        range_string: None = None,
        *,
        start_row: ConvertibleToInt,
        start_column: ConvertibleToInt,
        end_row: ConvertibleToInt,
        end_column: ConvertibleToInt,
    ) -> None: ...
    @overload
    def merge_cells(
        self,
        range_string: None,
        start_row: ConvertibleToInt,
        start_column: ConvertibleToInt,
        end_row: ConvertibleToInt,
        end_column: ConvertibleToInt,
    ) -> None: ...
    # Will always raise: TypeError: 'set' object is not subscriptable
    @property
    @deprecated("Use ws.merged_cells.ranges")
    def merged_cell_ranges(self) -> NoReturn: ...
    def unmerge_cells(
        self,
        range_string: str | None = None,
        start_row: int | None = None,
        start_column: int | None = None,
        end_row: int | None = None,
        end_column: int | None = None,
    ) -> None: ...
    def append(
        self,
        iterable: (
            list[_AnyCellValue]
            | tuple[_CellOrMergedCell | _CellGetValue, ...]
            | range
            | GeneratorType[_CellOrMergedCell | _CellGetValue, object, object]
            | dict[int | str, _AnyCellValue]
        ),
    ) -> None: ...
    def insert_rows(self, idx: int, amount: int = 1) -> None: ...
    def insert_cols(self, idx: int, amount: int = 1) -> None: ...
    def delete_rows(self, idx: int, amount: int = 1) -> None: ...
    def delete_cols(self, idx: int, amount: int = 1) -> None: ...
    def move_range(self, cell_range: CellRange | str, rows: int = 0, cols: int = 0, translate: bool = False) -> None: ...
    @property
    def print_title_rows(self) -> str | None: ...
    @print_title_rows.setter
    def print_title_rows(self, rows: str | None) -> None: ...
    @property
    def print_title_cols(self) -> str | None: ...
    @print_title_cols.setter
    def print_title_cols(self, cols: str | None) -> None: ...
    @property
    def print_titles(self) -> str: ...
    @property
    def print_area(self) -> str: ...
    @print_area.setter
    def print_area(self, value: str | Iterable[str] | None) -> None: ...
