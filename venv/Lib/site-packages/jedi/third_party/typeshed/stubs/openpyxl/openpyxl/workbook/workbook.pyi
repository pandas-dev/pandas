from _typeshed import Incomplete, Unused
from collections.abc import Iterator
from datetime import datetime
from typing import Any, Final, type_check_only
from typing_extensions import TypeAlias, deprecated
from zipfile import ZipFile

from openpyxl import _Decodable, _ZipFileFileWriteProtocol
from openpyxl.chartsheet.chartsheet import Chartsheet
from openpyxl.styles.named_styles import NamedStyle
from openpyxl.utils.indexed_list import IndexedList
from openpyxl.workbook.child import _WorkbookChild
from openpyxl.worksheet._read_only import ReadOnlyWorksheet
from openpyxl.worksheet._write_only import WriteOnlyWorksheet
from openpyxl.worksheet.worksheet import Worksheet

_WorkbookWorksheet: TypeAlias = Worksheet | WriteOnlyWorksheet | ReadOnlyWorksheet
_WorkbookSheet: TypeAlias = _WorkbookWorksheet | Chartsheet

# The type of worksheets in a workbook are the same as the aliases above.
# However, because Worksheet adds a lots of attributes that other _WorkbookChild subclasses
# don't have (ReadOnlyWorksheet doesn't even inherit from it), this ends up being too
# disruptive to the typical usage of openpyxl where sheets are just Worksheets.
# Using Any may just lose too much type information and duck-typing
# from Worksheet works great here. Allowing instance type check, even if direct
# type comparison might be wrong.
@type_check_only
class _WorksheetLike(  # type: ignore[misc] # Incompatible definitions, favor Worksheet
    Worksheet, WriteOnlyWorksheet, ReadOnlyWorksheet
): ...

@type_check_only
class _WorksheetOrChartsheetLike(  # type: ignore[misc] # Incompatible definitions, favor Worksheet
    Chartsheet, _WorksheetLike
): ...

INTEGER_TYPES: Final[tuple[type[int]]]

class Workbook:
    template: bool
    path: str
    defined_names: Incomplete
    properties: Incomplete
    security: Incomplete
    shared_strings: IndexedList[str]
    loaded_theme: Incomplete
    vba_archive: ZipFile | None
    is_template: bool
    code_name: Incomplete
    encoding: str
    iso_dates: Incomplete
    rels: Incomplete
    calculation: Incomplete
    views: Incomplete
    # Useful as a reference of what "sheets" can be for other types
    # ExcelReader can add ReadOnlyWorksheet in read_only mode.
    # _sheets: list[_WorksheetOrChartsheetLike]
    def __init__(self, write_only: bool = False, iso_dates: bool = False) -> None: ...
    @property
    def epoch(self) -> datetime: ...
    @epoch.setter
    def epoch(self, value: datetime) -> None: ...
    @property
    def read_only(self) -> bool: ...
    @property
    def data_only(self) -> bool: ...
    @property
    def write_only(self) -> bool: ...
    @property
    def excel_base_date(self) -> datetime: ...
    @property
    def active(self) -> _WorksheetOrChartsheetLike | None: ...
    @active.setter
    def active(self, value: Worksheet | Chartsheet | int) -> None: ...
    # read_only workbook cannot call this method
    # Could be generic based on write_only
    def create_sheet(
        self, title: str | _Decodable | None = None, index: int | None = None
    ) -> Any: ...  # AnyOf[WriteOnlyWorksheet, Worksheet]
    def move_sheet(self, sheet: Worksheet | str, offset: int = 0) -> None: ...
    def remove(self, worksheet: _WorkbookSheet) -> None: ...
    @deprecated("Use wb.remove(worksheet) or del wb[sheetname]")
    def remove_sheet(self, worksheet: _WorkbookSheet) -> None: ...
    def create_chartsheet(self, title: str | _Decodable | None = None, index: int | None = None) -> Chartsheet: ...
    @deprecated("Use wb[sheetname]")
    def get_sheet_by_name(self, name: str) -> _WorksheetOrChartsheetLike: ...
    def __contains__(self, key: str) -> bool: ...
    def index(self, worksheet: _WorkbookWorksheet) -> int: ...
    @deprecated("Use wb.index(worksheet)")
    def get_index(self, worksheet: _WorkbookWorksheet) -> int: ...
    def __getitem__(self, key: str) -> _WorksheetOrChartsheetLike: ...
    def __delitem__(self, key: str) -> None: ...
    def __iter__(self) -> Iterator[_WorksheetLike]: ...
    @deprecated("Use wb.sheetnames")
    def get_sheet_names(self) -> list[str]: ...
    @property
    def worksheets(self) -> list[_WorksheetLike]: ...
    @property
    def chartsheets(self) -> list[Chartsheet]: ...
    @property
    def sheetnames(self) -> list[str]: ...
    @deprecated("Assign scoped named ranges directly to worksheets or global ones to the workbook. Deprecated in 3.1")
    def create_named_range(
        self,
        name: str,
        worksheet: _WorkbookChild | ReadOnlyWorksheet | None = None,
        value: str | Incomplete | None = None,
        scope: Unused = None,
    ) -> None: ...
    def add_named_style(self, style: NamedStyle) -> None: ...
    @property
    def named_styles(self) -> list[str]: ...
    @property
    def mime_type(self) -> str: ...
    def save(self, filename: _ZipFileFileWriteProtocol) -> None: ...
    @property
    def style_names(self) -> list[str]: ...
    # A write_only and read_only workbooks can't use this method as it requires both reading and writing.
    # On an implementation level, a WorksheetCopy is created from the call to self.create_sheet,
    # but WorksheetCopy only works with Worksheet.
    def copy_worksheet(self, from_worksheet: Worksheet) -> Worksheet: ...
    def close(self) -> None: ...
