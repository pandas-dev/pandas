# Some documentations from upstream under MIT License. See authors in https://github.com/tafia/calamine
from __future__ import annotations

import datetime
import enum
import os
import types
import typing

@typing.type_check_only
class ReadBuffer(typing.Protocol):
    def seek(self, __offset: int, __whence: int = ...) -> int: ...
    def read(self, __size: int = ...) -> bytes | None: ...

@typing.final
class SheetTypeEnum(enum.Enum):
    WorkSheet = ...
    DialogSheet = ...
    MacroSheet = ...
    ChartSheet = ...
    Vba = ...

@typing.final
class SheetVisibleEnum(enum.Enum):
    Visible = ...
    """Visible."""
    Hidden = ...
    """Hidden."""
    VeryHidden = ...
    """The sheet is hidden and cannot be displayed using the user interface. It is supported only by Excel formats."""

@typing.final
class SheetMetadata:
    name: str
    """Name of sheet."""
    typ: SheetTypeEnum
    """Type of sheet.

    Only Excel formats support this. Default value for ODS is `WorkSheet`.
    """
    visible: SheetVisibleEnum
    """Visible of sheet."""

    def __new__(
        cls, name: str, typ: SheetTypeEnum, visible: SheetVisibleEnum
    ) -> SheetMetadata: ...

@typing.final
class CalamineSheet:
    name: str
    @property
    def height(self) -> int:
        """Get the row height of a sheet data.

        The height is defined as the number of rows between the start and end positions.
        """

    @property
    def width(self) -> int:
        """Get the column width of a sheet data.

        The width is defined as the number of columns between the start and end positions.
        """

    @property
    def total_height(self) -> int: ...
    @property
    def total_width(self) -> int: ...
    @property
    def start(self) -> tuple[int, int] | None:
        """Get top left cell position of a sheet data."""

    @property
    def end(self) -> tuple[int, int] | None:
        """Get bottom right cell position of a sheet data."""

    def to_python(
        self, skip_empty_area: bool = True, nrows: int | None = None
    ) -> list[
        list[
            int
            | float
            | str
            | bool
            | datetime.time
            | datetime.date
            | datetime.datetime
            | datetime.timedelta
        ]
    ]:
        """Retunrning data from sheet as list of lists.

        Args:
            skip_empty_area (bool):
                By default, calamine skips empty rows/cols before data.
                For suppress this behaviour, set `skip_empty_area` to `False`.
        """

    def iter_rows(
        self,
    ) -> typing.Iterator[
        list[
            int
            | float
            | str
            | bool
            | datetime.time
            | datetime.date
            | datetime.datetime
            | datetime.timedelta
        ]
    ]:
        """Retunrning data from sheet as iterator of lists."""

    @property
    def merged_cell_ranges(
        self,
    ) -> list[tuple[tuple[int, int], tuple[int, int]]] | None:
        """Return a copy of merged cell ranges.

        Support only for xlsx/xls.

        Returns:
            list of merged cell ranges (tuple[start coordinate, end coordinate]) or None for unsuported format
        """

@typing.final
class CalamineTable:
    name: str
    """Get the name of the table."""
    sheet: str
    """Get the name of the parent worksheet for a table."""
    columns: list[str]
    """Get the header names of the table columns.

    In Excel table headers can be hidden but the table will still have
    column header names.
    """
    @property
    def height(self) -> int:
        """Get the row height of a table data.

        The height is defined as the number of rows between the start and end positions.
        """

    @property
    def width(self) -> int:
        """Get the column width of a table data.

        The width is defined as the number of columns between the start and end positions.
        """

    @property
    def start(self) -> tuple[int, int] | None:
        """Get top left cell position of a table data."""

    @property
    def end(self) -> tuple[int, int] | None:
        """Get bottom right cell position of a table data."""

    def to_python(
        self,
    ) -> list[
        list[
            int
            | float
            | str
            | bool
            | datetime.time
            | datetime.date
            | datetime.datetime
            | datetime.timedelta
        ]
    ]:
        """Retunrning data from table as list of lists."""

@typing.final
class CalamineWorkbook:
    path: str | None
    """Path to file. `None` if bytes was loaded."""
    sheet_names: list[str]
    """All sheet names of this workbook, in workbook order."""
    sheets_metadata: list[SheetMetadata]
    """All sheets metadata of this workbook, in workbook order."""
    table_names: list[str] | None
    """All table names of this workbook."""
    @classmethod
    def from_object(
        cls, path_or_filelike: str | os.PathLike | ReadBuffer, load_tables: bool = False
    ) -> "CalamineWorkbook":
        """Determining type of pyobject and reading from it.

        Args:
            path_or_filelike (str | os.PathLike | ReadBuffer): path to file or IO (must imlpement read/seek methods).
            load_tables (bool): load Excel tables (supported for XLSX only).
        """

    @classmethod
    def from_path(
        cls, path: str | os.PathLike, load_tables: bool = False
    ) -> "CalamineWorkbook":
        """Reading file from path.

        Args:
            path (str | os.PathLike): path to file.
            load_tables (bool): load Excel tables (supported for XLSX only).
        """

    @classmethod
    def from_filelike(
        cls, filelike: ReadBuffer, load_tables: bool = False
    ) -> "CalamineWorkbook":
        """Reading file from IO.

        Args:
            filelike : IO (must imlpement read/seek methods).
            load_tables (bool): load Excel tables (supported for XLSX only).
        """

    def close(self) -> None:
        """Close the workbook.

        Drop internal rust structure from workbook (and close the file under the hood).
        `get_sheet_by_name`/`get_sheet_by_index` will raise WorkbookClosed after calling that method.

        Raises:
            WorkbookClosed: If workbook already closed.
        """

    def __enter__(self) -> "CalamineWorkbook": ...
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: types.TracebackType | None,
    ) -> None: ...
    def get_sheet_by_name(self, name: str) -> CalamineSheet:
        """Get worksheet by name.

        Args:
            name(str): name of worksheet

        Returns:
            CalamineSheet

        Raises:
            WorkbookClosed: If workbook already closed.
            WorksheetNotFound: If worksheet not found in workbook.
        """

    def get_sheet_by_index(self, index: int) -> CalamineSheet:
        """Get worksheet by index.

        Args:
            index(int): index of worksheet

        Returns:
            CalamineSheet

        Raises:
            WorkbookClosed: If workbook already closed.
            WorksheetNotFound: If worksheet not found in workbook.
        """

    def get_table_by_name(self, name: str) -> CalamineTable:
        """Get table by name.

        Args:
            name(str): name of table

        Returns:
            CalamineTable

        Raises:
            WorkbookClosed: If workbook already closed.
            WorksheetNotFound: If worksheet not found in workbook.
        """

class CalamineError(Exception): ...
class PasswordError(CalamineError): ...
class WorksheetNotFound(CalamineError): ...
class XmlError(CalamineError): ...
class ZipError(CalamineError): ...
class WorkbookClosed(CalamineError): ...
class TablesNotLoaded(CalamineError): ...
class TablesNotSupported(CalamineError): ...
class TableNotFound(CalamineError): ...

def load_workbook(
    path_or_filelike: str | os.PathLike | ReadBuffer, load_tables: bool = False
) -> CalamineWorkbook:
    """Determining type of pyobject and reading from it.

    Args:
        path_or_filelike (str | os.PathLike | ReadBuffer): path to file or IO (must imlpement read/seek methods).
        load_tables (bool): load Excel tables (supported for XLSX only).
    """

__all__ = [
    "CalamineError",
    "CalamineSheet",
    "CalamineTable",
    "CalamineWorkbook",
    "PasswordError",
    "SheetMetadata",
    "SheetTypeEnum",
    "SheetVisibleEnum",
    "TableNotFound",
    "TablesNotLoaded",
    "TablesNotSupported",
    "WorkbookClosed",
    "WorksheetNotFound",
    "XmlError",
    "ZipError",
    "load_workbook",
]
