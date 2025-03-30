from __future__ import annotations

import contextlib
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
    Hidden = ...
    VeryHidden = ...

@typing.final
class SheetMetadata:
    name: str
    typ: SheetTypeEnum
    visible: SheetVisibleEnum

    def __init__(
        self, name: str, typ: SheetTypeEnum, visible: SheetVisibleEnum
    ) -> None: ...

@typing.final
class CalamineSheet:
    name: str
    @property
    def height(self) -> int: ...
    @property
    def width(self) -> int: ...
    @property
    def total_height(self) -> int: ...
    @property
    def total_width(self) -> int: ...
    @property
    def start(self) -> tuple[int, int] | None: ...
    @property
    def end(self) -> tuple[int, int] | None: ...
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

@typing.final
class CalamineWorkbook(contextlib.AbstractContextManager):
    path: str | None
    sheet_names: list[str]
    sheets_metadata: list[SheetMetadata]
    @classmethod
    def from_object(
        cls, path_or_filelike: str | os.PathLike | ReadBuffer
    ) -> "CalamineWorkbook":
        """Determining type of pyobject and reading from it.

        Args:
            path_or_filelike (str | os.PathLike | ReadBuffer): path to file or IO (must imlpement read/seek methods).
        """

    @classmethod
    def from_path(cls, path: str | os.PathLike) -> "CalamineWorkbook":
        """Reading file from path.

        Args:
            path (str | os.PathLike): path to file.
        """

    @classmethod
    def from_filelike(cls, filelike: ReadBuffer) -> "CalamineWorkbook":
        """Reading file from IO.

        Args:
            filelike : IO (must imlpement read/seek methods).
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

class CalamineError(Exception): ...
class PasswordError(CalamineError): ...
class WorksheetNotFound(CalamineError): ...
class XmlError(CalamineError): ...
class ZipError(CalamineError): ...
class WorkbookClosed(CalamineError): ...

def load_workbook(
    path_or_filelike: str | os.PathLike | ReadBuffer,
) -> CalamineWorkbook:
    """Determining type of pyobject and reading from it.

    Args:
        path_or_filelike (str | os.PathLike | ReadBuffer): path to file or IO (must imlpement read/seek methods).
    """
