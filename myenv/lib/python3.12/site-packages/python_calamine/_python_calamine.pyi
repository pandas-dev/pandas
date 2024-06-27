from __future__ import annotations

import datetime
import enum
import os
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

        Parameters
        ----------
        skip_empty_area : bool
            By default, calamine skips empty rows/cols before data.
            For suppress this behaviour, set `skip_empty_area` to `False`.
        """

@typing.final
class CalamineWorkbook:
    path: str | None
    sheet_names: list[str]
    sheets_metadata: list[SheetMetadata]
    @classmethod
    def from_object(
        cls, path_or_filelike: str | os.PathLike | ReadBuffer
    ) -> "CalamineWorkbook":
        """Determining type of pyobject and reading from it.

        Parameters
        ----------
        path_or_filelike :
            - path (string),
            - pathlike (pathlib.Path),
            - IO (must imlpement read/seek methods).
        """

    @classmethod
    def from_path(cls, path: str | os.PathLike) -> "CalamineWorkbook":
        """Reading file from path.

        Parameters
        ----------
        path : path (string)."""

    @classmethod
    def from_filelike(cls, filelike: ReadBuffer) -> "CalamineWorkbook":
        """Reading file from IO.

        Parameters
        ----------
        filelike : IO (must imlpement read/seek methods).
        """

    def get_sheet_by_name(self, name: str) -> CalamineSheet: ...
    def get_sheet_by_index(self, index: int) -> CalamineSheet: ...

class CalamineError(Exception): ...
class PasswordError(CalamineError): ...
class WorksheetNotFound(CalamineError): ...
class XmlError(CalamineError): ...
class ZipError(CalamineError): ...

def load_workbook(
    path_or_filelike: str | os.PathLike | ReadBuffer,
) -> CalamineWorkbook:
    """Determining type of pyobject and reading from it.

    Parameters
    ----------
    path_or_filelike :
        - path (string),
        - pathlike (pathlib.Path),
        - IO (must imlpement read/seek methods).
    """
