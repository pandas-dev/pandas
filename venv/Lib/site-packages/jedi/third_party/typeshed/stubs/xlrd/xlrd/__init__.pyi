import sys
from _typeshed import SupportsWrite
from typing import Final

from . import timemachine as timemachine
from .biffh import (
    XL_CELL_BLANK as XL_CELL_BLANK,
    XL_CELL_BOOLEAN as XL_CELL_BOOLEAN,
    XL_CELL_DATE as XL_CELL_DATE,
    XL_CELL_EMPTY as XL_CELL_EMPTY,
    XL_CELL_ERROR as XL_CELL_ERROR,
    XL_CELL_NUMBER as XL_CELL_NUMBER,
    XL_CELL_TEXT as XL_CELL_TEXT,
    biff_text_from_num as biff_text_from_num,
    error_text_from_code as error_text_from_code,
)
from .book import Book as Book, colname as colname, open_workbook_xls as open_workbook_xls
from .formula import *
from .info import __VERSION__ as __VERSION__, __version__ as __version__
from .sheet import empty_cell as empty_cell
from .xldate import XLDateError as XLDateError, xldate_as_datetime as xldate_as_datetime, xldate_as_tuple as xldate_as_tuple

FILE_FORMAT_DESCRIPTIONS: Final[dict[str, str]]
ZIP_SIGNATURE: Final[bytes]
PEEK_SIZE: Final[int]

def inspect_format(path: str | None = None, content: bytes | None = None) -> str | None: ...
def open_workbook(
    filename: str | None = None,
    logfile: SupportsWrite[str] = sys.stdout,
    verbosity: int = 0,
    use_mmap: bool = True,
    file_contents: bytes | None = None,
    encoding_override: str | None = None,
    formatting_info: bool = False,
    on_demand: bool = False,
    ragged_rows: bool = False,
    ignore_workbook_corruption: bool = False,
) -> Book: ...
def dump(filename: str, outfile: SupportsWrite[str] = sys.stdout, unnumbered: bool = False) -> None: ...
def count_records(filename: str, outfile: SupportsWrite[str] = sys.stdout) -> None: ...
