from typing import Final, Literal
from typing_extensions import TypeAlias
from zipfile import ZipFile

from openpyxl import _ZipFileFileProtocol
from openpyxl.chartsheet.chartsheet import Chartsheet
from openpyxl.packaging.manifest import Manifest
from openpyxl.packaging.relationship import Relationship
from openpyxl.reader.workbook import WorkbookParser
from openpyxl.workbook import Workbook

_SupportedFormats: TypeAlias = Literal[".xlsx", ".xlsm", ".xltx", ".xltm"]
SUPPORTED_FORMATS: Final[tuple[_SupportedFormats, ...]]

class ExcelReader:
    archive: ZipFile
    valid_files: list[str]
    read_only: bool
    keep_vba: bool
    data_only: bool
    keep_links: bool
    rich_text: bool
    shared_strings: list[str]
    package: Manifest  # defined after call to read_manifest()
    parser: WorkbookParser  # defined after call to read_workbook()
    wb: Workbook  # defined after call to read_workbook()

    def __init__(
        self,
        fn: _ZipFileFileProtocol,
        read_only: bool = False,
        keep_vba: bool = False,
        data_only: bool = False,
        keep_links: bool = True,
        rich_text: bool = False,
    ) -> None: ...
    def read_manifest(self) -> None: ...
    def read_strings(self) -> None: ...
    def read_workbook(self) -> None: ...
    def read_properties(self) -> None: ...
    def read_custom(self) -> None: ...
    def read_theme(self) -> None: ...
    def read_chartsheet(self, sheet: Chartsheet, rel: Relationship) -> None: ...
    def read_worksheets(self) -> None: ...
    def read(self) -> None: ...

def load_workbook(
    filename: _ZipFileFileProtocol,
    read_only: bool = False,
    keep_vba: bool = False,
    data_only: bool = False,
    keep_links: bool = True,
    rich_text: bool = False,
) -> Workbook: ...
