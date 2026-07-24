from typing import Literal
from zipfile import ZipFile

from openpyxl import _ZipFileFileProtocol
from openpyxl.packaging.manifest import Manifest
from openpyxl.workbook.workbook import Workbook
from openpyxl.worksheet.worksheet import Worksheet

class ExcelWriter:
    workbook: Workbook
    manifest: Manifest
    vba_modified: set[str | None]
    def __init__(self, workbook: Workbook, archive: ZipFile) -> None: ...
    def write_data(self) -> None: ...
    def write_worksheet(self, ws: Worksheet) -> None: ...
    def save(self) -> None: ...

def save_workbook(workbook: Workbook, filename: _ZipFileFileProtocol) -> Literal[True]: ...
