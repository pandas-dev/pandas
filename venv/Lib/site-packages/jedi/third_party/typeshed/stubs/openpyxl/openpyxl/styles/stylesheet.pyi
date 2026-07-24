from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal, TypeVar
from typing_extensions import Self
from zipfile import ZipFile

from openpyxl.descriptors.base import Typed
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable, _ChildSerialisableTreeElement
from openpyxl.styles.cell_style import CellStyleList
from openpyxl.styles.colors import ColorList
from openpyxl.styles.named_styles import _NamedCellStyleList
from openpyxl.styles.numbers import NumberFormatList
from openpyxl.styles.table import TableStyleList
from openpyxl.workbook.workbook import Workbook
from openpyxl.xml.functions import Element

_WorkbookT = TypeVar("_WorkbookT", bound=Workbook)

class Stylesheet(Serialisable):
    tagname: ClassVar[str]
    numFmts: Typed[NumberFormatList, Literal[False]]
    fonts: Incomplete
    fills: Incomplete
    borders: Incomplete
    cellStyleXfs: Typed[CellStyleList, Literal[False]]
    cellXfs: Typed[CellStyleList, Literal[False]]
    cellStyles: Typed[_NamedCellStyleList, Literal[False]]
    dxfs: Incomplete
    tableStyles: Typed[TableStyleList, Literal[True]]
    colors: Typed[ColorList, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    number_formats: Incomplete
    cell_styles: Incomplete
    alignments: Incomplete
    protections: Incomplete
    named_styles: Incomplete
    def __init__(
        self,
        numFmts: NumberFormatList | None = None,
        fonts=(),
        fills=(),
        borders=(),
        cellStyleXfs: CellStyleList | None = None,
        cellXfs: CellStyleList | None = None,
        cellStyles: _NamedCellStyleList | None = None,
        dxfs=(),
        tableStyles: TableStyleList | None = None,
        colors: ColorList | None = None,
        extLst: Unused = None,
    ) -> None: ...
    @classmethod
    def from_tree(cls, node: _ChildSerialisableTreeElement) -> Self: ...
    @property
    def custom_formats(self) -> dict[int, str]: ...
    def to_tree(self, tagname: str | None = None, idx: Unused = None, namespace: str | None = None) -> Element: ...

def apply_stylesheet(archive: ZipFile, wb: _WorkbookT) -> _WorkbookT | None: ...
def write_stylesheet(wb: Workbook): ...
