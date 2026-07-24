from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias
from zipfile import ZipFile

from openpyxl.descriptors.base import Bool, Integer, NoneSet, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.nested import NestedText
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.packaging.relationship import Relationship
from openpyxl.xml.functions import Element

_ExternalCellType: TypeAlias = Literal["b", "d", "n", "e", "s", "str", "inlineStr"]

class ExternalCell(Serialisable):
    r: String[Literal[False]]
    t: NoneSet[_ExternalCellType]
    vm: Integer[Literal[True]]
    v: NestedText[str, Literal[True]]
    def __init__(
        self, r: str, t: _ExternalCellType | Literal["none"] | None = None, vm: ConvertibleToInt | None = None, v: object = None
    ) -> None: ...

class ExternalRow(Serialisable):
    r: Integer[Literal[False]]
    cell: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, r: ConvertibleToInt, cell=None) -> None: ...

class ExternalSheetData(Serialisable):
    sheetId: Integer[Literal[False]]
    refreshError: Bool[Literal[True]]
    row: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, sheetId: ConvertibleToInt, refreshError: _ConvertibleToBool | None = None, row=()) -> None: ...

class ExternalSheetDataSet(Serialisable):
    sheetData: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, sheetData=None) -> None: ...

class ExternalSheetNames(Serialisable):
    sheetName: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, sheetName=()) -> None: ...

class ExternalDefinedName(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    refersTo: String[Literal[True]]
    sheetId: Integer[Literal[True]]
    def __init__(self, name: str, refersTo: str | None = None, sheetId: ConvertibleToInt | None = None) -> None: ...

class ExternalBook(Serialisable):
    tagname: ClassVar[str]
    sheetNames: Typed[ExternalSheetNames, Literal[True]]
    definedNames: Incomplete
    sheetDataSet: Typed[ExternalSheetDataSet, Literal[True]]
    id: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        sheetNames: ExternalSheetNames | None = None,
        definedNames=(),
        sheetDataSet: ExternalSheetDataSet | None = None,
        id=None,
    ) -> None: ...

class ExternalLink(Serialisable):
    tagname: ClassVar[str]
    mime_type: str
    externalBook: Typed[ExternalBook, Literal[True]]
    file_link: Typed[Relationship, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self, externalBook: ExternalBook | None = None, ddeLink: Unused = None, oleLink: Unused = None, extLst: Unused = None
    ) -> None: ...
    def to_tree(self) -> Element: ...  # type: ignore[override]
    @property
    def path(self) -> str: ...

def read_external_link(archive: ZipFile, book_path: str) -> ExternalLink: ...
