from _typeshed import ConvertibleToInt, Incomplete, Unused
from collections.abc import Iterator
from typing import ClassVar, Final, Literal, overload
from typing_extensions import Self, TypeAlias

from openpyxl.descriptors import Strict, String
from openpyxl.descriptors.base import Alias, Bool, Integer, NoneSet, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable, _ChildSerialisableTreeElement
from openpyxl.worksheet.filters import AutoFilter, SortState
from openpyxl.xml.functions import Element

_TableColumnTotalsRowFunction: TypeAlias = Literal[
    "sum", "min", "max", "average", "count", "countNums", "stdDev", "var", "custom"
]
_TableTableType: TypeAlias = Literal["worksheet", "xml", "queryTable"]

TABLESTYLES: Final[tuple[str, ...]]
PIVOTSTYLES: Final[tuple[str, ...]]

class TableStyleInfo(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[True]]
    showFirstColumn: Bool[Literal[True]]
    showLastColumn: Bool[Literal[True]]
    showRowStripes: Bool[Literal[True]]
    showColumnStripes: Bool[Literal[True]]
    def __init__(
        self,
        name: str | None = None,
        showFirstColumn: _ConvertibleToBool | None = None,
        showLastColumn: _ConvertibleToBool | None = None,
        showRowStripes: _ConvertibleToBool | None = None,
        showColumnStripes: _ConvertibleToBool | None = None,
    ) -> None: ...

class XMLColumnProps(Serialisable):
    tagname: ClassVar[str]
    mapId: Integer[Literal[False]]
    xpath: String[Literal[False]]
    denormalized: Bool[Literal[True]]
    xmlDataType: String[Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        mapId: ConvertibleToInt,
        xpath: str,
        denormalized: _ConvertibleToBool | None = None,
        *,
        xmlDataType: str,
        extLst: Unused = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        mapId: ConvertibleToInt,
        xpath: str,
        denormalized: _ConvertibleToBool | None,
        xmlDataType: str,
        extLst: Unused = None,
    ) -> None: ...

class TableFormula(Serialisable):
    tagname: ClassVar[str]
    array: Bool[Literal[True]]
    attr_text: Incomplete
    text: Alias
    def __init__(self, array: _ConvertibleToBool | None = None, attr_text=None) -> None: ...

class TableColumn(Serialisable):
    tagname: ClassVar[str]
    id: Integer[Literal[False]]
    uniqueName: String[Literal[True]]
    name: String[Literal[False]]
    totalsRowFunction: NoneSet[_TableColumnTotalsRowFunction]
    totalsRowLabel: String[Literal[True]]
    queryTableFieldId: Integer[Literal[True]]
    headerRowDxfId: Integer[Literal[True]]
    dataDxfId: Integer[Literal[True]]
    totalsRowDxfId: Integer[Literal[True]]
    headerRowCellStyle: String[Literal[True]]
    dataCellStyle: String[Literal[True]]
    totalsRowCellStyle: String[Literal[True]]
    calculatedColumnFormula: Typed[TableFormula, Literal[True]]
    totalsRowFormula: Typed[TableFormula, Literal[True]]
    xmlColumnPr: Typed[XMLColumnProps, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        id: ConvertibleToInt,
        uniqueName: str | None = None,
        *,
        name: str,
        totalsRowFunction: _TableColumnTotalsRowFunction | Literal["none"] | None = None,
        totalsRowLabel: str | None = None,
        queryTableFieldId: ConvertibleToInt | None = None,
        headerRowDxfId: ConvertibleToInt | None = None,
        dataDxfId: ConvertibleToInt | None = None,
        totalsRowDxfId: ConvertibleToInt | None = None,
        headerRowCellStyle: str | None = None,
        dataCellStyle: str | None = None,
        totalsRowCellStyle: str | None = None,
        calculatedColumnFormula: TableFormula | None = None,
        totalsRowFormula: TableFormula | None = None,
        xmlColumnPr: XMLColumnProps | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        id: ConvertibleToInt,
        uniqueName: str | None,
        name: str,
        totalsRowFunction: _TableColumnTotalsRowFunction | Literal["none"] | None = None,
        totalsRowLabel: str | None = None,
        queryTableFieldId: ConvertibleToInt | None = None,
        headerRowDxfId: ConvertibleToInt | None = None,
        dataDxfId: ConvertibleToInt | None = None,
        totalsRowDxfId: ConvertibleToInt | None = None,
        headerRowCellStyle: str | None = None,
        dataCellStyle: str | None = None,
        totalsRowCellStyle: str | None = None,
        calculatedColumnFormula: TableFormula | None = None,
        totalsRowFormula: TableFormula | None = None,
        xmlColumnPr: XMLColumnProps | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...
    def __iter__(self) -> Iterator[tuple[str, str]]: ...
    @classmethod
    def from_tree(cls, node: _ChildSerialisableTreeElement) -> Self: ...

class TableNameDescriptor(String[Incomplete]):
    def __set__(self, instance: Serialisable | Strict, value) -> None: ...

class Table(Serialisable):
    mime_type: str
    tagname: ClassVar[str]
    id: Integer[Literal[False]]
    name: String[Literal[True]]
    displayName: Incomplete
    comment: String[Literal[True]]
    ref: Incomplete
    tableType: NoneSet[_TableTableType]
    headerRowCount: Integer[Literal[True]]
    insertRow: Bool[Literal[True]]
    insertRowShift: Bool[Literal[True]]
    totalsRowCount: Integer[Literal[True]]
    totalsRowShown: Bool[Literal[True]]
    published: Bool[Literal[True]]
    headerRowDxfId: Integer[Literal[True]]
    dataDxfId: Integer[Literal[True]]
    totalsRowDxfId: Integer[Literal[True]]
    headerRowBorderDxfId: Integer[Literal[True]]
    tableBorderDxfId: Integer[Literal[True]]
    totalsRowBorderDxfId: Integer[Literal[True]]
    headerRowCellStyle: String[Literal[True]]
    dataCellStyle: String[Literal[True]]
    totalsRowCellStyle: String[Literal[True]]
    connectionId: Integer[Literal[True]]
    autoFilter: Typed[AutoFilter, Literal[True]]
    sortState: Typed[SortState, Literal[True]]
    tableColumns: Incomplete
    tableStyleInfo: Typed[TableStyleInfo, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        id: ConvertibleToInt = 1,
        displayName=None,
        ref=None,
        name: str | None = None,
        comment: str | None = None,
        tableType: _TableTableType | Literal["none"] | None = None,
        headerRowCount: ConvertibleToInt | None = 1,
        insertRow: _ConvertibleToBool | None = None,
        insertRowShift: _ConvertibleToBool | None = None,
        totalsRowCount: ConvertibleToInt | None = None,
        totalsRowShown: _ConvertibleToBool | None = None,
        published: _ConvertibleToBool | None = None,
        headerRowDxfId: ConvertibleToInt | None = None,
        dataDxfId: ConvertibleToInt | None = None,
        totalsRowDxfId: ConvertibleToInt | None = None,
        headerRowBorderDxfId: ConvertibleToInt | None = None,
        tableBorderDxfId: ConvertibleToInt | None = None,
        totalsRowBorderDxfId: ConvertibleToInt | None = None,
        headerRowCellStyle: str | None = None,
        dataCellStyle: str | None = None,
        totalsRowCellStyle: str | None = None,
        connectionId: ConvertibleToInt | None = None,
        autoFilter: AutoFilter | None = None,
        sortState: SortState | None = None,
        tableColumns=(),
        tableStyleInfo: TableStyleInfo | None = None,
        extLst: Unused = None,
    ) -> None: ...
    def to_tree(self) -> Element: ...  # type: ignore[override]
    @property
    def path(self) -> str: ...
    @property
    def column_names(self) -> list[str]: ...

class TablePartList(Serialisable):
    tagname: ClassVar[str]
    # Overwritten by property below
    # count: Integer
    tablePart: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, count: Unused = None, tablePart=()) -> None: ...
    def append(self, part) -> None: ...
    @property
    def count(self) -> int: ...
    def __bool__(self) -> bool: ...

class TableList(dict[Incomplete, Incomplete]):
    def add(self, table) -> None: ...
    def get(self, name=None, table_range=None): ...
    def items(self): ...
