from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, Integer, Set, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

_TableStyleElementType: TypeAlias = Literal[
    "wholeTable",
    "headerRow",
    "totalRow",
    "firstColumn",
    "lastColumn",
    "firstRowStripe",
    "secondRowStripe",
    "firstColumnStripe",
    "secondColumnStripe",
    "firstHeaderCell",
    "lastHeaderCell",
    "firstTotalCell",
    "lastTotalCell",
    "firstSubtotalColumn",
    "secondSubtotalColumn",
    "thirdSubtotalColumn",
    "firstSubtotalRow",
    "secondSubtotalRow",
    "thirdSubtotalRow",
    "blankRow",
    "firstColumnSubheading",
    "secondColumnSubheading",
    "thirdColumnSubheading",
    "firstRowSubheading",
    "secondRowSubheading",
    "thirdRowSubheading",
    "pageFieldLabels",
    "pageFieldValues",
]

class TableStyleElement(Serialisable):
    tagname: ClassVar[str]
    type: Set[_TableStyleElementType]
    size: Integer[Literal[True]]
    dxfId: Integer[Literal[True]]
    def __init__(
        self, type: _TableStyleElementType, size: ConvertibleToInt | None = None, dxfId: ConvertibleToInt | None = None
    ) -> None: ...

class TableStyle(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    pivot: Bool[Literal[True]]
    table: Bool[Literal[True]]
    count: Integer[Literal[True]]
    tableStyleElement: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        name: str,
        pivot: _ConvertibleToBool | None = None,
        table: _ConvertibleToBool | None = None,
        count: ConvertibleToInt | None = None,
        tableStyleElement=(),
    ) -> None: ...

class TableStyleList(Serialisable):
    tagname: ClassVar[str]
    defaultTableStyle: String[Literal[True]]
    defaultPivotStyle: String[Literal[True]]
    tableStyle: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        count: Unused = None,
        defaultTableStyle: str | None = "TableStyleMedium9",
        defaultPivotStyle: str | None = "PivotStyleLight16",
        tableStyle=(),
    ) -> None: ...
    @property
    def count(self) -> int: ...
