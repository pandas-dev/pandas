from _typeshed import ConvertibleToInt, Incomplete
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, Integer, Set, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

_WebPublishItemSourceType: TypeAlias = Literal[
    "sheet", "printArea", "autoFilter", "range", "chart", "pivotTable", "query", "label"
]

class WebPublishItem(Serialisable):
    tagname: ClassVar[str]
    id: Integer[Literal[False]]
    divId: String[Literal[False]]
    sourceType: Set[_WebPublishItemSourceType]
    sourceRef: String[Literal[False]]
    sourceObject: String[Literal[True]]
    destinationFile: String[Literal[False]]
    title: String[Literal[True]]
    autoRepublish: Bool[Literal[True]]
    @overload
    def __init__(
        self,
        id: ConvertibleToInt,
        divId: str,
        sourceType: _WebPublishItemSourceType,
        sourceRef: str,
        sourceObject: str | None = None,
        *,
        destinationFile: str,
        title: str | None = None,
        autoRepublish: _ConvertibleToBool | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        id: ConvertibleToInt,
        divId: str,
        sourceType: _WebPublishItemSourceType,
        sourceRef: str,
        sourceObject: str | None,
        destinationFile: str,
        title: str | None = None,
        autoRepublish: _ConvertibleToBool | None = None,
    ) -> None: ...

class WebPublishItems(Serialisable):
    tagname: ClassVar[str]
    count: Integer[Literal[True]]
    webPublishItem: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, count: ConvertibleToInt | None = None, webPublishItem=None) -> None: ...
