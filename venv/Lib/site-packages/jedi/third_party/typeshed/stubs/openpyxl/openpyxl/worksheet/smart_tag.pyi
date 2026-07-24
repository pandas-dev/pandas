from _typeshed import ConvertibleToInt, Incomplete
from typing import ClassVar, Literal, overload

from openpyxl.descriptors.base import Bool, Integer, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

class CellSmartTagPr(Serialisable):
    tagname: ClassVar[str]
    key: String[Literal[False]]
    val: String[Literal[False]]
    def __init__(self, key: str, val: str) -> None: ...

class CellSmartTag(Serialisable):
    tagname: ClassVar[str]
    cellSmartTagPr: Incomplete
    type: Integer[Literal[False]]
    deleted: Bool[Literal[True]]
    xmlBased: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        cellSmartTagPr=(),
        *,
        type: ConvertibleToInt,
        deleted: _ConvertibleToBool | None = False,
        xmlBased: _ConvertibleToBool | None = False,
    ) -> None: ...
    @overload
    def __init__(
        self,
        cellSmartTagPr,
        type: ConvertibleToInt,
        deleted: _ConvertibleToBool | None = False,
        xmlBased: _ConvertibleToBool | None = False,
    ) -> None: ...

class CellSmartTags(Serialisable):
    tagname: ClassVar[str]
    cellSmartTag: Incomplete
    r: String[Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, cellSmartTag, r: str) -> None: ...

class SmartTags(Serialisable):
    tagname: ClassVar[str]
    cellSmartTags: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, cellSmartTags=()) -> None: ...
