from _typeshed import ConvertibleToInt
from re import Pattern
from typing import ClassVar, Final, Literal
from typing_extensions import Self

from openpyxl.descriptors import Strict
from openpyxl.descriptors.base import Alias, Bool, Integer, MatchPattern, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.xml.functions import Element

from ..xml._functions_overloads import _HasText

FONT_PATTERN: Final = '&"(?P<font>.+)"'
COLOR_PATTERN: Final = "&K(?P<color>[A-F0-9]{6})"
SIZE_REGEX: Final = r"&(?P<size>\d+\s?)"
FORMAT_REGEX: Final[Pattern[str]]

class _HeaderFooterPart(Strict):
    text: String[Literal[True]]
    font: String[Literal[True]]
    size: Integer[Literal[True]]
    RGB: ClassVar[str]
    color: MatchPattern[str, Literal[True]]
    def __init__(
        self, text: str | None = None, font: str | None = None, size: ConvertibleToInt | None = None, color: str | None = None
    ) -> None: ...
    def __bool__(self) -> bool: ...
    @classmethod
    def from_str(cls, text): ...

class HeaderFooterItem(Strict):
    left: Typed[_HeaderFooterPart, Literal[False]]
    center: Typed[_HeaderFooterPart, Literal[False]]
    centre: Alias
    right: Typed[_HeaderFooterPart, Literal[False]]
    def __init__(
        self,
        left: _HeaderFooterPart | None = None,
        right: _HeaderFooterPart | None = None,
        center: _HeaderFooterPart | None = None,
    ) -> None: ...
    def __bool__(self) -> bool: ...
    def to_tree(self, tagname: str) -> Element: ...
    @classmethod
    def from_tree(cls, node: _HasText) -> Self: ...

class HeaderFooter(Serialisable):
    tagname: ClassVar[str]
    differentOddEven: Bool[Literal[True]]
    differentFirst: Bool[Literal[True]]
    scaleWithDoc: Bool[Literal[True]]
    alignWithMargins: Bool[Literal[True]]
    oddHeader: Typed[HeaderFooterItem, Literal[True]]
    oddFooter: Typed[HeaderFooterItem, Literal[True]]
    evenHeader: Typed[HeaderFooterItem, Literal[True]]
    evenFooter: Typed[HeaderFooterItem, Literal[True]]
    firstHeader: Typed[HeaderFooterItem, Literal[True]]
    firstFooter: Typed[HeaderFooterItem, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        differentOddEven: _ConvertibleToBool | None = None,
        differentFirst: _ConvertibleToBool | None = None,
        scaleWithDoc: _ConvertibleToBool | None = None,
        alignWithMargins: _ConvertibleToBool | None = None,
        oddHeader: HeaderFooterItem | None = None,
        oddFooter: HeaderFooterItem | None = None,
        evenHeader: HeaderFooterItem | None = None,
        evenFooter: HeaderFooterItem | None = None,
        firstHeader: HeaderFooterItem | None = None,
        firstFooter: HeaderFooterItem | None = None,
    ) -> None: ...
    def __bool__(self) -> bool: ...
