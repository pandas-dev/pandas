from _typeshed import ConvertibleToFloat, ConvertibleToInt
from typing import ClassVar, Final, Literal
from typing_extensions import Self, TypeAlias

from openpyxl.descriptors.base import Alias, _ConvertibleToBool
from openpyxl.descriptors.nested import (
    NestedBool,
    NestedFloat,
    NestedInteger,
    NestedMinMax,
    NestedNoneSet,
    NestedString,
    _NestedNoneSetParam,
)
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.styles.colors import Color, ColorDescriptor

from ..xml._functions_overloads import _HasTagAndGet, _SupportsFindAndIterAndAttribAndText

_FontU: TypeAlias = Literal["single", "double", "singleAccounting", "doubleAccounting"]
_FontVertAlign: TypeAlias = Literal["superscript", "subscript", "baseline"]
_FontScheme: TypeAlias = Literal["major", "minor"]

class Font(Serialisable):
    UNDERLINE_DOUBLE: Final = "double"
    UNDERLINE_DOUBLE_ACCOUNTING: Final = "doubleAccounting"
    UNDERLINE_SINGLE: Final = "single"
    UNDERLINE_SINGLE_ACCOUNTING: Final = "singleAccounting"
    name: NestedString[Literal[True]]
    charset: NestedInteger[Literal[True]]
    family: NestedMinMax[float, Literal[True]]
    sz: NestedFloat[Literal[True]]
    size: Alias
    b: NestedBool[Literal[False]]
    bold: Alias
    i: NestedBool[Literal[False]]
    italic: Alias
    strike: NestedBool[Literal[True]]
    strikethrough: Alias
    outline: NestedBool[Literal[True]]
    shadow: NestedBool[Literal[True]]
    condense: NestedBool[Literal[True]]
    extend: NestedBool[Literal[True]]
    u: NestedNoneSet[_FontU]
    underline: Alias
    vertAlign: NestedNoneSet[_FontVertAlign]
    color: ColorDescriptor[Literal[True]]
    scheme: NestedNoneSet[_FontScheme]
    tagname: ClassVar[str]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        name: object = None,
        sz: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        b: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        i: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        charset: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        u: _NestedNoneSetParam[_FontU] = None,
        strike: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        color: str | Color | None = None,
        scheme: _NestedNoneSetParam[_FontScheme] = None,
        family: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        size: _HasTagAndGet[ConvertibleToFloat] | ConvertibleToFloat | None = None,
        bold: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool | None = None,
        italic: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool | None = None,
        strikethrough: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool | None = None,
        underline: _NestedNoneSetParam[_FontU] = None,
        vertAlign: _NestedNoneSetParam[_FontVertAlign] = None,
        outline: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        shadow: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        condense: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        extend: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
    ) -> None: ...
    @classmethod
    def from_tree(cls, node: _SupportsFindAndIterAndAttribAndText) -> Self: ...

DEFAULT_FONT: Final[Font]
