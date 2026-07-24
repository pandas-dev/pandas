from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, Integer, NoneSet, Typed, _ConvertibleToBool
from openpyxl.descriptors.nested import NestedString, NestedText, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.styles.colors import Color
from openpyxl.styles.fonts import Font, _FontScheme, _FontU, _FontVertAlign

from ..xml._functions_overloads import _HasTagAndGet

_PhoneticPropertiesType: TypeAlias = Literal["halfwidthKatakana", "fullwidthKatakana", "Hiragana", "noConversion"]
_PhoneticPropertiesAlignment: TypeAlias = Literal["noControl", "left", "center", "distributed"]

class PhoneticProperties(Serialisable):
    tagname: ClassVar[str]
    fontId: Integer[Literal[False]]
    type: NoneSet[_PhoneticPropertiesType]
    alignment: NoneSet[_PhoneticPropertiesAlignment]
    def __init__(
        self,
        fontId: ConvertibleToInt,
        type: _PhoneticPropertiesType | Literal["none"] | None = None,
        alignment: _PhoneticPropertiesAlignment | Literal["none"] | None = None,
    ) -> None: ...

_PhoneticProperties: TypeAlias = PhoneticProperties

class PhoneticText(Serialisable):
    tagname: ClassVar[str]
    sb: Integer[Literal[False]]
    eb: Integer[Literal[False]]
    t: NestedText[str, Literal[False]]
    text: Alias
    def __init__(self, sb: ConvertibleToInt, eb: ConvertibleToInt, t: object = None) -> None: ...

class InlineFont(Font):
    tagname: ClassVar[str]
    rFont: NestedString[Literal[True]]
    # Same as parent
    # charset = Font.charset
    # family = Font.family
    # b = Font.b
    # i = Font.i
    # strike = Font.strike
    # outline = Font.outline
    # shadow = Font.shadow
    # condense = Font.condense
    # extend = Font.extend
    # color = Font.color
    # sz = Font.sz
    # u = Font.u
    # vertAlign = Font.vertAlign
    # scheme = Font.scheme
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        rFont: object = None,
        charset: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        family: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        b: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        i: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        strike: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        outline: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        shadow: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        condense: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        extend: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        color: Color | None = None,
        sz: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        u: _NestedNoneSetParam[_FontU] = None,
        vertAlign: _NestedNoneSetParam[_FontVertAlign] = None,
        scheme: _NestedNoneSetParam[_FontScheme] = None,
    ) -> None: ...

class RichText(Serialisable):
    tagname: ClassVar[str]
    rPr: Typed[InlineFont, Literal[True]]
    font: Alias
    t: NestedText[str, Literal[True]]
    text: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, rPr: InlineFont | None = None, t: object = None) -> None: ...

class Text(Serialisable):
    tagname: ClassVar[str]
    t: NestedText[str, Literal[True]]
    plain: Alias
    r: Incomplete
    formatted: Alias
    rPh: Incomplete
    phonetic: Alias
    phoneticPr: Typed[_PhoneticProperties, Literal[True]]
    PhoneticProperties: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, t: object = None, r=(), rPh=(), phoneticPr: _PhoneticProperties | None = None) -> None: ...
    @property
    def content(self) -> str: ...
