from _typeshed import Self
from collections.abc import Iterator, Mapping, Sequence, Set as AbstractSet
from typing import Any, ClassVar, TypedDict, type_check_only

from pygments.token import _TokenType

ansicolors: AbstractSet[str]  # not intended to be mutable

@type_check_only
class _StyleDict(TypedDict):
    color: str | None
    bold: bool
    italic: bool
    underline: bool
    bgcolor: str | None
    border: str | None
    roman: bool | None  # lol yes, can be True or False or None
    sans: bool | None
    mono: bool | None
    ansicolor: str | None
    bgansicolor: str | None

class StyleMeta(type):
    def __new__(cls: type[Self], name: str, bases: tuple[type[Any], ...], dct: dict[str, Any]) -> Self: ...
    def style_for_token(cls, token: _TokenType) -> _StyleDict: ...
    def styles_token(cls, ttype: _TokenType) -> bool: ...
    def list_styles(cls) -> list[tuple[_TokenType, _StyleDict]]: ...
    def __iter__(cls) -> Iterator[tuple[_TokenType, _StyleDict]]: ...
    def __len__(cls) -> int: ...

class Style(metaclass=StyleMeta):
    background_color: ClassVar[str]
    highlight_color: ClassVar[str]
    line_number_color: ClassVar[str]
    line_number_background_color: ClassVar[str]
    line_number_special_color: ClassVar[str]
    line_number_special_background_color: ClassVar[str]
    styles: ClassVar[Mapping[_TokenType, str]]  # not intended to be mutable
    name: ClassVar[str]
    aliases: ClassVar[Sequence[str]]  # not intended to be mutable
    web_style_gallery_exclude: ClassVar[bool]
