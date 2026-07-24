from _typeshed import ConvertibleToInt, Incomplete, Unused
from collections.abc import Generator
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.cell import _CellOrMergedCell
from openpyxl.cell.text import Text
from openpyxl.comments.author import AuthorList
from openpyxl.comments.comments import Comment
from openpyxl.descriptors.base import Bool, Integer, Set, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.worksheet.ole import ObjectAnchor
from openpyxl.xml.functions import Element

_PropertiesTextHAlign: TypeAlias = Literal["left", "center", "right", "justify", "distributed"]
_PropertiesTextVAlign: TypeAlias = Literal["top", "center", "bottom", "justify", "distributed"]

class Properties(Serialisable):
    locked: Bool[Literal[True]]
    defaultSize: Bool[Literal[True]]
    _print: Bool[Literal[True]]  # Not private. Avoids name clash
    disabled: Bool[Literal[True]]
    uiObject: Bool[Literal[True]]
    autoFill: Bool[Literal[True]]
    autoLine: Bool[Literal[True]]
    altText: String[Literal[True]]
    textHAlign: Set[_PropertiesTextHAlign]
    textVAlign: Set[_PropertiesTextVAlign]
    lockText: Bool[Literal[True]]
    justLastX: Bool[Literal[True]]
    autoScale: Bool[Literal[True]]
    rowHidden: Bool[Literal[True]]
    colHidden: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    anchor: ObjectAnchor | None
    @overload
    def __init__(
        self,
        locked: _ConvertibleToBool | None = None,
        defaultSize: _ConvertibleToBool | None = None,
        _print: _ConvertibleToBool | None = None,
        disabled: _ConvertibleToBool | None = None,
        uiObject: _ConvertibleToBool | None = None,
        autoFill: _ConvertibleToBool | None = None,
        autoLine: _ConvertibleToBool | None = None,
        altText: str | None = None,
        *,
        textHAlign: _PropertiesTextHAlign,
        textVAlign: _PropertiesTextVAlign,
        lockText: _ConvertibleToBool | None = None,
        justLastX: _ConvertibleToBool | None = None,
        autoScale: _ConvertibleToBool | None = None,
        rowHidden: _ConvertibleToBool | None = None,
        colHidden: _ConvertibleToBool | None = None,
        anchor: ObjectAnchor | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        locked: _ConvertibleToBool | None,
        defaultSize: _ConvertibleToBool | None,
        _print: _ConvertibleToBool | None,
        disabled: _ConvertibleToBool | None,
        uiObject: _ConvertibleToBool | None,
        autoFill: _ConvertibleToBool | None,
        autoLine: _ConvertibleToBool | None,
        altText: str | None,
        textHAlign: _PropertiesTextHAlign,
        textVAlign: _PropertiesTextVAlign,
        lockText: _ConvertibleToBool | None = None,
        justLastX: _ConvertibleToBool | None = None,
        autoScale: _ConvertibleToBool | None = None,
        rowHidden: _ConvertibleToBool | None = None,
        colHidden: _ConvertibleToBool | None = None,
        anchor: ObjectAnchor | None = None,
    ) -> None: ...

class CommentRecord(Serialisable):
    tagname: ClassVar[str]
    ref: String[Literal[False]]
    authorId: Integer[Literal[False]]
    guid: Incomplete
    shapeId: Integer[Literal[True]]
    text: Typed[Text, Literal[False]]
    commentPr: Typed[Properties, Literal[True]]
    author: String[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    height: Incomplete
    width: Incomplete
    def __init__(
        self,
        ref: str = "",
        authorId: ConvertibleToInt = 0,
        guid=None,
        shapeId: ConvertibleToInt | None = 0,
        text: Text | None = None,
        commentPr: Properties | None = None,
        author: str | None = None,
        height: int = 79,
        width: int = 144,
    ) -> None: ...
    @classmethod
    def from_cell(cls, cell: _CellOrMergedCell): ...
    @property
    def content(self) -> str: ...

class CommentSheet(Serialisable):
    tagname: ClassVar[str]
    authors: Typed[AuthorList, Literal[False]]
    commentList: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    mime_type: str
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, authors: AuthorList, commentList=None, extLst: Unused = None) -> None: ...
    def to_tree(self) -> Element: ...  # type: ignore[override]
    @property
    def comments(self) -> Generator[tuple[str, Comment]]: ...
    @classmethod
    def from_comments(cls, comments): ...
    def write_shapes(self, vml=None): ...
    @property
    def path(self) -> str: ...
