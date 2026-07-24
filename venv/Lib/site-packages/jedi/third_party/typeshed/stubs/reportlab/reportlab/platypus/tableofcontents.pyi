from _typeshed import Unused
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Final, Literal, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias, Unpack

from reportlab.lib.styles import ParagraphStyle, PropertySet
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus.doctemplate import IndexingFlowable, _CanvasMaker
from reportlab.platypus.tables import TableStyle

_T = TypeVar("_T")
_Entry: TypeAlias = tuple[int, str, int] | tuple[int, str, int, str | None] | Sequence[int | str | None]
_SequencerFormat: TypeAlias = Literal["I", "i", "123", "ABC", "abc"]

@type_check_only
class _TableOfContentsKwargs(TypedDict, total=False):
    rightColumnWidth: float
    levelStyles: list[PropertySet]  # should be ParagraphStyle
    tableStyle: TableStyle
    dotsMinLevel: int
    formatter: Callable[[int], str] | None

@type_check_only
class _SimpleIndexKwargs(TypedDict, total=False):
    style: Iterable[PropertySet] | PropertySet | None  # should be ParagraphStyle
    dot: str | None
    tableStyle: TableStyle | None
    headers: bool
    name: str | None
    format: _SequencerFormat
    offset: int

__version__: Final[str]

def unquote(txt: str) -> str: ...
def drawPageNumbers(
    canvas: Canvas,
    style: PropertySet,  # should be ParagraphStyle
    pages: Iterable[tuple[int | str, Unused]],
    availWidth: float,
    availHeight: float,
    dot: str = " . ",
    formatter: Unused | None = None,
) -> None: ...

delta: float
epsilon: float
defaultLevelStyles: list[ParagraphStyle]
defaultTableStyle: TableStyle

class TableOfContents(IndexingFlowable):
    rightColumnWidth: float
    levelStyles: list[ParagraphStyle]
    tableStyle: TableStyle
    dotsMinLevel: int
    formatter: Callable[[int], str] | None
    def __init__(self, **kwds: Unpack[_TableOfContentsKwargs]) -> None: ...
    def isIndexing(self) -> Literal[1]: ...
    def isSatisfied(self) -> bool: ...
    def clearEntries(self) -> None: ...
    def getLevelStyle(self, n: int) -> ParagraphStyle: ...
    def addEntry(self, level: int, text: str, pageNum: int, key: str | None = None) -> None: ...
    def addEntries(self, listOfEntries: Iterable[_Entry]) -> None: ...

@overload
def makeTuple(x: tuple[_T, ...]) -> tuple[_T, ...]: ...
@overload
def makeTuple(x: list[_T]) -> tuple[_T, ...]: ...
@overload
def makeTuple(x: _T) -> tuple[_T, ...]: ...

class SimpleIndex(IndexingFlowable):
    # NOTE: Will be a list after getLevelStyle is called
    textStyle: ParagraphStyle | Iterable[ParagraphStyle] | list[ParagraphStyle]
    tableStyle: TableStyle
    dot: str | None
    headers: bool
    name: str
    formatFunc: Callable[[int], str]
    offset: float
    def __init__(self, **kwargs: Unpack[_SimpleIndexKwargs]) -> None: ...
    def getFormatFunc(self, formatName): ...
    def setup(
        self,
        style: PropertySet | None = None,  # should be ParagraphStyle
        dot: str | None = None,
        tableStyle: TableStyle | None = None,
        headers: bool = True,
        name: str | None = None,
        format: _SequencerFormat = "123",
        offset: float = 0,
    ) -> None: ...
    def __call__(self, canv: Canvas, kind: str | None, label: str) -> None: ...
    def getCanvasMaker(self, canvasmaker: _CanvasMaker = ...) -> _CanvasMaker: ...
    def isIndexing(self) -> Literal[1]: ...
    def isSatisfied(self) -> bool: ...
    def clearEntries(self) -> None: ...
    def addEntry(self, text: str, pageNum: tuple[int, str], key: str | None = None) -> None: ...
    def draw(self) -> None: ...
    def getLevelStyle(self, n: int) -> ParagraphStyle: ...

AlphabeticIndex = SimpleIndex

def listdiff(l1: list[Any], l2: list[_T]) -> tuple[int, list[_T]]: ...

class ReferenceText(IndexingFlowable):
    textPattern: str
    target: str
    paraStyle: ParagraphStyle
    def __init__(self, textPattern: str, targetKey: str) -> None: ...
