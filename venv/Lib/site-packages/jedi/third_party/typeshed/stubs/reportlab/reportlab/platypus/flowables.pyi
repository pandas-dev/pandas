from _typeshed import Incomplete, SupportsRead, Unused
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Literal, NoReturn, Protocol, type_check_only
from typing_extensions import Self, TypeAlias

from reportlab.lib.colors import Color
from reportlab.lib.styles import ListStyle, ParagraphStyle, PropertySet
from reportlab.pdfgen.canvas import Canvas
from reportlab.pdfgen.textobject import _Color
from reportlab.platypus.paragraph import Paragraph

__all__ = [
    "AnchorFlowable",
    "BalancedColumns",
    "BulletDrawer",
    "CallerMacro",
    "CondPageBreak",
    "DDIndenter",
    "DocAssert",
    "DocAssign",
    "DocExec",
    "DocIf",
    "DocPara",
    "DocWhile",
    "FailOnDraw",
    "FailOnWrap",
    "Flowable",
    "FrameBG",
    "FrameSplitter",
    "HRFlowable",
    "Image",
    "ImageAndFlowables",
    "KeepInFrame",
    "KeepTogether",
    "LIIndenter",
    "ListFlowable",
    "ListItem",
    "Macro",
    "NullDraw",
    "PTOContainer",
    "PageBreak",
    "PageBreakIfNotEmpty",
    "ParagraphAndImage",
    "Preformatted",
    "SetPageTopFlowables",
    "SetTopFlowables",
    "SlowPageBreak",
    "Spacer",
    "TopPadder",
    "TraceInfo",
    "UseUpSpace",
    "XBox",
    "splitLine",
    "splitLines",
    "PlacedStory",
]

_HAlignment: TypeAlias = Literal["LEFT", "CENTER", "CENTRE", "RIGHT", 0, 1, 2]
_VAlignment: TypeAlias = Literal["BOTTOM", "MIDDLE", "TOP"]
# FIXME: Consider using Sequence[Flowable] for covariance on list, even though
#        that will give false negatives for non list or tuple sequences, it also
#        would reduce type safety, since flowables don't copy the list
_FlowableSublist: TypeAlias = Flowable | list[Flowable] | tuple[Flowable, ...]
# NOTE: Technically can only be list or tuple, but would be annoying for variance
_NestedFlowable: TypeAlias = Flowable | Sequence[_NestedFlowable]

@type_check_only
class _StyledFlowableFactory(Protocol):
    # NOTE: We leave style at Any so people can specify a specifc property set
    def __call__(self, value: str, /, *, style: Any) -> Flowable: ...

class TraceInfo:
    srcFile: str
    startLineNo: int
    startLinePos: int
    endLineNo: int
    endLinePos: int
    def __init__(self) -> None: ...

class Flowable:
    width: float
    height: float
    wrapped: int
    hAlign: _HAlignment
    vAlign: _VAlignment
    encoding: str | None
    # NOTE: this only exists during drawing, splitting and wrapping
    canv: Canvas
    # NOTE: The following attributes will not exist on all flowables, but
    #       they need to be settable on individual instances
    keepWithNext: Incomplete
    spaceAfter: float
    spaceBefore: float
    def __init__(self) -> None: ...
    # NOTE: We pretend the optional internal _sW argument does not exist
    #       since not all flowables support it and we'd have to deal with
    #       a bunch of LSP errors. Conversely we will get type errors in
    #       subclasses that rely on the argument existing when called through
    #       super() inside their own implementation, so we can't really
    #       make everyone happy here, sigh...
    def drawOn(self, canvas: Canvas, x: float, y: float) -> None: ...
    def wrapOn(self, canv: Canvas, aW: float, aH: float) -> tuple[float, float]: ...
    def wrap(self, aW: float, aH: float) -> tuple[float, float]: ...
    def minWidth(self) -> float: ...
    def splitOn(self, canv: Canvas, aW: float, aH: float) -> list[Flowable]: ...
    def split(self, aW: float, aH: float, /) -> list[Flowable]: ...
    def getKeepWithNext(self): ...
    def getSpaceAfter(self) -> float: ...
    def getSpaceBefore(self) -> float: ...
    def isIndexing(self) -> int: ...
    def identity(self, maxLen: int | None = None) -> str: ...

class XBox(Flowable):
    text: str
    def __init__(self, width: float, height: float, text: str = "A Box") -> None: ...
    def draw(self) -> None: ...

def splitLines(lines, maximum_length, split_characters, new_line_characters): ...
def splitLine(line_to_split, lines_splitted, maximum_length, split_characters, new_line_characters) -> None: ...

class Preformatted(Flowable):
    style: ParagraphStyle
    bulletText: str | None
    lines: list[str]
    def __init__(
        self,
        text: str,
        # NOTE: Technically has to be a ParagraphStyle, but that would
        #       conflict with stylesheet["Style"] usage
        style: PropertySet,
        bulletText: str | None = None,
        dedent: int = 0,
        maxLineLength: int | None = None,
        splitChars: str | None = None,
        newLineChars: str = "",
    ) -> None: ...
    def draw(self) -> None: ...

class Image(Flowable):
    filename: str
    # these are lazy, but __getattr__ ensures the image gets loaded
    # as soon as these attributes are accessed
    imageWidth: int
    imageHeight: int
    drawWidth: float
    drawHeight: float
    def __init__(
        self,
        # TODO: I think this might also accept a PIL.Image and other
        #       kinds of path represenations, should be kept in sync
        #       with reportlab.lib.utils.ImageReader, except for the
        #       potential PIL.Image shortcut
        filename: str | SupportsRead[bytes] | Incomplete,
        width: float | None = None,
        height: float | None = None,
        kind: str = "direct",
        mask: str = "auto",
        lazy: int = 1,
        hAlign: _HAlignment = "CENTER",
        useDPI: bool = False,
    ) -> None: ...
    def draw(self) -> None: ...

class NullDraw(Flowable):
    def draw(self) -> None: ...

class Spacer(NullDraw):
    # NOTE: This may actually be a bug, it seems likely that Spacer is meant
    #       to set spaceBefore in the isGlue case.
    spacebefore: float
    def __init__(self, width: float, height: float, isGlue: bool = False) -> None: ...

class UseUpSpace(NullDraw): ...

class PageBreak(UseUpSpace):
    locChanger: int
    nextTemplate: str | None
    def __init__(self, nextTemplate: str | None = None) -> None: ...

class SlowPageBreak(PageBreak): ...
class PageBreakIfNotEmpty(PageBreak): ...

class CondPageBreak(Spacer):
    locChanger: int
    def __init__(self, height: float) -> None: ...

class _ContainerSpace:
    def getSpaceBefore(self) -> float: ...
    def getSpaceAfter(self) -> float: ...

class KeepTogether(_ContainerSpace, Flowable):
    splitAtTop: bool
    # TODO: Consider using Sequence[Flowable] for covariance, even if reportlab
    #       only supports list/tuple
    def __init__(self, flowables: _FlowableSublist | None, maxHeight=None) -> None: ...

class KeepTogetherSplitAtTop(KeepTogether):
    splitAtTop: bool

class Macro(Flowable):
    command: str
    def __init__(self, command: str) -> None: ...
    def draw(self) -> None: ...

class CallerMacro(Flowable):
    def __init__(
        self,
        drawCallable: Callable[[CallerMacro, float, float], object] | None = None,
        wrapCallable: Callable[[CallerMacro, float, float], object] | None = None,
    ) -> None: ...
    def draw(self) -> None: ...

class ParagraphAndImage(Flowable):
    P: Paragraph
    I: Image
    xpad: float
    ypad: float
    def __init__(self, P: Paragraph, I: Image, xpad: float = 3, ypad: float = 3, side: str = "right") -> None: ...
    def draw(self) -> None: ...

class FailOnWrap(NullDraw):
    def wrap(self, aW: float, aH: float) -> NoReturn: ...

class FailOnDraw(Flowable):
    def draw(self) -> NoReturn: ...

class HRFlowable(Flowable):
    width: float | str  # type: ignore[assignment]
    lineWidth: float
    lineCap: str
    color: _Color
    dash: Incomplete | None
    def __init__(
        self,
        width: float | str = "80%",
        thickness: float = 1,
        lineCap: str = "round",
        color: _Color = ...,
        spaceBefore: float = 1,
        spaceAfter: float = 1,
        hAlign: _HAlignment = "CENTER",
        vAlign: _VAlignment = "BOTTOM",
        dash=None,
    ) -> None: ...
    def draw(self) -> None: ...

class _Container(_ContainerSpace):
    def drawOn(self, canv: Canvas, x: float, y: float) -> None: ...
    def copyContent(self, content: _FlowableSublist | None = None) -> None: ...

class PTOContainer(_Container, Flowable):
    def __init__(
        self, content: _FlowableSublist | None, trailer: _FlowableSublist | None = None, header: _FlowableSublist | None = None
    ) -> None: ...

class KeepInFrame(_Container, Flowable):
    name: str
    maxWidth: float
    maxHeight: float
    mode: Literal["error", "continue", "shrink", "truncate"]
    mergespace: Incomplete | None
    fakeWidth: bool | None
    def __init__(
        self,
        maxWidth: float,
        maxHeight: float,
        content: list[Flowable] = [],
        mergeSpace: Incomplete | None = 1,
        mode: Literal["error", "continue", "shrink", "truncate"] = "shrink",
        name: str = "",
        hAlign: str = "LEFT",
        vAlign: str = "BOTTOM",
        fakeWidth: bool | None = None,
    ) -> None: ...

class PlacedStory(Flowable):
    def __init__(
        self,
        x,
        y,
        maxWidth: float,
        maxHeight: float,
        content: list[Flowable] = [],
        mergeSpace: Incomplete | None = 1,
        mode: Literal["error", "continue", "shrink", "truncate"] = "shrink",
        name: str = "",
        anchor: str = "sw",
        fakeWidth: bool | None = None,
        hAlign: str = "LEFT",
        vAlign: str = "BOTTOM",
        showBoundary=None,
        origin="page",
    ) -> None: ...
    def wrap(self, _aW: Unused, _aH: Unused) -> tuple[Literal[0], Literal[0]]: ...
    def drawOn(self, canv: Canvas, lx: float, ly: float, _sW=0) -> None: ...

class _FindSplitterMixin: ...

class ImageAndFlowables(_Container, _FindSplitterMixin, Flowable):
    imageHref: str | None
    def __init__(
        self,
        I: Image,
        F: _FlowableSublist | None,
        imageLeftPadding: float = 0,
        imageRightPadding: float = 3,
        imageTopPadding: float = 0,
        imageBottomPadding: float = 3,
        imageSide: str = "right",
        imageHref: str | None = None,
    ) -> None: ...
    def deepcopy(self) -> Self: ...

class BalancedColumns(_FindSplitterMixin, NullDraw):
    name: str
    showBoundary: Incomplete | None
    endSlack: float
    def __init__(
        self,
        F: _FlowableSublist | None,
        nCols: int = 2,
        needed: float = 72,
        spaceBefore: float = 0,
        spaceAfter: float = 0,
        showBoundary=None,
        leftPadding: float | None = None,
        innerPadding: float | None = None,
        rightPadding: float | None = None,
        topPadding: float | None = None,
        bottomPadding: float | None = None,
        name: str = "",
        endSlack: float = 0.1,
        boxStrokeColor: Color | None = None,
        boxStrokeWidth: float = 0,
        boxFillColor: Color | None = None,
        boxMargin: tuple[int, int, int, int] | tuple[int, int, int] | tuple[int, int] | tuple[int] | None = None,
        vLinesStrokeColor: Color | None = None,
        vLinesStrokeWidth: float | None = None,
    ) -> None: ...

class AnchorFlowable(Spacer):
    def __init__(self, name: str) -> None: ...

class FrameBG(AnchorFlowable):
    start: bool
    left: float
    right: float
    color: Color
    strokeWidth: float
    strokeColor: Color
    strokeDashArray: list[float] | tuple[float, ...] | None
    def __init__(
        self,
        color: Color | None = None,
        left: float | str = 0,
        right: float | str = 0,
        start: bool = True,
        strokeWidth: float | None = None,
        strokeColor: Color | None = None,
        strokeDashArray: list[float] | tuple[float, ...] | None = None,
    ) -> None: ...

class FrameSplitter(NullDraw):
    nextTemplate: str
    nextFrames: list[str]
    gap: float
    required: float
    adjustHeight: bool
    def __init__(
        self,
        nextTemplate: str,
        nextFrames: list[str] | None = [],
        gap: float = 10,
        required: float = 72,
        adjustHeight: bool = True,
    ) -> None: ...

class BulletDrawer:
    value: str
    def __init__(
        self,
        value: str = "0",
        bulletAlign: str = "left",
        bulletType: str = "1",
        bulletColor: str = "black",
        bulletFontName: str = "Helvetica",
        bulletFontSize: int = 12,
        bulletOffsetY: int = 0,
        bulletDedent: int = 0,
        bulletDir: str = "ltr",
        bulletFormat=None,
    ) -> None: ...
    def drawOn(self, indenter: DDIndenter, canv: Canvas, x: float, y: float) -> None: ...

class DDIndenter(Flowable):
    def __init__(self, flowable: Flowable, leftIndent: float = 0, rightIndent: float = 0) -> None: ...

class LIIndenter(DDIndenter):
    def __init__(
        self,
        flowable: Flowable,
        leftIndent: float = 0,
        rightIndent: float = 0,
        bullet=None,
        spaceBefore: float | None = None,
        spaceAfter: float | None = None,
    ) -> None: ...

class ListItem:
    # NOTE: style has to be a ListStyle, but this will be annoying with sheet["ul"]
    # TODO: Use Unpack for kwds with the ListStyle properties + value/spaceBefore/spaceAfter
    def __init__(self, flowables: _FlowableSublist, style: PropertySet | None = None, **kwds) -> None: ...

class ListFlowable(_Container, Flowable, _FindSplitterMixin):
    style: ListStyle
    # NOTE: style has to be a ListStyle, but this will be annoying with sheet["ul"]
    # TODO: Use Unpack for kwds with the ListStyle properties + spaceBefore/spaceAfter
    def __init__(self, flowables: Iterable[_NestedFlowable], start=None, style: PropertySet | None = None, **kwds) -> None: ...

class TopPadder(Flowable):
    # NOTE: TopPadder is mostly a transparent wrapper, we may consider trying
    #       something using __new__ in the future
    def __init__(self, f: Flowable) -> None: ...
    def __setattr__(self, a: str, v: Any) -> None: ...
    def __getattr__(self, a: str) -> Any: ...
    def __delattr__(self, a: str) -> None: ...

class DocAssign(NullDraw):
    args: tuple[Any, ...]
    def __init__(self, var: str, expr: object, life: str = "forever") -> None: ...
    def funcWrap(self, aW: float, aH: float) -> None: ...
    def func(self) -> None: ...

class DocExec(DocAssign):
    def __init__(self, stmt: str, lifetime: str = "forever") -> None: ...

class DocPara(DocAssign):
    expr: object
    format: str | None
    style: PropertySet | None
    klass: Incomplete
    escape: bool
    def __init__(
        self,
        expr: object,
        format: str | None = None,
        style: PropertySet | None = None,
        klass: _StyledFlowableFactory | None = None,
        escape: bool = True,
    ) -> None: ...
    def funcWrap(self, aW: float, aH: float) -> Any: ...
    def func(self) -> Any: ...
    def add_content(self, *args: Flowable) -> None: ...
    def get_value(self, aW: float, aH: float) -> str: ...

class DocAssert(DocPara):
    def __init__(self, cond: object, format: str | None = None) -> None: ...

class DocIf(DocPara):
    blocks: tuple[_FlowableSublist, _FlowableSublist]
    def __init__(self, cond: object, thenBlock: _FlowableSublist, elseBlock: _FlowableSublist = []) -> None: ...
    def checkBlock(self, block: _FlowableSublist) -> list[Flowable] | tuple[Flowable, ...]: ...

class DocWhile(DocIf):
    block: _FlowableSublist
    def __init__(self, cond: object, whileBlock: _FlowableSublist) -> None: ...

class SetTopFlowables(NullDraw):
    def __init__(self, F: list[Flowable], show: bool = False) -> None: ...

class SetPageTopFlowables(NullDraw):
    def __init__(self, F: list[Flowable], show: bool = False) -> None: ...
