from _typeshed import Incomplete
from abc import abstractmethod
from collections.abc import Callable
from typing import IO, Any, Literal, Protocol, type_check_only
from typing_extensions import Self, TypeAlias

from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus.flowables import Flowable
from reportlab.platypus.frames import Frame

__all__ = (
    "ActionFlowable",
    "BaseDocTemplate",
    "CurrentFrameFlowable",
    "FrameActionFlowable",
    "FrameBreak",
    "Indenter",
    "IndexingFlowable",
    "LayoutError",
    "LCActionFlowable",
    "NextFrameFlowable",
    "NextPageTemplate",
    "NotAtTopPageBreak",
    "NullActionFlowable",
    "PageAccumulator",
    "PageBegin",
    "PageTemplate",
    "SimpleDocTemplate",
)

# NOTE: Since we don't know what kind of DocTemplate we will use in a PageTemplate
#       we'll leave the second argument at Any, since the workaround with unbound
#       type vars didn't seem to work for this one
_PageCallback: TypeAlias = Callable[[Canvas, Any], object]

@type_check_only
class _CanvasMaker(Protocol):
    # NOTE: This matches a subset of Canvas.__init__
    def __call__(
        self,
        filename: str | IO[bytes],
        /,
        *,
        pagesize=None,
        pageCompression=None,
        invariant=None,
        enforceColorSpace=None,
        initialFontName=None,
        initialFontSize=None,
        initialLeading=None,
        cropBox=None,
        artBox=None,
        trimBox=None,
        bleedBox=None,
        lang=None,
    ) -> Canvas: ...

class LayoutError(Exception): ...

# NOTE: While internal, this is used as sentinel value in PageTemplate
#       and SimpleDocTemplate so in subclasses you may need to use this
def _doNothing(canvas: Canvas, doc: BaseDocTemplate) -> None: ...

class PTCycle(list[PageTemplate]):
    @property
    def next_value(self) -> PageTemplate: ...
    @property
    def peek(self) -> PageTemplate: ...

class IndexingFlowable(Flowable):
    def isIndexing(self) -> Literal[1]: ...
    def isSatisfied(self) -> int: ...
    def notify(self, kind: str, stuff: Any) -> None: ...
    def beforeBuild(self) -> None: ...
    def afterBuild(self) -> None: ...

class ActionFlowable(Flowable):
    # NOTE: Technically action always has to contain a string referencing
    #       a handle_ method on the DocTemplate, while the rest are the args
    #       that should be passed to that method, but since the default arg
    #       on __init__ violates that we might as well keep things simple
    action: tuple[Any, ...]
    def __init__(self, action: list[Any] | tuple[Any, ...] = ()) -> None: ...
    def apply(self, doc: BaseDocTemplate) -> None: ...
    def __call__(self) -> Self: ...

class NullActionFlowable(ActionFlowable): ...

class LCActionFlowable(ActionFlowable):
    locChanger: int
    def draw(self) -> None: ...

class NextFrameFlowable(ActionFlowable):
    locChanger: int
    def __init__(self, ix: int | str, resume: int = 0) -> None: ...

class CurrentFrameFlowable(LCActionFlowable):
    def __init__(self, ix: int | str, resume: int = 0) -> None: ...

class _FrameBreak(LCActionFlowable):
    def __call__(self, ix: int | str | None = None, resume: int = 0) -> Self: ...
    def apply(self, doc: BaseDocTemplate) -> None: ...

FrameBreak: _FrameBreak
PageBegin: LCActionFlowable

class FrameActionFlowable(Flowable):
    @abstractmethod
    def __init__(self, *arg: Any, **kw: Any) -> None: ...
    @abstractmethod
    def frameAction(self, frame: Frame) -> None: ...

class Indenter(FrameActionFlowable):
    width: float
    height: float
    left: float
    right: float
    def __init__(self, left: float | str = 0, right: float | str = 0) -> None: ...
    def frameAction(self, frame: Frame) -> None: ...

class NotAtTopPageBreak(FrameActionFlowable):
    locChanger: int
    nextTemplate: Incomplete
    def __init__(self, nextTemplate=None) -> None: ...
    def frameAction(self, frame: Frame) -> None: ...

class NextPageTemplate(ActionFlowable):
    locChanger: int
    def __init__(self, pt: str | int | list[str] | tuple[str, ...]) -> None: ...

class PageTemplate:
    id: str | None
    frames: list[Frame]
    onPage: _PageCallback
    onPageEnd: _PageCallback
    pagesize: tuple[float, float]
    autoNextPageTemplate: Incomplete
    cropBox: Incomplete
    artBox: Incomplete
    trimBox: Incomplete
    bleedBox: Incomplete
    def __init__(
        self,
        id: str | None = None,
        frames: list[Frame] | Frame = [],
        onPage: _PageCallback = ...,
        onPageEnd: _PageCallback = ...,
        pagesize: tuple[float, float] | None = None,
        autoNextPageTemplate=None,
        cropBox=None,
        artBox=None,
        trimBox=None,
        bleedBox=None,
    ) -> None: ...
    def beforeDrawPage(self, canv: Canvas, doc: BaseDocTemplate) -> None: ...
    def checkPageSize(self, canv: Canvas, doc: BaseDocTemplate) -> None: ...
    def afterDrawPage(self, canv: Canvas, doc: BaseDocTemplate) -> None: ...

class onDrawStr(str):
    onDraw: Callable[[Canvas, str | None, str], object]
    kind: str | None
    label: str
    def __new__(
        cls, value: object, onDraw: Callable[[Canvas, str | None, str], object], label: str, kind: str | None = None
    ) -> Self: ...
    def __getnewargs__(self) -> tuple[str, Callable[[Canvas, str | None, str], object], str, str | None]: ...  # type: ignore[override]

class PageAccumulator:
    name: str
    data: list[tuple[Any, ...]]
    def __init__(self, name: str | None = None) -> None: ...
    def reset(self) -> None: ...
    def add(self, *args) -> None: ...
    def onDrawText(self, *args) -> str: ...
    def __call__(self, canv: Canvas, kind: str | None, label: str) -> None: ...
    def attachToPageTemplate(self, pt: PageTemplate) -> None: ...
    def onPage(self, canv: Canvas, doc: BaseDocTemplate) -> None: ...
    def onPageEnd(self, canv: Canvas, doc: BaseDocTemplate) -> None: ...
    def pageEndAction(self, canv: Canvas, doc: BaseDocTemplate) -> None: ...
    def onDrawStr(self, value: object, *args) -> onDrawStr: ...

class BaseDocTemplate:
    filename: Incomplete
    pagesize: Incomplete
    pageTemplates: list[PageTemplate]
    showBoundary: Incomplete
    width: float
    height: float
    leftMargin: float
    rightMargin: float
    topMargin: float
    bottomMargin: float
    allowSplitting: Incomplete
    title: Incomplete | None
    author: Incomplete | None
    subject: Incomplete | None
    creator: Incomplete | None
    producer: Incomplete | None
    keywords: list[Incomplete]
    invariant: Incomplete | None
    pageCompression: Incomplete | None
    rotation: Incomplete
    encrypt: Incomplete | None
    cropMarks: Incomplete | None
    enforceColorSpace: Incomplete | None
    displayDocTitle: Incomplete | None
    lang: Incomplete | None
    initialFontName: Incomplete | None
    initialFontSize: Incomplete | None
    initialLeading: Incomplete | None
    cropBox: Incomplete | None
    artBox: Incomplete | None
    trimBox: Incomplete | None
    bleedBox: Incomplete | None
    keepTogetherClass: type[Flowable]
    hideToolbar: Incomplete | None
    hideMenubar: Incomplete | None
    hideWindowUI: Incomplete | None
    fitWindow: Incomplete | None
    centerWindow: Incomplete | None
    nonFullScreenPageMode: Incomplete | None
    direction: Incomplete | None
    viewArea: Incomplete | None
    viewClip: Incomplete | None
    printArea: Incomplete | None
    printClip: Incomplete | None
    printScaling: Incomplete | None
    duplex: Incomplete | None
    # NOTE: The following attributes only exist while/after pages are rendered
    pageTemplate: PageTemplate
    page: int
    frame: Frame
    canv: Canvas
    # TODO: Use TypedDict with Unpack for **kw
    def __init__(self, filename: str | IO[bytes], **kw) -> None: ...
    def setPageCallBack(self, func: Callable[[int], object] | None) -> None: ...
    def setProgressCallBack(self, func: Callable[[str, int], object] | None) -> None: ...
    def clean_hanging(self) -> None: ...
    def addPageTemplates(self, pageTemplates: list[PageTemplate] | tuple[PageTemplate, ...] | PageTemplate) -> None: ...
    def handle_documentBegin(self) -> None: ...
    def handle_pageBegin(self) -> None: ...
    def handle_pageEnd(self) -> None: ...
    def handle_pageBreak(self, slow: bool | None = None) -> None: ...
    def handle_frameBegin(self, resume: int = 0, pageTopFlowables=None) -> None: ...
    def handle_frameEnd(self, resume: int = 0) -> None: ...
    def handle_nextPageTemplate(self, pt: str | int | list[str] | tuple[str, ...]) -> None: ...
    def handle_nextFrame(self, fx: str | int, resume: int = 0) -> None: ...
    def handle_currentFrame(self, fx: str | int, resume: int = 0) -> None: ...
    def handle_breakBefore(self, flowables: list[Flowable]) -> None: ...
    def handle_keepWithNext(self, flowables: list[Flowable]) -> None: ...
    def handle_flowable(self, flowables: list[Flowable]) -> None: ...
    def build(
        self, flowables: list[Flowable], filename: str | IO[bytes] | None = None, canvasmaker: _CanvasMaker = ...
    ) -> None: ...
    def notify(self, kind: str, stuff: Any) -> None: ...
    def pageRef(self, label: str) -> None: ...
    def multiBuild(self, story: list[Flowable], maxPasses: int = 10, **buildKwds: Any) -> int: ...
    def afterInit(self) -> None: ...
    def beforeDocument(self) -> None: ...
    def beforePage(self) -> None: ...
    def afterPage(self) -> None: ...
    def filterFlowables(self, flowables: list[Flowable]) -> None: ...
    def afterFlowable(self, flowable: Flowable) -> None: ...
    def docAssign(self, var: str, expr: object, lifetime: str) -> None: ...
    def docExec(self, stmt: str, lifetime: str) -> None: ...
    def docEval(self, expr: str) -> Any: ...

class SimpleDocTemplate(BaseDocTemplate):
    def handle_pageBegin(self) -> None: ...
    def build(  # type: ignore[override]
        self,
        flowables: list[Flowable],
        onFirstPage: _PageCallback = ...,
        onLaterPages: _PageCallback = ...,
        canvasmaker: _CanvasMaker = ...,
    ) -> None: ...
