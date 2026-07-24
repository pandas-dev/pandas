from _typeshed import Incomplete, SupportsItems
from abc import abstractmethod
from collections.abc import Iterable, Sequence
from typing import Any, Final, Literal, NoReturn, TypedDict, type_check_only
from typing_extensions import Self, TypeAlias, Unpack

from reportlab.lib.colors import Color
from reportlab.lib.validators import NoneOr, Validator
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import Flowable
from reportlab.platypus.flowables import _HAlignment, _VAlignment

_IntBool: TypeAlias = Literal[0, 1]
_BoolLike: TypeAlias = _IntBool | bool
_PathOp: TypeAlias = (
    tuple[Literal["moveTo"], float, float]
    | tuple[Literal["lineTo"], float, float]
    | tuple[Literal["curveTo"], float, float, float, float, float, float]
    # close path may either be a tuple or just the string
    | Literal["closePath"]
    | tuple[Literal["closePath"]]
    # fallback for list that is not type safe
    | list[Any]
)

# NOTE: These are derived from _attrMap and can optionally be
#       verified at runtime
@type_check_only
class _GroupKwArgs(TypedDict, total=False):
    transform: tuple[float, float, float, float, float, float] | list[float] | list[int]
    # NOTE: This should be used with care, since it will replace elements
    #       it's mostly useful for circumventing validation logic and
    #       reusing the list, rather than populating a new list
    contents: list[Shape]
    strokeOverprint: _BoolLike
    fillOverprint: _BoolLike
    overprintMask: _BoolLike

@type_check_only
class _DrawingKwArgs(_GroupKwArgs, total=False):
    # TODO: Restrict to supported formats?
    formats: list[str] | tuple[str, ...]
    # NOTE: This looks like an implementation detail, so we may not
    #       want to include this in KwArgs
    canv: Canvas
    background: Shape | UserNode | None
    # NOTE: The runtime validation for alignments is incorrect, so
    #       we assume it is turned off and allow all valid values
    hAlign: _HAlignment
    vAlign: _VAlignment
    renderScale: float
    initialFontName: str | None
    initialFontSize: float | None

@type_check_only
class _LineShapeKwArgs(TypedDict, total=False):
    strokeColor: Color | None
    strokeWidth: float
    strokeLineCap: Literal[0, 1, 2]
    strokeLineJoin: Literal[0, 1, 2]
    strokeMiterLimit: float
    strokeDashArray: Sequence[float] | tuple[float, Sequence[float]]
    strokeOpacity: float | None
    strokeOverprint: _BoolLike
    overprintMask: _BoolLike

@type_check_only
class _PathKwArgs(_LineShapeKwArgs, total=False):
    fillColor: Color | None
    fillOpacity: float
    fillOverprint: _BoolLike

@type_check_only
class _AllPathKwArgs(_PathKwArgs, total=False):
    points: list[float] | None
    operators: list[float] | None
    isClipPath: _BoolLike
    autoclose: Literal["svg", "pdf"] | None
    fillMode: Literal[0, 1]

@type_check_only
class _SolidShapeKwArgs(_PathKwArgs, total=False):
    fillMode: Literal[0, 1]

@type_check_only
class _DefinePathKwArgs(_SolidShapeKwArgs, total=False):
    autoclose: Literal["svg", "pdf"] | None
    bbox: tuple[float, float, float, float] | None

@type_check_only
class _WedgeKwArgs(_SolidShapeKwArgs, total=False):
    radius1: float | None
    yradius1: float | None

@type_check_only
class _StringKwArgs(TypedDict, total=False):
    fontName: str
    fontSize: float
    fillColor: Color | None
    textAnchor: Literal["start", "middle", "end", "numeric"]
    encoding: str
    textRenderMode: Literal[0, 1, 2, 3, 4, 5, 6, 7]

__version__: Final[str]
isOpacity: NoneOr
NON_ZERO_WINDING: Final[str]
EVEN_ODD: Final[str]
STATE_DEFAULTS: Final[Incomplete]

class _DrawTimeResizeable: ...

class _SetKeyWordArgs:
    def __init__(self, keywords: SupportsItems[str, Any] = {}) -> None: ...

def getRectsBounds(rectList): ...
def getPathBounds(points): ...
def getPointsBounds(pointList): ...

class Shape(_SetKeyWordArgs, _DrawTimeResizeable):
    @abstractmethod
    def copy(self) -> Self: ...
    def getProperties(self, recur: int = 1) -> dict[str, Any]: ...
    def setProperties(self, props) -> None: ...
    def dumpProperties(self, prefix: str = "") -> None: ...
    def verify(self) -> None: ...
    @abstractmethod
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Group(Shape):
    contents: list[Shape]
    transform: tuple[float, float, float, float, float, float] | list[float] | list[int]
    def __init__(self, *elements: Shape | UserNode, **keywords: Unpack[_GroupKwArgs]) -> None: ...
    def add(self, node: Shape | UserNode, name: str | None = None) -> None: ...
    def insert(self, i: int, n: Shape | UserNode, name: str | None = None) -> None: ...
    def expandUserNodes(self) -> Group: ...
    def copy(self) -> Self: ...
    def rotate(self, theta: float, cx: float = 0, cy: float = 0) -> None: ...
    def translate(self, dx: float, dy: float = 0) -> None: ...
    def scale(self, sx: float, sy: float = 1) -> None: ...
    def skew(self, kx: float, ky: float = 0) -> None: ...
    def shift(self, x: float, y: float = 0) -> None: ...
    # NOTE: This changes the object to a Drawing, rather than returning
    #       a new one, which is not ideal...
    def asDrawing(self, width: float, height: float) -> None: ...
    def getContents(self) -> list[Shape | UserNode]: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Drawing(Group, Flowable):
    background: Shape | UserNode | None
    renderScale: float
    def __init__(
        self, width: float = 400, height: float = 200, *nodes: Shape | UserNode, **keywords: Unpack[_DrawingKwArgs]
    ) -> None: ...
    def draw(self, showBoundary=...) -> None: ...
    def expandUserNodes(self) -> Drawing: ...
    def asGroup(self, *args: Shape | UserNode, **kw: Unpack[_GroupKwArgs]) -> Group: ...
    def save(
        self,
        formats: Iterable[str] | None = None,
        verbose: bool | None = None,
        fnRoot: str | None = None,
        outDir: str | None = None,
        title: str = "",
        **kw,
    ): ...
    def asString(self, format: str, verbose: bool | None = None, preview: int = 0, **kw) -> str: ...
    def resized(
        self, kind: Literal["fit", "fitx", "fity"] = "fit", lpad: float = 0, rpad: float = 0, bpad: float = 0, tpad: float = 0
    ) -> Drawing: ...

class _DrawingEditorMixin: ...

@type_check_only
class _isStrokeDashArray(Validator):
    def test(self, x): ...

isStrokeDashArray: _isStrokeDashArray

class LineShape(Shape):
    strokeColor: Color | None
    strokeWidth: float
    strokeLineCap: Literal[0, 1, 2]
    strokeLineJoin: Literal[0, 1, 2]
    strokeMiterLimit: float
    strokeDashArray: Sequence[float] | tuple[float, Sequence[float]]
    strokeOpacity: float | None
    def __init__(self, kw: _LineShapeKwArgs) -> None: ...
    @abstractmethod
    def copy(self) -> Self: ...
    @abstractmethod
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Line(LineShape):
    x1: float
    y1: float
    x2: float
    y2: float
    def __init__(self, x1: float, y1: float, x2: float, y2: float, **kw: Unpack[_LineShapeKwArgs]) -> None: ...
    # NOTE: For some reason Line doesn't implement copy
    def copy(self) -> NoReturn: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class SolidShape(LineShape):
    fillColor: Color | None
    fillOpacity: float | None
    def __init__(self, kw: _SolidShapeKwArgs) -> None: ...
    @abstractmethod
    def copy(self) -> Self: ...
    @abstractmethod
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Path(SolidShape):
    points: list[float]
    operators: list[float]
    isClipPath: _BoolLike
    autoclose: Literal["svg", "pdf"] | None
    fillMode: Literal[0, 1]
    def __init__(
        self,
        points: list[float] | None = None,
        operators: list[float] | None = None,
        isClipPath: _BoolLike = 0,
        autoclose: Literal["svg", "pdf"] | None = None,
        fillMode: Literal[0, 1] = 0,
        **kw: Unpack[_PathKwArgs],
    ) -> None: ...
    def copy(self) -> Self: ...
    def moveTo(self, x: float, y: float) -> None: ...
    def lineTo(self, x: float, y: float) -> None: ...
    def curveTo(self, x1: float, y1: float, x2: float, y2: float, x3: float, y3: float) -> None: ...
    def closePath(self) -> None: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

EmptyClipPath: Final[Path]

def getArcPoints(
    centerx: float,
    centery: float,
    radius: float,
    startangledegrees: float,
    endangledegrees: float,
    yradius: float | None = None,
    degreedelta: float | None = None,
    reverse: _BoolLike | None = None,
) -> list[float]: ...

class ArcPath(Path):
    def addArc(
        self,
        centerx: float,
        centery: float,
        radius: float,
        startangledegrees: float,
        endangledegrees: float,
        yradius: float | None = None,
        degreedelta: float | None = None,
        moveTo: _BoolLike | None = None,
        reverse: _BoolLike | None = None,
    ) -> None: ...

def definePath(
    pathSegs: Iterable[_PathOp] = [], isClipPath: _BoolLike = 0, dx: float = 0, dy: float = 0, **kw: Unpack[_DefinePathKwArgs]
) -> Path: ...

class Rect(SolidShape):
    x: float
    y: float
    width: float
    height: float
    rx: float
    ry: float
    def __init__(
        self, x: float, y: float, width: float, height: float, rx: float = 0, ry: float = 0, **kw: _SolidShapeKwArgs
    ) -> None: ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Image(SolidShape):
    x: float
    y: float
    width: float
    height: float
    path: Incomplete
    def __init__(self, x: float, y: float, width: float, height: float, path, **kw: Unpack[_SolidShapeKwArgs]) -> None: ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Circle(SolidShape):
    cx: float
    cy: float
    r: float
    def __init__(self, cx: float, cy: float, r: float, **kw: Unpack[_SolidShapeKwArgs]) -> None: ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Ellipse(SolidShape):
    cx: float
    cy: float
    rx: float
    ry: float
    def __init__(self, cx: float, cy: float, rx: float, ry: float, **kw: Unpack[_SolidShapeKwArgs]) -> None: ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Wedge(SolidShape):
    centerx: float
    centery: float
    radius: float
    startangledegrees: float
    endangledegrees: float
    yradius: float | None
    annular: bool
    # NOTE: This one is not actually settable on the instance if runtime validation
    #       is turned on, but it seems bad to disallow it anyways
    degreedelta: float
    def __init__(
        self,
        centerx: float,
        centery: float,
        radius: float,
        startangledegrees: float,
        endangledegrees: float,
        yradius: float | None = None,
        annular: bool = False,
        **kw: Unpack[_WedgeKwArgs],
    ) -> None: ...
    def asPolygon(self) -> Path | Polygon: ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Polygon(SolidShape):
    points: list[float]
    def __init__(self, points: list[float] = [], **kw: Unpack[_SolidShapeKwArgs]) -> None: ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class PolyLine(LineShape):
    points: list[float]
    def __init__(self, points: list[float] = [], **kw: Unpack[_SolidShapeKwArgs]) -> None: ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class Hatching(Path):
    xyLists: Sequence[tuple[float, float]]
    angles: Sequence[float]
    spacings: Sequence[float]
    def __init__(
        self,
        spacings: float | Sequence[float] = 2,
        angles: float | Sequence[float] = 45,
        xyLists: Sequence[tuple[float, float] | list[float]] = [],
        **kwds: Unpack[_AllPathKwArgs],
    ) -> None: ...

def numericXShift(
    tA, text: str, w: float, fontName: str, fontSize: float, encoding: str | None = None, pivotCharacter: str = "."
) -> float: ...

class String(Shape):
    encoding: str
    x: float
    y: float
    text: str
    textAnchor: Literal["start", "middle", "end", "numeric"]
    fontName: str
    fontSize: float
    fillColor: Color | None
    def __init__(self, x: float, y: float, text: str, **kw: Unpack[_StringKwArgs]) -> None: ...
    def getEast(self): ...
    def copy(self) -> Self: ...
    def getBounds(self) -> tuple[float, float, float, float]: ...

class UserNode(_DrawTimeResizeable):
    @abstractmethod
    def provideNode(self) -> Shape: ...

class DirectDraw(Shape):
    @abstractmethod
    def drawDirectly(self, canvas: Canvas) -> None: ...

def test() -> None: ...
