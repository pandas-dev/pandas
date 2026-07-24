from _typeshed import Incomplete
from typing import Final, NoReturn

from reportlab.graphics.shapes import LineShape
from reportlab.graphics.widgetbase import Widget

__version__: Final[str]

def frange(start, end=None, inc=None): ...
def makeDistancesList(list): ...

class Grid(Widget):
    x: int
    y: int
    width: int
    height: int
    orientation: str
    useLines: int
    useRects: int
    delta: int
    delta0: int
    deltaSteps: Incomplete
    fillColor: Incomplete
    stripeColors: Incomplete
    strokeColor: Incomplete
    strokeWidth: int
    def __init__(self) -> None: ...
    def demo(self): ...
    def makeOuterRect(self): ...
    def makeLinePosList(self, start, isX: int = 0): ...
    def makeInnerLines(self): ...
    def makeInnerTiles(self): ...
    def draw(self): ...

class DoubleGrid(Widget):
    x: int
    y: int
    width: int
    height: int
    grid0: Incomplete
    grid1: Incomplete
    def __init__(self) -> None: ...
    def demo(self): ...
    def draw(self): ...

class ShadedRect(Widget):
    x: int
    y: int
    width: int
    height: int
    orientation: str
    numShades: int
    fillColorStart: Incomplete
    fillColorEnd: Incomplete
    strokeColor: Incomplete
    strokeWidth: int
    cylinderMode: int
    def __init__(self, **kw) -> None: ...
    def demo(self): ...
    def draw(self): ...

def colorRange(c0, c1, n): ...
def centroid(P): ...
def rotatedEnclosingRect(P, angle, rect): ...

class ShadedPolygon(Widget, LineShape):
    angle: int
    fillColorStart: Incomplete
    fillColorEnd: Incomplete
    cylinderMode: int
    numShades: int
    points: Incomplete
    def __init__(self, **kw) -> None: ...
    def draw(self): ...
    # NOTE: widgets don't implement this, only actual shapes
    def copy(self) -> NoReturn: ...
