from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.charts.linecharts import AbstractLineChart
from reportlab.graphics.charts.utils import *
from reportlab.graphics.shapes import Polygon, _SetKeyWordArgs
from reportlab.graphics.widgetbase import PropHolder
from reportlab.graphics.widgets.grids import ShadedPolygon
from reportlab.lib.attrmap import *
from reportlab.lib.validators import *

__version__: Final[str]

class LinePlotProperties(PropHolder): ...

class InFillValue(int):
    yValue: Incomplete
    def __new__(cls, v, yValue=None): ...

class Shader(_SetKeyWordArgs):
    def shade(self, lp, g, rowNo, rowColor, row) -> None: ...

class NoFiller:
    def fill(self, lp, g, rowNo, rowColor, points) -> None: ...

class Filler:
    __dict__: Incomplete
    def __init__(self, **kw) -> None: ...
    def fill(self, lp, g, rowNo, rowColor, points) -> None: ...

class ShadedPolyFiller(Filler, ShadedPolygon): ...
class PolyFiller(Filler, Polygon): ...

class LinePlot(AbstractLineChart):
    reversePlotOrder: int
    xValueAxis: Incomplete
    yValueAxis: Incomplete
    data: Incomplete
    lines: Incomplete
    lineLabels: Incomplete
    lineLabelFormat: Incomplete
    lineLabelArray: Incomplete
    lineLabelNudge: int
    annotations: Incomplete
    behindAxes: int
    gridFirst: int
    def __init__(self) -> None: ...
    @property
    def joinedLines(self): ...
    @joinedLines.setter
    def joinedLines(self, v) -> None: ...
    def demo(self): ...
    def calcPositions(self) -> None: ...
    def drawLabel(self, G, rowNo, colNo, x, y) -> None: ...
    def makeLines(self): ...
    def draw(self): ...
    def addCrossHair(self, name, xv, yv, strokeColor=..., strokeWidth: int = 1, beforeLines: bool = True): ...

class LinePlot3D(LinePlot):
    theta_x: float
    theta_y: float
    zDepth: int
    zSpace: int
    def calcPositions(self) -> None: ...
    def makeLines(self): ...

class SimpleTimeSeriesPlot(LinePlot):
    xValueAxis: Incomplete
    yValueAxis: Incomplete
    data: Incomplete
    def __init__(self) -> None: ...

class GridLinePlot(SimpleTimeSeriesPlot):
    scaleFactor: Incomplete
    background: Incomplete
    def __init__(self) -> None: ...
    def demo(self, drawing=None): ...
    def draw(self): ...

class AreaLinePlot(LinePlot):
    reversePlotOrder: int
    data: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class SplitLinePlot(AreaLinePlot):
    xValueAxis: Incomplete
    yValueAxis: Incomplete
    data: Incomplete
    def __init__(self) -> None: ...

class ScatterPlot(LinePlot):
    width: int
    height: int
    outerBorderOn: int
    outerBorderColor: Incomplete
    background: Incomplete
    xLabel: str
    yLabel: str
    data: Incomplete
    joinedLines: int
    leftPadding: int
    rightPadding: int
    topPadding: int
    bottomPadding: int
    x: Incomplete
    y: Incomplete
    lineLabelFormat: str
    lineLabelNudge: int
    def __init__(self) -> None: ...
    def demo(self, drawing=None): ...
    def draw(self): ...

def sample1a(): ...
def sample1b(): ...
def sample1c(): ...
def preprocessData(series): ...
def sample2(): ...
def sampleFillPairedData(): ...
