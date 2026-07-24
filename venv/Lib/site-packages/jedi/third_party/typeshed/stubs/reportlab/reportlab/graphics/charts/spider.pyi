from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.charts.areas import PlotArea
from reportlab.graphics.charts.piecharts import WedgeLabel
from reportlab.graphics.widgetbase import PropHolder
from reportlab.lib.attrmap import *

__version__: Final[str]

class StrandProperty(PropHolder):
    strokeWidth: int
    fillColor: Incomplete
    strokeColor: Incomplete
    strokeDashArray: Incomplete
    symbol: Incomplete
    symbolSize: int
    name: Incomplete
    def __init__(self) -> None: ...

class SpokeProperty(PropHolder):
    strokeWidth: float
    fillColor: Incomplete
    strokeColor: Incomplete
    strokeDashArray: Incomplete
    visible: int
    labelRadius: float
    def __init__(self, **kw) -> None: ...

class SpokeLabel(WedgeLabel):
    def __init__(self, **kw) -> None: ...

class StrandLabel(SpokeLabel):
    format: str
    dR: int
    def __init__(self, **kw) -> None: ...

class SpiderChart(PlotArea):
    def makeSwatchSample(self, rowNo, x, y, width, height): ...
    def getSeriesName(self, i, default=None): ...
    data: Incomplete
    labels: Incomplete
    startAngle: int
    direction: str
    strands: Incomplete
    spokes: Incomplete
    spokeLabels: Incomplete
    strandLabels: Incomplete
    x: int
    y: int
    width: int
    height: int
    def __init__(self) -> None: ...
    def demo(self): ...
    def normalizeData(self, outer: float = 0.0): ...
    def labelClass(self, kind): ...
    def draw(self): ...

def sample1(): ...
def sample2(): ...
