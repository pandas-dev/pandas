from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.charts.piecharts import AbstractPieChart, WedgeProperties
from reportlab.lib.attrmap import *

__version__: Final[str]

class SectorProperties(WedgeProperties): ...

class Doughnut(AbstractPieChart):
    x: int
    y: int
    width: int
    height: int
    data: Incomplete
    labels: Incomplete
    startAngle: int
    direction: str
    simpleLabels: int
    checkLabelOverlap: int
    sideLabels: int
    innerRadiusFraction: Incomplete
    slices: Incomplete
    angleRange: int
    def __init__(self, *, angleRange: int = 360, **kwds) -> None: ...
    def demo(self): ...
    def normalizeData(self, data=None): ...
    def makeSectors(self): ...
    def draw(self): ...

def sample1(): ...
def sample2(): ...
def sample3(): ...
def sample4(): ...
