from _typeshed import Incomplete

from reportlab.graphics.widgetbase import Widget
from reportlab.lib.attrmap import *
from reportlab.lib.validators import *

class SlideBox(Widget):
    labelFontName: str
    labelFontSize: int
    labelStrokeColor: Incomplete
    labelFillColor: Incomplete
    startColor: Incomplete
    endColor: Incomplete
    numberOfBoxes: int
    trianglePosition: int
    triangleHeight: Incomplete
    triangleWidth: Incomplete
    triangleFillColor: Incomplete
    triangleStrokeColor: Incomplete
    triangleStrokeWidth: float
    boxHeight: Incomplete
    boxWidth: Incomplete
    boxSpacing: Incomplete
    boxOutlineColor: Incomplete
    boxOutlineWidth: float
    leftPadding: int
    rightPadding: int
    topPadding: int
    bottomPadding: int
    background: Incomplete
    sourceLabelText: str
    sourceLabelOffset: Incomplete
    sourceLabelFontName: str
    sourceLabelFontSize: int
    sourceLabelFillColor: Incomplete
    def __init__(self) -> None: ...
    def demo(self, drawing=None): ...
    def draw(self): ...
