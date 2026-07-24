from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.widgetbase import Widget
from reportlab.lib.attrmap import *
from reportlab.lib.validators import *

__version__: Final[str]

class TableWidget(Widget):
    x: Incomplete
    y: Incomplete
    width: int
    height: int
    borderStrokeColor: Incomplete
    fillColor: Incomplete
    borderStrokeWidth: float
    horizontalDividerStrokeColor: Incomplete
    verticalDividerStrokeColor: Incomplete
    horizontalDividerStrokeWidth: float
    verticalDividerStrokeWidth: float
    dividerDashArray: Incomplete
    data: Incomplete
    boxAnchor: str
    fontSize: int
    fontColor: Incomplete
    alignment: str
    textAnchor: str
    def __init__(self, x: int = 10, y: int = 10, **kw) -> None: ...
    def demo(self): ...
    def draw(self): ...
    def preProcessData(self, data): ...
