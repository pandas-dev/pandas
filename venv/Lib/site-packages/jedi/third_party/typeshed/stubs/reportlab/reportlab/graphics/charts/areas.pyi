from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.widgetbase import Widget

__version__: Final[str]

class PlotArea(Widget):
    x: int
    y: int
    height: int
    width: int
    strokeColor: Incomplete
    strokeWidth: int
    fillColor: Incomplete
    background: Incomplete
    debug: int
    def __init__(self) -> None: ...
    def makeBackground(self): ...
