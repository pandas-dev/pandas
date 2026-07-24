from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.charts.utils import CustomDrawChanger
from reportlab.graphics.shapes import Drawing, Group
from reportlab.graphics.widgetbase import PropHolder, Widget
from reportlab.lib.attrmap import *

__version__: Final[str]

class Label(Widget):
    # TODO: This has more attributes.
    x: Incomplete
    y: Incomplete
    def __init__(self, **kw) -> None: ...
    @property
    def padding(self): ...
    @padding.setter
    def padding(self, p) -> None: ...
    def setText(self, text) -> None: ...
    def setOrigin(self, x, y) -> None: ...
    def demo(self) -> Drawing: ...
    def computeSize(self) -> None: ...
    def draw(self) -> Group: ...

class LabelDecorator:
    textAnchor: str
    boxAnchor: str
    def __init__(self) -> None: ...
    def decorate(self, l, L) -> None: ...
    def __call__(self, l) -> None: ...

isOffsetMode: Incomplete

class LabelOffset(PropHolder):
    posMode: str
    pos: int
    def __init__(self) -> None: ...

NoneOrInstanceOfLabelOffset: Incomplete

class PMVLabel(Label):
    def __init__(self, **kwds) -> None: ...

class BarChartLabel(PMVLabel):
    lineStrokeWidth: int
    lineStrokeColor: Incomplete
    fixedStart: Incomplete
    nudge: int
    def __init__(self, **kwds) -> None: ...

class NA_Label(BarChartLabel):
    text: str
    def __init__(self) -> None: ...

NoneOrInstanceOfNA_Label: Incomplete

class RedNegativeChanger(CustomDrawChanger):
    fillColor: Incomplete
    def __init__(self, fillColor=...) -> None: ...

class XLabel(Label):
    ddfKlass: Incomplete
    ddf: Incomplete
    def __init__(self, *args, **kwds) -> None: ...
    def computeSize(self) -> None: ...
