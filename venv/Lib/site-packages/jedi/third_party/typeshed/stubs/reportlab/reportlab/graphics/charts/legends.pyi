from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.widgetbase import PropHolder, Widget
from reportlab.lib.attrmap import *

__version__: Final[str]

class SubColProperty(PropHolder):
    dividerLines: int

class LegendCallout:
    def __call__(self, legend, g, thisx, y, colName) -> None: ...

class LegendSwatchCallout(LegendCallout):
    def __call__(self, legend, g, thisx, y, i, colName, swatch) -> None: ...  # type: ignore[override]

class LegendColEndCallout(LegendCallout):
    def __call__(self, legend, g, x, xt, y, width, lWidth) -> None: ...  # type: ignore[override]

class Legend(Widget):
    x: int
    y: int
    alignment: str
    deltax: int
    deltay: int
    autoXPadding: int
    autoYPadding: int
    dx: int
    dy: int
    swdx: int
    swdy: int
    dxTextSpace: int
    columnMaximum: int
    colorNamePairs: Incomplete
    fontName: Incomplete
    fontSize: Incomplete
    leading: Incomplete
    fillColor: Incomplete
    strokeColor: Incomplete
    strokeWidth: Incomplete
    swatchMarker: Incomplete
    boxAnchor: str
    yGap: int
    variColumn: int
    dividerLines: int
    dividerWidth: float
    dividerDashArray: Incomplete
    dividerColor: Incomplete
    dividerOffsX: Incomplete
    dividerOffsY: int
    colEndCallout: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...
    def demo(self): ...

class TotalAnnotator(LegendColEndCallout):
    lText: Incomplete
    rText: Incomplete
    fontName: Incomplete
    fontSize: Incomplete
    fillColor: Incomplete
    dy: Incomplete
    dx: Incomplete
    dly: Incomplete
    dlx: Incomplete
    strokeWidth: Incomplete
    strokeColor: Incomplete
    strokeDashArray: Incomplete
    def __init__(
        self,
        lText: str = "Total",
        rText: str = "0.0",
        fontName="Times-Roman",
        fontSize: int = 10,
        fillColor=...,
        strokeWidth: float = 0.5,
        strokeColor=...,
        strokeDashArray=None,
        dx: int = 0,
        dy: int = 0,
        dly: int = 0,
        dlx=(0, 0),
    ) -> None: ...
    def __call__(self, legend, g, x, xt, y, width, lWidth) -> None: ...  # type: ignore[override]

class LineSwatch(Widget):
    x: int
    y: int
    width: int
    height: int
    strokeColor: Incomplete
    strokeDashArray: Incomplete
    strokeWidth: int
    def __init__(self) -> None: ...
    def draw(self): ...

class LineLegend(Legend):
    dx: int
    dy: int
    def __init__(self) -> None: ...
