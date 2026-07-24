from _typeshed import Incomplete
from typing import Final

from reportlab.platypus import Flowable

__version__: Final[str]
captionStyle: Incomplete

class Figure(Flowable):
    width: Incomplete
    figureHeight: Incomplete
    caption: Incomplete
    captionFont: Incomplete
    captionSize: Incomplete
    captionTextColor: Incomplete
    captionBackColor: Incomplete
    captionGap: Incomplete
    captionAlign: Incomplete
    captionPosition: Incomplete
    captionHeight: int
    background: Incomplete
    border: Incomplete
    spaceBefore: Incomplete
    spaceAfter: Incomplete
    hAlign: Incomplete
    def __init__(
        self,
        width,
        height,
        caption: str = "",
        captionFont="Helvetica-Oblique",
        captionSize: int = 12,
        background=None,
        captionTextColor=...,
        captionBackColor=None,
        border=None,
        spaceBefore: int = 12,
        spaceAfter: int = 12,
        captionGap=None,
        captionAlign: str = "centre",
        captionPosition: str = "bottom",
        hAlign: str = "CENTER",
    ) -> None: ...
    height: Incomplete
    dx: Incomplete
    def wrap(self, availWidth, availHeight): ...
    def draw(self) -> None: ...
    def drawBorder(self) -> None: ...
    def drawBackground(self) -> None: ...
    def drawCaption(self) -> None: ...
    def drawFigure(self) -> None: ...

def drawPage(canvas, x, y, width, height) -> None: ...

class PageFigure(Figure):
    caption: str
    captionStyle: Incomplete
    background: Incomplete
    def __init__(self, background=None) -> None: ...
    def drawVirtualPage(self) -> None: ...
    def drawFigure(self) -> None: ...

class PlatPropFigure1(PageFigure):
    caption: str
    def __init__(self) -> None: ...
    def drawVirtualPage(self) -> None: ...

class FlexFigure(Figure):
    shrinkToFit: Incomplete
    growToFit: Incomplete
    scaleFactor: Incomplete
    background: Incomplete
    def __init__(
        self,
        width,
        height,
        caption,
        background=None,
        captionFont: str = "Helvetica-Oblique",
        captionSize: int = 8,
        captionTextColor=...,
        shrinkToFit: int = 1,
        growToFit: int = 1,
        spaceBefore: int = 12,
        spaceAfter: int = 12,
        captionGap: int = 9,
        captionAlign: str = "centre",
        captionPosition: str = "top",
        scaleFactor=None,
        hAlign: str = "CENTER",
        border: int = 1,
    ) -> None: ...
    def wrap(self, availWidth, availHeight): ...
    def split(self, availWidth, availHeight): ...

class ImageFigure(FlexFigure):
    filename: Incomplete
    def __init__(self, filename, caption, background=None, scaleFactor=None, hAlign: str = "CENTER", border=None) -> None: ...
    def drawFigure(self) -> None: ...

class DrawingFigure(FlexFigure):
    drawing: Incomplete
    growToFit: int
    def __init__(self, modulename, classname, caption, baseDir=None, background=None) -> None: ...
    def drawFigure(self) -> None: ...

def demo1(canvas) -> None: ...
def test1() -> None: ...
