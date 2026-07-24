from typing import Literal

from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus.flowables import Flowable

class Frame:
    id: str | None
    x1: float
    y1: float
    width: float
    height: float
    leftPadding: float
    bottomPadding: float
    rightPadding: float
    topPadding: float
    showBoundary: int
    def __init__(
        self,
        x1: float,
        y1: float,
        width: float,
        height: float,
        leftPadding: float = 6,
        bottomPadding: float = 6,
        rightPadding: float = 6,
        topPadding: float = 6,
        id: str | None = None,
        showBoundary: int = 0,
        overlapAttachedSpace=None,
        _debug=None,
    ) -> None: ...
    def add(self, flowable: Flowable, canv: Canvas, trySplit: int = 0) -> Literal[0, 1]: ...
    def split(self, flowable: Flowable, canv: Canvas) -> list[Flowable]: ...
    def drawBoundary(self, canv: Canvas) -> None: ...
    def addFromList(self, drawlist: list[Flowable], canv: Canvas) -> None: ...
    def add_generated_content(self, *C: Flowable) -> None: ...

__all__ = ("Frame",)
