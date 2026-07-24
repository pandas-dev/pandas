from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.shapes import Group
from reportlab.graphics.widgetbase import Widget
from reportlab.platypus import Flowable

__version__: Final[str]
adobe2codec: dict[str, str]

class CodeChartBase(Flowable):
    rows: Incomplete
    width: Incomplete
    height: Incomplete
    ylist: Incomplete
    xlist: Incomplete
    def calcLayout(self) -> None: ...
    def formatByte(self, byt) -> str: ...
    def drawChars(self, charList) -> None: ...
    def drawLabels(self, topLeft: str = "") -> None: ...

class SingleByteEncodingChart(CodeChartBase):
    codePoints: int
    faceName: Incomplete
    encodingName: Incomplete
    fontName: Incomplete
    charsPerRow: Incomplete
    boxSize: Incomplete
    hex: Incomplete
    rowLabels: Incomplete
    def __init__(
        self,
        faceName: str = "Helvetica",
        encodingName: str = "WinAnsiEncoding",
        charsPerRow: int = 16,
        boxSize: int = 14,
        hex: int = 1,
    ) -> None: ...
    def draw(self) -> None: ...

class KutenRowCodeChart(CodeChartBase):
    row: Incomplete
    codePoints: int
    boxSize: int
    charsPerRow: int
    rows: int
    rowLabels: Incomplete
    hex: int
    faceName: Incomplete
    encodingName: Incomplete
    fontName: Incomplete
    def __init__(self, row, faceName, encodingName) -> None: ...
    def makeRow(self, row) -> list[bytes | list[None]]: ...
    def draw(self) -> None: ...

class Big5CodeChart(CodeChartBase):
    row: Incomplete
    codePoints: int
    boxSize: int
    charsPerRow: int
    rows: int
    hex: int
    faceName: Incomplete
    encodingName: Incomplete
    rowLabels: Incomplete
    fontName: Incomplete
    def __init__(self, row, faceName, encodingName) -> None: ...
    def makeRow(self, row) -> list[bytes | list[None]]: ...
    def draw(self) -> None: ...

def hBoxText(msg, canvas, x, y, fontName) -> None: ...

class CodeWidget(Widget):
    x: int
    y: int
    width: int
    height: int
    def __init__(self) -> None: ...
    def draw(self) -> Group: ...

def test() -> None: ...
