from _typeshed import Incomplete
from typing import type_check_only

from reportlab.graphics.shapes import Rect
from reportlab.graphics.widgetbase import Widget
from reportlab.lib.validators import Validator
from reportlab.platypus.flowables import Flowable

__all__ = ["QrCodeWidget"]

@type_check_only
class _isLevel(Validator):
    def test(self, x): ...

isLevel: _isLevel

@type_check_only
class _isUnicodeOrQRList(Validator):
    def test(self, x): ...
    def normalize(self, x): ...

isUnicodeOrQRList: _isUnicodeOrQRList

class SRect(Rect):
    def __init__(self, x, y, width, height, fillColor=...) -> None: ...

class QrCodeWidget(Widget):
    codeName: str
    x: int
    y: int
    barFillColor: Incomplete
    barStrokeColor: Incomplete
    barStrokeWidth: int
    barHeight: Incomplete
    barWidth: Incomplete
    barBorder: int
    barLevel: str
    qrVersion: Incomplete
    value: Incomplete
    def __init__(self, value: str = "Hello World", **kw) -> None: ...
    def addData(self, value) -> None: ...
    def draw(self): ...

class QrCode(Flowable):
    height: Incomplete
    width: Incomplete
    qrBorder: int
    qrLevel: str
    qrVersion: Incomplete
    value: Incomplete
    qr: Incomplete
    def __init__(self, value=None, **kw) -> None: ...
    def addData(self, value) -> None: ...
    def draw(self) -> None: ...
    def rect(self, x, y, w, h) -> None: ...
