from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.charts.areas import PlotArea
from reportlab.lib.colors import black

class _BarcodeWidget(PlotArea):
    textColor = black
    barFillColor = black
    barStrokeColor: Incomplete
    barStrokeWidth: int
    x: int
    def __init__(self, _value: str = "", **kw) -> None: ...
    def rect(self, x, y, w, h, **kw) -> None: ...
    canv: Incomplete
    def draw(self): ...
    def annotate(self, x, y, text, fontName, fontSize, anchor: str = "middle") -> None: ...

class BarcodeI2of5(_BarcodeWidget):
    codeName: Final = "I2of5"
    def __init__(self, **kw) -> None: ...

class BarcodeCode128(_BarcodeWidget):
    codeName: Final = "Code128"
    def __init__(self, **kw) -> None: ...

class BarcodeStandard93(_BarcodeWidget):
    codeName: Final = "Standard93"
    def __init__(self, **kw) -> None: ...

class BarcodeExtended93(_BarcodeWidget):
    codeName: Final = "Extended93"
    def __init__(self, **kw) -> None: ...

class BarcodeStandard39(_BarcodeWidget):
    codeName: Final = "Standard39"
    def __init__(self, **kw) -> None: ...

class BarcodeExtended39(_BarcodeWidget):
    codeName: Final = "Extended39"
    def __init__(self, **kw) -> None: ...

class BarcodeMSI(_BarcodeWidget):
    codeName: Final = "MSI"
    def __init__(self, **kw) -> None: ...

class BarcodeCodabar(_BarcodeWidget):
    codeName: Final = "Codabar"
    def __init__(self, **kw) -> None: ...

class BarcodeCode11(_BarcodeWidget):
    codeName: Final = "Code11"
    def __init__(self, **kw) -> None: ...

class BarcodeFIM(_BarcodeWidget):
    codeName: Final = "FIM"
    def __init__(self, **kw) -> None: ...

class BarcodePOSTNET(_BarcodeWidget):
    codeName: Final = "POSTNET"
    def __init__(self, **kw) -> None: ...

class BarcodeUSPS_4State(_BarcodeWidget):
    codeName: Final = "USPS_4State"
    def __init__(self, **kw) -> None: ...

__all__ = (
    "BarcodeI2of5",
    "BarcodeCode128",
    "BarcodeStandard93",
    "BarcodeExtended93",
    "BarcodeStandard39",
    "BarcodeExtended39",
    "BarcodeMSI",
    "BarcodeCodabar",
    "BarcodeCode11",
    "BarcodeFIM",
    "BarcodePOSTNET",
    "BarcodeUSPS_4State",
)
