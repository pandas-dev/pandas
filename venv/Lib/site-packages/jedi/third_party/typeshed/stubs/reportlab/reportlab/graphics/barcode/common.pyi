from _typeshed import Incomplete

from reportlab.platypus.flowables import Flowable

class Barcode(Flowable):
    fontName: str
    fontSize: int
    humanReadable: int
    value: Incomplete
    gap: Incomplete
    def __init__(self, value: str = "", **kwd) -> None: ...
    valid: int
    validated: Incomplete
    def validate(self) -> None: ...
    encoded: Incomplete
    def encode(self) -> None: ...
    decomposed: Incomplete
    def decompose(self) -> None: ...
    barHeight: Incomplete
    def computeSize(self, *args) -> None: ...
    @property
    def width(self): ...
    @width.setter
    def width(self, v) -> None: ...
    @property
    def height(self): ...
    @height.setter
    def height(self, v) -> None: ...
    def draw(self) -> None: ...
    def drawHumanReadable(self) -> None: ...
    def rect(self, x, y, w, h) -> None: ...
    def annotate(self, x, y, text, fontName, fontSize, anchor: str = "middle") -> None: ...

class MultiWidthBarcode(Barcode):
    barHeight: Incomplete
    def computeSize(self, *args) -> None: ...
    def draw(self) -> None: ...

class I2of5(Barcode):
    patterns: Incomplete
    barHeight: Incomplete
    barWidth: Incomplete
    ratio: float
    checksum: int
    bearers: float
    bearerBox: bool
    quiet: int
    lquiet: Incomplete
    rquiet: Incomplete
    stop: int
    def __init__(self, value: str = "", **args) -> None: ...
    valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: Incomplete
    def encode(self) -> None: ...
    decomposed: Incomplete
    def decompose(self): ...

class MSI(Barcode):
    patterns: Incomplete
    stop: int
    barHeight: Incomplete
    barWidth: Incomplete
    ratio: float
    checksum: int
    bearers: float
    quiet: int
    lquiet: Incomplete
    rquiet: Incomplete
    def __init__(self, value: str = "", **args) -> None: ...
    valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: Incomplete
    def encode(self) -> None: ...
    decomposed: Incomplete
    def decompose(self): ...

class Codabar(Barcode):
    patterns: Incomplete
    values: Incomplete
    chars: Incomplete
    stop: int
    barHeight: Incomplete
    barWidth: Incomplete
    ratio: float
    checksum: int
    bearers: float
    quiet: int
    lquiet: Incomplete
    rquiet: Incomplete
    def __init__(self, value: str = "", **args) -> None: ...
    valid: int
    Valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: Incomplete
    def encode(self) -> None: ...
    decomposed: Incomplete
    def decompose(self): ...

class Code11(Barcode):
    chars: str
    patterns: Incomplete
    values: Incomplete
    stop: int
    barHeight: Incomplete
    barWidth: Incomplete
    ratio: float
    checksum: int
    bearers: float
    quiet: int
    lquiet: Incomplete
    rquiet: Incomplete
    def __init__(self, value: str = "", **args) -> None: ...
    valid: int
    Valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: Incomplete
    def encode(self) -> None: ...
    decomposed: Incomplete
    def decompose(self): ...
