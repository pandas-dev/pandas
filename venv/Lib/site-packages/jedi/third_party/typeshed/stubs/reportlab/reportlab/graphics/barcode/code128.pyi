from _typeshed import Incomplete

from reportlab.graphics.barcode.common import MultiWidthBarcode

starta: Incomplete
startb: Incomplete
startc: Incomplete
stop: Incomplete
seta: Incomplete
setb: Incomplete
setc: Incomplete
setmap: Incomplete
cStarts: Incomplete
tos: Incomplete

class Code128(MultiWidthBarcode):
    barWidth: Incomplete
    lquiet: Incomplete
    rquiet: Incomplete
    quiet: int
    barHeight: Incomplete
    def __init__(self, value: str = "", **args) -> None: ...
    valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: Incomplete
    def encode(self): ...
    decomposed: Incomplete
    def decompose(self): ...

class Code128Auto(Code128):
    encoded: Incomplete
    def encode(self): ...
