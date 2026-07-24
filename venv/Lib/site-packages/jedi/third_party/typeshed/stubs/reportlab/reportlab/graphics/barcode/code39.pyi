from _typeshed import Incomplete

from reportlab.graphics.barcode.common import Barcode

class _Code39Base(Barcode):
    barWidth: Incomplete
    lquiet: Incomplete
    rquiet: Incomplete
    quiet: int
    gap: Incomplete
    barHeight: Incomplete
    ratio: float
    checksum: int
    bearers: float
    stop: int
    def __init__(self, value: str = "", **args) -> None: ...
    decomposed: Incomplete
    def decompose(self): ...

class Standard39(_Code39Base):
    valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: Incomplete
    def encode(self): ...

class Extended39(_Code39Base):
    valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: str
    def encode(self): ...
