from _typeshed import Incomplete

from reportlab.graphics.barcode.common import MultiWidthBarcode

class _Code93Base(MultiWidthBarcode):
    barWidth: Incomplete
    lquiet: Incomplete
    rquiet: Incomplete
    quiet: int
    barHeight: Incomplete
    stop: int
    def __init__(self, value: str = "", **args) -> None: ...
    decomposed: Incomplete
    def decompose(self): ...

class Standard93(_Code93Base):
    valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: Incomplete
    def encode(self): ...

class Extended93(_Code93Base):
    valid: int
    validated: Incomplete
    def validate(self): ...
    encoded: str
    def encode(self): ...
