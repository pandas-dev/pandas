from _typeshed import Incomplete

from reportlab.graphics.barcode.common import Barcode

class ECC200DataMatrix(Barcode):
    barWidth: int
    row_modules: int
    col_modules: int
    row_regions: int
    col_regions: int
    cw_data: int
    cw_ecc: int
    row_usable_modules: Incomplete
    col_usable_modules: Incomplete
    def __init__(self, *args, **kwargs) -> None: ...
    valid: int
    validated: Incomplete
    def validate(self) -> None: ...
    encoded: Incomplete
    def encode(self): ...
    def computeSize(self, *args) -> None: ...
    def draw(self) -> None: ...

__all__ = ("ECC200DataMatrix",)
