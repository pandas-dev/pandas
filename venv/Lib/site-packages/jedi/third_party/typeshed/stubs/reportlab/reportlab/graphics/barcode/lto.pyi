from _typeshed import Incomplete

from reportlab.graphics.barcode.code39 import Standard39

class BaseLTOLabel(Standard39):
    LABELWIDTH: Incomplete
    LABELHEIGHT: Incomplete
    LABELROUND: Incomplete
    CODERATIO: float
    CODENOMINALWIDTH: Incomplete
    CODEBARHEIGHT: Incomplete
    CODEBARWIDTH: Incomplete
    CODEGAP = CODEBARWIDTH
    CODELQUIET: Incomplete
    CODERQUIET: Incomplete
    height: Incomplete
    border: Incomplete
    label: Incomplete
    def __init__(
        self, prefix: str = "", number=None, subtype: str = "1", border=None, checksum: bool = False, availheight=None
    ) -> None: ...
    def drawOn(self, canvas, x, y) -> None: ...

class VerticalLTOLabel(BaseLTOLabel):
    LABELFONT: Incomplete
    BLOCKWIDTH: Incomplete
    BLOCKHEIGHT: Incomplete
    LINEWIDTH: float
    NBBLOCKS: int
    COLORSCHEME: Incomplete
    colored: Incomplete
    def __init__(self, *args, **kwargs) -> None: ...
    def drawOn(self, canvas, x, y) -> None: ...

def test() -> None: ...
