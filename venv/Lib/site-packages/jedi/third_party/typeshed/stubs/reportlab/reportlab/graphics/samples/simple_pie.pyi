from _typeshed import Incomplete

from reportlab.graphics.samples.excelcolors import *
from reportlab.graphics.shapes import Drawing, _DrawingEditorMixin

class SimplePie(_DrawingEditorMixin, Drawing):
    background: Incomplete
    def __init__(self, width: int = 200, height: int = 150, *args, **kw) -> None: ...
