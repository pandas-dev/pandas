from reportlab.graphics.shapes import Drawing, _DrawingEditorMixin
from reportlab.graphics.widgetbase import Widget
from reportlab.lib.attrmap import *
from reportlab.lib.validators import *

class AdjustableArrow(Widget):
    def __init__(self, **kwds) -> None: ...
    def draw(self): ...

class AdjustableArrowDrawing(_DrawingEditorMixin, Drawing):
    def __init__(self, width: int = 100, height: int = 63, *args, **kw) -> None: ...
