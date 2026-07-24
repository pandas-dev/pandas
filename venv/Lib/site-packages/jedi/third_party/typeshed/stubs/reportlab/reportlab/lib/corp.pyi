from _typeshed import Incomplete
from collections.abc import Iterable
from typing import Final

from reportlab.graphics.shapes import Drawing, Group
from reportlab.graphics.widgetbase import Widget
from reportlab.lib.attrmap import *
from reportlab.lib.validators import *

__version__: Final[str]

class RL_CorpLogo(Widget):
    fillColor: Incomplete
    strokeColor: Incomplete
    strokeWidth: float
    background: Incomplete
    border: Incomplete
    borderWidth: int
    shadow: float
    height: int
    width: int
    x: int
    skewX: int
    showPage: int
    oColors: Incomplete
    pageColors: Incomplete
    prec: Incomplete
    def __init__(self) -> None: ...
    def demo(self) -> Drawing: ...
    @staticmethod
    def applyPrec(P, prec): ...
    def draw(self) -> Group: ...

class RL_CorpLogoReversed(RL_CorpLogo):
    background: Incomplete
    fillColor: Incomplete
    def __init__(self) -> None: ...

class RL_CorpLogoThin(Widget):
    fillColor: Incomplete
    strokeColor: Incomplete
    x: int
    y: int
    height: Incomplete
    width: Incomplete
    def __init__(self) -> None: ...
    def demo(self) -> Drawing: ...
    def draw(self) -> Group: ...

class ReportLabLogo:
    origin: Incomplete
    dimensions: Incomplete
    powered_by: Incomplete
    def __init__(self, atx: int = 0, aty: int = 0, width=180.0, height=108.0, powered_by: int = 0) -> None: ...
    def draw(self, canvas) -> None: ...

class RL_BusinessCard(Widget):
    fillColor: Incomplete
    strokeColor: Incomplete
    altStrokeColor: Incomplete
    x: int
    y: int
    height: Incomplete
    width: Incomplete
    borderWidth: Incomplete
    bleed: Incomplete
    cropMarks: int
    border: int
    name: str
    position: str
    telephone: str
    mobile: str
    fax: str
    email: str
    web: str
    rh_blurb_top: Incomplete
    def __init__(self) -> None: ...
    def demo(self) -> Drawing: ...
    def draw(self) -> Group: ...

def test(formats: Iterable[str] = ["pdf", "eps", "jpg", "gif", "svg"]) -> None: ...
