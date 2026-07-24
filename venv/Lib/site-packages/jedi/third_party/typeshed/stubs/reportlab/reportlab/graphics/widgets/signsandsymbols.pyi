from _typeshed import Incomplete
from typing import Final

from reportlab.graphics.widgetbase import Widget
from reportlab.lib.attrmap import *
from reportlab.lib.validators import *

__version__: Final[str]

class _Symbol(Widget):
    x: int
    size: int
    fillColor: Incomplete
    strokeColor: Incomplete
    strokeWidth: float
    def __init__(self) -> None: ...
    def demo(self): ...

class ETriangle(_Symbol):
    def __init__(self) -> None: ...
    def draw(self): ...

class RTriangle(_Symbol):
    x: int
    y: int
    size: int
    fillColor: Incomplete
    strokeColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class Octagon(_Symbol):
    x: int
    y: int
    size: int
    fillColor: Incomplete
    strokeColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class Crossbox(_Symbol):
    x: int
    y: int
    size: int
    fillColor: Incomplete
    crossColor: Incomplete
    strokeColor: Incomplete
    crosswidth: int
    def __init__(self) -> None: ...
    def draw(self): ...

class Tickbox(_Symbol):
    x: int
    y: int
    size: int
    tickColor: Incomplete
    strokeColor: Incomplete
    fillColor: Incomplete
    tickwidth: int
    def __init__(self) -> None: ...
    def draw(self): ...

class SmileyFace(_Symbol):
    x: int
    y: int
    size: int
    fillColor: Incomplete
    strokeColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class StopSign(_Symbol):
    x: int
    y: int
    size: int
    strokeColor: Incomplete
    fillColor: Incomplete
    stopColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class NoEntry(_Symbol):
    x: int
    y: int
    size: int
    strokeColor: Incomplete
    fillColor: Incomplete
    innerBarColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class NotAllowed(_Symbol):
    x: int
    y: int
    size: int
    strokeColor: Incomplete
    fillColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class NoSmoking(NotAllowed):
    def __init__(self) -> None: ...
    def draw(self): ...

class DangerSign(_Symbol):
    x: int
    y: int
    size: int
    strokeColor: Incomplete
    fillColor: Incomplete
    strokeWidth: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class YesNo(_Symbol):
    x: int
    y: int
    size: int
    tickcolor: Incomplete
    crosscolor: Incomplete
    testValue: int
    def __init__(self) -> None: ...
    def draw(self): ...
    def demo(self): ...

class FloppyDisk(_Symbol):
    x: int
    y: int
    size: int
    diskColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class ArrowOne(_Symbol):
    x: int
    y: int
    size: int
    fillColor: Incomplete
    strokeWidth: int
    strokeColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class ArrowTwo(ArrowOne):
    x: int
    y: int
    size: int
    fillColor: Incomplete
    strokeWidth: int
    strokeColor: Incomplete
    def __init__(self) -> None: ...
    def draw(self): ...

class CrossHair(_Symbol):
    x: int
    size: int
    fillColor: Incomplete
    strokeColor: Incomplete
    strokeWidth: float
    innerGap: str
    def __init__(self) -> None: ...
    def draw(self): ...

def test() -> None: ...
