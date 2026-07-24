from _typeshed import Incomplete
from typing import Literal
from typing_extensions import LiteralString

class Formatter:
    pattern: str
    def __init__(self, pattern: str) -> None: ...
    def format(self, obj: object) -> str: ...
    def __call__(self, x: object) -> str: ...

class DecimalFormatter(Formatter):
    calcPlaces: Incomplete
    places: int
    dot: Incomplete
    comma: Incomplete
    prefix: Incomplete
    suffix: Incomplete
    def __init__(
        self, places: int | Literal["auto"] = 2, decimalSep: str = ".", thousandSep=None, prefix=None, suffix=None
    ) -> None: ...
    def format(self, num) -> LiteralString: ...

__all__ = ("Formatter", "DecimalFormatter")
