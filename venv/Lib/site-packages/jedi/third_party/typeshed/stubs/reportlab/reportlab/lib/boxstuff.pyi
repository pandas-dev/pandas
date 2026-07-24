from _typeshed import Incomplete
from typing import Final, Literal, overload

__version__: Final[str]

@overload
def rectCorner(
    x, y, width, height, anchor: str = "sw", dims: Literal[True] = ...
) -> tuple[Incomplete, Incomplete, Incomplete, Incomplete]: ...
@overload
def rectCorner(x, y, width, height, anchor: str = "sw", dims: Literal[False] | None = False) -> tuple[Incomplete, Incomplete]: ...
def aspectRatioFix(
    preserve, anchor, x, y, width, height, imWidth, imHeight, anchorAtXY: bool = False
) -> tuple[Incomplete, Incomplete, Incomplete, Incomplete, Incomplete]: ...
