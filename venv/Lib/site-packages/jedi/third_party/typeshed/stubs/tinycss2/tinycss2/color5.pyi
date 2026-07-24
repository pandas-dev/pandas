from collections.abc import Iterable
from typing import Literal

from . import color4
from .ast import Node

COLOR_SCHEMES: set[str]
COLOR_SPACES: set[str]
D50: tuple[float, float, float]
D65: tuple[float, float, float]

class Color(color4.Color):
    COLOR_SPACES: set[str] | None

def parse_color(
    input: str | Iterable[Node], color_schemes: Literal["normal"] | Iterable[str] | None = None
) -> color4.Color | Color | Literal["currentcolor"] | None: ...
