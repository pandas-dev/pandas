from collections.abc import Iterable
from typing import NamedTuple

from .ast import Node

class RGBA(NamedTuple):
    red: float
    green: float
    blue: float
    alpha: float

def parse_color(input: str | Iterable[Node]) -> str | RGBA | None: ...
