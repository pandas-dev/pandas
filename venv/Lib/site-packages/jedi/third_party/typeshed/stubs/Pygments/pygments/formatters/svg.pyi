from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable
from typing import TypeVar

from pygments.formatter import Formatter
from pygments.token import _TokenType

_T = TypeVar("_T", str, bytes)

__all__ = ["SvgFormatter"]

class SvgFormatter(Formatter[_T]):
    nowrap: Incomplete
    fontfamily: Incomplete
    fontsize: Incomplete
    xoffset: Incomplete
    yoffset: Incomplete
    ystep: Incomplete
    spacehack: Incomplete
    linenos: Incomplete
    linenostart: Incomplete
    linenostep: Incomplete
    linenowidth: Incomplete
    def format_unencoded(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[str]) -> None: ...
