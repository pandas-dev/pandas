from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable
from typing import TypeVar

from pygments.formatter import Formatter
from pygments.token import _TokenType

_T = TypeVar("_T", str, bytes)

__all__ = ["RtfFormatter"]

class RtfFormatter(Formatter[_T]):
    fontface: Incomplete
    fontsize: Incomplete
    @staticmethod
    def hex_to_rtf_color(hex_color: str) -> str: ...
    def format_unencoded(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[str]) -> None: ...
