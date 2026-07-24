from _typeshed import SupportsWrite
from collections.abc import Iterable
from typing import TypeVar

from pygments.formatter import Formatter
from pygments.token import _TokenType

_T = TypeVar("_T", str, bytes)

__all__ = ["BBCodeFormatter"]

class BBCodeFormatter(Formatter[_T]):
    styles: dict[_TokenType, tuple[str, str]]
    def format_unencoded(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[str]) -> None: ...
