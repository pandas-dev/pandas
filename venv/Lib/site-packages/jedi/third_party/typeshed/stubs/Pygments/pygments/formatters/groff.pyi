from _typeshed import SupportsWrite
from collections.abc import Iterable
from typing import TypeVar

from pygments.formatter import Formatter
from pygments.token import _TokenType

__all__ = ["GroffFormatter"]

_T = TypeVar("_T", str, bytes)

class GroffFormatter(Formatter[_T]):
    monospaced: bool
    linenos: bool
    wrap: int
    styles: dict[_TokenType, tuple[str, str]]
    def __init__(self, **options) -> None: ...
    def format_unencoded(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[str]) -> None: ...
