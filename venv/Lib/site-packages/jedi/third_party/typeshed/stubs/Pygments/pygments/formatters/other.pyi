from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable
from typing import TypeVar

from pygments.formatter import Formatter
from pygments.token import _TokenType

_T = TypeVar("_T", str, bytes)

__all__ = ["NullFormatter", "RawTokenFormatter", "TestcaseFormatter"]

class NullFormatter(Formatter[_T]):
    def format(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[_T]) -> None: ...

class RawTokenFormatter(Formatter[bytes]):
    encoding: str
    compress: Incomplete
    error_color: Incomplete
    def format(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[bytes]) -> None: ...

class TestcaseFormatter(Formatter[_T]):
    def format(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[_T]) -> None: ...
