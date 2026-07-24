from _typeshed import Incomplete, SupportsWrite
from collections.abc import Iterable
from typing import TypeVar

from pygments.formatter import Formatter
from pygments.token import _TokenType

_T = TypeVar("_T", str, bytes)

__all__ = ["Terminal256Formatter", "TerminalTrueColorFormatter"]

class EscapeSequence:
    fg: Incomplete
    bg: Incomplete
    bold: Incomplete
    underline: Incomplete
    italic: Incomplete
    def __init__(self, fg=None, bg=None, bold: bool = False, underline: bool = False, italic: bool = False) -> None: ...
    def escape(self, attrs): ...
    def color_string(self): ...
    def true_color_string(self): ...
    def reset_string(self): ...

class Terminal256Formatter(Formatter[_T]):
    xterm_colors: Incomplete
    best_match: Incomplete
    style_string: Incomplete
    usebold: Incomplete
    useunderline: Incomplete
    useitalic: Incomplete
    linenos: Incomplete

    def format(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[_T]) -> None: ...
    def format_unencoded(self, tokensource: Iterable[tuple[_TokenType, str]], outfile: SupportsWrite[str]) -> None: ...

class TerminalTrueColorFormatter(Terminal256Formatter[_T]): ...
