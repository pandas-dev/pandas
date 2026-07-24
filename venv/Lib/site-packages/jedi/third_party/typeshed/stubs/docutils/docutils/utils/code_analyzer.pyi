from collections.abc import Generator, Iterable
from typing import Final, Literal
from typing_extensions import TypeAlias

from docutils import ApplicationError

_TokenNames: TypeAlias = Literal["long", "short", "none"]

__docformat__: Final = "reStructuredText"
with_pygments: bool
unstyled_tokens: list[str]

class LexerError(ApplicationError): ...

class Lexer:
    code: str
    language: str
    tokennames: _TokenNames
    lexer: Lexer | None
    def __init__(self, code: str, language: str, tokennames: _TokenNames = "short") -> None: ...
    def merge(self, tokens: Iterable[tuple[_TokenNames, str]]) -> Generator[tuple[_TokenNames, str]]: ...
    def __iter__(self) -> Generator[tuple[list[str], str]]: ...

class NumberLines:
    tokens: Iterable[tuple[list[str], str]]
    startline: int
    fmt_str: str
    def __init__(self, tokens: Iterable[tuple[list[str], str]], startline: int, endline: int) -> None: ...
    def __iter__(self) -> Generator[tuple[list[str], str]]: ...
