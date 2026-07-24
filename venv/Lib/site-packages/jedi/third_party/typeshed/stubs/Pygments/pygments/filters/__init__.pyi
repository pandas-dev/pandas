from _typeshed import ConvertibleToInt
from collections.abc import Callable, Generator, Iterable, Iterator
from re import Pattern
from typing import Any, ClassVar, Final, Literal

from pygments.filter import Filter
from pygments.lexer import Lexer
from pygments.token import _TokenType

def find_filter_class(filtername: str) -> type[Filter] | None: ...

# Keyword arguments are forwarded to the filter class.
def get_filter_by_name(filtername: str, **options: Any) -> Filter: ...
def get_all_filters() -> Generator[str]: ...

class CodeTagFilter(Filter):
    tag_re: Pattern[str]
    # Arbitrary additional keyword arguments are permitted and are stored in self.options.
    def __init__(
        self, *, codetags: str | list[str] | tuple[str, ...] = ["XXX", "TODO", "FIXME", "BUG", "NOTE"], **options: Any
    ) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class SymbolFilter(Filter):
    latex_symbols: ClassVar[dict[str, str]]
    isabelle_symbols: ClassVar[dict[str, str]]
    lang_map: ClassVar[dict[Literal["isabelle", "latex"], dict[str, str]]]
    symbols: dict[str, str]  # One of latex_symbols or isabelle_symbols.
    # Arbitrary additional keyword arguments are permitted and are stored in self.options.
    def __init__(self, *, lang: Literal["isabelle", "latex"] = "isabelle", **options: Any) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class KeywordCaseFilter(Filter):
    convert: Callable[[str], str]
    # Arbitrary additional keyword arguments are permitted and are stored in self.options.
    def __init__(self, *, case: Literal["lower", "upper", "capitalize"] = "lower", **options: Any) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class NameHighlightFilter(Filter):
    names: set[str]
    tokentype: _TokenType
    # Arbitrary additional keyword arguments are permitted and are stored in self.options.
    def __init__(
        self, *, names: str | list[str] | tuple[str, ...] = [], tokentype: str | _TokenType | None = None, **options: Any
    ) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class ErrorToken(Exception): ...

class RaiseOnErrorTokenFilter(Filter):
    exception: type[Exception]
    # Arbitrary additional keyword arguments are permitted and are stored in self.options.
    def __init__(self, *, excclass: type[Exception] = ..., **options: Any) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class VisibleWhitespaceFilter(Filter):
    spaces: str
    tabs: str
    newlines: str
    wstt: bool
    def __init__(
        self,
        *,
        spaces: str | bool = False,
        tabs: str | bool = False,
        newlines: str | bool = False,
        tabsize: ConvertibleToInt = 8,
        wstokentype: bool | int | str = True,  # Any value accepted by get_bool_opt.
        # Arbitrary additional keyword arguments are permitted and are stored in self.options.
        **options: Any,
    ) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class GobbleFilter(Filter):
    n: int
    # Arbitrary additional keyword arguments are permitted and are stored in self.options.
    def __init__(self, *, n: ConvertibleToInt = 0, **options: Any) -> None: ...
    def gobble(self, value: str, left: int) -> tuple[str, int]: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class TokenMergeFilter(Filter):
    # Arbitrary additional keyword arguments are permitted and are stored in self.options.
    def __init__(self, **options: Any) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

FILTERS: Final[dict[str, type[Filter]]]
