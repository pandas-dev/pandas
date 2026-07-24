from collections.abc import Iterable, Iterator
from typing import Any, ClassVar, Protocol, type_check_only

from pygments.lexer import Lexer
from pygments.token import _TokenType

@type_check_only
class _SimpleFilterFunction(Protocol):
    # Function that can looked up as a method on a FunctionFilter subclass.
    def __call__(
        self, self_: FunctionFilter, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]], options: dict[str, Any], /
    ) -> Iterator[tuple[_TokenType, str]]: ...

def apply_filters(
    stream: Iterable[tuple[_TokenType, str]], filters: Iterable[Filter], lexer: Lexer | None = None
) -> Iterator[tuple[_TokenType, str]]: ...
def simplefilter(f: _SimpleFilterFunction) -> type[FunctionFilter]: ...

class Filter:
    options: dict[str, Any]  # Arbitrary values used by subclasses.
    def __init__(self, **options: Any) -> None: ...  # ditto.
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...

class FunctionFilter(Filter):
    # Set to None in class, but overridden with a non-None value in the subclasses created by @simplefilter.
    function: ClassVar[_SimpleFilterFunction]
    # 'options' gets passed as a dict to 'function'; valid types depends on the wrapped function's signature.
    def __init__(self, **options: Any) -> None: ...
    def filter(self, lexer: Lexer | None, stream: Iterable[tuple[_TokenType, str]]) -> Iterator[tuple[_TokenType, str]]: ...
