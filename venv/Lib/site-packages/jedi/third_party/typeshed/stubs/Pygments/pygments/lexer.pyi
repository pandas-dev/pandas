from _typeshed import Incomplete
from collections.abc import Iterable, Iterator, Sequence
from re import Pattern, RegexFlag
from typing import ClassVar, Final

from pygments.token import _TokenType
from pygments.util import Future

__all__ = [
    "Lexer",
    "RegexLexer",
    "ExtendedRegexLexer",
    "DelegatingLexer",
    "LexerContext",
    "include",
    "inherit",
    "bygroups",
    "using",
    "this",
    "default",
    "words",
    "line_re",
]

line_re: Final[Pattern[str]]

class LexerMeta(type):
    def __new__(cls, name, bases, d): ...

class Lexer(metaclass=LexerMeta):
    name: ClassVar[str]  # Set to None, but always overridden with a non-None value in subclasses.
    aliases: ClassVar[Sequence[str]]  # not intended to be mutable
    filenames: ClassVar[Sequence[str]]
    alias_filenames: ClassVar[Sequence[str]]
    mimetypes: ClassVar[Sequence[str]]
    priority: ClassVar[float]
    url: ClassVar[str]  # Set to None, but always overridden with a non-None value in subclasses.
    version_added: ClassVar[str]  # Set to None, but always overridden with a non-None value in subclasses.
    options: Incomplete
    stripnl: Incomplete
    stripall: Incomplete
    ensurenl: Incomplete
    tabsize: Incomplete
    encoding: Incomplete
    filters: Incomplete
    def __init__(self, **options) -> None: ...
    def add_filter(self, filter_, **options) -> None: ...
    @staticmethod  # @staticmethod added by special handling in metaclass
    def analyse_text(text: str) -> float: ...
    def get_tokens(self, text: str, unfiltered: bool = False) -> Iterator[tuple[_TokenType, str]]: ...
    def get_tokens_unprocessed(self, text: str) -> Iterator[tuple[int, _TokenType, str]]: ...

class DelegatingLexer(Lexer):
    root_lexer: Incomplete
    language_lexer: Incomplete
    needle: Incomplete
    def __init__(self, _root_lexer, _language_lexer, _needle=..., **options) -> None: ...
    def get_tokens_unprocessed(self, text: str) -> Iterator[tuple[int, _TokenType, str]]: ...

class include(str): ...
class _inherit: ...

inherit: Incomplete

class combined(tuple[Incomplete, ...]):
    def __new__(cls, *args): ...
    def __init__(self, *args) -> None: ...

class _PseudoMatch:
    def __init__(self, start, text) -> None: ...
    def start(self, arg=None): ...
    def end(self, arg=None): ...
    def group(self, arg=None): ...
    def groups(self): ...
    def groupdict(self): ...

def bygroups(*args): ...

class _This: ...

this: Incomplete

def using(_other, **kwargs): ...

class default:
    state: Incomplete
    def __init__(self, state) -> None: ...

class words(Future):
    words: Incomplete
    prefix: Incomplete
    suffix: Incomplete
    def __init__(self, words, prefix: str = "", suffix: str = "") -> None: ...
    def get(self): ...

class RegexLexerMeta(LexerMeta):
    def process_tokendef(cls, name, tokendefs=None): ...
    def get_tokendefs(cls): ...
    def __call__(cls, *args, **kwds): ...

class RegexLexer(Lexer, metaclass=RegexLexerMeta):
    flags: ClassVar[RegexFlag]
    tokens: ClassVar[dict[str, list[Incomplete]]]
    def get_tokens_unprocessed(self, text: str, stack: Iterable[str] = ("root",)) -> Iterator[tuple[int, _TokenType, str]]: ...

class LexerContext:
    text: Incomplete
    pos: Incomplete
    end: Incomplete
    stack: Incomplete
    def __init__(self, text, pos, stack=None, end=None) -> None: ...

class ExtendedRegexLexer(RegexLexer):
    def get_tokens_unprocessed(  # type: ignore[override]
        self, text: str | None = None, context: LexerContext | None = None
    ) -> Iterator[tuple[int, _TokenType, str]]: ...

class ProfilingRegexLexerMeta(RegexLexerMeta): ...

class ProfilingRegexLexer(RegexLexer, metaclass=ProfilingRegexLexerMeta):
    def get_tokens_unprocessed(self, text: str, stack: Iterable[str] = ("root",)) -> Iterator[tuple[int, _TokenType, str]]: ...
