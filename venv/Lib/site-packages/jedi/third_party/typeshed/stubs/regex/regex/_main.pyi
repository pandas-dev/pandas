from _typeshed import ReadableBuffer, Unused
from collections.abc import Callable, Mapping
from types import GenericAlias
from typing import Any, AnyStr, Generic, Literal, TypeVar, final, overload
from typing_extensions import Self

from . import _regex
from ._regex_core import *

_T = TypeVar("_T")

__version__: str

# Sync with regex.__init__.__all__
__all__ = [
    "cache_all",
    "compile",
    "DEFAULT_VERSION",
    "escape",
    "findall",
    "finditer",
    "fullmatch",
    "match",
    "prefixmatch",
    "purge",
    "search",
    "split",
    "splititer",
    "sub",
    "subf",
    "subfn",
    "subn",
    "template",
    "Scanner",
    "A",
    "ASCII",
    "B",
    "BESTMATCH",
    "D",
    "DEBUG",
    "E",
    "ENHANCEMATCH",
    "S",
    "DOTALL",
    "F",
    "FULLCASE",
    "I",
    "IGNORECASE",
    "L",
    "LOCALE",
    "M",
    "MULTILINE",
    "P",
    "POSIX",
    "R",
    "REVERSE",
    "T",
    "TEMPLATE",
    "U",
    "UNICODE",
    "V0",
    "VERSION0",
    "V1",
    "VERSION1",
    "X",
    "VERBOSE",
    "W",
    "WORD",
    "error",
    "Regex",
    "__version__",
    "__doc__",
    "RegexFlag",
    "Pattern",
    "Match",
]

def compile(
    pattern: AnyStr | Pattern[AnyStr],
    flags: int = 0,
    ignore_unused: bool = False,
    cache_pattern: bool | None = None,
    **kwargs: Any,
) -> Pattern[AnyStr]: ...
@overload
def search(
    pattern: str | Pattern[str],
    string: str,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> Match[str] | None: ...
@overload
def search(
    pattern: bytes | Pattern[bytes],
    string: ReadableBuffer,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> Match[bytes] | None: ...
@overload
def match(
    pattern: str | Pattern[str],
    string: str,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> Match[str] | None: ...
@overload
def match(
    pattern: bytes | Pattern[bytes],
    string: ReadableBuffer,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> Match[bytes] | None: ...

prefixmatch = match

@overload
def fullmatch(
    pattern: str | Pattern[str],
    string: str,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> Match[str] | None: ...
@overload
def fullmatch(
    pattern: bytes | Pattern[bytes],
    string: ReadableBuffer,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> Match[bytes] | None: ...
@overload
def split(
    pattern: str | Pattern[str],
    string: str,
    maxsplit: int = 0,
    flags: int = 0,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> list[str | Any]: ...
@overload
def split(
    pattern: bytes | Pattern[bytes],
    string: ReadableBuffer,
    maxsplit: int = 0,
    flags: int = 0,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> list[bytes | Any]: ...
@overload
def splititer(
    pattern: str | Pattern[str],
    string: str,
    maxsplit: int = 0,
    flags: int = 0,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> _regex.Splitter[str]: ...
@overload
def splititer(
    pattern: bytes | Pattern[bytes],
    string: ReadableBuffer,
    maxsplit: int = 0,
    flags: int = 0,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> _regex.Splitter[bytes]: ...
@overload
def findall(
    pattern: str | Pattern[str],
    string: str,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    overlapped: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> list[Any]: ...
@overload
def findall(
    pattern: bytes | Pattern[bytes],
    string: ReadableBuffer,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    overlapped: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> list[Any]: ...
@overload
def finditer(
    pattern: str | Pattern[str],
    string: str,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    overlapped: bool = False,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> _regex.Scanner[str]: ...
@overload
def finditer(
    pattern: bytes | Pattern[bytes],
    string: ReadableBuffer,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    overlapped: bool = False,
    partial: bool = False,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> _regex.Scanner[bytes]: ...
@overload
def sub(
    pattern: str | Pattern[str],
    repl: str | Callable[[Match[str]], str],
    string: str,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> str: ...
@overload
def sub(
    pattern: bytes | Pattern[bytes],
    repl: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
    string: ReadableBuffer,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> bytes: ...
@overload
def subf(
    pattern: str | Pattern[str],
    format: str | Callable[[Match[str]], str],
    string: str,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> str: ...
@overload
def subf(
    pattern: bytes | Pattern[bytes],
    format: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
    string: ReadableBuffer,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> bytes: ...
@overload
def subn(
    pattern: str | Pattern[str],
    repl: str | Callable[[Match[str]], str],
    string: str,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> tuple[str, int]: ...
@overload
def subn(
    pattern: bytes | Pattern[bytes],
    repl: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
    string: ReadableBuffer,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> tuple[bytes, int]: ...
@overload
def subfn(
    pattern: str | Pattern[str],
    format: str | Callable[[Match[str]], str],
    string: str,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> tuple[str, int]: ...
@overload
def subfn(
    pattern: bytes | Pattern[bytes],
    format: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
    string: ReadableBuffer,
    count: int = 0,
    flags: int = 0,
    pos: int | None = None,
    endpos: int | None = None,
    concurrent: bool | None = None,
    timeout: float | None = None,
    ignore_unused: bool = False,
    **kwargs: Any,
) -> tuple[bytes, int]: ...
def purge() -> None: ...
@overload
def cache_all(value: bool = True) -> None: ...
@overload
def cache_all(value: None) -> bool: ...
def escape(pattern: AnyStr, special_only: bool = True, literal_spaces: bool = False) -> AnyStr: ...

DEFAULT_VERSION = RegexFlag.VERSION0

def template(pattern: AnyStr | Pattern[AnyStr], flags: int = 0) -> Pattern[AnyStr]: ...

Regex = compile

@final
class Pattern(Generic[AnyStr]):
    @property
    def flags(self) -> int: ...
    @property
    def groupindex(self) -> Mapping[str, int]: ...
    @property
    def groups(self) -> int: ...
    @property
    def pattern(self) -> AnyStr: ...
    @property
    def named_lists(self) -> Mapping[str, frozenset[AnyStr]]: ...
    @overload
    def search(
        self: Pattern[str],
        string: str,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> Match[str] | None: ...
    @overload
    def search(
        self: Pattern[bytes],
        string: ReadableBuffer,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> Match[bytes] | None: ...
    @overload
    def match(
        self: Pattern[str],
        string: str,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> Match[str] | None: ...
    @overload
    def match(
        self: Pattern[bytes],
        string: ReadableBuffer,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> Match[bytes] | None: ...
    prefixmatch = match
    @overload
    def fullmatch(
        self: Pattern[str],
        string: str,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> Match[str] | None: ...
    @overload
    def fullmatch(
        self: Pattern[bytes],
        string: ReadableBuffer,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> Match[bytes] | None: ...
    @overload
    def split(
        self: Pattern[str], string: str, maxsplit: int = 0, concurrent: bool | None = None, timeout: float | None = None
    ) -> list[str | Any]: ...
    @overload
    def split(
        self: Pattern[bytes],
        string: ReadableBuffer,
        maxsplit: int = 0,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> list[bytes | Any]: ...
    @overload
    def splititer(
        self: Pattern[str], string: str, maxsplit: int = 0, concurrent: bool | None = None, timeout: float | None = None
    ) -> _regex.Splitter[str]: ...
    @overload
    def splititer(
        self: Pattern[bytes],
        string: ReadableBuffer,
        maxsplit: int = 0,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> _regex.Splitter[bytes]: ...
    @overload
    def findall(
        self: Pattern[str],
        string: str,
        pos: int | None = None,
        endpos: int | None = None,
        overlapped: bool = False,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> list[Any]: ...
    @overload
    def findall(
        self: Pattern[bytes],
        string: ReadableBuffer,
        pos: int | None = None,
        endpos: int | None = None,
        overlapped: bool = False,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> list[Any]: ...
    @overload
    def finditer(
        self: Pattern[str],
        string: str,
        pos: int | None = None,
        endpos: int | None = None,
        overlapped: bool = False,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> _regex.Scanner[str]: ...
    @overload
    def finditer(
        self: Pattern[bytes],
        string: ReadableBuffer,
        pos: int | None = None,
        endpos: int | None = None,
        overlapped: bool = False,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> _regex.Scanner[bytes]: ...
    @overload
    def sub(
        self: Pattern[str],
        repl: str | Callable[[Match[str]], str],
        string: str,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> str: ...
    @overload
    def sub(
        self: Pattern[bytes],
        repl: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
        string: ReadableBuffer,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> bytes: ...
    @overload
    def subf(
        self: Pattern[str],
        format: str | Callable[[Match[str]], str],
        string: str,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> str: ...
    @overload
    def subf(
        self: Pattern[bytes],
        format: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
        string: ReadableBuffer,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> bytes: ...
    @overload
    def subn(
        self: Pattern[str],
        repl: str | Callable[[Match[str]], str],
        string: str,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> tuple[str, int]: ...
    @overload
    def subn(
        self: Pattern[bytes],
        repl: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
        string: ReadableBuffer,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> tuple[bytes, int]: ...
    @overload
    def subfn(
        self: Pattern[str],
        format: str | Callable[[Match[str]], str],
        string: str,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> tuple[str, int]: ...
    @overload
    def subfn(
        self: Pattern[bytes],
        format: ReadableBuffer | Callable[[Match[bytes]], ReadableBuffer],
        string: ReadableBuffer,
        count: int = 0,
        pos: int | None = None,
        endpos: int | None = None,
        concurrent: bool | None = None,
        timeout: float | None = None,
    ) -> tuple[bytes, int]: ...
    @overload
    def scanner(
        self: Pattern[str],
        string: str,
        pos: int | None = None,
        endpos: int | None = None,
        overlapped: bool = False,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> _regex.Scanner[str]: ...
    @overload
    def scanner(
        self: Pattern[bytes],
        string: bytes,
        pos: int | None = None,
        endpos: int | None = None,
        overlapped: bool = False,
        concurrent: bool | None = None,
        partial: bool = False,
        timeout: float | None = None,
    ) -> _regex.Scanner[bytes]: ...
    def __copy__(self) -> Self: ...
    def __deepcopy__(self, memo: Unused, /) -> Self: ...
    def __class_getitem__(cls, item: Any, /) -> GenericAlias: ...

@final
class Match(Generic[AnyStr]):
    @property
    def pos(self) -> int: ...
    @property
    def endpos(self) -> int: ...
    @property
    def lastindex(self) -> int | None: ...
    @property
    def lastgroup(self) -> str | None: ...
    @property
    def string(self) -> AnyStr: ...
    @property
    def re(self) -> Pattern[AnyStr]: ...
    @property
    def partial(self) -> bool: ...
    @property
    def regs(self) -> tuple[tuple[int, int], ...]: ...
    @property
    def fuzzy_counts(self) -> tuple[int, int, int]: ...
    @property
    def fuzzy_changes(self) -> tuple[list[int], list[int], list[int]]: ...
    @overload
    def group(self, group: Literal[0] = 0, /) -> AnyStr: ...
    @overload
    def group(self, group: int | str = ..., /) -> AnyStr | Any: ...
    @overload
    def group(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[AnyStr | Any, ...]: ...
    @overload
    def groups(self, default: None = None) -> tuple[AnyStr | Any, ...]: ...
    @overload
    def groups(self, default: _T) -> tuple[AnyStr | _T, ...]: ...
    @overload
    def groupdict(self, default: None = None) -> dict[str, AnyStr | Any]: ...
    @overload
    def groupdict(self, default: _T) -> dict[str, AnyStr | _T]: ...
    @overload
    def span(self, group: int | str = ..., /) -> tuple[int, int]: ...
    @overload
    def span(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[tuple[int, int], ...]: ...
    @overload
    def spans(self, group: int | str = ..., /) -> list[tuple[int, int]]: ...
    @overload
    def spans(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[list[tuple[int, int]], ...]: ...
    @overload
    def start(self, group: int | str = ..., /) -> int: ...
    @overload
    def start(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[int, ...]: ...
    @overload
    def starts(self, group: int | str = ..., /) -> list[int]: ...
    @overload
    def starts(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[list[int], ...]: ...
    @overload
    def end(self, group: int | str = ..., /) -> int: ...
    @overload
    def end(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[int, ...]: ...
    @overload
    def ends(self, group: int | str = ..., /) -> list[int]: ...
    @overload
    def ends(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[list[int], ...]: ...
    def expand(self, template: AnyStr, /) -> AnyStr: ...
    def expandf(self, format: AnyStr, /) -> AnyStr: ...
    @overload
    def captures(self, group: int | str = ..., /) -> list[AnyStr]: ...
    @overload
    def captures(self, group1: int | str, group2: int | str, /, *groups: int | str) -> tuple[list[AnyStr], ...]: ...
    def capturesdict(self) -> dict[str, list[AnyStr]]: ...
    def detach_string(self) -> None: ...
    def allcaptures(self) -> tuple[list[AnyStr]]: ...
    def allspans(self) -> tuple[list[tuple[int, int]]]: ...
    @overload
    def __getitem__(self, key: Literal[0], /) -> AnyStr: ...
    @overload
    def __getitem__(self, key: int | str, /) -> AnyStr | Any: ...
    def __copy__(self) -> Self: ...
    def __deepcopy__(self, memo: Unused, /) -> Self: ...
    def __class_getitem__(cls, item: Any, /) -> GenericAlias: ...
