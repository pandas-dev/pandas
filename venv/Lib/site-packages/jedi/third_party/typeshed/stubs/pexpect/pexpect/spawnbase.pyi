from asyncio import ReadTransport
from collections.abc import Awaitable, Callable, Iterable
from re import Match, Pattern
from typing import IO, AnyStr, Generic, Literal, Protocol, TextIO, overload, type_check_only
from typing_extensions import TypeAlias

from ._async import PatternWaiter
from .exceptions import EOF, TIMEOUT
from .expect import searcher_re, searcher_string

PY3: bool
text_type: type

class _NullCoder:
    @staticmethod
    def encode(b: str, final: bool = False): ...
    @staticmethod
    def decode(b: str, final: bool = False): ...

@type_check_only
class _Logfile(Protocol):
    def write(self, s, /) -> object: ...
    def flush(self) -> object: ...

_ErrorPattern: TypeAlias = type[EOF | TIMEOUT]
_InputStringPattern: TypeAlias = str | bytes | _ErrorPattern
_InputRePattern: TypeAlias = Pattern[str] | Pattern[bytes] | _InputStringPattern
_CompiledStringPattern: TypeAlias = AnyStr | _ErrorPattern
_CompiledRePattern: TypeAlias = Pattern[AnyStr] | _ErrorPattern
_Searcher: TypeAlias = searcher_string[AnyStr] | searcher_re[AnyStr]

class SpawnBase(Generic[AnyStr]):
    encoding: str | None
    pid: int | None
    flag_eof: bool
    stdin: TextIO
    stdout: TextIO
    stderr: TextIO
    searcher: None
    ignorecase: bool
    before: AnyStr | None
    after: _CompiledStringPattern[AnyStr] | None
    match: AnyStr | Match[AnyStr] | _ErrorPattern | None
    match_index: int | None
    terminated: bool
    exitstatus: int | None
    signalstatus: int | None
    status: int | None
    child_fd: int
    timeout: float | None
    delimiter: type[EOF]
    logfile: _Logfile | None
    logfile_read: _Logfile | None
    logfile_send: _Logfile | None
    maxread: int
    searchwindowsize: int | None
    delaybeforesend: float | None
    delayafterclose: float
    delayafterterminate: float
    delayafterread: float | None
    softspace: bool
    name: str
    closed: bool
    codec_errors: str
    string_type: type[AnyStr]
    buffer_type: IO[AnyStr]
    crlf: AnyStr
    allowed_string_types: tuple[type, ...]
    linesep: AnyStr
    write_to_stdout: Callable[[AnyStr], int]
    async_pw_transport: tuple[PatternWaiter[AnyStr], ReadTransport] | None
    def __init__(
        self,
        timeout: float | None = 30,
        maxread: int = 2000,
        searchwindowsize: int | None = None,
        logfile: _Logfile | None = None,
        encoding: str | None = None,
        codec_errors: str = "strict",
    ) -> None: ...
    @property
    def buffer(self) -> AnyStr: ...
    @buffer.setter
    def buffer(self, value: AnyStr) -> None: ...
    def read_nonblocking(self, size: int = 1, timeout: float | None = None) -> AnyStr: ...
    def compile_pattern_list(self, patterns: _InputRePattern | list[_InputRePattern]) -> list[_CompiledRePattern[AnyStr]]: ...
    @overload
    def expect(
        self,
        pattern: _InputRePattern | list[_InputRePattern],
        timeout: float | None = -1,
        searchwindowsize: int | None = -1,
        async_: Literal[False] = False,
    ) -> int: ...
    @overload
    def expect(
        self,
        pattern: _InputRePattern | list[_InputRePattern],
        timeout: float | None = -1,
        searchwindowsize: int | None = -1,
        *,
        async_: Literal[True],
    ) -> Awaitable[int]: ...
    @overload
    def expect_list(
        self,
        pattern_list: list[_CompiledRePattern[AnyStr]],
        timeout: float | None = -1,
        searchwindowsize: int | None = -1,
        async_: Literal[False] = False,
    ) -> int: ...
    @overload
    def expect_list(
        self,
        pattern_list: list[_CompiledRePattern[AnyStr]],
        timeout: float | None = -1,
        searchwindowsize: int | None = -1,
        *,
        async_: Literal[True],
    ) -> Awaitable[int]: ...
    @overload
    def expect_exact(
        self,
        pattern_list: _InputStringPattern | Iterable[_InputStringPattern],
        timeout: float | None = -1,
        searchwindowsize: int | None = -1,
        async_: Literal[False] = False,
    ) -> int: ...
    @overload
    def expect_exact(
        self,
        pattern_list: _InputStringPattern | Iterable[_InputStringPattern],
        timeout: float | None = -1,
        searchwindowsize: int | None = -1,
        *,
        async_: Literal[True],
    ) -> Awaitable[int]: ...
    def expect_loop(self, searcher: _Searcher[AnyStr], timeout: float | None = -1, searchwindowsize: int | None = -1) -> int: ...
    def read(self, size: int = -1) -> AnyStr: ...
    def readline(self, size: int = -1) -> AnyStr: ...
    def __iter__(self): ...
    def readlines(self, sizehint: int = -1) -> list[AnyStr]: ...
    def fileno(self) -> int: ...
    def flush(self) -> None: ...
    def isatty(self) -> bool: ...
    def __enter__(self): ...
    def __exit__(self, etype, evalue, tb) -> None: ...
