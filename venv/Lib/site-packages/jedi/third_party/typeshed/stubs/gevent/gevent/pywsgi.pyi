from _typeshed import OptExcInfo, StrOrBytesPath, SupportsWrite
from _typeshed.wsgi import WSGIApplication, WSGIEnvironment
from collections.abc import Callable, Container, Iterable, Iterator
from http.client import HTTPMessage
from io import BufferedIOBase, BufferedReader
from logging import Logger
from types import TracebackType
from typing import Any, ClassVar, Literal, Protocol, TypeVar, overload, type_check_only
from typing_extensions import Self

from gevent.baseserver import _Spawner
from gevent.server import StreamServer
from gevent.socket import socket as _GeventSocket
from gevent.ssl import SSLContext

__all__ = ["WSGIServer", "WSGIHandler", "LoggingLogAdapter", "Environ", "SecureEnviron", "WSGISecureEnviron"]

_T = TypeVar("_T")

@type_check_only
class _LogOutputStream(SupportsWrite[str], Protocol):
    def writelines(self, lines: Iterable[str], /) -> None: ...
    def flush(self) -> None: ...

class Input:
    __slots__ = (
        "rfile",
        "content_length",
        "socket",
        "position",
        "chunked_input",
        "chunk_length",
        "_chunked_input_error",
        "send_100_continue_enabled",
    )
    rfile: BufferedReader
    content_length: int | None
    socket: _GeventSocket | None
    position: int
    chunked_input: bool
    chunk_length: int
    send_100_continue_enabled: bool
    def __init__(
        self, rfile: BufferedReader, content_length: int | None, socket: _GeventSocket | None = None, chunked_input: bool = False
    ) -> None: ...
    def read(self, length: int | None = None) -> bytes: ...
    def readline(self, size: int | None = None) -> bytes: ...
    def readlines(self, hint: object | None = None) -> list[bytes]: ...
    def __iter__(self) -> Self: ...
    def next(self) -> bytes: ...
    __next__ = next

class OldMessage(HTTPMessage):
    status: str
    def __init__(self) -> None: ...
    @overload
    def getheader(self, name: str, default: None = None) -> str | None: ...
    @overload
    def getheader(self, name: str, default: _T) -> str | _T: ...
    @property
    def headers(self) -> Iterator[str]: ...
    @property
    def typeheader(self) -> str | None: ...

class WSGIHandler:
    protocol_version: str
    def MessageClass(self, fp: BufferedIOBase) -> OldMessage: ...
    status: str | None
    response_headers: list[tuple[str, str]] | None
    code: int | None
    provided_date: str | None
    provided_content_length: str | None
    close_connection: bool
    time_start: float
    time_finish: float
    headers_sent: bool
    response_use_chunked: bool
    connection_upgraded: bool
    environ: WSGIEnvironment | None
    application: WSGIApplication | None
    requestline: str | None
    response_length: int
    result: Iterable[bytes] | None
    wsgi_input: Input | None
    content_length: int
    headers: OldMessage
    request_version: str | None
    command: str | None
    path: str | None
    socket: _GeventSocket
    client_address: str
    server: WSGIServer
    rfile: BufferedReader
    def __init__(self, sock: _GeventSocket, address: str, server: WSGIServer) -> None: ...
    def handle(self) -> None: ...
    def read_request(self, raw_requestline: str) -> OldMessage: ...
    def log_error(self, msg: str, *args: object) -> None: ...
    def read_requestline(self) -> str: ...
    def handle_one_request(self) -> tuple[str, bytes] | Literal[True] | None: ...
    def finalize_headers(self) -> None: ...
    ApplicationError: type[AssertionError]
    def write(self, data: bytes) -> None: ...
    def start_response(
        self, status: str, headers: list[tuple[str, str]], exc_info: OptExcInfo | None = None
    ) -> Callable[[bytes], None]: ...
    def log_request(self) -> None: ...
    def format_request(self) -> str: ...
    def process_result(self) -> None: ...
    def run_application(self) -> None: ...
    ignored_socket_errors: tuple[int, ...]
    def handle_one_response(self) -> None: ...
    def handle_error(self, t: type[BaseException] | None, v: BaseException | None, tb: TracebackType | None) -> None: ...
    def get_environ(self) -> WSGIEnvironment: ...

class LoggingLogAdapter:
    __slots__ = ("_logger", "_level")
    def __init__(self, logger: Logger, level: int = 20) -> None: ...
    def write(self, msg: str) -> None: ...
    def flush(self) -> None: ...
    def writelines(self, lines: Iterable[str]) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def __setattr__(self, name: str, value: object) -> None: ...
    def __delattr__(self, name: str) -> None: ...

class Environ(WSGIEnvironment):
    __slots__ = ()

class SecureEnviron(Environ):
    __slots__ = ("secure_repr", "whitelist_keys", "print_masked_keys")
    default_secure_repr: ClassVar[bool]
    default_whitelist_keys: ClassVar[Container[str]]
    default_print_masked_keys: ClassVar[bool]
    secure_repr: bool
    whitelist_keys: Container[str]
    print_masked_keys: bool

class WSGISecureEnviron(SecureEnviron): ...

class WSGIServer(StreamServer):
    handler_class: type[WSGIHandler]
    log: _LogOutputStream
    error_log: _LogOutputStream
    environ_class: type[WSGIEnvironment]
    secure_environ_class: type[SecureEnviron]
    base_env: WSGIEnvironment
    application: WSGIApplication
    @overload
    def __init__(
        self,
        listener: _GeventSocket | tuple[str, int] | str,
        application: WSGIApplication | None = None,
        backlog: int | None = None,
        spawn: _Spawner = "default",
        log: str | Logger | _LogOutputStream | None = "default",
        error_log: str | Logger | _LogOutputStream | None = "default",
        handler_class: type[WSGIHandler] | None = None,
        environ: WSGIEnvironment | None = None,
        *,
        ssl_context: SSLContext,
        server_side: bool = True,
        do_handshake_on_connect: bool = True,
        suppress_ragged_eofs: bool = True,
    ) -> None: ...
    @overload
    def __init__(
        self,
        listener: _GeventSocket | tuple[str, int] | str,
        application: WSGIApplication | None = None,
        backlog: int | None = None,
        spawn: _Spawner = "default",
        log: str | Logger | _LogOutputStream | None = "default",
        error_log: str | Logger | _LogOutputStream | None = "default",
        handler_class: type[WSGIHandler] | None = None,
        environ: WSGIEnvironment | None = None,
        *,
        keyfile: StrOrBytesPath = ...,
        certfile: StrOrBytesPath = ...,
        server_side: bool = True,
        cert_reqs: int = ...,
        ssl_version: int = ...,
        ca_certs: str = ...,
        do_handshake_on_connect: bool = True,
        suppress_ragged_eofs: bool = True,
        ciphers: str = ...,
    ) -> None: ...
    environ: WSGIEnvironment
    def set_environ(self, environ: WSGIEnvironment | None = None) -> None: ...
    max_accept: int
    def set_max_accept(self) -> None: ...
    def get_environ(self) -> WSGIEnvironment: ...
    def init_socket(self) -> None: ...
    def update_environ(self) -> None: ...
    def handle(self, sock: _GeventSocket, address: str) -> None: ...
