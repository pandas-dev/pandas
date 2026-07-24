from _socket import _Address as _StrictAddress
from _typeshed import ReadableBuffer, StrOrBytesPath
from collections.abc import Callable
from typing import Any, ClassVar, TypedDict, overload, type_check_only
from typing_extensions import TypeAlias

from gevent.baseserver import BaseServer, _Spawner
from gevent.socket import socket as _GeventSocket
from gevent.ssl import SSLContext, wrap_socket as ssl_wrap_socket

# For simplicity we treat _Address as Any, we could be more strict and use the definition
# from the stdlib _socket.pyi. But that would exclude some potentially valid handlers.
_Address: TypeAlias = Any

@type_check_only
class _SSLArguments(TypedDict, total=False):
    keyfile: StrOrBytesPath
    certfile: StrOrBytesPath
    server_side: bool
    cert_reqs: int
    ssl_version: int
    ca_certs: str
    suppress_ragged_eofs: bool
    do_handshake_on_connect: bool
    ciphers: str

class StreamServer(BaseServer[_GeventSocket, _Address]):
    backlog: int
    reuse_addr: ClassVar[int | None]
    wrap_socket = ssl_wrap_socket
    ssl_args: _SSLArguments | None
    @overload
    def __init__(
        self,
        listener: _GeventSocket | tuple[str, int] | str,
        handle: Callable[[_GeventSocket, _Address], object] | None = None,
        backlog: int | None = None,
        spawn: _Spawner = "default",
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
        handle: Callable[[_GeventSocket, _Address], object] | None = None,
        backlog: int | None = None,
        spawn: _Spawner = "default",
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
    @property
    def ssl_enabled(self) -> bool: ...
    @classmethod
    def get_listener(cls, address: _StrictAddress, backlog: int | None = None, family: int | None = None) -> _GeventSocket: ...
    def do_read(self) -> tuple[_GeventSocket, _Address]: ...
    def do_close(self, sock: _GeventSocket, address: _Address) -> None: ...
    def wrap_socket_and_handle(self, client_socket: _GeventSocket, address: _StrictAddress) -> Any: ...

class DatagramServer(BaseServer[_GeventSocket, _Address]):
    reuse_addr: ClassVar[int | None]
    def __init__(
        self,
        listener: _GeventSocket | tuple[str, int] | str,
        handle: Callable[[_GeventSocket, _Address], object] | None = None,
        spawn: _Spawner = "default",
    ) -> None: ...
    @classmethod
    def get_listener(cls, address: _StrictAddress, family: int | None = None) -> _GeventSocket: ...
    def do_read(self) -> tuple[_GeventSocket, _Address]: ...
    @overload
    def sendto(self, data: ReadableBuffer, address: _StrictAddress, /) -> int: ...
    @overload
    def sendto(self, data: ReadableBuffer, flags: int, address: _StrictAddress, /) -> int: ...

__all__ = ["StreamServer", "DatagramServer"]
