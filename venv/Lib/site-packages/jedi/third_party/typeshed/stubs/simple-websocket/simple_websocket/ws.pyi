import socket
import threading
from _typeshed import FileDescriptorLike
from _typeshed.wsgi import WSGIEnvironment
from collections.abc import Callable
from selectors import SelectorKey, _EventMask
from ssl import SSLContext
from typing import Any, Protocol, type_check_only

from wsproto import ConnectionType, WSConnection
from wsproto.events import Request
from wsproto.frame_protocol import CloseReason

@type_check_only
class _ThreadClassProtocol(Protocol):
    name: str
    # this accepts any callable as the target, like `threading.Thread`
    def __init__(self, target: Callable[..., Any]) -> None: ...
    def start(self) -> None: ...

@type_check_only
class _EventClassProtocol(Protocol):
    def clear(self) -> None: ...
    def set(self) -> None: ...
    def wait(self, timeout: float | None = None) -> bool: ...

@type_check_only
class _SelectorClassProtocol(Protocol):
    # the signature of `register` here is the same as `selectors._BaseSelectorImpl` from the stdlib
    def register(self, fileobj: FileDescriptorLike, events: _EventMask, data: Any = None) -> SelectorKey: ...
    # the signature of `select` here is the same as `selectors.DefaultSelector` from the stdlib
    def select(self, timeout: float | None = None) -> list[tuple[SelectorKey, _EventMask]]: ...
    def close(self) -> None: ...

class Base:
    subprotocol: str | None
    sock: socket.socket | None
    receive_bytes: int
    ping_interval: float | None
    max_message_size: int | None
    pong_received: bool
    input_buffer: list[bytes | str]
    incoming_message: bytes | str | None
    incoming_message_len: int
    connected: bool
    is_server: bool
    close_reason: CloseReason
    close_message: str | None
    selector_class: type[_SelectorClassProtocol]
    event: _EventClassProtocol | threading.Event
    ws: WSConnection
    thread: _ThreadClassProtocol | threading.Thread
    def __init__(
        self,
        sock: socket.socket | None = None,
        connection_type: ConnectionType | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
        thread_class: type[_ThreadClassProtocol] | None = None,
        event_class: type[_EventClassProtocol] | None = None,
        selector_class: type[_SelectorClassProtocol] | None = None,
    ) -> None: ...
    def handshake(self) -> None: ...
    # data can be antyhing. a special case is made for `bytes`, anything else is converted to `str`.
    def send(self, data: bytes | Any) -> None: ...
    def receive(self, timeout: float | None = None) -> bytes | str | None: ...
    def close(self, reason: CloseReason | None = None, message: str | None = None) -> None: ...
    def choose_subprotocol(self, request: Request) -> str | None: ...

class Server(Base):
    environ: WSGIEnvironment
    subprotocols: list[str]
    mode: str
    connected: bool
    def __init__(
        self,
        environ: WSGIEnvironment,
        subprotocols: list[str] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
        thread_class: type[_ThreadClassProtocol] | None = None,
        event_class: type[_EventClassProtocol] | None = None,
        selector_class: type[_SelectorClassProtocol] | None = None,
    ) -> None: ...
    @classmethod
    def accept(
        cls,
        environ: WSGIEnvironment,
        subprotocols: list[str] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
        thread_class: type[_ThreadClassProtocol] | None = None,
        event_class: type[_EventClassProtocol] | None = None,
        selector_class: type[_SelectorClassProtocol] | None = None,
    ) -> Server: ...
    def handshake(self) -> None: ...
    def choose_subprotocol(self, request: Request) -> str | None: ...

class Client(Base):
    host: str
    port: int
    path: str
    subprotocols: list[str]
    extra_headeers: list[tuple[bytes, bytes]]
    subprotocol: str | None
    connected: bool
    def __init__(
        self,
        url: str,
        subprotocols: list[str] | None = None,
        headers: dict[bytes, bytes] | list[tuple[bytes, bytes]] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
        ssl_context: SSLContext | None = None,
        thread_class: type[_ThreadClassProtocol] | None = None,
        event_class: type[_EventClassProtocol] | None = None,
    ) -> None: ...
    @classmethod
    def connect(
        cls,
        url: str,
        subprotocols: list[str] | None = None,
        headers: dict[bytes, bytes] | list[tuple[bytes, bytes]] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
        ssl_context: SSLContext | None = None,
        thread_class: type[_ThreadClassProtocol] | None = None,
        event_class: type[_EventClassProtocol] | None = None,
    ) -> Client: ...
    def handshake(self) -> None: ...
    def close(self, reason: CloseReason | None = None, message: str | None = None) -> None: ...
