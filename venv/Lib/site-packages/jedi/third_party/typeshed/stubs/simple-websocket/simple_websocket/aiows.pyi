import asyncio
import socket
from _typeshed import Incomplete, Unused
from _typeshed.wsgi import WSGIEnvironment
from collections.abc import Awaitable, Callable
from ssl import SSLContext
from typing import Any, Literal, TypedDict, type_check_only

from wsproto import ConnectionType, WSConnection
from wsproto.events import Request
from wsproto.frame_protocol import CloseReason

from .asgi import WebSocketASGI, _SocketDataBase, _SocketDataBytes, _SocketDataProtocol, _SocketDataStr

class AioBase:
    subprotocol: str | None
    connection_type: ConnectionType
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
    close_message: str
    rsock: asyncio.StreamReader
    wsock: asyncio.StreamWriter
    event: asyncio.Event
    ws: WSConnection | None
    task: asyncio.Task[None]
    def __init__(
        self,
        connection_type: ConnectionType | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
    ) -> None: ...
    async def connect(self) -> None: ...
    async def handshake(self) -> None: ...
    # data can be antyhing. a special case is made for `bytes`, anything else is converted to `str`.
    async def send(self, data: bytes | Any) -> None: ...
    async def receive(self, timeout: float | None = None) -> bytes | str | None: ...
    async def close(self, reason: CloseReason | None = None, message: str | None = None) -> None: ...
    def choose_subprotocol(self, request: Request) -> str | None: ...

@type_check_only
class _AioServerRequest(TypedDict):
    # this is `aiohttp.web.Request`
    aiohttp: Incomplete
    sock: None
    headers: None

class AioServer(AioBase):
    request: _AioServerRequest
    headers: dict[str, Any]
    subprotocols: list[str]
    is_server: Literal[True]
    mode: str
    connected: bool
    def __init__(
        self,
        request: _AioServerRequest,
        subprotocols: list[str] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
    ) -> None: ...
    @classmethod
    async def accept(
        cls,
        # this is `aiohttp.web.Request`
        aiohttp=None,
        asgi: (
            tuple[
                WSGIEnvironment,
                Callable[[], Awaitable[_SocketDataBytes | _SocketDataStr]],
                Callable[[_SocketDataBase | _SocketDataProtocol | _SocketDataBytes | _SocketDataStr], Awaitable[None]],
            ]
            | None
        ) = None,
        sock: socket.socket | None = None,
        headers: dict[str, Any] | None = None,
        subprotocols: list[str] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
    ) -> WebSocketASGI | AioServer: ...
    async def handshake(self) -> None: ...
    def choose_subprotocol(self, request: Request) -> str | None: ...

class AioClient(AioBase):
    url: str
    ssl_context: SSLContext | None
    is_secure: bool
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
        headers: dict[str, Any] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
        ssl_context: SSLContext | None = None,
    ) -> None: ...
    # the source code itself has this override
    @classmethod
    async def connect(  # type: ignore[override]
        cls,
        url: str,
        subprotocols: list[str] | None = None,
        headers: dict[str, Any] | None = None,
        receive_bytes: int = 4096,
        ping_interval: float | None = None,
        max_message_size: int | None = None,
        ssl_context: SSLContext | None = None,
        thread_class: Unused = None,
        event_class: Unused = None,
    ) -> AioClient: ...
    async def handshake(self) -> None: ...
    async def close(self, reason: CloseReason | None = None, message: str | None = None) -> None: ...
