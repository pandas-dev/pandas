from _typeshed.wsgi import WSGIEnvironment
from collections.abc import Awaitable, Callable
from typing import TypedDict, type_check_only

@type_check_only
class _SocketDataBase(TypedDict):
    type: str

@type_check_only
class _SocketDataProtocol(_SocketDataBase):
    subprotocol: str | None

@type_check_only
class _SocketDataStr(_SocketDataBase):
    text: str

@type_check_only
class _SocketDataBytes(_SocketDataBase):
    bytes: bytes

class WebSocketASGI:
    subprotocols: list[str]
    subprotocol: str
    connected: bool
    # this is set in `close` to `False`
    conncted: bool
    def __init__(
        self,
        scope: WSGIEnvironment,
        receive: Callable[[], Awaitable[_SocketDataBytes | _SocketDataStr]],
        send: Callable[[_SocketDataBase | _SocketDataProtocol | _SocketDataBytes | _SocketDataStr], Awaitable[None]],
        subprotocols: list[str] | None = None,
    ) -> None: ...
    @classmethod
    async def accept(
        cls,
        scope: WSGIEnvironment,
        receive: Callable[[], Awaitable[_SocketDataBytes | _SocketDataStr]],
        send: Callable[[_SocketDataBase | _SocketDataProtocol | _SocketDataBytes | _SocketDataStr], Awaitable[None]],
        subprotocols: list[str] | None = None,
    ) -> WebSocketASGI: ...
    async def receive(self) -> bytes | str: ...
    async def send(self, data: bytes | str) -> None: ...
    async def close(self) -> None: ...
