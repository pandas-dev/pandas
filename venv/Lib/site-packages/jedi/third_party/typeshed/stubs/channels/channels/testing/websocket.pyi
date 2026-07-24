from collections.abc import Iterable
from typing import Any, Literal, TypedDict, overload, type_check_only
from typing_extensions import NotRequired, TypeAlias

from asgiref.typing import ASGIVersions
from channels.testing.application import ApplicationCommunicator
from channels.utils import _ChannelApplication

@type_check_only
class _WebsocketTestScope(TypedDict, total=False):
    spec_version: int
    type: Literal["websocket"]
    asgi: ASGIVersions
    http_version: str
    scheme: str
    path: str
    raw_path: bytes
    query_string: bytes
    root_path: str
    headers: Iterable[tuple[bytes, bytes]] | None
    client: tuple[str, int] | None
    server: tuple[str, int | None] | None
    subprotocols: Iterable[str] | None
    state: NotRequired[dict[str, Any]]
    extensions: dict[str, dict[object, object]] | None

_Connected: TypeAlias = bool
_CloseCodeOrAcceptSubProtocol: TypeAlias = int | str | None
_WebsocketConnectResponse: TypeAlias = tuple[_Connected, _CloseCodeOrAcceptSubProtocol]

class WebsocketCommunicator(ApplicationCommunicator):
    scope: _WebsocketTestScope
    response_headers: list[tuple[bytes, bytes]] | None

    def __init__(
        self,
        application: _ChannelApplication,
        path: str,
        headers: Iterable[tuple[bytes, bytes]] | None = None,
        subprotocols: Iterable[str] | None = None,
        spec_version: int | None = None,
    ) -> None: ...
    async def connect(self, timeout: float = 1) -> _WebsocketConnectResponse: ...
    async def send_to(self, text_data: str | None = None, bytes_data: bytes | None = None) -> None: ...
    async def receive_from(self, timeout: float = 1) -> str | bytes: ...

    # These overloads reflect common usage, where users typically send and receive `dict[str, Any]`.
    # The base case allows `Any` to support broader `json.dumps` / `json.loads` compatibility.
    @overload
    async def send_json_to(self, data: dict[str, Any]) -> None: ...
    @overload
    async def send_json_to(self, data: Any) -> None: ...
    async def receive_json_from(self, timeout: float = 1) -> Any: ...
    async def disconnect(self, code: int = 1000, timeout: float = 1) -> None: ...
