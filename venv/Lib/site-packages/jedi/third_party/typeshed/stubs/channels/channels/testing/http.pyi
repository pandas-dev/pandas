from collections.abc import Iterable
from typing import Literal, TypedDict, type_check_only

from channels.testing.application import ApplicationCommunicator
from channels.utils import _ChannelApplication

# HTTP test-specific response type
@type_check_only
class _HTTPTestResponse(TypedDict, total=False):
    status: int
    headers: Iterable[tuple[bytes, bytes]]
    body: bytes

@type_check_only
class _HTTPTestScope(TypedDict, total=False):
    type: Literal["http"]
    http_version: str
    method: str
    scheme: str
    path: str
    raw_path: bytes
    query_string: bytes
    root_path: str
    headers: Iterable[tuple[bytes, bytes]] | None
    client: tuple[str, int] | None
    server: tuple[str, int | None] | None

class HttpCommunicator(ApplicationCommunicator):
    scope: _HTTPTestScope
    body: bytes
    sent_request: bool

    def __init__(
        self,
        application: _ChannelApplication,
        method: str,
        path: str,
        body: bytes = b"",
        headers: Iterable[tuple[bytes, bytes]] | None = None,
    ) -> None: ...
    async def get_response(self, timeout: float = 1) -> _HTTPTestResponse: ...
