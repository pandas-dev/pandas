from typing import Literal
from typing_extensions import Self

from gunicorn.asgi.unreader import AsyncUnreader
from gunicorn.config import Config
from gunicorn.uwsgi.message import UWSGIRequest

from .._types import _AddressType

class AsyncUWSGIRequest(UWSGIRequest):
    cfg: Config
    unreader: AsyncUnreader  # type: ignore[assignment]
    peer_addr: _AddressType
    remote_addr: _AddressType
    req_number: int
    method: str | None
    uri: str | None
    path: str | None
    query: str | None
    fragment: str | None
    version: tuple[int, int]
    headers: list[tuple[str, str]]
    trailers: list[tuple[str, str]]
    scheme: Literal["https", "http"]
    must_close: bool
    uwsgi_vars: dict[str, str]
    modifier1: int
    modifier2: int
    proxy_protocol_info: dict[str, str | int | None] | None  # TODO: Use TypedDict
    content_length: int
    chunked: bool

    def __init__(self, cfg: Config, unreader: AsyncUnreader, peer_addr: _AddressType, req_number: int = 1) -> None: ...
    @classmethod
    async def parse(cls, cfg: Config, unreader: AsyncUnreader, peer_addr: _AddressType, req_number: int = 1) -> Self: ...  # type: ignore[override]
    async def read_body(self, size: int = 8192) -> bytes: ...
    async def drain_body(self) -> None: ...
    def get_header(self, name: str) -> str | None: ...
