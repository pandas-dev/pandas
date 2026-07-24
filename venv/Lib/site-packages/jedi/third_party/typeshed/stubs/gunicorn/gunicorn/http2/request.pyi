from _typeshed import Incomplete, ReadableBuffer
from collections.abc import Iterator
from typing import Literal

from gunicorn.config import Config
from gunicorn.http2.stream import HTTP2Stream

from .._types import _AddressType

class HTTP2Body:
    def __init__(self, data: ReadableBuffer) -> None: ...
    def read(self, size: int | None = None) -> bytes: ...
    def readline(self, size: int | None = None) -> bytes: ...
    def readlines(self, hint: int | None = None) -> list[bytes]: ...
    def __iter__(self) -> Iterator[bytes]: ...
    def __len__(self) -> int: ...
    def close(self) -> None: ...

class HTTP2Request:
    stream: HTTP2Stream
    cfg: Config
    peer_addr: _AddressType
    remote_addr: _AddressType
    version: tuple[int, int]
    method: str
    scheme: Literal["https", "http"]
    uri: str
    path: str
    query: str
    fragment: str
    headers: list[tuple[str, str]]
    trailers: list[tuple[str, Incomplete]]
    body: HTTP2Body
    must_close: bool
    req_number: int
    proxy_protocol_info: dict[str, str | int | None] | None  # TODO: Use TypedDict
    priority_weight: int
    priority_depends_on: int

    def __init__(self, stream: HTTP2Stream, cfg: Config, peer_addr: _AddressType) -> None: ...
    def force_close(self) -> None: ...
    def should_close(self) -> bool: ...
    def get_header(self, name: str) -> str | None: ...
    @property
    def content_length(self) -> int | None: ...
    @property
    def content_type(self) -> str | None: ...

__all__ = ["HTTP2Request", "HTTP2Body"]
