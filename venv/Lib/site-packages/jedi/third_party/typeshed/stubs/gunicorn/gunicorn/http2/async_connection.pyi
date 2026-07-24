from _typeshed import Incomplete
from asyncio import StreamReader, StreamWriter
from collections.abc import Iterable
from typing import ClassVar

from gunicorn.config import Config
from gunicorn.http2.connection import _H2Connection
from gunicorn.http2.request import HTTP2Request
from gunicorn.http2.stream import HTTP2Stream

from .._types import _AddressType

class AsyncHTTP2Connection:
    READ_BUFFER_SIZE: ClassVar[int]
    cfg: Config
    reader: StreamReader
    writer: StreamWriter
    client_addr: _AddressType
    streams: dict[int, HTTP2Stream]
    initial_window_size: int
    max_concurrent_streams: int
    max_frame_size: int
    max_header_list_size: int
    h2_conn: _H2Connection

    def __init__(self, cfg: Config, reader: StreamReader, writer: StreamWriter, client_addr: _AddressType) -> None: ...
    async def initiate_connection(self) -> None: ...
    async def receive_data(self, timeout: float | None = None) -> list[HTTP2Request]: ...
    async def send_informational(self, stream_id: int, status: int, headers: Iterable[tuple[str, Incomplete]]) -> None: ...
    async def send_response(
        self, stream_id: int, status: int, headers: Iterable[tuple[str, Incomplete]], body: bytes | None = None
    ) -> bool: ...
    async def send_data(self, stream_id: int, data: bytes, end_stream: bool = False) -> bool: ...
    async def send_trailers(self, stream_id: int, trailers: Iterable[tuple[str, Incomplete]]) -> bool: ...
    async def send_error(self, stream_id: int, status_code: int, message: str | None = None) -> None: ...
    async def reset_stream(self, stream_id: int, error_code: int = 0x8) -> None: ...
    async def close(self, error_code: int = 0x0, last_stream_id: int | None = None) -> None: ...
    @property
    def is_closed(self) -> bool: ...
    def cleanup_stream(self, stream_id: int) -> None: ...

__all__ = ["AsyncHTTP2Connection"]
