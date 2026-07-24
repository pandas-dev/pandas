from _typeshed import Incomplete
from collections.abc import Iterable
from ssl import SSLSocket
from typing import Any, ClassVar
from typing_extensions import TypeAlias

from gunicorn.config import Config
from gunicorn.http2.request import HTTP2Request
from gunicorn.http2.stream import HTTP2Stream

from .._types import _AddressType

_H2Connection: TypeAlias = Any  # h2.connection.H2Connection class

class HTTP2ServerConnection:
    READ_BUFFER_SIZE: ClassVar[int]
    cfg: Config
    sock: SSLSocket
    client_addr: _AddressType
    streams: dict[int, HTTP2Stream]
    initial_window_size: int
    max_concurrent_streams: int
    max_frame_size: int
    max_header_list_size: int
    h2_conn: _H2Connection

    def __init__(self, cfg: Config, sock: SSLSocket, client_addr: _AddressType) -> None: ...
    def initiate_connection(self) -> None: ...
    def receive_data(self, data: bytes | None = None) -> list[HTTP2Request]: ...
    def send_informational(self, stream_id: int, status: int, headers: Iterable[tuple[str, Incomplete]]) -> None: ...
    def send_response(
        self, stream_id: int, status: int, headers: Iterable[tuple[str, Incomplete]], body: bytes | None = None
    ) -> bool: ...
    def send_data(self, stream_id: int, data: bytes, end_stream: bool = False) -> bool: ...
    def send_trailers(self, stream_id: int, trailers: Iterable[tuple[str, Incomplete]]) -> bool: ...
    def send_error(self, stream_id: int, status_code: int, message: str | None = None) -> None: ...
    def reset_stream(self, stream_id: int, error_code: int = 0x8) -> None: ...
    def close(self, error_code: int = 0x0, last_stream_id: int | None = None) -> None: ...
    @property
    def is_closed(self) -> bool: ...
    def cleanup_stream(self, stream_id: int) -> None: ...

__all__ = ["HTTP2ServerConnection"]
