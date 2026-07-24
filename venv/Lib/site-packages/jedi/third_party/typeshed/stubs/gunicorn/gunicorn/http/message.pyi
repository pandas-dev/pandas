import io
import re
from enum import IntEnum
from typing import Final

from gunicorn.config import Config
from gunicorn.http.body import Body
from gunicorn.http.unreader import Unreader

from .._types import _AddressType

PP_V2_SIGNATURE: Final = b"\x0d\x0a\x0d\x0a\x00\x0d\x0a\x51\x55\x49\x54\x0a"

class PPCommand(IntEnum):
    LOCAL = 0x0
    PROXY = 0x1

class PPFamily(IntEnum):
    UNSPEC = 0x0
    INET = 0x1
    INET6 = 0x2
    UNIX = 0x3

class PPProtocol(IntEnum):
    UNSPEC = 0x0
    STREAM = 0x1
    DGRAM = 0x2

MAX_REQUEST_LINE: Final = 8190
MAX_HEADERS: Final = 32768
DEFAULT_MAX_HEADERFIELD_SIZE: Final = 8190
RFC9110_5_6_2_TOKEN_SPECIALS: Final = r"!#$%&'*+-.^_`|~"
TOKEN_RE: Final[re.Pattern[str]]
METHOD_BADCHAR_RE: Final[re.Pattern[str]]
VERSION_RE: Final[re.Pattern[str]]
RFC9110_5_5_INVALID_AND_DANGEROUS: Final[re.Pattern[str]]

class Message:
    cfg: Config
    unreader: Unreader
    peer_addr: _AddressType
    remote_addr: _AddressType
    version: tuple[int, int] | None
    headers: list[tuple[str, str]]
    trailers: list[tuple[str, str]]
    body: Body | None
    scheme: str
    must_close: bool
    limit_request_fields: int
    limit_request_field_size: int
    max_buffer_headers: int

    def __init__(self, cfg: Config, unreader: Unreader, peer_addr: _AddressType) -> None: ...
    def force_close(self) -> None: ...
    def parse(self, unreader: Unreader) -> bytes: ...
    def parse_headers(self, data: bytes, from_trailer: bool = False) -> list[tuple[str, str]]: ...
    def set_body_reader(self) -> None: ...
    def should_close(self) -> bool: ...

class Request(Message):
    method: str | None
    uri: str | None
    path: str | None
    query: str | None
    fragment: str | None
    limit_request_line: int
    req_number: int
    proxy_protocol_info: dict[str, str | int | None] | None  # TODO: Use TypedDict

    def __init__(self, cfg: Config, unreader: Unreader, peer_addr: _AddressType, req_number: int = 1) -> None: ...
    def get_data(self, unreader: Unreader, buf: io.BytesIO, stop: bool = False) -> None: ...
    def parse(self, unreader: Unreader) -> bytes: ...
    def read_into(self, unreader: Unreader, buf: bytearray, stop: bool | None = False) -> None: ...
    def read_line(self, unreader: Unreader, buf: bytearray, limit: int = 0) -> tuple[bytes, bytearray]: ...
    def read_bytes(self, unreader: Unreader, buf: bytearray, count: int) -> tuple[bytes, bytearray]: ...
    def proxy_protocol_access_check(self) -> None: ...
    def parse_request_line(self, line_bytes: bytes) -> None: ...
    def set_body_reader(self) -> None: ...
