from typing import Final, Literal

from gunicorn.config import Config
from gunicorn.http.body import Body
from gunicorn.http.unreader import Unreader

from .._types import _AddressType

MAX_UWSGI_VARS: Final = 1000

class UWSGIRequest:
    cfg: Config
    unreader: Unreader
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
    body: Body | None
    scheme: Literal["https", "http"]
    must_close: bool
    uwsgi_vars: dict[str, str]
    modifier1: int
    modifier2: int
    proxy_protocol_info: dict[str, str | int | None] | None  # TODO: Use TypedDict

    def __init__(self, cfg: Config, unreader: Unreader, peer_addr: _AddressType, req_number: int = 1) -> None: ...
    def force_close(self) -> None: ...
    def parse(self, unreader: Unreader) -> bytes: ...
    def set_body_reader(self) -> None: ...
    def should_close(self) -> bool: ...
