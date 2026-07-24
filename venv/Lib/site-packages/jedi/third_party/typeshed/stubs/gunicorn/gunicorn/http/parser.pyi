import socket
from collections.abc import Callable, Iterable, Iterator
from typing import ClassVar

from gunicorn.config import Config
from gunicorn.http.message import Request
from gunicorn.http.unreader import Unreader

from .._types import _AddressType

class Parser:
    # TODO: Use Protocol instead of Request class
    mesg_class: ClassVar[type[Request] | None]
    cfg: Config
    unreader: Unreader
    mesg: Request | None
    source_addr: _AddressType
    req_count: int

    def __init__(self, cfg: Config, source: socket.socket | Iterable[bytes], source_addr: _AddressType) -> None: ...
    def __iter__(self) -> Iterator[Request]: ...
    def finish_body(self) -> None: ...
    def __next__(self) -> Request: ...

    next: Callable[[Parser], Request]

class RequestParser(Parser):
    mesg_class: ClassVar[type[Request]]
