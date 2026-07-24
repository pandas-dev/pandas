import socket
from collections.abc import Iterable
from typing import Literal, overload

from gunicorn.config import Config
from gunicorn.http.message import Message as Message, Request as Request
from gunicorn.http.parser import RequestParser as RequestParser
from gunicorn.http2.connection import HTTP2ServerConnection
from gunicorn.uwsgi.parser import UWSGIParser

from .._types import _AddressType

@overload
def get_parser(
    cfg: Config,
    source: socket.socket | Iterable[bytes],
    source_addr: _AddressType,
    http2_connection: Literal[False] | None = False,
) -> UWSGIParser | RequestParser: ...
@overload
def get_parser(
    cfg: Config, source: socket.socket | Iterable[bytes], source_addr: _AddressType, http2_connection: Literal[True] = ...
) -> HTTP2ServerConnection: ...

__all__ = ["Message", "Request", "RequestParser", "get_parser"]
