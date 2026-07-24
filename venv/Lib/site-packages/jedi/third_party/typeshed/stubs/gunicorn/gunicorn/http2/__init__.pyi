from typing import Final

from .async_connection import AsyncHTTP2Connection
from .connection import HTTP2ServerConnection

H2_MIN_VERSION: Final[tuple[int, int, int]]

def is_http2_available() -> bool: ...
def get_h2_version() -> tuple[int, int, int]: ...
def get_http2_connection_class() -> type[HTTP2ServerConnection]: ...
def get_async_http2_connection_class() -> type[AsyncHTTP2Connection]: ...

__all__ = [
    "is_http2_available",
    "get_h2_version",
    "get_http2_connection_class",
    "get_async_http2_connection_class",
    "H2_MIN_VERSION",
]
