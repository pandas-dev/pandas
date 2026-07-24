from _typeshed import StrOrBytesPath
from typing import Any, overload

from gevent.baseserver import _Spawner
from gevent.server import StreamServer, _Address
from gevent.socket import socket as _GeventSocket
from gevent.ssl import SSLContext

class BackdoorServer(StreamServer):
    locals: dict[str, Any]
    banner: str | None
    @overload
    def __init__(
        self,
        listener: _GeventSocket | tuple[str, int] | str,
        locals: dict[str, Any] | None = None,
        banner: str | None = None,
        *,
        backlog: int | None = None,
        spawn: _Spawner = "default",
        ssl_context: SSLContext,
        server_side: bool = True,
        do_handshake_on_connect: bool = True,
        suppress_ragged_eofs: bool = True,
    ) -> None: ...
    @overload
    def __init__(
        self,
        listener: _GeventSocket | tuple[str, int] | str,
        locals: dict[str, Any] | None = None,
        banner: str | None = None,
        *,
        backlog: int | None = None,
        spawn: _Spawner = "default",
        keyfile: StrOrBytesPath = ...,
        certfile: StrOrBytesPath = ...,
        server_side: bool = True,
        cert_reqs: int = ...,
        ssl_version: int = ...,
        ca_certs: str = ...,
        do_handshake_on_connect: bool = True,
        suppress_ragged_eofs: bool = True,
        ciphers: str = ...,
    ) -> None: ...
    def handle(self, conn: _GeventSocket, _address: _Address) -> None: ...

__all__ = ["BackdoorServer"]
