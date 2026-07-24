import socket

from gunicorn.http import Request
from gunicorn.http2.connection import HTTP2ServerConnection
from gunicorn.workers import base

from .._types import _AddressType

ALREADY_HANDLED: object

class AsyncWorker(base.Worker):
    worker_connections: int
    alive: bool

    def timeout_ctx(self) -> None: ...
    def is_already_handled(self, respiter: object) -> bool: ...
    def handle(self, listener: socket.socket, client: socket.socket, addr: _AddressType) -> None: ...
    def handle_http2(self, listener: socket.socket, client: socket.socket, addr: _AddressType) -> None: ...
    def handle_http2_request(
        self, listener_name: _AddressType, req: Request, sock: socket.socket, addr: _AddressType, h2_conn: HTTP2ServerConnection
    ) -> None: ...
    def handle_request(self, listener_name: _AddressType, req: Request, sock: socket.socket, addr: _AddressType) -> bool: ...
