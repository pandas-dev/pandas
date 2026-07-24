from collections.abc import Sequence
from socket import _RetAddress, socket
from threading import Condition, Lock

from waitress.adjustments import Adjustments
from waitress.buffers import OverflowableBuffer
from waitress.parser import HTTPRequestParser
from waitress.server import BaseWSGIServer
from waitress.task import ErrorTask, WSGITask

from . import wasyncore
from .wasyncore import _SocketMap

class ClientDisconnected(Exception): ...

class HTTPChannel(wasyncore.dispatcher):
    task_class: type[WSGITask]
    error_task_class: type[ErrorTask]
    parser_class: type[HTTPRequestParser]
    request: HTTPRequestParser | None
    last_activity: float
    will_close: bool
    close_when_flushed: bool
    requests: Sequence[HTTPRequestParser]
    sent_continue: bool
    total_outbufs_len: int
    current_outbuf_count: int
    server: BaseWSGIServer
    adj: Adjustments
    outbufs: Sequence[OverflowableBuffer]
    creation_time: float
    sendbuf_len: int
    task_lock: Lock
    outbuf_lock: Condition
    connected: bool
    addr: _RetAddress
    def __init__(
        self, server: BaseWSGIServer, sock: socket, addr: _RetAddress, adj: Adjustments, map: _SocketMap | None = None
    ) -> None: ...
    def check_client_disconnected(self) -> None: ...
    def writable(self) -> bool: ...
    def handle_write(self) -> None: ...
    def readable(self) -> bool: ...
    def handle_read(self) -> None: ...
    def send_continue(self) -> None: ...
    def received(self, data: bytes) -> bool: ...
    def handle_close(self) -> None: ...
    def add_channel(self, map: _SocketMap | None = None) -> None: ...
    def del_channel(self, map: _SocketMap | None = None) -> None: ...
    def write_soon(self, data: bytes) -> int: ...
    def service(self) -> None: ...
    def cancel(self) -> None: ...
