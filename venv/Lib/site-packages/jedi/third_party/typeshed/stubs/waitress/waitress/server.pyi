import sys
from _typeshed import Unused
from _typeshed.wsgi import WSGIApplication
from collections.abc import Callable, Sequence
from socket import _RetAddress, socket
from typing import Literal

from waitress import wasyncore
from waitress.adjustments import Adjustments, _AdjustmentsParams
from waitress.channel import HTTPChannel
from waitress.task import Task, ThreadedTaskDispatcher
from waitress.wasyncore import _SocketMap

def create_server(
    application: WSGIApplication,
    map: _SocketMap | None = None,
    _start: bool = True,
    _sock: socket | None = None,
    _dispatcher: ThreadedTaskDispatcher | None = None,
    **kw: _AdjustmentsParams,
) -> MultiSocketServer | BaseWSGIServer: ...

class MultiSocketServer:
    asyncore = wasyncore
    adj: Adjustments
    map: _SocketMap | None
    effective_listen: Sequence[tuple[str, int]]
    task_dispatcher: ThreadedTaskDispatcher
    def __init__(
        self,
        map: _SocketMap | None = None,
        adj: Adjustments | None = None,
        effective_listen: Sequence[tuple[str, int]] | None = None,
        dispatcher: ThreadedTaskDispatcher | None = None,
        # Can be None, but print_listen will fail
        log_info: Callable[[str], Unused] | None = None,
    ) -> None: ...
    def print_listen(self, format_str: str) -> None: ...
    def run(self) -> None: ...
    def close(self) -> None: ...

class BaseWSGIServer(wasyncore.dispatcher):
    channel_class: type[HTTPChannel]
    next_channel_cleanup: int
    socketmod: socket
    asyncore = wasyncore
    in_connection_overflow: bool
    sockinfo: tuple[int, int, int | None, _RetAddress]
    family: int
    socktype: int
    application: WSGIApplication
    adj: Adjustments
    trigger: int
    task_dispatcher: ThreadedTaskDispatcher
    server_name: str
    active_channels: HTTPChannel
    def __init__(
        self,
        application: WSGIApplication,
        map: _SocketMap | None = None,
        _start: bool = True,
        _sock: socket | None = None,
        dispatcher: ThreadedTaskDispatcher | None = None,
        adj: Adjustments | None = None,
        sockinfo: tuple[int, int, int | None, _RetAddress] | None = None,
        bind_socket: bool = True,
        **kw: _AdjustmentsParams,
    ) -> None: ...
    def bind_server_socket(self) -> None: ...
    def getsockname(self) -> tuple[str, str]: ...
    accepting: bool
    def accept_connections(self) -> None: ...
    def add_task(self, task: Task) -> None: ...
    def readable(self) -> bool: ...
    def writable(self) -> bool: ...
    def handle_read(self) -> None: ...
    def handle_connect(self) -> None: ...
    def handle_accept(self) -> None: ...
    def run(self) -> None: ...
    def pull_trigger(self) -> None: ...
    def set_socket_options(self, conn: socket) -> None: ...
    def fix_addr(self, addr: _RetAddress) -> _RetAddress: ...
    def maintenance(self, now: int) -> None: ...
    def print_listen(self, format_str: str) -> None: ...
    def close(self) -> None: ...

class TcpWSGIServer(BaseWSGIServer):
    def bind_server_socket(self) -> None: ...
    def getsockname(self) -> tuple[str, str]: ...
    def set_socket_options(self, conn: socket) -> None: ...

if sys.platform != "win32":
    class UnixWSGIServer(BaseWSGIServer):
        def __init__(
            self,
            application: WSGIApplication,
            map: _SocketMap | None = None,
            _start: bool = True,
            _sock: socket | None = None,
            dispatcher: ThreadedTaskDispatcher | None = None,
            adj: Adjustments | None = None,
            sockinfo: tuple[int, int, int | None, _RetAddress] | None = None,
            **kw: _AdjustmentsParams,
        ) -> None: ...
        def bind_server_socket(self) -> None: ...
        def getsockname(self) -> tuple[str, str]: ...
        def fix_addr(self, addr: _RetAddress) -> tuple[Literal["localhost"], None]: ...

WSGIServer = TcpWSGIServer
