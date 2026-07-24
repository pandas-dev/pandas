import socket
from types import FrameType
from typing import ClassVar

from gunicorn.app.base import BaseApplication
from gunicorn.config import Config
from gunicorn.glogging import Logger as GLogger
from gunicorn.http import Request
from gunicorn.workers.workertmp import WorkerTmp

from .._types import _AddressType, _WSGIAppType
from ..reloader import _ReloaderType

class Worker:
    SIGNALS: ClassVar[list[int]]
    PIPE: ClassVar[list[int]]
    age: int
    pid: str
    ppid: int
    sockets: list[socket.socket]
    app: BaseApplication
    timeout: int
    cfg: Config
    booted: bool
    aborted: bool
    reloader: _ReloaderType | None
    nr: int
    max_requests: int
    alive: bool
    log: GLogger
    tmp: WorkerTmp
    wait_fds: list[socket.socket | int]
    wsgi: _WSGIAppType

    def __init__(
        self, age: int, ppid: int, sockets: list[socket.socket], app: BaseApplication, timeout: int, cfg: Config, log: GLogger
    ) -> None: ...
    def notify(self) -> None: ...
    def run(self) -> None: ...
    def init_process(self) -> None: ...
    def load_wsgi(self) -> None: ...
    def init_signals(self) -> None: ...
    def handle_usr1(self, sig: int, frame: FrameType | None) -> None: ...
    def handle_exit(self, sig: int, frame: FrameType | None) -> None: ...
    def handle_quit(self, sig: int, frame: FrameType | None) -> None: ...
    def handle_abort(self, sig: int, frame: FrameType | None) -> None: ...
    def handle_error(self, req: Request | None, client: socket.socket, addr: _AddressType, exc: BaseException) -> None: ...
    def handle_winch(self, sig: int, fname: str | None) -> None: ...
