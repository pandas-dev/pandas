from types import FrameType
from typing import Any
from typing_extensions import TypeAlias

from gunicorn.workers.base import Worker

IOLoop: TypeAlias = Any  # tornado IOLoop class
PeriodicCallback: TypeAlias = Any  # tornado PeriodicCallback class
_HTTPServer: TypeAlias = Any  # tornado httpserver.HTTPServer class

class TornadoWorker(Worker):
    alive: bool
    server_alive: bool
    ioloop: IOLoop
    callbacks: list[PeriodicCallback]
    server: _HTTPServer

    @classmethod
    def setup(cls) -> None: ...
    def handle_exit(self, sig: int, frame: FrameType | None) -> None: ...
    def handle_request(self) -> None: ...
    def watchdog(self) -> None: ...
    def heartbeat(self) -> None: ...
    def init_process(self) -> None: ...
    def run(self) -> None: ...
