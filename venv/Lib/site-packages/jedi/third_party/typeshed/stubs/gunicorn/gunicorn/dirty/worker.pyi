from _typeshed import Incomplete
from asyncio import StreamReader, StreamWriter
from collections.abc import Iterable, Mapping
from signal import Signals
from typing import Any, ClassVar, Literal

from gunicorn.config import Config
from gunicorn.glogging import Logger as GLogger
from gunicorn.workers.workertmp import WorkerTmp

class DirtyWorker:
    SIGNALS: ClassVar[list[Signals]]
    age: int
    pid: int | Literal["[booting]"]
    ppid: int
    app_paths: list[str]
    cfg: Config
    log: GLogger
    socket_path: str
    booted: bool
    aborted: bool
    alive: bool
    tmp: WorkerTmp
    apps: dict[str, Incomplete]

    def __init__(self, age: int, ppid: int, app_paths: list[str], cfg: Config, log: GLogger, socket_path: str) -> None: ...
    def notify(self) -> None: ...
    def init_process(self) -> None: ...
    def init_signals(self) -> None: ...
    def load_apps(self) -> None: ...
    def run(self) -> None: ...
    async def handle_connection(self, reader: StreamReader, writer: StreamWriter) -> None: ...
    async def handle_request(self, message: dict[str, Incomplete], writer: StreamWriter) -> None: ...
    async def execute(
        self, app_path: str, action: str, args: Iterable[Any], kwargs: Mapping[str, Any]
    ) -> Any: ...  # Arguments and result depend on method name passed to action
