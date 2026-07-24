from _typeshed import Unused
from _typeshed.wsgi import WSGIEnvironment
from collections import deque
from collections.abc import Callable, Mapping, Sequence
from logging import Logger
from threading import Condition, Lock

from .channel import HTTPChannel
from .utilities import Error

rename_headers: Mapping[str, str]
hop_by_hop: frozenset[str]

class ThreadedTaskDispatcher:
    stop_count: int
    active_count: int
    logger: Logger
    queue_logger: Logger
    threads: set[int]
    queue: deque[Task]
    lock: Lock
    queue_cv: Condition
    thread_exit_cv: Condition
    def start_new_thread(self, target: Callable[[int], Unused], thread_no: int) -> None: ...
    def handler_thread(self, thread_no: int) -> None: ...
    def set_thread_count(self, count: int) -> None: ...
    def add_task(self, task: Task) -> None: ...
    def shutdown(self, cancel_pending: bool = True, timeout: int = 5) -> bool: ...

class Task:
    close_on_finish: bool
    status: str
    wrote_header: bool
    start_time: int
    content_length: int | None
    content_bytes_written: int
    logged_write_excess: bool
    logged_write_no_body: bool
    complete: bool
    chunked_response: bool
    logger: Logger
    channel: HTTPChannel
    request: Error
    response_headers: Sequence[tuple[str, str]]
    version: str
    def __init__(self, channel: HTTPChannel, request: Error) -> None: ...
    def service(self) -> None: ...
    @property
    def has_body(self) -> bool: ...
    def set_close_on_finish(self) -> None: ...
    def build_response_header(self) -> bytes: ...
    def remove_content_length_header(self) -> None: ...
    def start(self) -> None: ...
    def finish(self) -> None: ...
    def write(self, data: bytes) -> None: ...

class ErrorTask(Task):
    complete: bool
    status: str
    close_on_finish: bool
    content_length: int | None
    def execute(self) -> None: ...

class WSGITask(Task):
    environ: WSGIEnvironment | None
    response_headers: Sequence[tuple[str, str]]
    complete: bool
    status: str
    content_length: int | None
    close_on_finish: bool
    def execute(self) -> None: ...
    def get_environment(self) -> WSGIEnvironment: ...
