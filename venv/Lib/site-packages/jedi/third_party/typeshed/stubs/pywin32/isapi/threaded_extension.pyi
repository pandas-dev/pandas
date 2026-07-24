import threading
from _typeshed import Unused
from collections.abc import Callable
from typing import Final

import isapi.simple

ISAPI_REQUEST: Final = 1
ISAPI_SHUTDOWN: Final = 2

class WorkerThread(threading.Thread):
    running: bool
    io_req_port: int
    extension: ThreadPoolExtension
    def __init__(self, extension: ThreadPoolExtension, io_req_port: int) -> None: ...
    def call_handler(self, cblock) -> None: ...

class ThreadPoolExtension(isapi.simple.SimpleExtension):
    max_workers: int
    worker_shutdown_wait: int
    workers: list[WorkerThread]
    dispatch_map: dict[int, Callable[..., Unused]]
    io_req_port: int
    def GetExtensionVersion(self, vi) -> None: ...
    def HttpExtensionProc(self, control_block) -> int: ...
    def TerminateExtension(self, status) -> None: ...
    def DispatchConnection(self, errCode, bytes, key, overlapped) -> None: ...
    def Dispatch(self, ecb) -> None: ...
    def HandleDispatchError(self, ecb) -> None: ...
