import enum
import io
import queue
import threading
from collections.abc import Callable, Iterable
from typing import Any, ClassVar, TypeVar
from typing_extensions import TypeAlias

import gdb

_T = TypeVar("_T")
_F = TypeVar("_F", bound=Callable[..., Any])

DAPQueue: TypeAlias = queue.SimpleQueue[_T]

class DAPException(Exception): ...

def parse_and_eval(expression: str, global_context: bool = False) -> gdb.Value: ...

# target and args are passed to gdb.Thread
def start_thread(name: str, target: Callable[..., object], args: Iterable[Any] = ()) -> gdb.Thread: ...
def start_dap(target: Callable[..., object]) -> None: ...
def in_gdb_thread(func: _F) -> _F: ...
def in_dap_thread(func: _F) -> _F: ...

class LogLevel(enum.IntEnum):
    DEFAULT = 1
    FULL = 2

class LogLevelParam(gdb.Parameter):
    def __init__(self) -> None: ...

class LoggingParam(gdb.Parameter):
    lock: ClassVar[threading.Lock]
    log_file: io.TextIOWrapper | None
    def __init__(self) -> None: ...
    def get_set_string(self) -> str: ...

dap_log: LoggingParam

def log(something: object, level: LogLevel = LogLevel.DEFAULT) -> None: ...
def thread_log(something: object, level: LogLevel = LogLevel.DEFAULT) -> None: ...
def log_stack(level: LogLevel = LogLevel.DEFAULT) -> None: ...
def exec_and_log(cmd: str) -> None: ...
