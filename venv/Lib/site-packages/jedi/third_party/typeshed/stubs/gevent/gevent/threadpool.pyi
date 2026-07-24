import concurrent.futures
from collections.abc import Callable
from typing import Any, Generic, TypeVar
from typing_extensions import ParamSpec, TypeAlias

from gevent._threading import Queue
from gevent._types import _AsyncWatcher, _Watcher
from gevent.event import AsyncResult, _OptExcInfo, _ValueSource
from gevent.greenlet import Greenlet
from gevent.hub import Hub
from gevent.pool import GroupMappingMixin

_T = TypeVar("_T")
_P = ParamSpec("_P")
_TaskItem: TypeAlias = tuple[Callable[..., Any], tuple[Any, ...], dict[str, Any], ThreadResult[Any]]
_Receiver: TypeAlias = Callable[[_ValueSource[_T]], object]

class ThreadPool(GroupMappingMixin):
    __slots__ = (
        "hub",
        "_maxsize",
        "manager",
        "pid",
        "fork_watcher",
        "_available_worker_threads_greenlet_sem",
        "_worker_greenlets",
        "task_queue",
        "_idle_task_timeout",
    )
    hub: Hub
    pid: int
    manager: Greenlet[..., Any] | None
    task_queue: Queue[_TaskItem]
    fork_watcher: _Watcher
    def __init__(self, maxsize: int, hub: Hub | None = None, idle_task_timeout: int = -1) -> None: ...
    @property
    def maxsize(self) -> int: ...
    @maxsize.setter
    def maxsize(self, value: int) -> None: ...
    @property
    def size(self) -> int: ...
    @size.setter
    def size(self, value: int) -> None: ...
    def __len__(self) -> int: ...
    def join(self) -> None: ...
    def kill(self) -> None: ...
    def adjust(self) -> None: ...
    def spawn(self, func: Callable[_P, _T], *args: _P.args, **kwargs: _P.kwargs) -> AsyncResult[_T]: ...  # type: ignore[override]

class ThreadResult(Generic[_T]):
    __slots__ = ("exc_info", "async_watcher", "_call_when_ready", "value", "context", "hub", "receiver")
    receiver: _Receiver[_T]
    hub: Hub
    context: object | None
    value: _T | None
    exc_info: _OptExcInfo | tuple[()]
    async_watcher: _AsyncWatcher
    def __init__(self, receiver: _Receiver[_T], hub: Hub, call_when_ready: Callable[[], object]) -> None: ...
    @property
    def exception(self) -> BaseException | None: ...
    def destroy_in_main_thread(self) -> None: ...
    def set(self, value: _T) -> None: ...
    def handle_error(self, context: object, exc_info: _OptExcInfo) -> None: ...
    def successful(self) -> bool: ...

class ThreadPoolExecutor(concurrent.futures.ThreadPoolExecutor):
    kill = concurrent.futures.ThreadPoolExecutor.shutdown

__all__ = ["ThreadPool", "ThreadResult", "ThreadPoolExecutor"]
