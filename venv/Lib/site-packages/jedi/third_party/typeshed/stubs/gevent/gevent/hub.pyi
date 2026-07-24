from collections.abc import Callable
from types import TracebackType
from typing import Any, Generic, Protocol, TextIO, TypeVar, overload, type_check_only
from typing_extensions import ParamSpec

import gevent._hub_local
import gevent._waiter
import greenlet
from gevent._hub_primitives import WaitOperationsGreenlet
from gevent._ident import IdentRegistry
from gevent._monitor import PeriodicMonitoringThread
from gevent._types import _Loop, _Watcher
from gevent._util import Lazy, readproperty
from gevent.greenlet import Greenlet
from gevent.resolver import AbstractResolver
from gevent.threadpool import ThreadPool

_T = TypeVar("_T")
_P = ParamSpec("_P")

GreenletExit = greenlet.GreenletExit
getcurrent = greenlet.getcurrent
get_hub = gevent._hub_local.get_hub
Waiter = gevent._waiter.Waiter

@type_check_only
class _DefaultReturnProperty(Protocol[_T]):
    @overload
    def __get__(self, obj: None, owner: type[object] | None = None) -> property: ...
    @overload
    def __get__(self, obj: object, owner: type[object] | None = None) -> _T: ...
    def __set__(self, obj: object, value: _T | None) -> None: ...
    def __del__(self) -> None: ...

def spawn_raw(function: Callable[..., object], *args: object, **kwargs: object) -> greenlet.greenlet: ...
def sleep(seconds: float = 0, ref: bool = True) -> None: ...
def idle(priority: int = 0) -> None: ...
def kill(greenlet: greenlet.greenlet, exception: type[BaseException] | BaseException = ...) -> None: ...

class signal(Generic[_P]):
    greenlet_class: type[Greenlet[..., Any]] | None
    hub: Hub
    watcher: _Watcher
    handler: Callable[_P, object]
    # we can't use _P.args/_P.kwargs here because pyright will complain
    # mypy doesn't seem to mind though
    args: tuple[Any, ...]
    kwargs: dict[str, Any]
    def __init__(self, signalnum: int, handler: Callable[_P, object], *args: _P.args, **kwargs: _P.kwargs) -> None: ...
    @property
    def ref(self) -> bool: ...
    @ref.setter
    def ref(self, value: bool) -> None: ...
    def cancel(self) -> None: ...
    def handle(self) -> None: ...

def reinit(hub: Hub | None = None) -> None: ...

class Hub(WaitOperationsGreenlet):
    SYSTEM_ERROR: tuple[type[BaseException], ...]
    NOT_ERROR: tuple[type[BaseException], ...]
    threadpool_size: int
    periodic_monitoring_thread: PeriodicMonitoringThread | None
    thread_ident: int
    name: str
    loop: _Loop
    format_context: Callable[[object], str]
    minimal_ident: int
    @overload
    def __init__(self, loop: _Loop, default: None = None) -> None: ...
    @overload
    def __init__(self, loop: None = None, default: bool | None = None) -> None: ...
    @Lazy
    def ident_registry(self) -> IdentRegistry: ...
    @property
    def loop_class(self) -> type[_Loop]: ...
    @property
    def backend(self) -> int | str: ...
    @property
    def main_hub(self) -> bool: ...
    def handle_error(
        self,
        context: object | None,
        type: type[BaseException] | None,
        value: BaseException | str | None,
        tb: TracebackType | None,
    ) -> None: ...
    def handle_system_error(
        self, type: type[BaseException], value: BaseException | None, tb: TracebackType | None = None
    ) -> None: ...
    @readproperty
    def exception_stream(self) -> TextIO | None: ...
    def print_exception(
        self, context: object | None, t: type[BaseException] | None, v: BaseException | str | None, tb: TracebackType | None
    ) -> None: ...
    def run(self) -> None: ...
    def start_periodic_monitoring_thread(self) -> PeriodicMonitoringThread: ...
    def join(self, timeout: float | None = None) -> bool: ...
    def destroy(self, destroy_loop: bool | None = None) -> None: ...
    @property
    def resolver_class(self) -> type[AbstractResolver]: ...
    resolver: _DefaultReturnProperty[AbstractResolver]
    @property
    def threadpool_class(self) -> type[ThreadPool]: ...
    threadpool: _DefaultReturnProperty[ThreadPool]

class linkproxy:
    __slots__ = ["callback", "obj"]
    callback: Callable[[object], object]
    obj: object
    def __init__(self, callback: Callable[[_T], object], obj: _T) -> None: ...
    def __call__(self, *args: object) -> None: ...

__all__ = ["getcurrent", "GreenletExit", "spawn_raw", "sleep", "kill", "signal", "reinit", "get_hub", "Hub", "Waiter"]
