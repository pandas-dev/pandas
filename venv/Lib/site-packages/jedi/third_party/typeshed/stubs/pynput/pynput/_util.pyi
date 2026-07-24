import threading
from _typeshed import OptExcInfo
from collections.abc import Callable
from queue import Queue
from types import ModuleType, TracebackType
from typing import Any, ClassVar, Generic, TypedDict, TypeVar, type_check_only
from typing_extensions import ParamSpec, Self

_T = TypeVar("_T")
_AbstractListenerT = TypeVar("_AbstractListenerT", bound=AbstractListener)
_P = ParamSpec("_P")

@type_check_only
class _RESOLUTIONS(TypedDict):
    darwin: str
    uinput: str
    xorg: str

RESOLUTIONS: _RESOLUTIONS

def backend(package: str) -> ModuleType: ...
def prefix(base: type | tuple[type | tuple[Any, ...], ...], cls: type) -> str | None: ...

class AbstractListener(threading.Thread):
    class StopException(Exception): ...
    _HANDLED_EXCEPTIONS: ClassVar[tuple[type | tuple[Any, ...], ...]]  # undocumented
    _suppress: bool  # undocumented
    _running: bool  # undocumented
    _thread: threading.Thread  # undocumented
    _condition: threading.Condition  # undocumented
    _ready: bool  # undocumented
    _queue: Queue[OptExcInfo | None]  # undocumented
    daemon: bool
    def __init__(self, suppress: bool = False, **kwargs: Callable[..., bool | None] | None) -> None: ...
    @property
    def suppress(self) -> bool: ...
    @property
    def running(self) -> bool: ...
    def stop(self) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...
    def wait(self) -> None: ...
    def run(self) -> None: ...
    @classmethod
    def _emitter(cls, f: Callable[_P, _T]) -> Callable[_P, _T]: ...  # undocumented
    def _mark_ready(self) -> None: ...  # undocumented
    def _run(self) -> None: ...  # undocumented
    def _stop_platform(self) -> None: ...  # undocumented
    def join(self, timeout: float | None = None, *args: Any) -> None: ...

class Events(Generic[_T, _AbstractListenerT]):
    _Listener: type[_AbstractListenerT] | None  # undocumented

    class Event:
        def __eq__(self, other: object) -> bool: ...

    _event_queue: Queue[_T]  # undocumented
    _sentinel: object  # undocumented
    _listener: _AbstractListenerT  # undocumented
    start: Callable[[], None]
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...
    def __iter__(self) -> Self: ...
    def __next__(self) -> _T: ...
    def get(self, timeout: float | None = None) -> _T | None: ...
    def _event_mapper(self, event: Callable[_P, object]) -> Callable[_P, None]: ...

class NotifierMixin: ...
