from collections.abc import Callable
from types import TracebackType
from typing import Any, Literal, Protocol, TypeVar, overload, type_check_only
from typing_extensions import ParamSpec, Self

from gevent._types import _TimerWatcher

_T = TypeVar("_T")
_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")
_TimeoutT = TypeVar("_TimeoutT", bound=Timeout)
_P = ParamSpec("_P")

@type_check_only
class _HasSeconds(Protocol):
    @property
    def seconds(self) -> float | int: ...

class Timeout(BaseException):
    seconds: float | None
    exception: type[BaseException] | BaseException | None
    timer: _TimerWatcher
    def __init__(
        self,
        seconds: float | None = None,
        exception: type[BaseException] | BaseException | None = None,
        ref: bool = True,
        priority: int = -1,
    ) -> None: ...
    def start(self) -> None: ...
    @overload
    @classmethod
    def start_new(
        cls, timeout: None | float = None, exception: type[BaseException] | BaseException | None = None, ref: bool = True
    ) -> Self: ...
    @overload
    @classmethod
    def start_new(cls, timeout: _TimeoutT) -> _TimeoutT: ...
    @property
    def pending(self) -> bool: ...
    def cancel(self) -> None: ...
    def close(self) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self, typ: type[BaseException] | None, value: BaseException | None, tb: TracebackType | None
    ) -> Literal[True] | None: ...
    def __lt__(self, other: _HasSeconds | float) -> bool: ...

# when timeout_value is provided we unfortunately get no type checking on *args, **kwargs, because
# ParamSpec does not allow mixing in additional keyword arguments
@overload
def with_timeout(
    seconds: float | None, function: Callable[..., _T1], *args: Any, timeout_value: _T2, **kwds: Any
) -> _T1 | _T2: ...
@overload
def with_timeout(seconds: float | None, function: Callable[_P, _T], *args: _P.args, **kwds: _P.kwargs) -> _T: ...

__all__ = ["Timeout", "with_timeout"]
