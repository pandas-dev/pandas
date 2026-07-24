from types import TracebackType
from typing import Generic, TypeVar, final, overload
from typing_extensions import TypeAlias

from gevent.event import _ValueSource
from gevent.hub import Hub
from greenlet import greenlet as greenlet_t

__all__ = ["Waiter"]

_T = TypeVar("_T")
# this is annoying, it's due to them using *throw args, rather than just storing them in standardized form
_ThrowArgs: TypeAlias = (
    tuple[()]
    | tuple[BaseException]
    | tuple[BaseException, None]
    | tuple[BaseException, None, TracebackType | None]
    | tuple[type[BaseException]]
    | tuple[type[BaseException], BaseException | object]
    | tuple[type[BaseException], BaseException | object, TracebackType | None]
)

class Waiter(Generic[_T]):
    __slots__ = ["hub", "greenlet", "value", "_exception"]
    @property
    def hub(self) -> Hub: ...  # readonly in Cython
    @property
    def greenlet(self) -> greenlet_t | None: ...  # readonly in Cython
    @property
    def value(self) -> _T | None: ...  # readonly in Cython
    def __init__(self, hub: Hub | None = None) -> None: ...
    def clear(self) -> None: ...
    def ready(self) -> bool: ...
    def successful(self) -> bool: ...
    @property
    def exc_info(self) -> _ThrowArgs | None: ...
    def switch(self, value: _T) -> None: ...
    @overload
    def throw(self, typ: type[BaseException], val: BaseException | object = None, tb: TracebackType | None = None, /) -> None: ...
    @overload
    def throw(self, typ: BaseException = ..., val: None = None, tb: TracebackType | None = None, /) -> None: ...
    def get(self) -> _T: ...
    def __call__(self, source: _ValueSource[_T]) -> None: ...

@final
class MultipleWaiter(Waiter[_T]):
    __slots__ = ["_values"]
