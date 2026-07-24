from collections.abc import Callable, Iterable
from typing import Any, TypeVar
from typing_extensions import ParamSpec, Self, disjoint_base

from gevent.greenlet import Greenlet
from gevent.queue import UnboundQueue

_T = TypeVar("_T")
_P = ParamSpec("_P")

# this matches builtins.map to some degree, but since it is an non-public API type that just gets
# returned by some public API functions, we don't bother adding a whole bunch of overloads to handle
# the case of 1-n Iterables being passed in and just go for the fully unsafe signature
# we do the crazy overloads instead in the functions that create these objects
@disjoint_base
class IMapUnordered(Greenlet[_P, _T]):
    finished: bool
    # it may contain an undocumented Failure object
    queue: UnboundQueue[_T | object]
    def __init__(self, func: Callable[_P, _T], iterable: Iterable[Any], spawn: Callable[_P, Greenlet[_P, _T]]) -> None: ...
    def __iter__(self) -> Self: ...
    def __next__(self) -> _T: ...

@disjoint_base
class IMap(IMapUnordered[_P, _T]):
    index: int

__all__ = ["IMapUnordered", "IMap"]
