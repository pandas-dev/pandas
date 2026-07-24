from collections.abc import Callable, Collection, Iterable, Iterator
from typing import Any, TypeVar, overload
from typing_extensions import ParamSpec

from gevent._imap import IMap, IMapUnordered
from gevent.greenlet import Greenlet
from gevent.queue import Full as QueueFull

__all__ = ["Group", "Pool", "PoolFull"]

_T = TypeVar("_T")
_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")
_T3 = TypeVar("_T3")
_T4 = TypeVar("_T4")
_T5 = TypeVar("_T5")
_S = TypeVar("_S")
_P = ParamSpec("_P")

class GroupMappingMixin:
    __slots__ = ()
    def spawn(self, func: Callable[_P, _T], *args: _P.args, **kwargs: _P.kwargs) -> Greenlet[_P, _T]: ...
    # we would like to use ParamSpec for these, but since args and kwds are passed in as is
    # pyright will complain if we use _P.args/_P.kwargs, it appears to work on mypy though
    # we can probably get away with Sequence and Mapping instead of tuple and dict, but for
    # now we will be strict, just to be safe
    def apply_cb(
        self,
        func: Callable[..., _T],
        args: tuple[Any, ...] | None = None,
        kwds: dict[str, Any] | None = None,
        callback: Callable[[_T], object] | None = None,
    ) -> _T: ...
    # The ParamSpec of the spawned greenlet can differ from the one being passed in, but the return type will match
    def apply_async(
        self,
        func: Callable[..., _T],
        args: tuple[Any, ...] | None = None,
        kwds: dict[str, Any] | None = None,
        callback: Callable[[_T], object] | None = None,
    ) -> Greenlet[..., _T]: ...
    def apply(self, func: Callable[..., _T], args: tuple[Any, ...] | None = None, kwds: dict[str, Any] | None = None) -> _T: ...
    def map(self, func: Callable[[_T], _S], iterable: Iterable[_T]) -> list[_S]: ...
    def map_cb(
        self, func: Callable[[_T], _S], iterable: Iterable[_T], callback: Callable[[list[_S]], object] | None = None
    ) -> list[_S]: ...
    def map_async(
        self, func: Callable[[_T], _S], iterable: Iterable[_T], callback: Callable[[list[_S]], object] | None = None
    ) -> Greenlet[..., list[_S]]: ...
    @overload
    def imap(self, func: Callable[[_T1], _S], iter1: Iterable[_T1], /, *, maxsize: int | None = None) -> IMap[[_T1], _S]: ...
    @overload
    def imap(
        self, func: Callable[[_T1, _T2], _S], iter1: Iterable[_T1], iter2: Iterable[_T2], /, *, maxsize: int | None = None
    ) -> IMap[[_T1, _T2], _S]: ...
    @overload
    def imap(
        self,
        func: Callable[[_T1, _T2, _T3], _S],
        iter1: Iterable[_T1],
        iter2: Iterable[_T2],
        iter3: Iterable[_T3],
        /,
        *,
        maxsize: int | None = None,
    ) -> IMap[[_T1, _T2, _T3], _S]: ...
    @overload
    def imap(
        self,
        func: Callable[[_T1, _T2, _T3, _T4], _S],
        iter1: Iterable[_T1],
        iter2: Iterable[_T2],
        iter3: Iterable[_T3],
        iter4: Iterable[_T4],
        /,
        *,
        maxsize: int | None = None,
    ) -> IMap[[_T1, _T2, _T3, _T4], _S]: ...
    @overload
    def imap(
        self,
        func: Callable[[_T1, _T2, _T3, _T4, _T5], _S],
        iter1: Iterable[_T1],
        iter2: Iterable[_T2],
        iter3: Iterable[_T3],
        iter4: Iterable[_T4],
        iter5: Iterable[_T5],
        /,
        *,
        maxsize: int | None = None,
    ) -> IMap[[_T1, _T2, _T3, _T4, _T5], _S]: ...
    @overload
    def imap(
        self,
        func: Callable[_P, _S],
        iter1: Iterable[Any],
        iter2: Iterable[Any],
        iter3: Iterable[Any],
        iter4: Iterable[Any],
        iter5: Iterable[Any],
        iter6: Iterable[Any],
        /,
        *iterables: Iterable[Any],
        maxsize: int | None = None,
    ) -> IMap[_P, _S]: ...
    @overload
    def imap_unordered(
        self, func: Callable[[_T1], _S], iter1: Iterable[_T1], /, *, maxsize: int | None = None
    ) -> IMapUnordered[[_T1], _S]: ...
    @overload
    def imap_unordered(
        self, func: Callable[[_T1, _T2], _S], iter1: Iterable[_T1], iter2: Iterable[_T2], /, *, maxsize: int | None = None
    ) -> IMapUnordered[[_T1, _T2], _S]: ...
    @overload
    def imap_unordered(
        self,
        func: Callable[[_T1, _T2, _T3], _S],
        iter1: Iterable[_T1],
        iter2: Iterable[_T2],
        iter3: Iterable[_T3],
        /,
        *,
        maxsize: int | None = None,
    ) -> IMapUnordered[[_T1, _T2, _T3], _S]: ...
    @overload
    def imap_unordered(
        self,
        func: Callable[[_T1, _T2, _T3, _T4], _S],
        iter1: Iterable[_T1],
        iter2: Iterable[_T2],
        iter3: Iterable[_T3],
        iter4: Iterable[_T4],
        /,
        *,
        maxsize: int | None = None,
    ) -> IMapUnordered[[_T1, _T2, _T3, _T4], _S]: ...
    @overload
    def imap_unordered(
        self,
        func: Callable[[_T1, _T2, _T3, _T4, _T5], _S],
        iter1: Iterable[_T1],
        iter2: Iterable[_T2],
        iter3: Iterable[_T3],
        iter4: Iterable[_T4],
        iter5: Iterable[_T5],
        /,
        *,
        maxsize: int | None = None,
    ) -> IMapUnordered[[_T1, _T2, _T3, _T4, _T5], _S]: ...
    @overload
    def imap_unordered(
        self,
        func: Callable[_P, _S],
        iter1: Iterable[Any],
        iter2: Iterable[Any],
        iter3: Iterable[Any],
        iter4: Iterable[Any],
        iter5: Iterable[Any],
        iter6: Iterable[Any],
        /,
        *iterables: Iterable[Any],
        maxsize: int | None = None,
    ) -> IMapUnordered[_P, _S]: ...

# TODO: Consider making these generic in Greenlet. The drawback would be, that it
#       wouldn't be possible to mix Greenlets with different return values/ParamSpecs
#       unless you bind Grenlet[..., object], but in that case all the spawn/apply/map
#       methods become less helpful, because the return types cannot be as specific...
#       We would need higher-kinded TypeVars if we wanted to give up neither
class Group(GroupMappingMixin):
    greenlet_class: type[Greenlet[..., Any]]
    greenlets: set[Greenlet[..., Any]]
    dying: set[Greenlet[..., Any]]
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, grenlets: Collection[Greenlet[..., object]], /) -> None: ...
    def __len__(self) -> int: ...
    def __contains__(self, item: Greenlet[..., object]) -> bool: ...
    def __iter__(self) -> Iterator[Greenlet[..., object]]: ...
    def add(self, greenlet: Greenlet[..., object]) -> None: ...
    def discard(self, greenlet: Greenlet[..., object]) -> None: ...
    def start(self, greenlet: Greenlet[..., object]) -> None: ...
    def join(self, timeout: float | None = None, raise_error: bool = False) -> bool: ...
    def kill(
        self, exception: type[BaseException] | BaseException = ..., block: bool = True, timeout: float | None = None
    ) -> None: ...
    def killone(
        self,
        greenlet: Greenlet[..., object],
        exception: type[BaseException] | BaseException = ...,
        block: bool = True,
        timeout: float | None = None,
    ) -> None: ...
    def full(self) -> bool: ...
    def wait_available(self, timeout: float | None = None) -> int | None: ...

class PoolFull(QueueFull): ...

class Pool(Group):
    size: int | None
    def __init__(self, size: int | None = None, greenlet_class: type[Greenlet[..., object]] | None = None) -> None: ...
    def wait_available(self, timeout: float | None = None) -> int: ...
    def free_count(self) -> int: ...
    def start(self, greenlet: Greenlet[..., object], blocking: bool = True, timeout: float | None = None) -> None: ...
    def add(self, greenlet: Greenlet[..., object], blocking: bool = True, timeout: float | None = None) -> None: ...
