import asyncio
import dataclasses
import sys
from asyncio.coroutines import _is_coroutine  # type: ignore[attr-defined]
from functools import _CacheInfo, _make_key, partial, partialmethod
from typing import (
    Any,
    Callable,
    Coroutine,
    Generic,
    Hashable,
    Optional,
    OrderedDict,
    Set,
    Type,
    TypedDict,
    TypeVar,
    Union,
    cast,
    final,
    overload,
)


if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self


__version__ = "2.0.5"

__all__ = ("alru_cache",)


_T = TypeVar("_T")
_R = TypeVar("_R")
_Coro = Coroutine[Any, Any, _R]
_CB = Callable[..., _Coro[_R]]
_CBP = Union[_CB[_R], "partial[_Coro[_R]]", "partialmethod[_Coro[_R]]"]


@final
class _CacheParameters(TypedDict):
    typed: bool
    maxsize: Optional[int]
    tasks: int
    closed: bool


@final
@dataclasses.dataclass
class _CacheItem(Generic[_R]):
    fut: "asyncio.Future[_R]"
    later_call: Optional[asyncio.Handle]

    def cancel(self) -> None:
        if self.later_call is not None:
            self.later_call.cancel()
            self.later_call = None


@final
class _LRUCacheWrapper(Generic[_R]):
    def __init__(
        self,
        fn: _CB[_R],
        maxsize: Optional[int],
        typed: bool,
        ttl: Optional[float],
    ) -> None:
        try:
            self.__module__ = fn.__module__
        except AttributeError:
            pass
        try:
            self.__name__ = fn.__name__
        except AttributeError:
            pass
        try:
            self.__qualname__ = fn.__qualname__
        except AttributeError:
            pass
        try:
            self.__doc__ = fn.__doc__
        except AttributeError:
            pass
        try:
            self.__annotations__ = fn.__annotations__
        except AttributeError:
            pass
        try:
            self.__dict__.update(fn.__dict__)
        except AttributeError:
            pass
        # set __wrapped__ last so we don't inadvertently copy it
        # from the wrapped function when updating __dict__
        self._is_coroutine = _is_coroutine
        self.__wrapped__ = fn
        self.__maxsize = maxsize
        self.__typed = typed
        self.__ttl = ttl
        self.__cache: OrderedDict[Hashable, _CacheItem[_R]] = OrderedDict()
        self.__closed = False
        self.__hits = 0
        self.__misses = 0
        self.__tasks: Set["asyncio.Task[_R]"] = set()

    def cache_invalidate(self, /, *args: Hashable, **kwargs: Any) -> bool:
        key = _make_key(args, kwargs, self.__typed)

        cache_item = self.__cache.pop(key, None)
        if cache_item is None:
            return False
        else:
            cache_item.cancel()
            return True

    def cache_clear(self) -> None:
        self.__hits = 0
        self.__misses = 0

        for c in self.__cache.values():
            if c.later_call:
                c.later_call.cancel()
        self.__cache.clear()
        self.__tasks.clear()

    async def cache_close(self, *, wait: bool = False) -> None:
        self.__closed = True

        tasks = list(self.__tasks)
        if not tasks:
            return

        if not wait:
            for task in tasks:
                if not task.done():
                    task.cancel()

        await asyncio.gather(*tasks, return_exceptions=True)

    def cache_info(self) -> _CacheInfo:
        return _CacheInfo(
            self.__hits,
            self.__misses,
            self.__maxsize,
            len(self.__cache),
        )

    def cache_parameters(self) -> _CacheParameters:
        return _CacheParameters(
            maxsize=self.__maxsize,
            typed=self.__typed,
            tasks=len(self.__tasks),
            closed=self.__closed,
        )

    def _cache_hit(self, key: Hashable) -> None:
        self.__hits += 1
        self.__cache.move_to_end(key)

    def _cache_miss(self, key: Hashable) -> None:
        self.__misses += 1

    def _task_done_callback(
        self, fut: "asyncio.Future[_R]", key: Hashable, task: "asyncio.Task[_R]"
    ) -> None:
        self.__tasks.discard(task)

        if task.cancelled():
            fut.cancel()
            self.__cache.pop(key, None)
            return

        exc = task.exception()
        if exc is not None:
            fut.set_exception(exc)
            self.__cache.pop(key, None)
            return

        cache_item = self.__cache.get(key)
        if self.__ttl is not None and cache_item is not None:
            loop = asyncio.get_running_loop()
            cache_item.later_call = loop.call_later(
                self.__ttl, self.__cache.pop, key, None
            )

        fut.set_result(task.result())

    async def __call__(self, /, *fn_args: Any, **fn_kwargs: Any) -> _R:
        if self.__closed:
            raise RuntimeError(f"alru_cache is closed for {self}")

        loop = asyncio.get_running_loop()

        key = _make_key(fn_args, fn_kwargs, self.__typed)

        cache_item = self.__cache.get(key)

        if cache_item is not None:
            self._cache_hit(key)
            if not cache_item.fut.done():
                return await asyncio.shield(cache_item.fut)

            return cache_item.fut.result()

        fut = loop.create_future()
        coro = self.__wrapped__(*fn_args, **fn_kwargs)
        task: asyncio.Task[_R] = loop.create_task(coro)
        self.__tasks.add(task)
        task.add_done_callback(partial(self._task_done_callback, fut, key))

        self.__cache[key] = _CacheItem(fut, None)

        if self.__maxsize is not None and len(self.__cache) > self.__maxsize:
            dropped_key, cache_item = self.__cache.popitem(last=False)
            cache_item.cancel()

        self._cache_miss(key)
        return await asyncio.shield(fut)

    def __get__(
        self, instance: _T, owner: Optional[Type[_T]]
    ) -> Union[Self, "_LRUCacheWrapperInstanceMethod[_R, _T]"]:
        if owner is None:
            return self
        else:
            return _LRUCacheWrapperInstanceMethod(self, instance)


@final
class _LRUCacheWrapperInstanceMethod(Generic[_R, _T]):
    def __init__(
        self,
        wrapper: _LRUCacheWrapper[_R],
        instance: _T,
    ) -> None:
        try:
            self.__module__ = wrapper.__module__
        except AttributeError:
            pass
        try:
            self.__name__ = wrapper.__name__
        except AttributeError:
            pass
        try:
            self.__qualname__ = wrapper.__qualname__
        except AttributeError:
            pass
        try:
            self.__doc__ = wrapper.__doc__
        except AttributeError:
            pass
        try:
            self.__annotations__ = wrapper.__annotations__
        except AttributeError:
            pass
        try:
            self.__dict__.update(wrapper.__dict__)
        except AttributeError:
            pass
        # set __wrapped__ last so we don't inadvertently copy it
        # from the wrapped function when updating __dict__
        self._is_coroutine = _is_coroutine
        self.__wrapped__ = wrapper.__wrapped__
        self.__instance = instance
        self.__wrapper = wrapper

    def cache_invalidate(self, /, *args: Hashable, **kwargs: Any) -> bool:
        return self.__wrapper.cache_invalidate(self.__instance, *args, **kwargs)

    def cache_clear(self) -> None:
        self.__wrapper.cache_clear()

    async def cache_close(
        self, *, cancel: bool = False, return_exceptions: bool = True
    ) -> None:
        await self.__wrapper.cache_close()

    def cache_info(self) -> _CacheInfo:
        return self.__wrapper.cache_info()

    def cache_parameters(self) -> _CacheParameters:
        return self.__wrapper.cache_parameters()

    async def __call__(self, /, *fn_args: Any, **fn_kwargs: Any) -> _R:
        return await self.__wrapper(self.__instance, *fn_args, **fn_kwargs)


def _make_wrapper(
    maxsize: Optional[int],
    typed: bool,
    ttl: Optional[float] = None,
) -> Callable[[_CBP[_R]], _LRUCacheWrapper[_R]]:
    def wrapper(fn: _CBP[_R]) -> _LRUCacheWrapper[_R]:
        origin = fn

        while isinstance(origin, (partial, partialmethod)):
            origin = origin.func

        if not asyncio.iscoroutinefunction(origin):
            raise RuntimeError(f"Coroutine function is required, got {fn!r}")

        # functools.partialmethod support
        if hasattr(fn, "_make_unbound_method"):
            fn = fn._make_unbound_method()

        return _LRUCacheWrapper(cast(_CB[_R], fn), maxsize, typed, ttl)

    return wrapper


@overload
def alru_cache(
    maxsize: Optional[int] = 128,
    typed: bool = False,
    *,
    ttl: Optional[float] = None,
) -> Callable[[_CBP[_R]], _LRUCacheWrapper[_R]]:
    ...


@overload
def alru_cache(
    maxsize: _CBP[_R],
    /,
) -> _LRUCacheWrapper[_R]:
    ...


def alru_cache(
    maxsize: Union[Optional[int], _CBP[_R]] = 128,
    typed: bool = False,
    *,
    ttl: Optional[float] = None,
) -> Union[Callable[[_CBP[_R]], _LRUCacheWrapper[_R]], _LRUCacheWrapper[_R]]:
    if maxsize is None or isinstance(maxsize, int):
        return _make_wrapper(maxsize, typed, ttl)
    else:
        fn = cast(_CB[_R], maxsize)

        if callable(fn) or hasattr(fn, "_make_unbound_method"):
            return _make_wrapper(128, False, None)(fn)

        raise NotImplementedError(f"{fn!r} decorating is not supported")
