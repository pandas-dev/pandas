from __future__ import annotations

__all__ = (
    "accumulate",
    "batched",
    "Chain",
    "combinations",
    "combinations_with_replacement",
    "compress",
    "count",
    "cycle",
    "dropwhile",
    "filterfalse",
    "groupby",
    "islice",
    "pairwise",
    "permutations",
    "product",
    "repeat",
    "starmap",
    "tee",
    "takewhile",
    "zip_longest",
)

import itertools
import operator
import sys
from collections.abc import (
    AsyncGenerator,
    AsyncIterable,
    AsyncIterator,
    Awaitable,
    Callable,
    Iterable,
    Iterator,
)
from dataclasses import dataclass, field
from typing import Any, Generic, TypeVar, cast, overload

from ._core._synchronization import Lock
from ._core._tasks import CancelScope
from .lowlevel import cancel_shielded_checkpoint, checkpoint, checkpoint_if_cancelled

T = TypeVar("T")
R = TypeVar("R")
_tee_end = object()


@dataclass(eq=False)
class _IterableAsyncIterator(AsyncIterator[T]):
    iterator: Iterator[T]

    async def __anext__(self) -> T:
        await checkpoint_if_cancelled()
        try:
            result = next(self.iterator)
        except StopIteration:
            await cancel_shielded_checkpoint()
            raise StopAsyncIteration from None

        await cancel_shielded_checkpoint()
        return result


def _iterate(iterable: Iterable[T] | AsyncIterable[T]) -> AsyncIterator[T]:
    if isinstance(iterable, AsyncIterator):
        return iterable

    if isinstance(iterable, AsyncIterable):
        return iterable.__aiter__()

    return _IterableAsyncIterator(iter(iterable))


@dataclass(eq=False)
class _TeeLink(Generic[T]):
    value: object | None = None
    next: _TeeLink[T] | None = None
    filled: bool = False


@dataclass(eq=False)
class _TeeState(Generic[T]):
    iterator: AsyncIterator[T]
    lock: Lock = field(default_factory=Lock)

    async def fill(self, link: _TeeLink[T]) -> bool:
        if link.filled:
            return False

        async with self.lock:
            if link.filled:
                return True

            link.value = await anext(self.iterator, _tee_end)
            if link.value is not _tee_end:
                link.next = _TeeLink()

            link.filled = True
            return True


class _TeeAsyncIterator(AsyncIterator[T]):
    _state: _TeeState[T]
    _link: _TeeLink[T]
    _element_yielded: bool

    def __init__(
        self, iterable: Iterable[T] | AsyncIterable[T] | _TeeAsyncIterator[T]
    ) -> None:
        if isinstance(iterable, _TeeAsyncIterator):
            self._state = iterable._state
            self._link = iterable._link
        else:
            self._state = _TeeState(_iterate(iterable))
            self._link = _TeeLink()

        self._element_yielded = False

    async def __anext__(self) -> T:
        had_yieldpoint = await self._state.fill(self._link)
        if self._link.value is _tee_end:
            if not self._element_yielded:
                await checkpoint()

            raise StopAsyncIteration

        if not had_yieldpoint:
            await checkpoint_if_cancelled()

        self._element_yielded = True
        value = cast(T, self._link.value)
        next_link = self._link.next
        assert next_link is not None
        self._link = next_link
        if not had_yieldpoint:
            await cancel_shielded_checkpoint()

        return value


async def _operator_add(x: T, y: T) -> T:
    return operator.add(x, y)


async def accumulate(
    iterable: Iterable[T] | AsyncIterable[T],
    function: Callable[[T, T], Awaitable[T]] = _operator_add,
    *,
    initial: T | None = None,
) -> AsyncGenerator[T, None]:
    iterator = _iterate(iterable)
    if initial is None:
        try:
            total = await anext(iterator)
        except StopAsyncIteration:
            await checkpoint()
            return
    else:
        await checkpoint_if_cancelled()
        total = initial
        await cancel_shielded_checkpoint()

    yield total

    async for element in iterator:
        total = await function(total, element)
        yield total


async def batched(
    iterable: Iterable[T] | AsyncIterable[T], n: int, *, strict: bool = False
) -> AsyncGenerator[tuple[T, ...], None]:
    if n < 1:
        raise ValueError("n must be at least one")

    iterator = _iterate(iterable)

    while True:
        batch: list[T] = []
        for _ in range(n):
            try:
                batch.append(await anext(iterator))
            except StopAsyncIteration:
                if not batch:
                    await checkpoint()
                    return
                if strict:
                    raise ValueError("batched(): incomplete batch") from None

                yield tuple(batch)
                return

        yield tuple(batch)


class Chain:
    def __call__(
        self, *iterables: Iterable[T] | AsyncIterable[T]
    ) -> AsyncGenerator[T, None]:
        return self.from_iterable(iterables)

    async def from_iterable(
        self,
        iterables: (
            Iterable[Iterable[T] | AsyncIterable[T]]
            | AsyncIterable[Iterable[T] | AsyncIterable[T]]
        ),
    ) -> AsyncGenerator[T, None]:
        element_yielded = False
        outer_iter = _iterate(iterables)

        try:
            async for iterable in outer_iter:
                async for element in _iterate(iterable):
                    element_yielded = True
                    yield element
        finally:
            aclose = getattr(outer_iter, "aclose", None)
            if aclose is not None:
                with CancelScope(shield=True):
                    await aclose()

        if not element_yielded:
            await checkpoint()


chain: Chain = Chain()


async def combinations(
    iterable: Iterable[T] | AsyncIterable[T], r: int
) -> AsyncGenerator[tuple[T, ...], None]:
    pool: list[T] = [element async for element in _iterate(iterable)]
    async for combination in _iterate(itertools.combinations(pool, r)):
        yield combination


async def combinations_with_replacement(
    iterable: Iterable[T] | AsyncIterable[T], r: int
) -> AsyncGenerator[tuple[T, ...], None]:
    pool: list[T] = [element async for element in _iterate(iterable)]
    async for combination in _iterate(itertools.combinations_with_replacement(pool, r)):
        yield combination


async def compress(
    data: Iterable[T] | AsyncIterable[T],
    selectors: Iterable[object] | AsyncIterable[object],
) -> AsyncGenerator[T, None]:
    data_iterator = _iterate(data)
    selector_iterator = _iterate(selectors)
    element_yielded = False

    while True:
        try:
            datum = await anext(data_iterator)
            selector = await anext(selector_iterator)
        except StopAsyncIteration:
            if not element_yielded:
                await checkpoint()

            return

        if selector:
            element_yielded = True
            yield datum


async def count(start: int = 0, step: int = 1) -> AsyncGenerator[int, None]:
    n = start
    while True:
        await checkpoint_if_cancelled()
        value = n
        n += step
        await cancel_shielded_checkpoint()
        yield value


async def cycle(
    iterable: Iterable[T] | AsyncIterable[T],
) -> AsyncGenerator[T, None]:
    saved: list[T] = []
    async for element in _iterate(iterable):
        saved.append(element)
        yield element

    if not saved:
        await checkpoint()
        return

    while True:
        for element in saved:
            await checkpoint()
            yield element


async def dropwhile(
    predicate: Callable[[T], Awaitable[object]],
    iterable: Iterable[T] | AsyncIterable[T],
) -> AsyncGenerator[T, None]:
    element_yielded = False
    dropping = True

    async for element in _iterate(iterable):
        if dropping and await predicate(element):
            continue

        dropping = False
        element_yielded = True
        yield element

    if not element_yielded:
        await checkpoint()


async def filterfalse(
    predicate: Callable[[T], Awaitable[object]],
    iterable: Iterable[T] | AsyncIterable[T],
) -> AsyncGenerator[T, None]:
    element_yielded = False

    async for element in _iterate(iterable):
        if not await predicate(element):
            element_yielded = True
            yield element

    if not element_yielded:
        await checkpoint()


@overload
def groupby(
    iterable: Iterable[T] | AsyncIterable[T],
) -> AsyncGenerator[tuple[T, list[T]], None]: ...


@overload
def groupby(
    iterable: Iterable[T] | AsyncIterable[T],
    key: Callable[[T], Awaitable[R]],
) -> AsyncGenerator[tuple[R, list[T]], None]: ...


async def groupby(
    iterable: Iterable[T] | AsyncIterable[T],
    key: Callable[[T], Awaitable[object]] | None = None,
) -> AsyncGenerator[tuple[object, list[T]], None]:
    iterator = _iterate(iterable)
    try:
        element = await anext(iterator)
    except StopAsyncIteration:
        await checkpoint()
        return

    group_key = element if key is None else await key(element)
    values = [element]

    async for element in iterator:
        next_key = element if key is None else await key(element)
        if next_key != group_key:
            completed_group = group_key, values
            group_key = next_key
            values = [element]
            yield completed_group
        else:
            values.append(element)

    yield group_key, values


@overload
def islice(
    iterable: Iterable[T] | AsyncIterable[T],
    stop: int | None,
    /,
) -> AsyncGenerator[T, None]: ...


@overload
def islice(
    iterable: Iterable[T] | AsyncIterable[T],
    start: int | None,
    stop: int | None,
    step: int | None = 1,
    /,
) -> AsyncGenerator[T, None]: ...


async def islice(
    iterable: Iterable[T] | AsyncIterable[T],
    *args: int | None,
) -> AsyncGenerator[T, None]:
    if not args:
        raise TypeError("islice expected at least 2 arguments, got 1")
    if len(args) > 3:
        raise TypeError(f"islice expected at most 4 arguments, got {len(args) + 1}")

    slice_args = slice(*args)

    start_message = (
        "Indices for islice() must be None or an integer: 0 <= x <= sys.maxsize."
    )
    stop_message = (
        "Stop argument for islice() must be None or an integer: 0 <= x <= sys.maxsize."
    )
    step_message = "Step for islice() must be a positive integer or None."

    def normalize_index(value: object, message: str) -> int:
        try:
            index = operator.index(cast(Any, value))
        except TypeError:
            raise ValueError(message) from None

        if index < 0 or index > sys.maxsize:
            raise ValueError(message)

        return index

    start = (
        0
        if slice_args.start is None
        else normalize_index(slice_args.start, start_message)
    )
    stop = (
        None
        if slice_args.stop is None
        else normalize_index(slice_args.stop, stop_message)
    )
    step = (
        1 if slice_args.step is None else normalize_index(slice_args.step, step_message)
    )

    if step <= 0:
        raise ValueError(step_message)

    if stop == 0 or start == stop:
        await checkpoint()
        return

    iterator = _iterate(iterable)
    index = 0
    element_yielded = False

    while stop is None or index < stop:
        try:
            element = await anext(iterator)
        except StopAsyncIteration:
            if not element_yielded:
                await checkpoint()

            return

        if index >= start and (index - start) % step == 0:
            index += 1
            element_yielded = True
            yield element
        else:
            index += 1

    if not element_yielded:
        await checkpoint()


async def pairwise(
    iterable: Iterable[T] | AsyncIterable[T],
) -> AsyncGenerator[tuple[T, T], None]:
    iterator = _iterate(iterable)
    try:
        previous = await anext(iterator)
    except StopAsyncIteration:
        await checkpoint()
        return

    element_yielded = False
    async for element in iterator:
        element_yielded = True
        pair = (previous, element)
        previous = element
        yield pair

    if not element_yielded:
        await checkpoint()


async def permutations(
    iterable: Iterable[T] | AsyncIterable[T], r: int | None = None
) -> AsyncGenerator[tuple[T, ...], None]:
    pool: list[T] = [element async for element in _iterate(iterable)]
    n = len(pool)
    if r is None:
        r = n
    elif not isinstance(r, int):
        raise TypeError("Expected int as r")
    elif r < 0:
        raise ValueError("r must be non-negative")

    async for permutation in _iterate(itertools.permutations(pool, r)):
        yield permutation


async def product(
    *iterables: Iterable[T] | AsyncIterable[T], repeat: int = 1
) -> AsyncGenerator[tuple[T, ...], None]:
    repeat = operator.index(repeat)
    if repeat < 0:
        raise ValueError("repeat argument cannot be negative")

    pools: list[tuple[T, ...]] = []
    for iterable in iterables:
        pool: list[T] = [element async for element in _iterate(iterable)]
        pools.append(tuple(pool))

    async for value in _iterate(itertools.product(*pools, repeat=repeat)):
        yield value


async def repeat(element: T, times: int | None = None) -> AsyncGenerator[T, None]:
    if times is None:
        while True:
            await checkpoint()
            yield element

    remaining = operator.index(cast(Any, times))
    if remaining <= 0:
        await checkpoint()
        return

    while remaining > 0:
        await checkpoint_if_cancelled()
        remaining -= 1
        await cancel_shielded_checkpoint()
        yield element


async def starmap(
    function: Callable[..., Awaitable[R]],
    iterable: (
        Iterable[Iterable[object] | AsyncIterable[object]]
        | AsyncIterable[Iterable[object] | AsyncIterable[object]]
    ),
) -> AsyncGenerator[R, None]:
    result_yielded = False

    async for args_iterable in _iterate(iterable):
        args = [element async for element in _iterate(args_iterable)]
        result_yielded = True
        yield await function(*args)

    if not result_yielded:
        await checkpoint()


def tee(
    iterable: Iterable[T] | AsyncIterable[T], n: int = 2
) -> tuple[AsyncIterator[T], ...]:
    n = operator.index(cast(Any, n))
    if n < 0:
        raise ValueError("n must be >= 0")
    if n == 0:
        return ()

    iterator = _TeeAsyncIterator(iterable)
    iterators: list[AsyncIterator[T]] = [iterator]
    iterators.extend(_TeeAsyncIterator(iterator) for _ in range(n - 1))
    return tuple(iterators)


async def takewhile(
    predicate: Callable[[T], Awaitable[object]],
    iterable: Iterable[T] | AsyncIterable[T],
) -> AsyncGenerator[T, None]:
    element_yielded = False

    async for element in _iterate(iterable):
        if not await predicate(element):
            if not element_yielded:
                await checkpoint()

            return

        element_yielded = True
        yield element

    if not element_yielded:
        await checkpoint()


async def zip_longest(
    *iterables: Iterable[object] | AsyncIterable[object],
    fillvalue: object = None,
) -> AsyncGenerator[tuple[object, ...], None]:
    iterators = [_iterate(iterable) for iterable in iterables]
    num_active = len(iterators)
    if not num_active:
        await checkpoint()
        return

    active = [True] * num_active
    tuple_yielded = False

    while True:
        values: list[object] = []
        for index, iterator in enumerate(iterators):
            if not active[index]:
                values.append(fillvalue)
                continue

            try:
                value = await anext(iterator)
            except StopAsyncIteration:
                active[index] = False
                num_active -= 1
                if not num_active:
                    if not tuple_yielded:
                        await checkpoint()

                    return

                value = fillvalue

            values.append(value)

        tuple_yielded = True
        yield tuple(values)
