"""Synchronization decorators and calling-convention markers/bridges.

Provides ``synchronized`` for thread and async locking, ``mark_as_sync``
and ``mark_as_async`` for declaring the effective calling convention of a
wrapped callable (without converting it), and ``async_to_sync`` /
``sync_to_async`` for bridging between the two.
"""

import asyncio
import sys
from functools import partial
from inspect import (
    CO_ASYNC_GENERATOR,
    CO_COROUTINE,
    CO_GENERATOR,
    CO_ITERABLE_COROUTINE,
    iscoroutinefunction,
)
from threading import Lock, RLock

from .__wrapt__ import BoundFunctionWrapper, CallableObjectProxy, FunctionWrapper
from .decorators import decorator

# Calling-convention marker wrappers. These manipulate __code__.co_flags
# so that inspect.iscoroutinefunction() reports the intended calling
# convention, which lets stdlib code and the synchronized() decorator
# auto-select the correct sync or async wrapping behaviour even when
# stacked decorators change the effective convention (for example an
# inner decorator that invokes an async def via asyncio.run()).


class _SyncCodeProxy(CallableObjectProxy):

    def __init__(self, wrapped, generator=None):
        super().__init__(wrapped)
        self._self_generator = generator

    @property
    def co_flags(self):
        original = self.__wrapped__.co_flags
        # Strip async-axis and iterable-coroutine bits; sync means neither
        # coroutine function nor async generator nor types.coroutine-style.
        flags = original & ~(CO_COROUTINE | CO_ASYNC_GENERATOR | CO_ITERABLE_COROUTINE)
        if self._self_generator is True:
            flags |= CO_GENERATOR
        elif self._self_generator is False:
            flags &= ~CO_GENERATOR
        else:
            # Auto: if input was an async generator, preserve generator-ness
            # on the sync side by setting CO_GENERATOR. Otherwise leave
            # CO_GENERATOR as-is (already copied from the wrapped flags).
            if original & CO_ASYNC_GENERATOR:
                flags |= CO_GENERATOR
        return flags


class _SyncFunctionSurrogate(CallableObjectProxy):

    def __init__(self, wrapped, generator=None):
        super().__init__(wrapped)
        self._self_generator = generator

    @property
    def __code__(self):
        return _SyncCodeProxy(self.__wrapped__.__code__, self._self_generator)


class _BoundSyncFunctionWrapper(BoundFunctionWrapper):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._self_is_not_coroutine = True

    @property
    def __func__(self):
        return _SyncFunctionSurrogate(
            self.__wrapped__.__func__, self._self_parent._self_generator
        )


class _SyncFunctionWrapper(FunctionWrapper):

    __bound_function_wrapper__ = _BoundSyncFunctionWrapper

    def __init__(self, wrapped, wrapper, generator=None):
        super().__init__(wrapped, wrapper)
        self._self_is_not_coroutine = True
        self._self_generator = generator

    @property
    def __code__(self):
        return _SyncCodeProxy(self.__wrapped__.__code__, self._self_generator)


class _AsyncCodeProxy(CallableObjectProxy):

    def __init__(self, wrapped, generator=None):
        super().__init__(wrapped)
        self._self_generator = generator

    @property
    def co_flags(self):
        original = self.__wrapped__.co_flags
        # Strip all four convention bits; we reassert the right ones below.
        flags = original & ~(
            CO_GENERATOR | CO_COROUTINE | CO_ITERABLE_COROUTINE | CO_ASYNC_GENERATOR
        )
        if self._self_generator is True:
            flags |= CO_ASYNC_GENERATOR
        elif self._self_generator is False:
            flags |= CO_COROUTINE
        else:
            # Auto: if input was a generator (sync or async), produce an
            # async generator; otherwise produce a coroutine function.
            if original & (CO_GENERATOR | CO_ASYNC_GENERATOR):
                flags |= CO_ASYNC_GENERATOR
            else:
                flags |= CO_COROUTINE
        return flags


class _AsyncFunctionSurrogate(CallableObjectProxy):

    def __init__(self, wrapped, generator=None):
        super().__init__(wrapped)
        self._self_generator = generator

    @property
    def __code__(self):
        return _AsyncCodeProxy(self.__wrapped__.__code__, self._self_generator)


class _BoundAsyncFunctionWrapper(BoundFunctionWrapper):

    @property
    def __func__(self):
        return _AsyncFunctionSurrogate(
            self.__wrapped__.__func__, self._self_parent._self_generator
        )


class _AsyncFunctionWrapper(FunctionWrapper):

    __bound_function_wrapper__ = _BoundAsyncFunctionWrapper

    def __init__(self, wrapped, wrapper, generator=None):
        super().__init__(wrapped, wrapper)
        self._self_generator = generator

    @property
    def __code__(self):
        return _AsyncCodeProxy(self.__wrapped__.__code__, self._self_generator)


def mark_as_sync(wrapped=None, /, *, generator=None):
    """Mark a callable as synchronous from the perspective of calling
    convention detection. The returned wrapper is a pass-through that
    reports `inspect.iscoroutinefunction()` as False regardless of
    whether the underlying callable is declared `async def`. Useful
    when a stacked decorator has already collapsed an async function
    into a synchronous one (for example by using `asyncio.run()`).

    The `generator` keyword toggles the sync generator bit
    (`CO_GENERATOR`) on the resulting wrapper. Tri-state:

    - `None` (default): auto. Preserve generator-ness from the input --
      if the input was an async generator, the wrapper reports as a sync
      generator; otherwise CO_GENERATOR is copied through unchanged.
    - `True`: force CO_GENERATOR on. Wrapper reports as a sync generator.
    - `False`: force CO_GENERATOR off. Wrapper reports as a plain sync
      function even if the input had CO_GENERATOR set.

    Regardless of `generator`, CO_COROUTINE, CO_ASYNC_GENERATOR, and
    CO_ITERABLE_COROUTINE are all cleared (sync means none of those).
    """

    def _decorator(wrapped):
        def _wrapper(wrapped, instance, args, kwargs):
            return wrapped(*args, **kwargs)

        return _SyncFunctionWrapper(wrapped, _wrapper, generator=generator)

    if wrapped is None:
        return _decorator
    return _decorator(wrapped)


def mark_as_async(wrapped=None, /, *, generator=None):
    """Mark a callable as asynchronous from the perspective of calling
    convention detection. The returned wrapper reports
    `inspect.iscoroutinefunction()` as True regardless of whether the
    underlying callable is declared `async def`. Useful when a stacked
    decorator returns a coroutine from a plain `def` wrapper.

    The `generator` keyword chooses between coroutine function and
    async generator reporting. Tri-state:

    - `None` (default): auto. If the input was a sync or async
      generator, the wrapper reports as an async generator
      (`CO_ASYNC_GENERATOR`); otherwise it reports as a coroutine
      function (`CO_COROUTINE`).
    - `True`: force async generator reporting (`CO_ASYNC_GENERATOR` set,
      `CO_COROUTINE` cleared). These two flags are mutually exclusive at
      the CPython code-object level.
    - `False`: force coroutine function reporting (`CO_COROUTINE` set,
      `CO_ASYNC_GENERATOR` cleared).

    CO_GENERATOR and CO_ITERABLE_COROUTINE are always cleared (the
    async path does not use either).
    """

    async def _wrapper(wrapped, instance, args, kwargs):
        return wrapped(*args, **kwargs)

    def _decorator(wrapped):
        return _AsyncFunctionWrapper(wrapped, _wrapper, generator=generator)

    if wrapped is None:
        return _decorator
    return _decorator(wrapped)


def async_to_sync(wrapped):
    """Adapt an async callable so it can be called synchronously. Each
    call runs the coroutine to completion via `asyncio.run()`. The
    returned wrapper reports as synchronous under
    `inspect.iscoroutinefunction()`. Naming follows the asgiref
    convention."""

    def wrapper(wrapped, instance, args, kwargs):
        return asyncio.run(wrapped(*args, **kwargs))

    return _SyncFunctionWrapper(wrapped, wrapper)


def sync_to_async(wrapped):
    """Adapt a sync callable so it can be awaited. Each call dispatches
    the synchronous work to the default executor via
    `loop.run_in_executor()`. The returned wrapper reports as
    asynchronous under `inspect.iscoroutinefunction()`. Naming follows
    the asgiref convention."""

    async def wrapper(wrapped, instance, args, kwargs):
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, partial(wrapped, *args, **kwargs))

    return _AsyncFunctionWrapper(wrapped, wrapper)


def _synchronized_is_async_lock(obj):
    return iscoroutinefunction(getattr(obj, "acquire", None))


def _synchronized_is_async_callable(obj):
    # Walk the __wrapped__ chain, returning True as soon as any layer
    # declares itself a coroutine function. A sync marker wrapper can
    # carry an authoritative `_self_is_not_coroutine` attribute that
    # short-circuits the walk before it descends into a genuinely
    # async inner layer. Cycle / runaway-chain protection modelled on
    # inspect.unwrap().

    memo = {id(obj): obj}
    recursion_limit = sys.getrecursionlimit()
    target = obj

    while True:
        if isinstance(target, (classmethod, staticmethod)):
            inner = getattr(target, "__wrapped__", None)
            if inner is None:
                inner = target.__func__
            target = inner
            id_target = id(target)
            if id_target in memo or len(memo) >= recursion_limit:
                raise ValueError("wrapper loop when unwrapping {!r}".format(obj))
            memo[id_target] = target
            continue

        if getattr(target, "_self_is_not_coroutine", False):
            return False

        if iscoroutinefunction(target):
            return True

        next_target = getattr(target, "__wrapped__", None)
        if next_target is None or next_target is target:
            return False
        target = next_target
        id_target = id(target)
        if id_target in memo or len(memo) >= recursion_limit:
            raise ValueError("wrapper loop when unwrapping {!r}".format(obj))
        memo[id_target] = target


# Decorator for implementing thread synchronization. It can be used as a
# decorator, in which case the synchronization context is determined by
# what type of function is wrapped, or it can also be used as a context
# manager, where the user needs to supply the correct synchronization
# context. It is also possible to supply an object which appears to be a
# synchronization primitive of some sort, by virtue of having release()
# and acquire() methods. In that case that will be used directly as the
# synchronization primitive without creating a separate lock against the
# derived or supplied context.


def synchronized(wrapped):
    """Depending on the nature of the `wrapped` object, will either return a
    decorator which can be used to wrap a function or method, or a context
    manager, both of which will act accordingly depending on how used, to
    synchronize access to calling of the wrapped function, or the block of
    code within the context manager. If it is an object which is a
    synchronization primitive, such as a threading Lock, RLock, Semaphore,
    Condition, or Event, then it is assumed that the object is to be used
    directly as the synchronization primitive, otherwise a lock is created
    automatically and attached to the wrapped object and used as the
    synchronization primitive.

    Async functions are supported: if the wrapped callable is an async
    function, an `asyncio.Lock` is created for the context and the wrapper
    awaits the lock. The returned object also exposes `__aenter__` and
    `__aexit__` so it can be used with `async with` to synchronise a block
    of code using an independent per-context `asyncio.Lock`. If an object
    with coroutine `acquire`/`release` methods (such as an `asyncio.Lock`)
    is supplied directly, the returned decorator and context manager will
    use it via the async protocol.
    """

    # Determine if being passed an object which is a synchronization
    # primitive. We can't check by type for Lock, RLock, Semaphore etc,
    # as the means of creating them isn't the type. Therefore use the
    # existence of acquire() and release() methods. This is more
    # extensible anyway as it allows custom synchronization mechanisms.

    if hasattr(wrapped, "acquire") and hasattr(wrapped, "release"):
        # We remember what the original lock is and then return a new
        # decorator which accesses and locks it. When returning the new
        # decorator we wrap it with an object proxy so we can override
        # the context manager methods in case it is being used to wrap
        # synchronized statements with a 'with' statement.

        lock = wrapped

        if _synchronized_is_async_lock(lock):

            @decorator
            async def _synchronized(wrapped, instance, args, kwargs):
                async with lock:
                    return await wrapped(*args, **kwargs)

            class _AsyncSynchronizedLockProxy(CallableObjectProxy):

                async def __aenter__(self):
                    await lock.acquire()
                    return lock

                async def __aexit__(self, *args):
                    lock.release()

            return _AsyncSynchronizedLockProxy(wrapped=_synchronized)

        @decorator
        def _synchronized(wrapped, instance, args, kwargs):
            # Execute the wrapped function while the original supplied
            # lock is held.

            with lock:
                return wrapped(*args, **kwargs)

        class _SynchronizedLockProxy(CallableObjectProxy):

            def __enter__(self):
                lock.acquire()
                return lock

            def __exit__(self, *args):
                lock.release()

        return _SynchronizedLockProxy(wrapped=_synchronized)

    # Following only apply when the lock is being created automatically
    # based on the context of what was supplied. In this case we supply
    # a final decorator, but need to use FunctionWrapper directly as we
    # want to derive from it to add context manager methods in case it is
    # being used to wrap synchronized statements with a 'with' statement.

    def _synchronized_lock(context):
        # Attempt to retrieve the lock for the specific context.

        lock = vars(context).get("_synchronized_lock", None)

        if lock is None:
            # There is no existing lock defined for the context we
            # are dealing with so we need to create one. This needs
            # to be done in a way to guarantee there is only one
            # created, even if multiple threads try and create it at
            # the same time. We can't always use the setdefault()
            # method on the __dict__ for the context. This is the
            # case where the context is a class, as __dict__ is
            # actually a dictproxy. What we therefore do is use a
            # meta lock on this wrapper itself, to control the
            # creation and assignment of the lock attribute against
            # the context.

            with synchronized._synchronized_meta_lock:
                # We need to check again for whether the lock we want
                # exists in case two threads were trying to create it
                # at the same time and were competing to create the
                # meta lock.

                lock = vars(context).get("_synchronized_lock", None)

                if lock is None:
                    lock = RLock()
                    setattr(context, "_synchronized_lock", lock)

        return lock

    def _synchronized_async_lock(context):
        # Per-context asyncio.Lock, created lazily on first use. Created
        # under the shared meta lock so creation is safe across threads;
        # the meta lock is never held across an await. asyncio.Lock is
        # not reentrant.

        lock = vars(context).get("_synchronized_async_lock", None)

        if lock is None:
            with synchronized._synchronized_meta_lock:
                lock = vars(context).get("_synchronized_async_lock", None)

                if lock is None:
                    lock = asyncio.Lock()
                    setattr(context, "_synchronized_async_lock", lock)

        return lock

    def _synchronized_wrapper(wrapped, instance, args, kwargs):
        # Execute the wrapped function while the lock for the
        # desired context is held. If instance is None then the
        # wrapped function is used as the context.

        with _synchronized_lock(instance if instance is not None else wrapped):
            return wrapped(*args, **kwargs)

    async def _synchronized_async_wrapper(wrapped, instance, args, kwargs):
        async with _synchronized_async_lock(
            instance if instance is not None else wrapped
        ):
            return await wrapped(*args, **kwargs)

    class _SynchronizedFunctionWrapper(FunctionWrapper):

        def __enter__(self):
            self._self_lock = _synchronized_lock(self.__wrapped__)
            self._self_lock.acquire()
            return self._self_lock

        def __exit__(self, *args):
            self._self_lock.release()

        async def __aenter__(self):
            self._self_async_lock = _synchronized_async_lock(self.__wrapped__)
            await self._self_async_lock.acquire()
            return self._self_async_lock

        async def __aexit__(self, *args):
            self._self_async_lock.release()

    if _synchronized_is_async_callable(wrapped):
        return _SynchronizedFunctionWrapper(
            wrapped=wrapped, wrapper=_synchronized_async_wrapper
        )

    return _SynchronizedFunctionWrapper(wrapped=wrapped, wrapper=_synchronized_wrapper)


synchronized._synchronized_meta_lock = Lock()  # type: ignore[attr-defined]
