"""Async wrapper around :class:`SoftReadWriteLock` for use with ``asyncio``."""

from __future__ import annotations

import asyncio
import functools
import os
from contextlib import asynccontextmanager
from typing import TYPE_CHECKING, ParamSpec, TypeVar

from ._sync import SoftReadWriteLock

if TYPE_CHECKING:
    from collections.abc import AsyncGenerator, Callable
    from concurrent import futures
    from types import TracebackType

_P = ParamSpec("_P")
_R = TypeVar("_R")


class AsyncSoftReadWriteLock:
    """
    Async wrapper around :class:`SoftReadWriteLock` for ``asyncio`` applications.

    The sync class's blocking filesystem operations run on a thread pool via ``loop.run_in_executor()``. The
    underlying :class:`SoftReadWriteLock` handles reentrancy, upgrade/downgrade rules, fork handling, heartbeat and
    TTL stale detection, and singleton behavior.

    :param lock_file: path to the lock file; sidecar state/write/readers live next to it
    :param timeout: maximum wait time in seconds; ``-1`` means block indefinitely
    :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately on contention
    :param is_singleton: if ``True``, reuse existing :class:`SoftReadWriteLock` instances per resolved path
    :param heartbeat_interval: seconds between heartbeat refreshes; default 30 s
    :param stale_threshold: seconds of mtime inactivity before a marker is stale; defaults to ``3 * heartbeat_interval``
    :param poll_interval: seconds between acquire retries under contention; default 0.25 s
    :param loop: event loop for ``run_in_executor``; ``None`` uses the running loop
    :param executor: executor for ``run_in_executor``; ``None`` uses the default executor

    .. versionadded:: 3.27.0

    """

    def __init__(  # ruff:ignore[too-many-arguments]  # public constructor: one parameter per documented lock option
        self,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        *,
        blocking: bool = True,
        is_singleton: bool = True,
        heartbeat_interval: float = 30.0,
        stale_threshold: float | None = None,
        poll_interval: float = 0.25,
        loop: asyncio.AbstractEventLoop | None = None,
        executor: futures.Executor | None = None,
    ) -> None:
        self._creator_pid = os.getpid()
        self._lock = SoftReadWriteLock(
            lock_file,
            timeout,
            blocking=blocking,
            is_singleton=is_singleton,
            heartbeat_interval=heartbeat_interval,
            stale_threshold=stale_threshold,
            poll_interval=poll_interval,
        )
        self._loop = loop
        self._executor = executor

    @property
    def lock_file(self) -> str:
        """The path to the lock file passed to the constructor."""
        return self._lock.lock_file

    @property
    def timeout(self) -> float:
        """The default timeout applied when ``acquire_read`` / ``acquire_write`` is called without one."""
        return self._lock.timeout

    @property
    def blocking(self) -> bool:
        """Whether ``acquire_*`` defaults to blocking; ``False`` makes contention raise immediately."""
        return self._lock.blocking

    @property
    def loop(self) -> asyncio.AbstractEventLoop | None:
        """The event loop used for ``run_in_executor``, or ``None`` for the running loop."""
        return self._loop

    @property
    def executor(self) -> futures.Executor | None:
        """The executor used for ``run_in_executor``, or ``None`` for the default executor."""
        return self._executor

    @asynccontextmanager
    async def read_lock(self, timeout: float | None = None, *, blocking: bool | None = None) -> AsyncGenerator[None]:
        """
        Async context manager that acquires and releases a shared read lock.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately; ``None`` uses the instance default

        :raises RuntimeError: if a write lock is already held on this instance
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        await self.acquire_read(timeout, blocking=blocking)
        try:
            yield
        finally:
            await self.release()

    @asynccontextmanager
    async def write_lock(self, timeout: float | None = None, *, blocking: bool | None = None) -> AsyncGenerator[None]:
        """
        Async context manager that acquires and releases an exclusive write lock.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately; ``None`` uses the instance default

        :raises RuntimeError: if a read lock is already held, or a write lock is held by a different thread
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        await self.acquire_write(timeout, blocking=blocking)
        try:
            yield
        finally:
            await self.release()

    async def acquire_read(
        self, timeout: float | None = None, *, blocking: bool | None = None
    ) -> AsyncAcquireSoftReadWriteReturnProxy:
        """
        Acquire a shared read lock.

        See :meth:`SoftReadWriteLock.acquire_read` for reentrancy / upgrade / fork semantics. The blocking work runs
        inside ``run_in_executor`` so other coroutines on the same loop keep progressing while this call waits.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately; ``None`` uses the instance default

        :returns: a proxy usable as an async context manager to release the lock

        :raises RuntimeError: if a write lock is already held, if this instance was invalidated by
            :func:`os.fork`, or if :meth:`close` was called
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        self._raise_if_inherited()
        await self._run(self._lock.acquire_read, timeout, blocking=blocking)
        return AsyncAcquireSoftReadWriteReturnProxy(lock=self)

    async def acquire_write(
        self, timeout: float | None = None, *, blocking: bool | None = None
    ) -> AsyncAcquireSoftReadWriteReturnProxy:
        """
        Acquire an exclusive write lock.

        See :meth:`SoftReadWriteLock.acquire_write` for the two-phase writer-preferring semantics. The blocking work
        runs inside ``run_in_executor``.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately; ``None`` uses the instance default

        :returns: a proxy usable as an async context manager to release the lock

        :raises RuntimeError: if a read lock is already held, if a write lock is held by a different thread, if
            this instance was invalidated by :func:`os.fork`, or if :meth:`close` was called
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        self._raise_if_inherited()
        await self._run(self._lock.acquire_write, timeout, blocking=blocking)
        return AsyncAcquireSoftReadWriteReturnProxy(lock=self)

    async def release(self, *, force: bool = False) -> None:
        """
        Release one level of the current lock.

        :param force: if ``True``, release the lock completely regardless of the current lock level

        :raises RuntimeError: if no lock is currently held and *force* is ``False``

        """
        if self._creator_pid == os.getpid():
            await self._run(self._lock.release, force=force)

    async def close(self) -> None:
        """Release any held lock and release the underlying filesystem resources. Idempotent."""
        if self._creator_pid == os.getpid():
            await self._run(self._lock.close)

    def _raise_if_inherited(self) -> None:
        if self._creator_pid != os.getpid():  # pragma: forked child
            msg = f"AsyncSoftReadWriteLock on {self.lock_file} was inherited across fork; construct a new instance"
            raise RuntimeError(msg)

    async def _run(self, func: Callable[_P, _R], *args: _P.args, **kwargs: _P.kwargs) -> _R:
        loop = self._loop or asyncio.get_running_loop()
        return await loop.run_in_executor(self._executor, functools.partial(func, *args, **kwargs))


class AsyncAcquireSoftReadWriteReturnProxy:
    """Async context-aware object that releases an :class:`AsyncSoftReadWriteLock` on exit."""

    def __init__(self, lock: AsyncSoftReadWriteLock) -> None:
        self.lock = lock

    async def __aenter__(self) -> AsyncSoftReadWriteLock:
        return self.lock

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        await self.lock.release()


__all__ = [
    "AsyncAcquireSoftReadWriteReturnProxy",
    "AsyncSoftReadWriteLock",
]
