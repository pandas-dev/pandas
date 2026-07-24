"""An asyncio-based implementation of the file lock."""

from __future__ import annotations

import asyncio
import contextlib
import logging
import os
import time
from dataclasses import dataclass
from inspect import iscoroutinefunction
from threading import local
from typing import TYPE_CHECKING, Final, NoReturn, TypeVar, cast

from ._api import (
    _UNSET_FILE_MODE,
    BaseFileLock,
    CloseErrorPolicy,
    ContextErrorPolicy,
    FileLockContext,
    FileLockMeta,
    _append_exception_context,
    _canonical,
    _ExtraValue,
    _fork_transition,
    _grouped_errors,
    _raise_body_and_release,
    _raise_chained_errors,
    _raise_cleanup_errors,
    _raise_grouped_errors,
    _register_fork_object,
)
from ._async import (
    _AsyncTransitionGate,
    _AsyncTransitionUnavailableError,
    _BackendOutcome,
    _capture_awaitable,
    _capture_call,
    _drain_future,
    _future_result,
    _wait_until_done,
)
from ._error import Timeout
from ._lease import SoftFileLease
from ._soft import SoftFileLock
from ._strict import StrictSoftFileLock
from ._unix import UnixFileLock
from ._windows import WindowsFileLock

if TYPE_CHECKING:
    import sys
    from collections.abc import Awaitable, Callable, Coroutine, Hashable
    from concurrent import futures
    from types import TracebackType

    if sys.version_info >= (3, 11):  # pragma: no cover (py311+)
        from typing import Self
    else:  # pragma: no cover (<py311)
        from typing_extensions import Self


_LOGGER: Final[logging.Logger] = logging.getLogger("filelock")
_ASYNC_RELEASE_CANCELLATION_ERRORS: Final[str] = "lock release cancellation and backend release both failed"
_ASYNC_CONTEXT_RELEASE_ERRORS: Final[str] = "context body, release cancellation, and backend release failed"
_ASYNC_RELEASE_CANCELLATION_MARKER_ATTR: Final[str] = "_filelock_async_release_cancellation"
_ASYNC_RELEASE_CANCELLATION_MARKER: Final[list[None]] = []

_AT = TypeVar("_AT", bound="BaseAsyncFileLock")


class AsyncFileLockMeta(FileLockMeta):
    def __call__(  # ruff:ignore[too-many-arguments]  # forwards the public constructor's documented parameters
        cls: type[_AT],  # ruff:ignore[invalid-first-argument-name-for-method]  # metaclass __call__ receives the class being constructed
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        mode: int = _UNSET_FILE_MODE,
        thread_local: bool = False,  # ruff:ignore[boolean-type-hint-positional-argument, boolean-default-value-positional-argument]  # public API: positional bool kept for backwards compatibility
        *,
        blocking: bool = True,
        is_singleton: bool = False,
        poll_interval: float = 0.05,
        lifetime: float | None = None,
        context_error_policy: ContextErrorPolicy = "chain",
        close_error_policy: CloseErrorPolicy = "default",
        fallback_to_soft: bool = True,
        preserve_lock_file: bool = False,
        on_acquired: Callable[[int], None] | None = None,
        loop: asyncio.AbstractEventLoop | None = None,
        run_in_executor: bool = True,
        executor: futures.Executor | None = None,
        **kwargs: _ExtraValue,
    ) -> _AT:
        if thread_local and run_in_executor:
            msg = "run_in_executor is not supported when thread_local is True"
            raise ValueError(msg)
        return super().__call__(  # a subclass may add options of its own, as AsyncSoftFileLease does
            **kwargs,
            lock_file=lock_file,
            timeout=timeout,
            mode=mode,
            thread_local=thread_local,
            blocking=blocking,
            is_singleton=is_singleton,
            poll_interval=poll_interval,
            lifetime=lifetime,
            context_error_policy=context_error_policy,
            close_error_policy=close_error_policy,
            fallback_to_soft=fallback_to_soft,
            preserve_lock_file=preserve_lock_file,
            on_acquired=on_acquired,
            loop=loop,
            run_in_executor=run_in_executor,
            executor=executor,
        )


class BaseAsyncFileLock(BaseFileLock, metaclass=AsyncFileLockMeta):
    """
    Base class for asynchronous file locks.

    .. versionadded:: 3.15.0

    """

    _deadlock_holder_desc: str = "BaseAsyncFileLock instance in this task"
    _constructor_lifetime_warning_stacklevel: int = 4

    @staticmethod
    def _deadlock_scope() -> Hashable | None:
        # One event loop thread runs every task, so a thread-scoped registry cannot tell a same-task reacquire
        # (a real deadlock: the polling task never reaches its own release) from another task queuing behind the
        # holder (no deadlock: each poll yields, so the holder runs on and releases). Only the first may fail
        # fast, so scope holders to the task.
        return asyncio.current_task()

    def __init__(  # ruff:ignore[too-many-arguments]  # public constructor: one parameter per documented lock option
        self,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        mode: int = _UNSET_FILE_MODE,
        thread_local: bool = False,  # ruff:ignore[boolean-type-hint-positional-argument, boolean-default-value-positional-argument]  # public API: positional bool kept for backwards compatibility
        *,
        blocking: bool = True,
        is_singleton: bool = False,
        poll_interval: float = 0.05,
        lifetime: float | None = None,
        context_error_policy: ContextErrorPolicy = "chain",
        close_error_policy: CloseErrorPolicy = "default",
        fallback_to_soft: bool = True,
        preserve_lock_file: bool = False,
        on_acquired: Callable[[int], None] | None = None,
        loop: asyncio.AbstractEventLoop | None = None,
        run_in_executor: bool = True,
        executor: futures.Executor | None = None,
    ) -> None:
        """
        Create a new lock object.

        :param lock_file: path to the file
        :param timeout: default timeout when acquiring the lock, in seconds. It will be used as fallback value in the
            acquire method, if no timeout value (``None``) is given. If you want to disable the timeout, set it to a
            negative value. A timeout of 0 means that there is exactly one attempt to acquire the file lock.
        :param mode: file permissions for the lockfile. When not specified, the OS controls permissions via umask and
            default ACLs, preserving POSIX default ACL inheritance in shared directories.
        :param thread_local: Whether this object's internal context should be thread local or not. If this is set to
            ``False`` then the lock will be reentrant across threads. When ``True`` (the default), **all fields of the
            lock's internal context are per-thread**, including the configuration values ``poll_interval``, ``timeout``,
            ``blocking``, ``mode``, and ``lifetime``. Setting one of these properties from one thread does not change
            the value seen by another thread; threads that did not perform the write continue to see the value supplied
            at construction time. If you need configuration values to be visible across threads, construct the lock
            with ``thread_local=False``.
        :param blocking: whether the lock should be blocking or not
        :param is_singleton: If this is set to ``True`` then only one instance of this class will be created per lock
            file. This is useful if you want to use the lock object for reentrant locking without needing to pass the
            same object around.
        :param poll_interval: default interval for polling the lock file, in seconds. It will be used as fallback value
            in the acquire method, if no poll_interval value (``None``) is given.
        :param lifetime: for :class:`AsyncSoftFileLock`, the age in seconds after which a waiting process may delete
            the marker, even while its holder remains alive. This legacy expiry mode does not provide strict mutual
            exclusion. ``None`` (the default) disables age-based expiry. Native OS locks (:class:`AsyncFileLock`)
            cannot be revoked by file age and ignore a non-``None`` ``lifetime`` with a warning.
        :param context_error_policy: how a context manager reconciles a failure in its body with a failure while
            releasing on exit. ``"chain"`` (the default) keeps Python's behavior: the release error propagates with the
            body error in its ``__context__``. ``"group"`` raises a :class:`BaseExceptionGroup` holding the body error
            first and the release error second, so neither hides the other.
        :param close_error_policy: for native locks (:class:`AsyncFileLock`), what to do with an ``os.close`` failure
            after the OS unlock has already committed. ``"default"`` keeps each platform's historical behavior,
            ``"raise"`` always propagates the ``OSError``, and ``"suppress"`` always ignores it.
        :param fallback_to_soft: for :class:`AsyncFileLock`, whether to fall back to soft existence locking when
            ``flock`` returns ``ENOSYS``. ``True`` (default) keeps the fallback; ``False`` propagates the error.
        :param preserve_lock_file: for native locks (:class:`AsyncFileLock`), whether filelock promises not to unlink
            the lock pathname on release. ``False`` (default) keeps each backend's cleanup; ``True`` keeps a stable file
            identity (Windows skips its unlink, Unix refuses the ``ENOSYS`` soft fallback). :class:`AsyncSoftFileLock`
            rejects ``True``.
        :param on_acquired: for native locks (:class:`AsyncFileLock`), a callable invoked with the borrowed lock
            descriptor once per physical acquisition, after the lock is held but before
            :meth:`~BaseAsyncFileLock.acquire` returns. With ``run_in_executor=True`` (the default) it runs in the
            backend executor. It must not close or unlock the descriptor; a raise rolls the acquisition back.
            :class:`AsyncSoftFileLock` rejects it.
        :param loop: The event loop to use. If not specified, the running event loop will be used.
        :param run_in_executor: If this is set to ``True`` then the lock will be acquired in an executor.
        :param executor: The executor to use. If not specified, the default executor will be used.

        """
        self._creator_pid = os.getpid()
        self._is_thread_local = thread_local
        self._is_singleton = is_singleton
        self._context_error_policy = context_error_policy  # already validated by the metaclass
        self._close_error_policy = close_error_policy  # already validated by the metaclass
        self._fallback_to_soft = fallback_to_soft
        self._preserve_lock_file = preserve_lock_file  # already validated by the metaclass
        self._on_acquired = on_acquired  # already validated by the metaclass
        self._transition_gate: Final[_AsyncTransitionGate] = _AsyncTransitionGate()

        self._context: AsyncFileLockContext = (AsyncThreadLocalFileContext if thread_local else AsyncFileLockContext)(
            lock_file=os.fspath(lock_file),
            timeout=timeout,
            mode=mode,
            blocking=blocking,
            poll_interval=poll_interval,
            lifetime=lifetime,
            loop=loop,
            run_in_executor=run_in_executor,
            executor=executor,
        )
        _register_fork_object(self)

    @property
    def run_in_executor(self) -> bool:
        """Whether run in executor."""
        return self._context.run_in_executor

    @property
    def executor(self) -> futures.Executor | None:
        """The executor."""
        return self._context.executor

    @executor.setter
    def executor(self, value: futures.Executor | None) -> None:  # pragma: no cover
        """
        Change the executor.

        :param futures.Executor | None value: the new executor or ``None``

        """
        self._context.executor = value

    @property
    def loop(self) -> asyncio.AbstractEventLoop | None:
        """The event loop."""
        return self._context.loop

    async def acquire(  # ty: ignore[invalid-method-override]
        self,
        timeout: float | None = None,
        poll_interval: float | None = None,
        *,
        blocking: bool | None = None,
        cancel_check: Callable[[], bool] | None = None,
    ) -> AsyncAcquireReturnProxy:
        """
        Try to acquire the file lock.

        :param timeout: maximum wait time for acquiring the lock, ``None`` means use the default
            :attr:`~BaseFileLock.timeout` is and if ``timeout < 0``, there is no timeout and this method will block
            until the lock could be acquired
        :param poll_interval: interval of trying to acquire the lock file, ``None`` means use the default
            :attr:`~BaseFileLock.poll_interval`
        :param blocking: defaults to True. If False, function will return immediately if it cannot obtain a lock on the
            first attempt. Otherwise, this method will block until the timeout expires or the lock is acquired.
        :param cancel_check: a callable returning ``True`` when the acquisition should be canceled. Checked on each poll
            iteration. When triggered, raises :class:`~Timeout` just like an expired timeout.

        :returns: a context object that will unlock the file when the context is exited

        :raises Timeout: if fails to acquire lock within the timeout period

        .. code-block:: python

            # You can use this method in the context manager (recommended)
            with lock.acquire():
                pass

            # Or use an equivalent try-finally construct:
            lock.acquire()
            try:
                pass
            finally:
                lock.release()

        """
        self._raise_if_inherited()
        if timeout is None:
            timeout = self._context.timeout

        if blocking is None:
            blocking = self._context.blocking

        if poll_interval is None:
            poll_interval = self._context.poll_interval

        start_time = time.perf_counter()
        try:
            return await self._acquire_with_admission(
                blocking=blocking,
                cancel_check=cancel_check,
                timeout=timeout,
                poll_interval=poll_interval,
                start_time=start_time,
            )
        except _AsyncTransitionUnavailableError:
            raise Timeout(self.lock_file) from None

    async def _acquire_with_admission(
        self,
        *,
        blocking: bool,
        cancel_check: Callable[[], bool] | None,
        timeout: float,
        poll_interval: float,
        start_time: float,
    ) -> AsyncAcquireReturnProxy:
        async with self._transition_gate.hold_for_acquire(
            blocking=blocking,
            cancel_check=cancel_check,
            deadline=None if timeout < 0 else start_time + timeout,
            poll_interval=poll_interval,
        ):
            # A canceled provisional acquire must finish rollback before another caller can claim its descriptor.
            canonical = _canonical(self.lock_file)
            self._context.lock_counter += 1
            self._raise_if_would_deadlock(canonical, timeout=timeout, blocking=blocking)
            self._context.claim_root = canonical
            try:
                await self._async_poll_until_acquired(
                    blocking=blocking,
                    cancel_check=cancel_check,
                    timeout=timeout,
                    poll_interval=poll_interval,
                    start_time=start_time,
                )
            except BaseException:
                self._reconcile_failed_acquire(canonical)
                raise
            finally:
                self._context.claim_root = None
            self._commit_acquire(canonical)
            return AsyncAcquireReturnProxy(lock=self)

    async def _async_poll_until_acquired(
        self,
        *,
        blocking: bool,
        cancel_check: Callable[[], bool] | None,
        timeout: float,
        poll_interval: float,
        start_time: float,
    ) -> None:
        lock_id = id(self)
        lock_filename = self.lock_file
        while True:
            self._raise_if_inherited()
            if not self.is_locked:
                self._try_break_expired_lock()
                _LOGGER.debug("Attempting to acquire lock %s on %s", lock_id, lock_filename)
                await self._run_acquire_attempt()
                self._raise_if_inherited()
            if self.is_locked:
                _LOGGER.debug("Lock %s acquired on %s", lock_id, lock_filename)
                return
            if self._check_give_up(
                blocking=blocking,
                cancel_check=cancel_check,
                timeout=timeout,
                start_time=start_time,
            ):
                raise Timeout(lock_filename)
            _LOGGER.debug("Lock %s not acquired on %s, waiting %s seconds ...", lock_id, lock_filename, poll_interval)
            await asyncio.sleep(poll_interval)

    async def _run_acquire_attempt(self) -> None:
        acquire_future = self._start_internal_method(
            self._acquire_with_fork_tracking_async
            if iscoroutinefunction(self._acquire)
            else self._acquire_with_fork_tracking
        )
        try:
            await _wait_until_done(acquire_future)
        except asyncio.CancelledError as cancellation:
            acquire_error: BaseException | None = None
            try:
                await _drain_future(acquire_future)
            except BaseException as error:  # ruff:ignore[blind-except]  # reported with the cancellation below
                acquire_error = error

            rollback_error: BaseException | None = None
            if self.is_locked:
                try:
                    await _drain_future(self._start_tracked_release())
                except BaseException as error:  # ruff:ignore[blind-except]  # reported with the cancellation below
                    rollback_error = error

            if acquire_error is not None:
                if rollback_error is not None:  # pragma: needs fcntl
                    self._raise_cancelled_errors(
                        "lock acquisition cancellation, backend attempt, and rollback failed",
                        cancellation,
                        acquire_error,
                        rollback_error,
                    )
                self._raise_cancelled_errors(
                    "lock acquisition cancellation and backend attempt both failed", cancellation, acquire_error
                )
            if rollback_error is not None:
                self._raise_cancelled_errors(
                    "lock acquisition cancellation and rollback both failed", cancellation, rollback_error
                )
            raise
        try:
            _future_result(acquire_future)
        except asyncio.CancelledError as acquire_error:
            await self._rollback_backend_cancelled_acquire(acquire_error)

    async def _rollback_backend_cancelled_acquire(self, acquire_error: asyncio.CancelledError) -> NoReturn:
        if self.is_locked:
            try:
                await _drain_future(self._start_tracked_release())
            except BaseException as rollback_error:  # ruff:ignore[blind-except]  # both backend errors must surface
                self._raise_acquire_rollback_errors(acquire_error, rollback_error)
        raise acquire_error

    def _raise_acquire_rollback_errors(self, acquire_error: BaseException, rollback_error: BaseException) -> NoReturn:
        if self._context_error_policy == "group":
            _raise_grouped_errors("lock acquisition backend and rollback both failed", acquire_error, rollback_error)
        _raise_chained_errors(acquire_error, rollback_error)

    def _raise_cancelled_errors(
        self,
        message: str,
        cancellation: asyncio.CancelledError,
        first_error: BaseException,
        second_error: BaseException | None = None,
    ) -> NoReturn:
        if self._context_error_policy == "group":
            marker = (
                (_ASYNC_RELEASE_CANCELLATION_MARKER_ATTR, _ASYNC_RELEASE_CANCELLATION_MARKER)
                if message == _ASYNC_RELEASE_CANCELLATION_ERRORS
                else None
            )
            if second_error is None:
                _raise_grouped_errors(message, cancellation, first_error, marker=marker)
            _raise_grouped_errors(message, cancellation, first_error, second_error, marker=marker)
        if (context := first_error.__context__) is not None and context is not cancellation:
            if (cancellation_context := cancellation.__context__) is not None:
                _append_exception_context(context, cancellation_context)
            cancellation.__context__ = context
        first_error.__context__ = cancellation
        _raise_chained_errors(first_error, second_error)

    async def release(self, force: bool = False) -> None:  # ty: ignore[invalid-method-override]  # ruff:ignore[boolean-type-hint-positional-argument, boolean-default-value-positional-argument]  # public API: positional bool kept for backwards compatibility
        """
        Release the file lock. The lock is only completely released when the lock counter reaches 0. The lock file
        itself may be deleted automatically, the behavior is platform-specific.

        :param force: If true, the lock counter is ignored and the lock is released in every case.

        """
        if self._creator_pid != os.getpid() or not self.is_locked:
            return
        async with self._transition_gate.hold():
            await self._release_serialized(force=force)

    async def _release_serialized(self, *, force: bool) -> None:
        if not self.is_locked:
            return
        if not force and self._context.lock_counter > 1:
            self._context.lock_counter -= 1
            return

        lock_id, lock_filename = id(self), self.lock_file
        _LOGGER.debug("Attempting to release lock %s on %s", lock_id, lock_filename)
        release_future = self._start_tracked_release()
        try:
            await _wait_until_done(release_future)
        except asyncio.CancelledError as cancellation:
            try:
                await _drain_future(release_future)
            except BaseException as release_error:  # ruff:ignore[blind-except]  # cancellation and backend failure must both surface
                self._commit_release_if_released()
                self._raise_cancelled_errors(_ASYNC_RELEASE_CANCELLATION_ERRORS, cancellation, release_error)
            self._commit_release()
            raise
        try:
            _future_result(release_future)
        except BaseException:  # state follows the backend for control-flow exceptions too
            self._commit_release_if_released()
            raise
        self._commit_release()
        _LOGGER.debug("Lock %s released on %s", lock_id, lock_filename)

    def _commit_release_if_released(self) -> None:
        # Commit only when the backend actually unlocked (close or unlink failed after the OS unlock). If the lock is
        # still held, keep the counter so a later release can retry.
        if not self.is_locked:
            self._commit_release()

    def _start_internal_method(
        self, method: Callable[[], None] | Callable[[], Coroutine[None, None, None]]
    ) -> asyncio.Future[_BackendOutcome[None]]:
        if iscoroutinefunction(method):
            return asyncio.create_task(_capture_awaitable(cast("Callable[[], Coroutine[None, None, None]]", method)()))
        loop = asyncio.get_running_loop()
        sync_method = cast("Callable[[], None]", method)
        if self.run_in_executor:
            return loop.run_in_executor(self.executor, _capture_call, sync_method)
        future: asyncio.Future[_BackendOutcome[None]] = loop.create_future()
        future.set_result(_capture_call(sync_method))
        return future

    def _start_tracked_release(self) -> asyncio.Future[_BackendOutcome[None]]:
        return self._start_internal_method(
            self._release_with_fork_tracking_async
            if iscoroutinefunction(self._release)
            else self._release_with_fork_tracking
        )

    async def _acquire_with_fork_tracking_async(self) -> None:
        with _fork_transition(self):
            try:
                await cast("Callable[[], Awaitable[None]]", self._acquire)()
            except BaseException as acquisition_error:
                await self._rollback_failed_acquire_async(acquisition_error)
                raise
            try:
                self._register_context_descriptor()
            except BaseException as registration_error:  # cancellation must roll back the descriptor
                await self._rollback_failed_registration_async(registration_error)
                raise
        if self.is_locked:
            await self._invoke_on_acquired_async()

    async def _rollback_failed_acquire_async(self, acquisition_error: BaseException) -> None:
        if not self.is_locked:
            return
        registration_error: BaseException | None = None
        tracking_error: BaseException | None = None
        try:
            self._register_context_descriptor()
        except BaseException as error:  # ruff:ignore[blind-except]  # preserve registration and acquisition failures
            registration_error = error
            try:
                # Rollback may fail too; retain the fd so a child can close it without another identity probe.
                self._register_unverified_context_descriptor()
            except BaseException as error:  # ruff:ignore[blind-except]  # pragma: no cover - allocation/control-flow during fallback
                tracking_error = error
        try:
            await _drain_future(self._start_tracked_release())
        except BaseException as rollback_error:  # ruff:ignore[blind-except]  # preserve rollback and acquisition failures
            if registration_error is None and tracking_error is None:
                self._raise_acquire_rollback_errors(acquisition_error, rollback_error)
            _raise_cleanup_errors(
                "lock acquisition cleanup failed",
                acquisition_error,
                registration_error,
                tracking_error,
                rollback_error,
            )
        if registration_error is not None:
            _raise_cleanup_errors(
                "lock acquisition cleanup failed", acquisition_error, registration_error, tracking_error
            )

    async def _rollback_failed_registration_async(self, registration_error: BaseException) -> None:
        tracking_error: BaseException | None = None
        try:
            # Rollback may fail too; retain the fd so a child can close it without another identity probe.
            self._register_unverified_context_descriptor()
        except BaseException as error:  # ruff:ignore[blind-except]  # pragma: no cover - allocation/control-flow during fallback
            tracking_error = error
        try:
            await _drain_future(self._start_tracked_release())
        except BaseException as rollback_error:  # ruff:ignore[blind-except]  # preserve rollback and registration failures
            _raise_cleanup_errors(
                "descriptor registration cleanup failed", registration_error, tracking_error, rollback_error
            )
        if tracking_error is not None:  # pragma: no cover - requires failed in-memory fallback
            _raise_cleanup_errors("descriptor registration cleanup failed", registration_error, tracking_error)

    async def _invoke_on_acquired_async(self) -> None:
        if self._on_acquired is None or self._context.lock_counter != 1:
            return
        try:
            self._on_acquired(cast("int", self._context.lock_file_fd))
        except BaseException as callback_error:  # caller control-flow errors must release the lock
            callback_context = callback_error.__context__
            try:
                await _drain_future(self._start_tracked_release())
            except BaseException as release_error:  # ruff:ignore[blind-except]  # both errors surface via the group below
                _raise_body_and_release(callback_error, release_error)
            callback_error.__context__ = callback_context
            raise

    async def _release_with_fork_tracking_async(self) -> None:
        with _fork_transition(self):
            try:
                await cast("Callable[[], Awaitable[None]]", self._release)()
            finally:
                self._unregister_released_descriptor()

    def __enter__(self) -> NoReturn:
        """Sync context manager entry is not supported because lock acquisition is a coroutine."""
        msg = "Use `async with`: acquire/release are coroutines and cannot be awaited in a sync context manager."
        raise NotImplementedError(msg)

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Sync context manager exit is not supported because lock release is a coroutine."""
        msg = "Use `async with`: acquire/release are coroutines and cannot be awaited in a sync context manager."
        raise NotImplementedError(msg)

    async def __aenter__(self) -> Self:
        """
        Acquire the lock.

        :returns: the lock object

        """
        await self.acquire()
        return self

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """
        Release the lock, reconciling a release failure with any body failure per :attr:`context_error_policy`.

        :param exc_type: the exception type if raised
        :param exc_value: the exception value if raised
        :param traceback: the exception traceback if raised

        """
        await self._release_in_context(exc_value)

    async def _release_in_context(  # ty: ignore[invalid-method-override]
        self, body_error: BaseException | None
    ) -> None:
        # The async counterpart of BaseFileLock._release_in_context: await release, then apply the same policy.
        try:
            await self.release()
        except BaseException as release_error:
            if body_error is None:
                raise
            if self._context_error_policy == "chain":
                _append_exception_context(release_error, body_error)
                raise
            if (
                errors := _grouped_errors(
                    release_error,
                    _ASYNC_RELEASE_CANCELLATION_ERRORS,
                    (_ASYNC_RELEASE_CANCELLATION_MARKER_ATTR, _ASYNC_RELEASE_CANCELLATION_MARKER),
                )
            ) is not None:
                match errors:
                    # The marker is only ever set on a (CancelledError, backend_error) pair, so this always matches.
                    case (asyncio.CancelledError() as cancellation, backend_error):  # pragma: no branch
                        _raise_grouped_errors(_ASYNC_CONTEXT_RELEASE_ERRORS, body_error, cancellation, backend_error)
            _raise_body_and_release(body_error, release_error)

    def __del__(self) -> None:
        """Release on deletion; safe to call during GC even when no event loop is running."""
        if vars(self).get("_creator_pid") != os.getpid():
            return  # pragma: forked child
        with contextlib.suppress(Exception):
            try:
                loop = asyncio.get_running_loop()
            except RuntimeError:
                loop = self._context.loop if self._context.loop and not self._context.loop.is_closed() else None
            if loop is None:
                return
            if not loop.is_running():  # pragma: no cover
                loop.run_until_complete(self.release(force=True))
            else:
                loop.create_task(self.release(force=True))


@dataclass
class AsyncFileLockContext(FileLockContext):
    """A dataclass which holds the context for a ``BaseAsyncFileLock`` object."""

    #: Whether run in executor
    run_in_executor: bool = True

    #: The executor
    executor: futures.Executor | None = None

    #: The loop
    loop: asyncio.AbstractEventLoop | None = None


class AsyncThreadLocalFileContext(AsyncFileLockContext, local):
    """A thread local version of the ``FileLockContext`` class."""


class AsyncAcquireReturnProxy:
    """A context-aware object that will release the lock file when exiting."""

    def __init__(self, lock: BaseAsyncFileLock) -> None:  # ruff:ignore[undocumented-public-init]  # trivial release-on-exit proxy
        self.lock = lock

    async def __aenter__(self) -> BaseAsyncFileLock:  # ruff:ignore[undocumented-magic-method]  # returns the wrapped lock
        return self.lock

    async def __aexit__(  # ruff:ignore[undocumented-magic-method]  # releases the wrapped lock on exit
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        await self.lock._release_in_context(exc_value)  # ruff:ignore[private-member-access]  # releases the wrapped lock's context


class AsyncSoftFileLock(SoftFileLock, BaseAsyncFileLock):
    """Simply watches the existence of the lock file."""

    _lifetime_replacements: tuple[str, str] | None = ("AsyncStrictSoftFileLock", "AsyncSoftFileLease")


class AsyncStrictSoftFileLock(StrictSoftFileLock, BaseAsyncFileLock):
    """Run strict owner-claim locking without blocking the event loop."""


class AsyncSoftFileLease(SoftFileLease, BaseAsyncFileLock):
    """Existence lock whose claim expires, so a peer may take it while the previous holder still runs."""


class AsyncUnixFileLock(UnixFileLock, BaseAsyncFileLock):
    """Uses the :func:`fcntl.flock` to hard lock the lock file on unix systems."""


class AsyncWindowsFileLock(WindowsFileLock, BaseAsyncFileLock):
    """Uses the :func:`msvcrt.locking` to hard lock the lock file on windows systems."""


__all__ = [
    "AsyncAcquireReturnProxy",
    "AsyncSoftFileLease",
    "AsyncSoftFileLock",
    "AsyncStrictSoftFileLock",
    "AsyncUnixFileLock",
    "AsyncWindowsFileLock",
    "BaseAsyncFileLock",
]
