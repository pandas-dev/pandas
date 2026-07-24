from __future__ import annotations

import math
import sys
from collections.abc import (
    Coroutine,
    Generator,
)
from contextlib import (
    contextmanager,
)
from enum import Enum, auto
from inspect import iscoroutine
from types import TracebackType
from typing import Any, Generic, final

from ..abc import TaskGroup, TaskStatus
from ._eventloop import get_async_backend, get_cancelled_exc_class
from ._exceptions import TaskCancelled, TaskFailed, TaskNotFinished

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar

if sys.version_info >= (3, 11):
    from typing import Never, TypeVarTuple
else:
    from typing_extensions import Never, TypeVarTuple

T = TypeVar("T")
T_co = TypeVar("T_co", covariant=True)
T_startval = TypeVar("T_startval", covariant=True, default=Never)
PosArgsT = TypeVarTuple("PosArgsT")


class _IgnoredTaskStatus(TaskStatus[object]):
    def started(self, value: object = None) -> None:
        pass


TASK_STATUS_IGNORED = _IgnoredTaskStatus()


class CancelScope:
    """
    Wraps a unit of work that can be made separately cancellable.

    :param deadline: The time (clock value) when this scope is cancelled automatically
    :param shield: ``True`` to shield the cancel scope from external cancellation
    :raises NoEventLoopError: if no supported asynchronous event loop is running in the
        current thread
    """

    __slots__ = ("__weakref__",)

    def __new__(
        cls, *, deadline: float = math.inf, shield: bool = False
    ) -> CancelScope:
        return get_async_backend().create_cancel_scope(shield=shield, deadline=deadline)

    def cancel(self, reason: str | None = None) -> None:
        """
        Cancel this scope immediately.

        :param reason: a message describing the reason for the cancellation

        """
        raise NotImplementedError

    @property
    def deadline(self) -> float:
        """
        The time (clock value) when this scope is cancelled automatically.

        Will be ``float('inf')`` if no timeout has been set.

        """
        raise NotImplementedError

    @deadline.setter
    def deadline(self, value: float) -> None:
        raise NotImplementedError

    @property
    def cancel_called(self) -> bool:
        """``True`` if :meth:`cancel` has been called."""
        raise NotImplementedError

    @property
    def cancelled_caught(self) -> bool:
        """
        ``True`` if this scope suppressed a cancellation exception it itself raised.

        This is typically used to check if any work was interrupted, or to see if the
        scope was cancelled due to its deadline being reached. The value will, however,
        only be ``True`` if the cancellation was triggered by the scope itself (and not
        an outer scope).

        """
        raise NotImplementedError

    @property
    def shield(self) -> bool:
        """
        ``True`` if this scope is shielded from external cancellation.

        While a scope is shielded, it will not receive cancellations from outside.

        """
        raise NotImplementedError

    @shield.setter
    def shield(self, value: bool) -> None:
        raise NotImplementedError

    def __enter__(self) -> CancelScope:
        raise NotImplementedError

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> bool:
        raise NotImplementedError


@contextmanager
def fail_after(
    delay: float | None, shield: bool = False
) -> Generator[CancelScope, None, None]:
    """
    Create a context manager which raises a :class:`TimeoutError` if does not finish in
    time.

    :param delay: maximum allowed time (in seconds) before raising the exception, or
        ``None`` to disable the timeout
    :param shield: ``True`` to shield the cancel scope from external cancellation
    :return: a context manager that yields a cancel scope
    :rtype: :class:`~typing.ContextManager`\\[:class:`~anyio.CancelScope`\\]
    :raises NoEventLoopError: if no supported asynchronous event loop is running in the
        current thread

    """
    current_time = get_async_backend().current_time
    deadline = (current_time() + delay) if delay is not None else math.inf
    with get_async_backend().create_cancel_scope(
        deadline=deadline, shield=shield
    ) as cancel_scope:
        yield cancel_scope

    if cancel_scope.cancelled_caught and current_time() >= cancel_scope.deadline:
        raise TimeoutError


def move_on_after(delay: float | None, shield: bool = False) -> CancelScope:
    """
    Create a cancel scope with a deadline that expires after the given delay.

    :param delay: maximum allowed time (in seconds) before exiting the context block, or
        ``None`` to disable the timeout
    :param shield: ``True`` to shield the cancel scope from external cancellation
    :return: a cancel scope
    :raises NoEventLoopError: if no supported asynchronous event loop is running in the
        current thread

    """
    deadline = (
        (get_async_backend().current_time() + delay) if delay is not None else math.inf
    )
    return get_async_backend().create_cancel_scope(deadline=deadline, shield=shield)


def current_effective_deadline() -> float:
    """
    Return the nearest deadline among all the cancel scopes effective for the current
    task.

    :return: a clock value from the event loop's internal clock (or ``float('inf')`` if
        there is no deadline in effect, or ``float('-inf')`` if the current scope has
        been cancelled)
    :rtype: float
    :raises NoEventLoopError: if no supported asynchronous event loop is running in the
        current thread

    """
    return get_async_backend().current_effective_deadline()


def create_task_group() -> TaskGroup:
    """
    Create a task group.

    :return: a task group
    :raises NoEventLoopError: if no supported asynchronous event loop is running in the
        current thread

    """
    return get_async_backend().create_task_group()


@final
class TaskHandle(Generic[T_co, T_startval]):
    """
    Returned from the task-spawning methods of :class:`TaskGroup`. Can be awaited on to
    get the return value of the task (or the raised exception). If the task was
    terminated by a :exc:`BaseException`, :exc:`TaskFailed` will be raised (or its
    subclass :exc:`TaskCancelled` if the task was cancelled).

    .. versionadded:: 4.14.0
    """

    class Status(Enum):
        """
        The status of a task handle.

        .. attribute:: PENDING

            The task has not finished yet.
        .. attribute:: FINISHED

            The task has finished with a return value.
        .. attribute:: CANCELLING

            The task has been cancelled but has not finished yet.
        .. attribute:: CANCELLED

            The task was cancelled and has finished since.
        .. attribute:: FAILED

            The task raised an exception.
        """

        PENDING = auto()
        FINISHED = auto()
        CANCELLING = auto()
        CANCELLED = auto()
        FAILED = auto()

    __slots__ = (
        "__weakref__",
        "_coro",
        "_name",
        "_cancel_scope",
        "_finished_event",
        "_return_value",
        "_start_value",
        "_exception",
    )

    _return_value: T_co
    _start_value: T_startval

    def __init__(self, coro: Coroutine[Any, Any, T_co], name: object) -> None:
        from ._synchronization import Event

        self._coro = coro
        self._cancel_scope = CancelScope()
        self._finished_event = Event()
        self._exception: BaseException | None = None

        if name is not None:
            self._name = str(name)
        elif iscoroutine(coro):
            self._name = coro.__qualname__
        else:
            self._name = str(coro)  # coroutine-like object (e.g. asend() objects)

    async def _run_coro(self) -> None:
        __tracebackhide__ = True

        with self._cancel_scope:
            try:
                retval = await self._coro
            except BaseException as exc:
                self._exception = exc
                raise
            else:
                self._return_value = retval
            finally:
                self._finished_event.set()
                del self  # Break the reference cycle

    def cancel(self) -> None:
        """
        Set the task to a cancelled state.

        This will interrupt any interruptible asynchronous operation, and will cause
        any further awaits on this task to get immediately cancelled, unless done in
        a shielded cancel scope.

        If the task has already finished, this method has no effect.
        """
        if not self._finished_event.is_set():
            self._cancel_scope.cancel()

    @property
    def coro(self) -> Coroutine[Any, Any, T_co]:
        """
        The coroutine object that was passed to one of the task-spawning methods in
        :class:`TaskGroup`.
        """
        return self._coro

    @property
    def status(self) -> TaskHandle.Status:
        """
        The current status of the task.

        Every task starts in the :attr:`~TaskHandle.Status.PENDING` state.
        If a task is cancelled while in this state, it will transition to the
        :attr:`~TaskHandle.Status.CANCELLING` state. When the task finishes, it will
        transition to one of the three final states (
        :attr:`~TaskHandle.Status.FINISHED`, :attr:`~TaskHandle.Status.FAILED`, or
        :attr:`~TaskHandle.Status.CANCELLING`) depending on the exception the task
        raised, if any. No other status transitions will happen.
        """
        if not self._finished_event.is_set():
            if self._cancel_scope.cancel_called:
                return TaskHandle.Status.CANCELLING
            else:
                return TaskHandle.Status.PENDING
        elif self._exception is not None:
            if isinstance(self._exception, get_cancelled_exc_class()):
                return TaskHandle.Status.CANCELLED
            else:
                return TaskHandle.Status.FAILED
        else:
            return TaskHandle.Status.FINISHED

    @property
    def name(self) -> str:
        """The name of the task."""
        return self._name

    @property
    def exception(self) -> BaseException | None:
        """
        The exception raised by the task, or ``None`` if it finished without raising.

        :raises TaskNotFinished: if the task has not finished yet
        :raises TaskCancelled: if the task was cancelled

        """
        match self.status:
            case TaskHandle.Status.PENDING:
                raise TaskNotFinished("the task has not finished yet")
            case TaskHandle.Status.FINISHED:
                return None
            case TaskHandle.Status.CANCELLING:
                raise TaskCancelled("the task was cancelled")
            case TaskHandle.Status.CANCELLED:
                raise TaskCancelled("the task was cancelled") from self._exception
            case TaskHandle.Status.FAILED:
                return self._exception

    @property
    def return_value(self) -> T_co:
        """
        The return value of the task.

        :raises TaskNotFinished: if the task has not finished yet
        :raises TaskCancelled: if the task was cancelled
        :raises TaskFailed: if the task raised an exception

        """
        match self.status:
            case TaskHandle.Status.PENDING:
                raise TaskNotFinished("the task has not finished yet")
            case TaskHandle.Status.FINISHED:
                return self._return_value
            case TaskHandle.Status.CANCELLING:
                raise TaskCancelled("the task was cancelled")
            case TaskHandle.Status.CANCELLED:
                raise TaskCancelled("the task was cancelled") from self._exception
            case TaskHandle.Status.FAILED:
                raise TaskFailed("the task raised an exception") from self._exception

    @property
    def start_value(self) -> T_startval:
        """
        The value passed to :meth:`task_status.started() <.abc.TaskStatus.started>`,

        :raises RuntimeError: if the task was not started with :meth:`TaskGroup.start()
            <.abc.TaskGroup.start>`
        """
        try:
            return self._start_value
        except AttributeError:
            raise RuntimeError(
                "the task was not started with TaskGroup.start()"
            ) from None

    async def wait(self) -> None:
        """
        Wait for the task to finish.

        This method will return as soon as the task has finished, no matter how it
        happened.
        """
        await self._finished_event.wait()

    def __await__(self) -> Generator[Any, Any, T_co]:
        yield from self._finished_event.wait().__await__()
        return self.return_value

    def __repr__(self) -> str:
        return (
            f"<{self.__class__.__name__} {self.status.name.lower()} "
            f"name={self._name!r} coro={self._coro!r}>"
        )
