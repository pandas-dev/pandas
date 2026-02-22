from __future__ import annotations

import contextlib
import inspect
import logging
import os
import pathlib
import sys
import time
import warnings
from abc import ABCMeta, abstractmethod
from dataclasses import dataclass
from threading import local
from typing import TYPE_CHECKING, Any, cast
from weakref import WeakValueDictionary

from ._error import Timeout

#: Sentinel indicating that no explicit file permission mode was passed.
#: When used, lock files are created with 0o666 (letting umask and default ACLs control the final permissions)
#: and fchmod is skipped so that POSIX default ACL inheritance is preserved.
_UNSET_FILE_MODE: int = -1

if TYPE_CHECKING:
    from collections.abc import Callable
    from types import TracebackType

    from ._read_write import ReadWriteLock

    if sys.version_info >= (3, 11):  # pragma: no cover (py311+)
        from typing import Self
    else:  # pragma: no cover (<py311)
        from typing_extensions import Self


_LOGGER = logging.getLogger("filelock")

# On Windows os.path.realpath calls CreateFileW with share_mode=0, which blocks concurrent DeleteFileW and causes
# livelocks under threaded contention with SoftFileLock. os.path.abspath is purely string-based and avoids this.
_canonical = os.path.abspath if sys.platform == "win32" else os.path.realpath


class _ThreadLocalRegistry(local):
    def __init__(self) -> None:
        super().__init__()
        self.held: dict[str, int] = {}


_registry = _ThreadLocalRegistry()


# This is a helper class which is returned by :meth:`BaseFileLock.acquire` and wraps the lock to make sure __enter__
# is not called twice when entering the with statement. If we would simply return *self*, the lock would be acquired
# again in the *__enter__* method of the BaseFileLock, but not released again automatically. issue #37 (memory leak)
class AcquireReturnProxy:
    """A context-aware object that will release the lock file when exiting."""

    def __init__(self, lock: BaseFileLock | ReadWriteLock) -> None:
        self.lock: BaseFileLock | ReadWriteLock = lock

    def __enter__(self) -> BaseFileLock | ReadWriteLock:
        return self.lock

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.lock.release()


@dataclass
class FileLockContext:
    """A dataclass which holds the context for a ``BaseFileLock`` object."""

    # The context is held in a separate class to allow optional use of thread local storage via the
    # ThreadLocalFileContext class.

    #: The path to the lock file.
    lock_file: str

    #: The default timeout value.
    timeout: float

    #: The mode for the lock files
    mode: int

    #: Whether the lock should be blocking or not
    blocking: bool

    #: The default polling interval value.
    poll_interval: float

    #: The lock lifetime in seconds; ``None`` means the lock never expires.
    lifetime: float | None = None

    #: The file descriptor for the *_lock_file* as it is returned by the os.open() function, not None when lock held
    lock_file_fd: int | None = None

    #: The lock counter is used for implementing the nested locking mechanism.
    lock_counter: int = 0  # When the lock is acquired is increased and the lock is only released, when this value is 0


class ThreadLocalFileContext(FileLockContext, local):
    """A thread local version of the ``FileLockContext`` class."""


class FileLockMeta(ABCMeta):
    _instances: WeakValueDictionary[str, BaseFileLock]

    def __call__(  # noqa: PLR0913
        cls,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        mode: int = _UNSET_FILE_MODE,
        thread_local: bool = True,  # noqa: FBT001, FBT002
        *,
        blocking: bool = True,
        is_singleton: bool = False,
        poll_interval: float = 0.05,
        lifetime: float | None = None,
        **kwargs: Any,  # capture remaining kwargs for subclasses  # noqa: ANN401
    ) -> BaseFileLock:
        if is_singleton:
            instance = cls._instances.get(str(lock_file))
            if instance:
                params_to_check = {
                    "thread_local": (thread_local, instance.is_thread_local()),
                    "timeout": (timeout, instance.timeout),
                    "mode": (mode, instance._context.mode),  # noqa: SLF001
                    "blocking": (blocking, instance.blocking),
                    "poll_interval": (poll_interval, instance.poll_interval),
                    "lifetime": (lifetime, instance.lifetime),
                }

                non_matching_params = {
                    name: (passed_param, set_param)
                    for name, (passed_param, set_param) in params_to_check.items()
                    if passed_param != set_param
                }
                if not non_matching_params:
                    return cast("BaseFileLock", instance)

                # parameters do not match; raise error
                msg = "Singleton lock instances cannot be initialized with differing arguments"
                msg += "\nNon-matching arguments: "
                for param_name, (passed_param, set_param) in non_matching_params.items():
                    msg += f"\n\t{param_name} (existing lock has {set_param} but {passed_param} was passed)"
                raise ValueError(msg)

        # Workaround to make `__init__`'s params optional in subclasses
        # E.g. virtualenv changes the signature of the `__init__` method in the `BaseFileLock` class descendant
        # (https://github.com/tox-dev/filelock/pull/340)

        all_params = {
            "timeout": timeout,
            "mode": mode,
            "thread_local": thread_local,
            "blocking": blocking,
            "is_singleton": is_singleton,
            "poll_interval": poll_interval,
            "lifetime": lifetime,
            **kwargs,
        }

        present_params = inspect.signature(cls.__init__).parameters
        init_params = {key: value for key, value in all_params.items() if key in present_params}

        instance = super().__call__(lock_file, **init_params)

        if is_singleton:
            cls._instances[str(lock_file)] = instance

        return cast("BaseFileLock", instance)


class BaseFileLock(contextlib.ContextDecorator, metaclass=FileLockMeta):
    """
    Abstract base class for a file lock object.

    Provides a reentrant, cross-process exclusive lock backed by OS-level primitives. Subclasses implement the actual
    locking mechanism (:class:`UnixFileLock <filelock.UnixFileLock>`, :class:`WindowsFileLock
    <filelock.WindowsFileLock>`, :class:`SoftFileLock <filelock.SoftFileLock>`).

    """

    _instances: WeakValueDictionary[str, BaseFileLock]

    def __init_subclass__(cls, **kwargs: dict[str, Any]) -> None:
        """Setup unique state for lock subclasses."""
        super().__init_subclass__(**kwargs)
        cls._instances = WeakValueDictionary()

    def __init__(  # noqa: PLR0913
        self,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        mode: int = _UNSET_FILE_MODE,
        thread_local: bool = True,  # noqa: FBT001, FBT002
        *,
        blocking: bool = True,
        is_singleton: bool = False,
        poll_interval: float = 0.05,
        lifetime: float | None = None,
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
            ``False`` then the lock will be reentrant across threads.
        :param blocking: whether the lock should be blocking or not
        :param is_singleton: If this is set to ``True`` then only one instance of this class will be created per lock
            file. This is useful if you want to use the lock object for reentrant locking without needing to pass the
            same object around.
        :param poll_interval: default interval for polling the lock file, in seconds. It will be used as fallback value
            in the acquire method, if no poll_interval value (``None``) is given.
        :param lifetime: maximum time in seconds a lock can be held before it is considered expired. When set, a waiting
            process will break a lock whose file modification time is older than ``lifetime`` seconds. ``None`` (the
            default) means locks never expire.

        """
        self._is_thread_local = thread_local
        self._is_singleton = is_singleton

        # Create the context. Note that external code should not work with the context directly and should instead use
        # properties of this class.
        kwargs: dict[str, Any] = {
            "lock_file": os.fspath(lock_file),
            "timeout": timeout,
            "mode": mode,
            "blocking": blocking,
            "poll_interval": poll_interval,
            "lifetime": lifetime,
        }
        self._context: FileLockContext = (ThreadLocalFileContext if thread_local else FileLockContext)(**kwargs)

    def is_thread_local(self) -> bool:
        """:returns: a flag indicating if this lock is thread local or not"""
        return self._is_thread_local

    @property
    def is_singleton(self) -> bool:
        """
        :returns: a flag indicating if this lock is singleton or not

        .. versionadded:: 3.13.0

        """
        return self._is_singleton

    @property
    def lock_file(self) -> str:
        """:returns: path to the lock file"""
        return self._context.lock_file

    @property
    def timeout(self) -> float:
        """
        :returns: the default timeout value, in seconds

        .. versionadded:: 2.0.0

        """
        return self._context.timeout

    @timeout.setter
    def timeout(self, value: float | str) -> None:
        """
        Change the default timeout value.

        :param value: the new value, in seconds

        """
        self._context.timeout = float(value)

    @property
    def blocking(self) -> bool:
        """
        :returns: whether the locking is blocking or not

        .. versionadded:: 3.14.0

        """
        return self._context.blocking

    @blocking.setter
    def blocking(self, value: bool) -> None:
        """
        Change the default blocking value.

        :param value: the new value as bool

        """
        self._context.blocking = value

    @property
    def poll_interval(self) -> float:
        """
        :returns: the default polling interval, in seconds

        .. versionadded:: 3.24.0

        """
        return self._context.poll_interval

    @poll_interval.setter
    def poll_interval(self, value: float) -> None:
        """
        Change the default polling interval.

        :param value: the new value, in seconds

        """
        self._context.poll_interval = value

    @property
    def lifetime(self) -> float | None:
        """
        :returns: the lock lifetime in seconds, or ``None`` if the lock never expires

        .. versionadded:: 3.24.0

        """
        return self._context.lifetime

    @lifetime.setter
    def lifetime(self, value: float | None) -> None:
        """
        Change the lock lifetime.

        :param value: the new value in seconds, or ``None`` to disable expiration

        """
        self._context.lifetime = value

    @property
    def mode(self) -> int:
        """:returns: the file permissions for the lockfile"""
        return 0o644 if self._context.mode == _UNSET_FILE_MODE else self._context.mode

    @property
    def has_explicit_mode(self) -> bool:
        """:returns: whether the file permissions were explicitly set"""
        return self._context.mode != _UNSET_FILE_MODE

    def _open_mode(self) -> int:
        """:returns: the mode for os.open() — 0o666 when unset (let umask/ACLs decide), else the explicit mode"""
        return 0o666 if self._context.mode == _UNSET_FILE_MODE else self._context.mode

    def _try_break_expired_lock(self) -> None:
        """Remove the lock file if its modification time exceeds the configured :attr:`lifetime`."""
        if (lifetime := self._context.lifetime) is None:
            return
        with contextlib.suppress(OSError):
            if time.time() - pathlib.Path(self.lock_file).stat().st_mtime < lifetime:
                return
            break_path = f"{self.lock_file}.break.{os.getpid()}"
            pathlib.Path(self.lock_file).rename(break_path)
            pathlib.Path(break_path).unlink()

    @abstractmethod
    def _acquire(self) -> None:
        """If the file lock could be acquired, self._context.lock_file_fd holds the file descriptor of the lock file."""
        raise NotImplementedError

    @abstractmethod
    def _release(self) -> None:
        """Releases the lock and sets self._context.lock_file_fd to None."""
        raise NotImplementedError

    @property
    def is_locked(self) -> bool:
        """
        :returns: A boolean indicating if the lock file is holding the lock currently.

        .. versionchanged:: 2.0.0

            This was previously a method and is now a property.

        """
        return self._context.lock_file_fd is not None

    @property
    def lock_counter(self) -> int:
        """:returns: The number of times this lock has been acquired (but not yet released)."""
        return self._context.lock_counter

    @staticmethod
    def _check_give_up(  # noqa: PLR0913
        lock_id: int,
        lock_filename: str,
        *,
        blocking: bool,
        cancel_check: Callable[[], bool] | None,
        timeout: float,
        start_time: float,
    ) -> bool:
        if blocking is False:
            _LOGGER.debug("Failed to immediately acquire lock %s on %s", lock_id, lock_filename)
            return True
        if cancel_check is not None and cancel_check():
            _LOGGER.debug("Cancellation requested for lock %s on %s", lock_id, lock_filename)
            return True
        if 0 <= timeout < time.perf_counter() - start_time:
            _LOGGER.debug("Timeout on acquiring lock %s on %s", lock_id, lock_filename)
            return True
        return False

    def acquire(  # noqa: C901
        self,
        timeout: float | None = None,
        poll_interval: float | None = None,
        *,
        poll_intervall: float | None = None,
        blocking: bool | None = None,
        cancel_check: Callable[[], bool] | None = None,
    ) -> AcquireReturnProxy:
        """
        Try to acquire the file lock.

        :param timeout: maximum wait time for acquiring the lock, ``None`` means use the default :attr:`~timeout` is and
            if ``timeout < 0``, there is no timeout and this method will block until the lock could be acquired
        :param poll_interval: interval of trying to acquire the lock file, ``None`` means use the default
            :attr:`~poll_interval`
        :param poll_intervall: deprecated, kept for backwards compatibility, use ``poll_interval`` instead
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

        .. versionchanged:: 2.0.0

            This method returns now a *proxy* object instead of *self*, so that it can be used in a with statement
            without side effects.

        """
        # Use the default timeout, if no timeout is provided.
        if timeout is None:
            timeout = self._context.timeout

        if blocking is None:
            blocking = self._context.blocking

        if poll_intervall is not None:
            msg = "use poll_interval instead of poll_intervall"
            warnings.warn(msg, DeprecationWarning, stacklevel=2)
            poll_interval = poll_intervall

        poll_interval = poll_interval if poll_interval is not None else self._context.poll_interval

        # Increment the number right at the beginning. We can still undo it, if something fails.
        self._context.lock_counter += 1

        lock_id = id(self)
        lock_filename = self.lock_file
        canonical = _canonical(lock_filename)

        would_block = self._context.lock_counter == 1 and not self.is_locked and timeout < 0 and blocking
        if would_block and (existing := _registry.held.get(canonical)) is not None and existing != lock_id:
            self._context.lock_counter -= 1
            msg = (
                f"Deadlock: lock '{lock_filename}' is already held by a different "
                f"FileLock instance in this thread. Use is_singleton=True to "
                f"enable reentrant locking across instances."
            )
            raise RuntimeError(msg)

        start_time = time.perf_counter()
        try:
            while True:
                if not self.is_locked:
                    self._try_break_expired_lock()
                    _LOGGER.debug("Attempting to acquire lock %s on %s", lock_id, lock_filename)
                    self._acquire()
                if self.is_locked:
                    _LOGGER.debug("Lock %s acquired on %s", lock_id, lock_filename)
                    break
                if self._check_give_up(
                    lock_id,
                    lock_filename,
                    blocking=blocking,
                    cancel_check=cancel_check,
                    timeout=timeout,
                    start_time=start_time,
                ):
                    raise Timeout(lock_filename)  # noqa: TRY301
                msg = "Lock %s not acquired on %s, waiting %s seconds ..."
                _LOGGER.debug(msg, lock_id, lock_filename, poll_interval)
                time.sleep(poll_interval)
        except BaseException:
            self._context.lock_counter = max(0, self._context.lock_counter - 1)
            if self._context.lock_counter == 0:
                _registry.held.pop(canonical, None)
            raise
        if self._context.lock_counter == 1:
            _registry.held[canonical] = lock_id
        return AcquireReturnProxy(lock=self)

    def release(self, force: bool = False) -> None:  # noqa: FBT001, FBT002
        """
        Release the file lock. The lock is only completely released when the lock counter reaches 0. The lock file
        itself is not automatically deleted.

        :param force: If true, the lock counter is ignored and the lock is released in every case.

        """
        if self.is_locked:
            self._context.lock_counter -= 1

            if self._context.lock_counter == 0 or force:
                lock_id, lock_filename = id(self), self.lock_file

                _LOGGER.debug("Attempting to release lock %s on %s", lock_id, lock_filename)
                self._release()
                self._context.lock_counter = 0
                _registry.held.pop(_canonical(lock_filename), None)
                _LOGGER.debug("Lock %s released on %s", lock_id, lock_filename)

    def __enter__(self) -> Self:
        """
        Acquire the lock.

        :returns: the lock object

        """
        self.acquire()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """
        Release the lock.

        :param exc_type: the exception type if raised
        :param exc_value: the exception value if raised
        :param traceback: the exception traceback if raised

        """
        self.release()

    def __del__(self) -> None:
        """Called when the lock object is deleted."""
        self.release(force=True)


__all__ = [
    "_UNSET_FILE_MODE",
    "AcquireReturnProxy",
    "BaseFileLock",
]
