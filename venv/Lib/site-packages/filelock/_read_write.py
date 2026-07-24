from __future__ import annotations

import logging
import os
import pathlib
import sqlite3
import sys
import threading
import time
from contextlib import contextmanager, suppress
from typing import TYPE_CHECKING, ClassVar, Final, Literal, TypeAlias, cast
from weakref import WeakValueDictionary

from ._api import (
    AcquireReturnProxy,
    _ensure_current_process,
    _fork_transition,
    _raise_chained_errors,
    _register_fork_class,
    _register_fork_object,
)
from ._error import Timeout

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self

_LOGGER: Final[logging.Logger] = logging.getLogger("filelock")
_GETPID: Final[Callable[[], int]] = os.getpid
_IS_PYPY: Final[bool] = sys.implementation.name == "pypy"
_NEEDS_CONNECTION_ESCROW: Final[bool] = (
    hasattr(os, "register_at_fork") and sys.implementation.name == "cpython" and sys.version_info < (3, 12)
)
_ConnectionParameter: TypeAlias = (
    str | bytes | os.PathLike[str] | os.PathLike[bytes] | float | int | type[sqlite3.Connection] | None
)
_DatabaseIdentity: TypeAlias = tuple[int, int]

# sqlite3_busy_timeout() accepts a C int, max 2_147_483_647 on 32-bit. Use a lower value to be safe (~23 days).
_MAX_SQLITE_TIMEOUT_MS: Final[int] = 2_000_000_000 - 1
_UNSAFE_FORK_EXIT_STATUS: Final[int] = 70


class _SQLiteTransitionContext(threading.local):
    depth: int = 0


_SQLITE_TRANSITION_CONTEXT: Final = _SQLiteTransitionContext()


class _ConnectionEscrow:
    def __init__(self) -> None:
        self._lock = threading.RLock()
        self._functions: tuple[Callable[[sqlite3.Connection], None], Callable[[sqlite3.Connection], None]] | None = None

    def functions(
        self,
    ) -> tuple[Callable[[sqlite3.Connection], None], Callable[[sqlite3.Connection], None]] | None:
        if not _NEEDS_CONNECTION_ESCROW:
            return None  # pragma: >=3.12 cover
        with self._lock:  # pragma: <3.12 cover  # pragma: needs fork
            if self._functions is None:
                import ctypes  # ruff:ignore[import-outside-top-level]  # keep optional ctypes and its audited dlsym out of ordinary imports

                function_type = ctypes.PYFUNCTYPE(None, ctypes.py_object)
                increment_address = ctypes.cast(ctypes.pythonapi.Py_IncRef, ctypes.c_void_p).value
                decrement_address = ctypes.cast(ctypes.pythonapi.Py_DecRef, ctypes.c_void_p).value
                if increment_address is None or decrement_address is None:  # pragma: no cover - resolved CPython API
                    msg = "CPython reference functions have no address"
                    raise RuntimeError(msg)
                self._functions = (
                    cast("Callable[[sqlite3.Connection], None]", function_type(increment_address)),
                    cast("Callable[[sqlite3.Connection], None]", function_type(decrement_address)),
                )
            return self._functions

    def _reset_after_fork_in_child(self) -> None:  # pragma: forked child
        self._lock = threading.RLock()


_CONNECTION_ESCROW: Final = _ConnectionEscrow()


class _ForkedDatabaseRegistry:
    def __init__(self) -> None:
        self._lock = threading.RLock()
        self._paths: set[pathlib.Path] = set()
        self._identities: set[_DatabaseIdentity] = set()
        self._sqlite_used = False
        self._all_paths_poisoned = False

    def raise_if_poisoned(self, path: pathlib.Path) -> None:
        identity = self.identity(path)
        with self._lock:
            all_paths_poisoned = self._all_paths_poisoned
            poisoned = (
                all_paths_poisoned or path in self._paths or (identity is not None and identity in self._identities)
            )
        if poisoned:  # pragma: needs fork
            msg = (
                "ReadWriteLock is unavailable in a PyPy fork child; exec or exit before using it"
                if all_paths_poisoned
                else f"SQLite database {path!s} was active across fork(); exec or exit before using it in the child"
            )
            raise RuntimeError(msg)

    def poison_after_fork(self, path: pathlib.Path, identity: _DatabaseIdentity | None) -> None:
        self._paths.add(path)
        if identity is not None:
            self._identities.add(identity)

    def note_sqlite_use(self) -> None:
        if _IS_PYPY:
            with self._lock:
                self._sqlite_used = True

    def _reset_after_fork_in_child(self) -> None:  # pragma: forked child
        self._lock = threading.RLock()
        self._all_paths_poisoned = self._all_paths_poisoned or (_IS_PYPY and self._sqlite_used)
        self._sqlite_used = False

    @staticmethod
    def identity(path: pathlib.Path) -> _DatabaseIdentity | None:
        try:
            stat_result = path.stat()
        except OSError:
            return None
        return stat_result.st_dev, stat_result.st_ino


_FORKED_DATABASES: Final = _ForkedDatabaseRegistry()


class _ForkSafeConnection(sqlite3.Connection):
    _creator_pid: int
    _decrement_escrow: Callable[[sqlite3.Connection], None] | None

    def __new__(
        cls,
        *_args: _ConnectionParameter,
        **_kwargs: _ConnectionParameter,
    ) -> Self:
        connection = super().__new__(cls)
        connection._creator_pid = _GETPID()
        connection._decrement_escrow = None
        return connection

    def close(self) -> None:
        with _sqlite_transition():
            if _GETPID() != self._creator_pid:  # pragma: needs fork
                return
            with _fork_transition():
                sqlite3.Connection.close(self)
                if (decrement := self._decrement_escrow) is not None:  # pragma: <3.12 cover  # pragma: needs fork
                    self._decrement_escrow = None
                    decrement(self)

    def acquire_escrow(  # pragma: <3.12 cover  # pragma: needs fork
        self,
        functions: tuple[Callable[[sqlite3.Connection], None], Callable[[sqlite3.Connection], None]] | None,
    ) -> None:
        # The caller only reaches here holding the escrow functions; it skips the call entirely without them.
        if functions is not None:  # pragma: no branch
            increment, decrement = functions
            increment(self)
            self._decrement_escrow = decrement

    def __del__(self) -> None:
        with suppress(sqlite3.Error, RuntimeError):
            self.close()


class _ReadWriteLockMeta(type):
    """
    Resolve singleton instances for ``is_singleton=True`` construction.

    This logic lives here rather than in ReadWriteLock.get_lock so ``ReadWriteLock(path)`` returns cached instances
    without a 2-arg ``super()`` call that type checkers cannot verify.

    """

    _instances: WeakValueDictionary[pathlib.Path, ReadWriteLock]
    _instances_lock: threading.RLock
    _instances_pid: int
    _instances_under_construction: set[pathlib.Path]

    def __call__(
        cls,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        *,
        blocking: bool = True,
        is_singleton: bool = True,
    ) -> ReadWriteLock:
        _ensure_current_process()
        if cls._instances_pid != _GETPID():
            cls._reset_class_after_fork()
        construction_pid = _GETPID()
        if not is_singleton:
            instance = super().__call__(lock_file, timeout, blocking=blocking, is_singleton=is_singleton)
            if _GETPID() != construction_pid:  # pragma: forked child
                msg = "ReadWriteLock construction cannot continue after fork"
                raise RuntimeError(msg)
            return instance

        normalized = pathlib.Path(lock_file).resolve()
        with cls._instances_lock:
            if normalized not in cls._instances:
                if normalized in cls._instances_under_construction:  # pragma: no cover - exercised in an audit callback
                    msg = f"Singleton lock construction is already active for {lock_file!s}"
                    raise RuntimeError(msg)
                construction_registry = cls._instances_under_construction
                construction_registry.add(normalized)
                try:
                    instance = super().__call__(lock_file, timeout, blocking=blocking, is_singleton=is_singleton)
                finally:
                    if _GETPID() == construction_pid:
                        construction_registry.discard(normalized)
                if _GETPID() != construction_pid:
                    msg = "ReadWriteLock construction cannot continue after fork"
                    raise RuntimeError(msg)
                cls._instances[normalized] = instance
            else:
                instance = cls._instances[normalized]

            if instance.timeout != timeout or instance.blocking != blocking:
                msg = (
                    f"Singleton lock created with timeout={instance.timeout}, blocking={instance.blocking},"
                    f" cannot be changed to timeout={timeout}, blocking={blocking}"
                )
                raise ValueError(msg)
            return instance

    def _reset_class_after_fork(cls) -> None:  # pragma: forked child
        cls._instances = WeakValueDictionary()
        cls._instances_lock = threading.RLock()
        cls._instances_pid = _GETPID()
        cls._instances_under_construction = set()


class ReadWriteLock(metaclass=_ReadWriteLockMeta):
    """
    Cross-process read-write lock backed by SQLite.

    Allows concurrent shared readers or a single exclusive writer. The lock is reentrant within the same mode (multiple
    ``acquire_read`` calls nest, as do multiple ``acquire_write`` calls from the same thread), but upgrading from read
    to write or downgrading from write to read raises :class:`RuntimeError`. Write locks are pinned to the thread that
    acquired them.

    By default, ``is_singleton=True``: calling ``ReadWriteLock(path)`` with the same resolved path returns the same
    instance. The path is handed to :func:`sqlite3.connect` as given, so a ``.db`` extension is a convention rather
    than a requirement; the filesystem must be one the active SQLite VFS supports.

    :param lock_file: path to the SQLite database file used as the lock
    :param timeout: maximum wait time in seconds; ``-1`` means block indefinitely
    :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately when the lock is unavailable
    :param is_singleton: if ``True``, reuse existing instances for the same resolved path

    .. versionadded:: 3.21.0

    """

    _instances: WeakValueDictionary[pathlib.Path, ReadWriteLock] = WeakValueDictionary()
    _instances_lock = threading.RLock()
    _instances_pid = _GETPID()
    _instances_under_construction: ClassVar[set[pathlib.Path]] = set()

    def __init_subclass__(cls) -> None:
        super().__init_subclass__()
        cls._instances = WeakValueDictionary()
        cls._instances_lock = threading.RLock()
        cls._instances_pid = _GETPID()
        cls._instances_under_construction = set()
        _register_fork_class(cls)

    @classmethod
    def get_lock(
        cls, lock_file: str | os.PathLike[str], timeout: float = -1, *, blocking: bool = True
    ) -> ReadWriteLock:
        """
        Return the singleton :class:`ReadWriteLock` for *lock_file*.

        :param lock_file: path to the SQLite database file used as the lock
        :param timeout: maximum wait time in seconds; ``-1`` means block indefinitely
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately when the lock is unavailable

        :returns: the singleton lock instance

        :raises ValueError: if an instance already exists for this path with different *timeout* or *blocking* values

        """
        return cls(lock_file, timeout, blocking=blocking)

    def __init__(
        self,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        *,
        blocking: bool = True,
        is_singleton: bool = True,  # ruff:ignore[unused-method-argument]  # consumed by _ReadWriteLockMeta.__call__
    ) -> None:
        self.lock_file = os.fspath(lock_file)
        self._canonical_path = pathlib.Path(lock_file).resolve()
        _FORKED_DATABASES.raise_if_poisoned(self._canonical_path)
        self.timeout = timeout
        self.blocking = blocking
        self._transaction_lock = threading.Lock()  # serializes the (possibly blocking) SQLite transaction work
        self._internal_lock = threading.Lock()  # protects _lock_level / _current_mode updates and rollback
        self._lock_level = 0
        self._current_mode: Literal["read", "write"] | None = None
        self._write_thread_id: int | None = None
        self._acquisition_thread_ids: set[int] = set()
        self._con: _ForkSafeConnection | None = None
        self._connection_transaction_released = True
        self._connection_identity: _DatabaseIdentity | None = None
        self._closed = False
        self._creator_pid = _GETPID()
        self._fork_invalidated = False
        _register_fork_object(self)
        with _fork_transition(), _sqlite_transition():
            validation_connection = self._open_connection(sqlite_timeout=5.0)
            validation_connection.close()

    def acquire_read(self, timeout: float = -1, *, blocking: bool = True) -> AcquireReturnProxy:
        """
        Acquire a shared read lock.

        If this instance already holds a read lock, the lock level is incremented (reentrant). Attempting to acquire a
        read lock while holding a write lock raises :class:`RuntimeError` (downgrade not allowed).

        :param timeout: maximum wait time in seconds; ``-1`` means block indefinitely
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately when the lock is unavailable

        :returns: a proxy that can be used as a context manager to release the lock

        :raises RuntimeError: if a write lock is already held on this instance
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        return self._acquire("read", timeout, blocking=blocking)

    def acquire_write(self, timeout: float = -1, *, blocking: bool = True) -> AcquireReturnProxy:
        """
        Acquire an exclusive write lock.

        If this instance already holds a write lock from the same thread, the lock level is incremented (reentrant).
        Attempting to acquire a write lock while holding a read lock raises :class:`RuntimeError` (upgrade not allowed).
        Write locks are pinned to the acquiring thread: a different thread trying to re-enter also raises
        :class:`RuntimeError`.

        :param timeout: maximum wait time in seconds; ``-1`` means block indefinitely
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately when the lock is unavailable

        :returns: a proxy that can be used as a context manager to release the lock

        :raises RuntimeError: if a read lock is already held, or a write lock is held by a different thread
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        return self._acquire("write", timeout, blocking=blocking)

    @contextmanager
    def read_lock(self, timeout: float | None = None, *, blocking: bool | None = None) -> Generator[None]:
        """
        Context manager that acquires and releases a shared read lock.

        Falls back to instance defaults for *timeout* and *blocking* when ``None``.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately; ``None`` uses the instance default

        """
        if timeout is None:
            timeout = self.timeout
        if blocking is None:
            blocking = self.blocking
        self.acquire_read(timeout, blocking=blocking)
        try:
            yield
        finally:
            self.release()

    @contextmanager
    def write_lock(self, timeout: float | None = None, *, blocking: bool | None = None) -> Generator[None]:
        """
        Context manager that acquires and releases an exclusive write lock.

        Falls back to instance defaults for *timeout* and *blocking* when ``None``.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately; ``None`` uses the instance default

        """
        if timeout is None:
            timeout = self.timeout
        if blocking is None:
            blocking = self.blocking
        self.acquire_write(timeout, blocking=blocking)
        try:
            yield
        finally:
            self.release()

    def release(self, *, force: bool = False) -> None:
        """
        Release one level of the current lock.

        When the lock level reaches zero the underlying SQLite transaction is rolled back, releasing the database lock.

        :param force: if ``True``, release the lock completely regardless of the current lock level

        :raises RuntimeError: if no lock is currently held and *force* is ``False``

        """
        with _fork_transition():
            _ensure_current_process()
            if self._inherited:  # pragma: needs fork
                return
            self._raise_if_acquiring("release")
            self._release(force=force, close=False)

    def close(self) -> None:
        """
        Release the lock (if held) and close the underlying SQLite connection.

        After calling this method, the lock instance is no longer usable.

        """
        with _fork_transition():
            _ensure_current_process()
            if self._inherited:  # pragma: needs fork
                return
            self._raise_if_acquiring("close")
            self._release(force=True, close=True)

    def _release(self, *, force: bool, close: bool) -> None:
        with self._transaction_lock, self._internal_lock:
            if self._lock_level == 0:
                if force and self._con is None:
                    if close:
                        self._closed = True
                    return
                if not force:
                    msg = f"Cannot release a lock on {self.lock_file} (lock id: {id(self)}) that is not held"
                    raise RuntimeError(msg)
            if not force and self._lock_level > 1:
                self._lock_level -= 1
                return
            try:
                self._finish_connection()
            except sqlite3.Error:
                if self._connection_transaction_released:
                    self._clear_lock_state()
                raise
            self._clear_lock_state()
            if close:
                self._closed = True

    def _clear_lock_state(self) -> None:
        self._lock_level = 0
        self._current_mode = None
        self._write_thread_id = None

    def __del__(self) -> None:
        if _GETPID() == getattr(self, "_creator_pid", None) and (connection := getattr(self, "_con", None)) is not None:
            with suppress(sqlite3.Error, RuntimeError):
                connection.close()

    def _reset_after_fork_in_child(self) -> None:  # pragma: forked child
        if self._con is not None:
            _FORKED_DATABASES.poison_after_fork(self._canonical_path, self._connection_identity)
        self._con = None
        self._connection_transaction_released = True
        self._connection_identity = None
        self._transaction_lock = threading.Lock()
        self._internal_lock = threading.Lock()
        self._clear_lock_state()
        self._acquisition_thread_ids = set()
        self._fork_invalidated = True

    @property
    def _inherited(self) -> bool:
        return self._fork_invalidated or _GETPID() != self._creator_pid

    def _raise_if_unusable(self) -> None:
        _ensure_current_process()
        if self._inherited:  # pragma: needs fork
            msg = f"ReadWriteLock on {self.lock_file} was invalidated by fork(); construct a new instance"
            raise RuntimeError(msg)
        if self._closed:
            msg = "Cannot operate on a closed database."
            raise sqlite3.ProgrammingError(msg)

    def _acquire(self, mode: Literal["read", "write"], timeout: float, *, blocking: bool) -> AcquireReturnProxy:
        with _fork_transition():
            self._raise_if_unusable()
            operation_pid = _GETPID()
            thread_id = threading.get_ident()
            with self._internal_lock:
                if self._lock_level > 0:
                    return self._validate_reentrant(mode)
                if thread_id in self._acquisition_thread_ids:  # pragma: no cover - exercised in an audit callback
                    msg = f"Cannot acquire ReadWriteLock on {self.lock_file} while acquisition is active in this thread"
                    raise RuntimeError(msg)
                self._acquisition_thread_ids.add(thread_id)
            try:
                start_time = time.perf_counter()
                self._acquire_transaction_lock(blocking=blocking, timeout=timeout)
                try:
                    self._raise_if_unusable()
                    return self._do_acquire_inner(
                        mode,
                        timeout,
                        blocking=blocking,
                        operation_pid=operation_pid,
                        start_time=start_time,
                    )
                finally:
                    self._transaction_lock.release()
            finally:
                with self._internal_lock:
                    self._acquisition_thread_ids.discard(thread_id)

    def _do_acquire_inner(
        self,
        mode: Literal["read", "write"],
        timeout: float,
        *,
        blocking: bool,
        operation_pid: int,
        start_time: float,
    ) -> AcquireReturnProxy:
        # Double-check: another thread may have acquired the lock while we waited on _transaction_lock.
        with self._internal_lock:
            if self._lock_level > 0:
                return self._validate_reentrant(mode)
        if self._con is not None:
            self._finish_connection()
        try:
            self._open_for_acquisition(
                timeout,
                blocking=blocking,
                operation_pid=operation_pid,
                start_time=start_time,
            )
            self._configure_and_begin(
                mode,
                timeout,
                blocking=blocking,
                operation=(operation_pid, start_time),
            )
            self._raise_if_process_changed(operation_pid)
        except BaseException as error:
            acquisition_error: BaseException
            if isinstance(error, sqlite3.OperationalError) and "database is locked" in str(error):
                acquisition_error = Timeout(self.lock_file)
            else:
                acquisition_error = error
            try:
                self._finish_connection()
            except sqlite3.Error as cleanup_error:
                _raise_chained_errors(acquisition_error, cleanup_error)
            if acquisition_error is not error:
                raise acquisition_error from None
            raise
        with self._internal_lock:
            self._raise_if_process_changed(operation_pid)
            self._current_mode = mode
            self._lock_level = 1
            if mode == "write":
                self._write_thread_id = threading.get_ident()
        return AcquireReturnProxy(lock=self)

    def _open_for_acquisition(self, timeout: float, *, blocking: bool, operation_pid: int, start_time: float) -> None:
        with _sqlite_transition():
            sqlite_timeout = (
                timeout_for_sqlite(
                    timeout,
                    blocking=blocking,
                    already_waited=time.perf_counter() - start_time,
                )
                / 1000
            )
            connection = self._open_connection(sqlite_timeout=sqlite_timeout)
            self._con, self._connection_transaction_released, self._connection_identity = (
                connection,
                False,
                _FORKED_DATABASES.identity(self._canonical_path),
            )
            self._raise_if_process_changed(operation_pid)

    def _configure_and_begin(
        self,
        mode: Literal["read", "write"],
        timeout: float,
        *,
        blocking: bool,
        operation: tuple[int, float],
    ) -> None:
        with _sqlite_transition():
            operation_pid, start_time = operation
            connection = cast("_ForkSafeConnection", self._con)
            waited = time.perf_counter() - start_time
            timeout_ms = timeout_for_sqlite(timeout, blocking=blocking, already_waited=waited)
            self._raise_if_process_changed(operation_pid)
            connection.executescript(f"PRAGMA busy_timeout={timeout_ms}; PRAGMA journal_mode=MEMORY;").close()
            # Use legacy journal mode (not WAL) because WAL does not block readers while a concurrent EXCLUSIVE
            # write transaction is active, which makes read-write locking impossible without modifying table data.
            # MEMORY is safe here since no writes happen, so a crash cannot corrupt the DB.
            # See https://sqlite.org/lang_transaction.html#deferred_immediate_and_exclusive_transactions
            #
            # Recompute the remaining timeout after the blocking journal_mode pragma.
            waited = time.perf_counter() - start_time
            recomputed = timeout_for_sqlite(timeout, blocking=blocking, already_waited=waited)
            self._raise_if_process_changed(operation_pid)
            statements = f"PRAGMA busy_timeout={recomputed}; " if recomputed != timeout_ms else ""
            statements += "BEGIN EXCLUSIVE TRANSACTION;" if mode == "write" else "BEGIN TRANSACTION;"
            if mode == "read":
                # SQLite takes the SHARED lock only when a statement reads; BEGIN alone stays deferred.
                # https://www.sqlite.org/lockingv3.html#transaction_control
                statements += " SELECT name FROM sqlite_schema LIMIT 1;"
            connection.executescript(statements).close()

    def _open_connection(self, *, sqlite_timeout: float) -> _ForkSafeConnection:
        with _sqlite_transition():
            creator_pid = _GETPID()
            functions = _CONNECTION_ESCROW.functions()
            if _GETPID() != creator_pid:  # pragma: forked child
                msg = "SQLite connection construction cannot continue after fork"
                raise RuntimeError(msg)
            connection = _connect(
                os.fspath(self._canonical_path),
                factory=_ForkSafeConnection,
                timeout=sqlite_timeout,
            )
            if functions is not None:  # pragma: <3.12 cover  # pragma: needs fork
                connection.acquire_escrow(functions)
            if _GETPID() != creator_pid:  # pragma: forked child
                _FORKED_DATABASES.poison_after_fork(
                    self._canonical_path,
                    _FORKED_DATABASES.identity(self._canonical_path),
                )
                msg = "SQLite connection construction cannot continue after fork"
                raise RuntimeError(msg)
            return connection

    def _finish_connection(self) -> None:
        with _sqlite_transition():
            if (connection := self._con) is None:
                return
            rollback_error: sqlite3.Error | None = None
            if not self._connection_transaction_released:
                if connection.in_transaction:
                    try:
                        connection.rollback()
                    except sqlite3.Error as error:
                        if connection.in_transaction:
                            raise
                        self._connection_transaction_released = True
                        rollback_error = error
                    else:
                        self._connection_transaction_released = True
                else:
                    self._connection_transaction_released = True
            try:
                connection.close()
            except sqlite3.Error as close_error:
                if rollback_error is not None:
                    _raise_chained_errors(rollback_error, close_error)
                raise
            self._con = None
            self._connection_transaction_released = True
            self._connection_identity = None
            if rollback_error is not None:
                raise rollback_error

    def _validate_reentrant(self, mode: Literal["read", "write"]) -> AcquireReturnProxy:
        if self._current_mode != mode:
            opposite = "write" if mode == "read" else "read"
            direction = "downgrade" if mode == "read" else "upgrade"
            msg = (
                f"Cannot acquire {mode} lock on {self.lock_file} (lock id: {id(self)}): "
                f"already holding a {opposite} lock ({direction} not allowed)"
            )
            raise RuntimeError(msg)
        if mode == "write" and (cur := threading.get_ident()) != self._write_thread_id:
            msg = (
                f"Cannot acquire write lock on {self.lock_file} (lock id: {id(self)}) "
                f"from thread {cur} while it is held by thread {self._write_thread_id}"
            )
            raise RuntimeError(msg)
        self._lock_level += 1
        return AcquireReturnProxy(lock=self)

    def _acquire_transaction_lock(self, *, blocking: bool, timeout: float) -> None:
        if not blocking:
            acquired = self._transaction_lock.acquire(blocking=False)
        elif timeout == -1:
            acquired = self._transaction_lock.acquire(blocking=True)
        else:
            acquired = self._transaction_lock.acquire(blocking=True, timeout=timeout)
        if not acquired:
            raise Timeout(self.lock_file) from None

    def _raise_if_acquiring(self, operation: Literal["acquire", "close", "release"]) -> None:
        with self._internal_lock:
            active_in_current_thread = threading.get_ident() in self._acquisition_thread_ids
        if active_in_current_thread:  # pragma: no cover - exercised in an audit callback
            msg = f"Cannot {operation} ReadWriteLock on {self.lock_file} while acquisition is active in this thread"
            raise RuntimeError(msg)

    def _raise_if_process_changed(self, operation_pid: int) -> None:
        if _GETPID() != operation_pid or self._inherited:  # pragma: forked child
            msg = f"ReadWriteLock on {self.lock_file} was invalidated by fork(); construct a new instance"
            raise RuntimeError(msg)


def _connect(database: str, *, factory: type[_ForkSafeConnection], timeout: float) -> _ForkSafeConnection:
    _FORKED_DATABASES.note_sqlite_use()
    return sqlite3.connect(
        database,
        check_same_thread=False,
        factory=factory,
        cached_statements=0,
        timeout=timeout,
    )


@contextmanager
def _sqlite_transition() -> Generator[None]:
    _SQLITE_TRANSITION_CONTEXT.depth += 1
    try:
        yield
    finally:
        _SQLITE_TRANSITION_CONTEXT.depth -= 1


def _abort_forked_sqlite_transition() -> None:  # pragma: forked child
    if _SQLITE_TRANSITION_CONTEXT.depth:
        os._exit(_UNSAFE_FORK_EXIT_STATUS)  # inherited SQLite handles cannot be used or closed safely


def _track_sqlite_use(event: str, _args: tuple[object, ...]) -> None:  # audit payloads are heterogeneous
    if event == "sqlite3.connect":
        _FORKED_DATABASES.note_sqlite_use()


def timeout_for_sqlite(timeout: float, *, blocking: bool, already_waited: float) -> int:
    if blocking is False:
        return 0

    if timeout == -1:
        return _MAX_SQLITE_TIMEOUT_MS

    if timeout < 0:
        msg = "timeout must be a non-negative number or -1"
        raise ValueError(msg)

    timeout_ms = int((max(timeout - already_waited, 0) if timeout > 0 else timeout) * 1000)
    if timeout_ms > _MAX_SQLITE_TIMEOUT_MS or timeout_ms < 0:
        _LOGGER.warning("timeout %s is too large for SQLite, using %s ms instead", timeout, _MAX_SQLITE_TIMEOUT_MS)
        return _MAX_SQLITE_TIMEOUT_MS
    return timeout_ms


_register_fork_object(_CONNECTION_ESCROW)
_register_fork_object(_FORKED_DATABASES)
_register_fork_class(ReadWriteLock)
if _IS_PYPY:
    sys.addaudithook(_track_sqlite_use)  # pragma: pypy cover
if hasattr(os, "register_at_fork"):  # pragma: needs fork
    os.register_at_fork(after_in_child=_abort_forked_sqlite_transition)

__all__ = [
    "ReadWriteLock",
    "timeout_for_sqlite",
]
