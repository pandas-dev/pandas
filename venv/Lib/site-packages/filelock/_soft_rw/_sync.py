"""Cross-process and cross-host reader/writer lock built on :class:`SoftFileLock` primitives."""

from __future__ import annotations

import atexit
import hmac
import os
import re
import secrets
import socket
import stat
import sys
import threading
import time
import uuid
from contextlib import closing, contextmanager, suppress
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Final, Literal
from weakref import WeakValueDictionary

from filelock._api import (
    AcquireReturnProxy,
    _ensure_current_process,
    _fork_transition,
    _raise_grouped_errors,
    _register_fork_class,
    _register_fork_object,
    _register_owned_descriptor,
    _unregister_owned_descriptor,
)
from filelock._error import Timeout
from filelock._soft import SoftFileLock
from filelock._util import ensure_directory_exists, touch, write_all

if TYPE_CHECKING:
    from collections.abc import Callable, Generator


_Mode = Literal["read", "write"]
_BREAK_SUFFIX: Final[str] = ".break"
_MAX_MARKER_SIZE: Final[int] = 1024
_O_NOFOLLOW: Final[int] = getattr(os, "O_NOFOLLOW", 0)
_O_NONBLOCK: Final[int] = getattr(os, "O_NONBLOCK", 0)
# dirfd-relative I/O is a Unix-only optimization; Windows cannot ``os.open()`` a directory at all, and
# its ``os`` module skips dir_fd support entirely. When disabled, callers fall back to full-path ops.
_SUPPORTS_DIR_FD: Final[bool] = sys.platform != "win32" and os.open in os.supports_dir_fd

_ALL_INSTANCES: Final[WeakValueDictionary[int, SoftReadWriteLock]] = WeakValueDictionary()
_ALL_INSTANCES_LOCK: threading.Lock = threading.Lock()
_SINGLETONS_UNDER_CONSTRUCTION: Final[set[Path]] = set()


class _SoftRWMeta(type):
    _instances: WeakValueDictionary[Path, SoftReadWriteLock]
    _instances_lock: threading.RLock

    def __call__(  # ruff:ignore[too-many-arguments]  # forwards the public constructor's documented parameters
        cls,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        *,
        blocking: bool = True,
        is_singleton: bool = True,
        heartbeat_interval: float = 30.0,
        stale_threshold: float | None = None,
        poll_interval: float = 0.25,
    ) -> SoftReadWriteLock:
        _ensure_current_process()
        if not is_singleton:
            return super().__call__(
                lock_file,
                timeout,
                blocking=blocking,
                is_singleton=is_singleton,
                heartbeat_interval=heartbeat_interval,
                stale_threshold=stale_threshold,
                poll_interval=poll_interval,
            )

        normalized = Path(lock_file).resolve()
        with cls._instances_lock:
            instance = cls._instances.get(normalized)
            if instance is None:
                if normalized in _SINGLETONS_UNDER_CONSTRUCTION:  # pragma: needs fork
                    msg = f"Singleton lock construction is already active for {lock_file!s}"
                    raise RuntimeError(msg)
                construction_pid = os.getpid()
                _SINGLETONS_UNDER_CONSTRUCTION.add(normalized)
                try:
                    instance = super().__call__(
                        lock_file,
                        timeout,
                        blocking=blocking,
                        is_singleton=is_singleton,
                        heartbeat_interval=heartbeat_interval,
                        stale_threshold=stale_threshold,
                        poll_interval=poll_interval,
                    )
                finally:
                    _SINGLETONS_UNDER_CONSTRUCTION.discard(normalized)
                if os.getpid() != construction_pid:  # pragma: needs fork
                    msg = "Lock construction cannot continue after fork; construct a new lock in the child"
                    raise RuntimeError(msg)
                cls._instances[normalized] = instance
            elif instance.timeout != timeout or instance.blocking != blocking:
                msg = (
                    f"Singleton lock created with timeout={instance.timeout}, blocking={instance.blocking},"
                    f" cannot be changed to timeout={timeout}, blocking={blocking}"
                )
                raise ValueError(msg)
            return instance


class SoftReadWriteLock(metaclass=_SoftRWMeta):
    """
    Cross-process and cross-host reader/writer lock built on :class:`SoftFileLock` primitives.

    Use this class instead of :class:`~filelock.ReadWriteLock` when the lock file lives on a network
    filesystem (NFS, Lustre with ``-o flock``, HPC cluster shared storage). ``ReadWriteLock`` is backed
    by SQLite and cannot run on NFS because SQLite's ``fcntl`` locking is unreliable there.

    Layout on disk for a lock at ``foo.lock``:

    - ``foo.lock.state`` — a :class:`SoftFileLock` taken only during state transitions (microseconds).
    - ``foo.lock.write`` — writer marker; its presence means a writer is claiming or holding the lock.
    - ``foo.lock.readers/<host>.<pid>.<uuid>`` — one file per reader.

    Each marker stores a random token (``secrets.token_hex(16)``), the holder's pid, and the holder's
    hostname. A daemon heartbeat thread refreshes ``mtime`` on every held marker. A marker whose mtime
    has not advanced in ``stale_threshold`` seconds may be evicted by any process on any host, giving
    correct behavior when a compute node crashes with a lock held.

    Writer acquire is two-phase and writer-preferring: phase 1 claims ``.write`` (blocking any new
    reader), phase 2 waits for existing readers to drain. Writer starvation is impossible.

    Reentrancy, upgrade/downgrade rules, thread pinning, and singleton caching by resolved path match
    :class:`~filelock.ReadWriteLock`.

    Forking invalidates the inherited instance in the child so the child cannot double-own the lock with its parent;
    ``release()`` on that instance is a no-op, and the child must construct a new instance if it needs a lock.

    Trust boundary: protects against same-UID non-cooperating processes (one host or cross-host) and
    same-host different-UID users via ``0o600`` / ``0o700`` permissions. Does not protect against root
    compromise, NTP tampering on same-UID cross-host nodes, or multi-tenant mounts where hostile
    co-tenants share the UID.

    :param lock_file: path to the lock file; sidecar state/write/readers live next to it
    :param timeout: maximum wait time in seconds; ``-1`` means block indefinitely
    :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately on contention
    :param is_singleton: if ``True``, reuse existing instances for the same resolved path
    :param heartbeat_interval: seconds between heartbeat refreshes; default 30 s
    :param stale_threshold: seconds of ``mtime`` inactivity before a marker is stale; defaults to
        ``3 * heartbeat_interval``, matching etcd's ``LeaseKeepAlive`` convention
    :param poll_interval: seconds between acquire retries under contention; default 0.25 s

    .. versionadded:: 3.27.0

    """

    _instances: WeakValueDictionary[Path, SoftReadWriteLock] = WeakValueDictionary()
    _instances_lock = threading.RLock()

    def __init__(  # ruff:ignore[too-many-arguments]  # public constructor: one parameter per documented lock option
        self,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        *,
        blocking: bool = True,
        is_singleton: bool = True,  # ruff:ignore[unused-method-argument]  # consumed by _SoftRWMeta.__call__
        heartbeat_interval: float = 30.0,
        stale_threshold: float | None = None,
        poll_interval: float = 0.25,
    ) -> None:
        self._creator_pid = os.getpid()
        if heartbeat_interval <= 0:
            msg = f"heartbeat_interval must be positive, got {heartbeat_interval}"
            raise ValueError(msg)
        if stale_threshold is None:
            stale_threshold = heartbeat_interval * 3
        if stale_threshold <= heartbeat_interval:
            msg = f"stale_threshold must exceed heartbeat_interval ({stale_threshold} <= {heartbeat_interval})"
            raise ValueError(msg)
        if poll_interval <= 0:
            msg = f"poll_interval must be positive, got {poll_interval}"
            raise ValueError(msg)

        self.lock_file: str = os.fspath(lock_file)
        self.timeout: float = timeout
        self.blocking: bool = blocking
        self.heartbeat_interval: float = heartbeat_interval
        self.stale_threshold: float = stale_threshold
        self.poll_interval: float = poll_interval

        self._paths = _Paths(
            state=f"{self.lock_file}.state",
            write=f"{self.lock_file}.write",
            readers=f"{self.lock_file}.readers",
        )
        ensure_directory_exists(self.lock_file)
        self._locks = _Locks(
            internal=threading.Lock(),
            transaction=threading.Lock(),
            state=SoftFileLock(self._paths.state, timeout=-1),
        )
        self._readers_dir_fd: int | None = None
        self._readers_dir_fd_token: int | None = None
        self._hold: _Hold | None = None
        self._closed: bool = False

        with _ALL_INSTANCES_LOCK:
            _ALL_INSTANCES[id(self)] = self
        _register_fork_object(self)

    @classmethod
    def _reset_class_after_fork(cls) -> None:  # pragma: forked child
        global _ALL_INSTANCES_LOCK  # ruff:ignore[global-statement]  # rebinds the module lock to a fresh one in the fork child
        _ALL_INSTANCES_LOCK = threading.Lock()
        cls._instances = WeakValueDictionary()
        cls._instances_lock = threading.RLock()
        _SINGLETONS_UNDER_CONSTRUCTION.clear()

    @contextmanager
    def read_lock(self, timeout: float | None = None, *, blocking: bool | None = None) -> Generator[None]:
        """
        Context manager that acquires and releases a shared read lock.

        Falls back to instance defaults for *timeout* and *blocking* when ``None``.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately; ``None`` uses the instance default

        :raises RuntimeError: if a write lock is already held on this instance
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
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

        :raises RuntimeError: if a read lock is already held, or a write lock is held by a different thread
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        self.acquire_write(timeout, blocking=blocking)
        try:
            yield
        finally:
            self.release()

    def acquire_read(self, timeout: float | None = None, *, blocking: bool | None = None) -> AcquireReturnProxy:
        """
        Acquire a shared read lock.

        If this instance already holds a read lock, the lock level is incremented (reentrant). Attempting to acquire a
        read lock while holding a write lock raises :class:`RuntimeError` (downgrade not allowed). On the 0→1
        transition a daemon heartbeat thread is started that refreshes the reader marker's ``mtime`` every
        ``heartbeat_interval`` seconds so peers on other hosts do not evict the marker as stale.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default; ``-1`` means block
            indefinitely
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately when the lock is unavailable;
            ``None`` uses the instance default

        :returns: a proxy that can be used as a context manager to release the lock

        :raises RuntimeError: if a write lock is already held on this instance, if this instance was invalidated by
            :func:`os.fork`, or if :meth:`close` was called
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        return self._acquire("read", timeout, blocking=blocking)

    def acquire_write(self, timeout: float | None = None, *, blocking: bool | None = None) -> AcquireReturnProxy:
        """
        Acquire an exclusive write lock.

        If this instance already holds a write lock from the same thread, the lock level is incremented (reentrant).
        Attempting to acquire a write lock while holding a read lock raises :class:`RuntimeError` (upgrade not
        allowed). Write locks are pinned to the acquiring thread: a different thread trying to re-enter also raises
        :class:`RuntimeError`.

        Writer acquisition runs in two phases. Phase 1 atomically claims ``<path>.write`` via ``O_CREAT | O_EXCL``,
        which immediately blocks any new reader on any host. Phase 2 waits for existing readers to drain. Writer
        starvation is impossible: new readers see ``<path>.write`` during phase 2 and wait behind the pending writer.

        :param timeout: maximum wait time in seconds, or ``None`` to use the instance default; ``-1`` means block
            indefinitely
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately when the lock is unavailable;
            ``None`` uses the instance default

        :returns: a proxy that can be used as a context manager to release the lock

        :raises RuntimeError: if a read lock is already held, if a write lock is held by a different thread, if this
            instance was invalidated by :func:`os.fork`, or if :meth:`close` was called
        :raises Timeout: if the lock cannot be acquired within *timeout* seconds

        """
        return self._acquire("write", timeout, blocking=blocking)

    @classmethod
    def get_lock(
        cls,
        lock_file: str | os.PathLike[str],
        timeout: float = -1,
        *,
        blocking: bool = True,
    ) -> SoftReadWriteLock:
        """
        Return the singleton :class:`SoftReadWriteLock` for *lock_file*.

        :param lock_file: path to the lock file; sidecar state/write/readers live next to it
        :param timeout: maximum wait time in seconds; ``-1`` means block indefinitely
        :param blocking: if ``False``, raise :class:`~filelock.Timeout` immediately when the lock is unavailable

        :returns: the singleton lock instance

        :raises ValueError: if an instance already exists for this path with different *timeout* or *blocking* values

        """
        return cls(lock_file, timeout, blocking=blocking)

    def close(self) -> None:
        """
        Release any held lock and release internal filesystem resources.

        Idempotent. After calling this method the instance can no longer acquire locks — subsequent acquires raise
        :class:`RuntimeError`. A fork-invalidated instance is closed without raising.
        """
        if self._creator_pid != os.getpid():  # pragma: forked child
            return
        self.release(force=True)
        with self._locks.internal:
            if self._closed:
                return
            self._closed = True
            if self._readers_dir_fd is not None:  # pragma: needs dir-fd
                with _fork_transition():
                    if self._readers_dir_fd_token is not None:  # pragma: needs dir-fd
                        _unregister_owned_descriptor(self._readers_dir_fd_token)
                    self._readers_dir_fd_token = None
                    fd, self._readers_dir_fd = self._readers_dir_fd, None
                    with suppress(OSError):  # pragma: needs dir-fd
                        os.close(fd)

    def release(self, *, force: bool = False) -> None:
        """
        Release one level of the current lock.

        When the lock level reaches zero the heartbeat thread is stopped and the held marker file is unlinked. On a
        fork-invalidated instance (that is, the child of a :func:`os.fork` call made while the parent held a lock)
        this method is a no-op so inherited ``with`` blocks can unwind cleanly in the child.

        :param force: if ``True``, release the lock completely regardless of the current lock level

        :raises RuntimeError: if no lock is currently held and *force* is ``False``

        """
        if self._creator_pid != os.getpid():  # pragma: forked child
            return
        with self._locks.internal:
            hold = self._hold
            if hold is None:
                if force:
                    return
                msg = f"Cannot release a lock on {self.lock_file} (lock id: {id(self)}) that is not held"
                raise RuntimeError(msg)
            if force:
                hold.level = 0
            else:
                hold.level -= 1
            if hold.level > 0:
                return
            self._hold = None

        # Order matters: signal → join → unlink. A late tick on a deleted marker is harmless and the
        # heartbeat's token check would catch a re-acquisition race, but joining first removes that race.
        hold.heartbeat_stop.set()
        hold.heartbeat_thread.join(timeout=self.heartbeat_interval + 1.0)
        if hold.is_reader:
            _unlink(hold.marker_name, dir_fd=self._readers_dir_fd)
        else:
            self._unlink_writer_marker_if_ours(hold.token)

    def _unlink_writer_marker_if_ours(self, token: str) -> None:
        # Remove the writer marker only while it still carries our token. If this holder was paused long
        # enough (a stop-the-world GC pause, SIGSTOP, a suspended VM) for a peer to evict the marker as
        # stale and claim the writer slot itself, the file now at <path>.write is the peer's live marker;
        # unlinking it by path would let a second writer through and break mutual exclusion. The state lock
        # serializes this against a concurrent break/claim, and the heartbeat is already stopped, so the
        # token we read is authoritative. Mirrors the token re-check the stale-break path already does.
        with self._locks.state:
            if (read := _read_marker(self._paths.write)) is None:
                return
            info, _ = read
            if info is None or not hmac.compare_digest(info.token, token):
                return
            _unlink(self._paths.write)

    def _acquire(
        self,
        mode: _Mode,
        timeout: float | None,
        *,
        blocking: bool | None,
    ) -> AcquireReturnProxy:
        if self._creator_pid != os.getpid():  # pragma: forked child
            msg = f"SoftReadWriteLock on {self.lock_file} was invalidated by fork(); construct a new instance"
            raise RuntimeError(msg)
        timeout = self.timeout if timeout is None else timeout
        blocking = self.blocking if blocking is None else blocking

        with self._locks.internal:
            if self._closed:
                msg = f"SoftReadWriteLock on {self.lock_file} has been closed"
                raise RuntimeError(msg)
            if self._hold is not None:
                return self._validate_reentrant(mode)

        start = time.perf_counter()
        if not blocking:
            acquired = self._locks.transaction.acquire(blocking=False)
        elif timeout == -1:
            acquired = self._locks.transaction.acquire(blocking=True)
        else:
            acquired = self._locks.transaction.acquire(blocking=True, timeout=timeout)
        if not acquired:
            raise Timeout(self.lock_file) from None
        try:
            return self._do_acquire_inner(mode, timeout, start, blocking=blocking)
        finally:
            self._locks.transaction.release()

    def _do_acquire_inner(
        self,
        mode: _Mode,
        effective_timeout: float,
        start: float,
        *,
        blocking: bool,
    ) -> AcquireReturnProxy:
        with self._locks.internal:
            if self._hold is not None:
                return self._validate_reentrant(mode)
        deadline = None if effective_timeout == -1 else start + effective_timeout
        token = secrets.token_hex(16)
        if mode == "write":
            marker_name, is_reader = self._acquire_writer_slot(token, deadline=deadline, blocking=blocking)
        else:
            marker_name, is_reader = self._acquire_reader_slot(token, deadline=deadline, blocking=blocking)
        stop_event = threading.Event()
        heartbeat = _HeartbeatThread(
            refresh=self._refresh_marker,
            interval=self.heartbeat_interval,
            stop_event=stop_event,
            name=f"filelock-heartbeat-{id(self):x}",
        )
        with self._locks.internal:
            self._hold = _Hold(
                level=1,
                mode=mode,
                write_thread_id=threading.get_ident() if mode == "write" else None,
                marker_name=marker_name,
                is_reader=is_reader,
                token=token,
                heartbeat_thread=heartbeat,
                heartbeat_stop=stop_event,
            )
        heartbeat.start()
        return AcquireReturnProxy(lock=self)

    def _validate_reentrant(self, mode: _Mode) -> AcquireReturnProxy:
        hold = self._hold
        assert hold is not None  # ruff:ignore[assert]  # callers dispatch here only inside the self._hold is not None branch
        if hold.mode != mode:
            opposite = "write" if mode == "read" else "read"
            direction = "downgrade" if mode == "read" else "upgrade"
            msg = (
                f"Cannot acquire {mode} lock on {self.lock_file} (lock id: {id(self)}): "
                f"already holding a {opposite} lock ({direction} not allowed)"
            )
            raise RuntimeError(msg)
        if mode == "write" and (cur := threading.get_ident()) != hold.write_thread_id:
            msg = (
                f"Cannot acquire write lock on {self.lock_file} (lock id: {id(self)}) "
                f"from thread {cur} while it is held by thread {hold.write_thread_id}"
            )
            raise RuntimeError(msg)
        hold.level += 1
        return AcquireReturnProxy(lock=self)

    def _acquire_writer_slot(
        self,
        token: str,
        *,
        deadline: float | None,
        blocking: bool,
    ) -> tuple[str, bool]:
        # Phase 2 scans readers/ via dirfd (where supported), so we need it open even though writers never
        # create files inside.
        self._open_readers_dir()

        def try_claim_writer() -> bool:
            with self._locks.state:
                return self._claim_writer_marker(token)

        def readers_drained_touching() -> bool:
            with self._locks.state:
                # A peer may replace an expired marker while this process pauses. Refresh only our token; touching a
                # successor's marker would let this acquisition proceed without owning the writer slot.
                if not self._touch_writer_marker_if_ours(token) and not self._claim_writer_marker(token):
                    return False
                self._break_stale_readers(time.time())
                return not self._any_readers()

        self._wait_for(try_claim_writer, deadline=deadline, blocking=blocking)
        try:
            self._wait_for(readers_drained_touching, deadline=deadline, blocking=blocking)
        except Timeout:
            # Give up our writer claim so readers can make progress again, but only while the marker is
            # still ours: a peer may have evicted it as stale and claimed the slot while phase 2 waited.
            self._unlink_writer_marker_if_ours(token)
            raise
        return self._paths.write, False

    def _claim_writer_marker(self, token: str) -> bool:
        # Claim the writer slot for ``token``. Must be called holding ``self._locks.state``. Evicts a
        # stale marker first, then refuses to claim while a live ``.write`` exists so a peer holding the
        # slot is waited out instead of overwritten.
        _break_stale_marker(self._paths.write, stale_threshold=self.stale_threshold, now=time.time())
        if _file_exists(self._paths.write):
            return False
        try:
            _atomic_create_marker(self._paths.write, token)
        except FileExistsError:
            return False
        return True

    def _touch_writer_marker_if_ours(self, token: str) -> bool:
        # Refresh the writer marker through a single O_NOFOLLOW fd, but only while it still carries our
        # token. Returns False when the marker is gone or now belongs to a peer that reclaimed the slot,
        # so the caller can re-claim rather than keep a stranger's marker alive. Mirrors _refresh_marker.
        fd = _open_marker(self._paths.write)
        if fd is None:
            return False
        try:
            try:
                data = os.read(fd, _MAX_MARKER_SIZE + 1)
            except OSError:  # pragma: no cover - e.g. EAGAIN from a hostile FIFO that has a writer attached
                return False
            info = _parse_marker_bytes(data)
            if info is None or not hmac.compare_digest(info.token, token):
                return False
            with suppress(OSError):
                touch(self._paths.write, fd=fd)
            return True
        finally:
            os.close(fd)

    def _acquire_reader_slot(
        self,
        token: str,
        *,
        deadline: float | None,
        blocking: bool,
    ) -> tuple[str, bool]:
        self._open_readers_dir()
        reader_name = f"{uuid.uuid4().hex}.{os.getpid()}"
        dir_fd = self._readers_dir_fd
        full_reader_path = str(Path(self._paths.readers) / reader_name)

        def try_claim_reader() -> bool:
            with self._locks.state:
                _break_stale_marker(self._paths.write, stale_threshold=self.stale_threshold, now=time.time())
                if _file_exists(self._paths.write):
                    return False
                if dir_fd is not None:  # pragma: needs dir-fd
                    _atomic_create_marker(reader_name, token, dir_fd=dir_fd)
                else:  # pragma: win32 cover
                    _atomic_create_marker(full_reader_path, token)
                return True

        self._wait_for(try_claim_reader, deadline=deadline, blocking=blocking)
        return (reader_name if dir_fd is not None else full_reader_path), True

    def _wait_for(
        self,
        predicate: Callable[[], bool],
        *,
        deadline: float | None,
        blocking: bool,
    ) -> None:
        while True:
            if predicate():
                return
            now = time.perf_counter()
            if not blocking:
                raise Timeout(self.lock_file)
            if deadline is not None and now >= deadline:
                raise Timeout(self.lock_file)
            sleep_for = self.poll_interval
            if deadline is not None:
                sleep_for = min(sleep_for, max(deadline - now, 0.0))
            time.sleep(sleep_for)

    def _open_readers_dir(self) -> None:
        readers_path = Path(self._paths.readers)
        with suppress(FileExistsError):
            readers_path.mkdir(mode=0o700)
        # mkdir has no O_NOFOLLOW, so verify via lstat that we did not land on an attacker-placed symlink
        # or a regular file before we open or scan inside.
        st = os.lstat(self._paths.readers)
        if stat.S_ISLNK(st.st_mode) or not stat.S_ISDIR(st.st_mode):
            msg = f"{self._paths.readers} exists but is not a directory or is a symlink; refusing to use it"
            raise RuntimeError(msg)
        if self._readers_dir_fd is None and _SUPPORTS_DIR_FD:  # pragma: needs dir-fd
            with _fork_transition():
                fd = os.open(self._paths.readers, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | _O_NOFOLLOW)
                try:
                    token = _register_owned_descriptor(fd)
                except BaseException as registration_error:
                    try:
                        os.close(fd)
                    except BaseException as close_error:  # ruff:ignore[blind-except]  # both errors surface via the group below
                        _raise_grouped_errors(
                            "reader directory registration and descriptor close both failed",
                            registration_error,
                            close_error,
                        )
                    raise
                self._readers_dir_fd = fd
                self._readers_dir_fd_token = token

    def _any_readers(self) -> bool:
        with closing(self._iter_reader_entries()) as entries:
            for _ in entries:
                return True
        return False

    def _iter_reader_entries(self) -> Generator[tuple[str, bool]]:
        """
        Yield ``(name, dirfd_relative)`` pairs for every live reader marker.

        ``dirfd_relative`` is ``True`` when *name* should be passed to ``dir_fd=``-aware syscalls; ``False``
        when *name* is a full path because dirfd-relative I/O is unavailable on this platform.

        A consumer that stops early must close this generator: while suspended it holds the ``scandir`` handle open,
        and leaving that to the collector surfaces as an unraisable exception inside whatever runs next.
        """
        if self._readers_dir_fd is not None:  # pragma: needs dir-fd
            with os.scandir(self._readers_dir_fd) as it:
                for entry in it:
                    if not _is_housekeeping_name(entry.name):
                        yield entry.name, True
            return
        readers_path = Path(self._paths.readers)  # pragma: win32 cover
        with os.scandir(readers_path) as it:  # pragma: win32 cover
            for entry in it:  # pragma: win32 cover
                if not _is_housekeeping_name(entry.name):  # pragma: win32 cover
                    yield str(readers_path / entry.name), False  # pragma: win32 cover

    def _break_stale_readers(self, now: float) -> None:
        names: list[tuple[str, int | None]] = []
        try:
            with closing(self._iter_reader_entries()) as entries:
                for name, dirfd_relative in entries:
                    names.append((name, self._readers_dir_fd if dirfd_relative else None))
        except OSError:  # pragma: no cover - transient NFS scandir hiccup
            return
        for name, fd in names:
            _break_stale_marker(name, stale_threshold=self.stale_threshold, now=now, dir_fd=fd)

    def _refresh_marker(self) -> bool:
        with self._locks.internal:
            hold = self._hold
            if hold is None:  # pragma: no cover - race between stop_event.set and join
                return False
            marker_name = hold.marker_name
            token = hold.token
            dir_fd = self._readers_dir_fd if hold.is_reader else None

        # Open once with O_NOFOLLOW and touch that exact descriptor. Refreshing through the verified fd
        # (instead of re-opening by name) closes the window where a peer unlinks our marker and drops a symlink
        # or a different file at the path between the read and the touch: utime then lands on the inode we
        # verified, or nowhere. Only an unambiguous loss stops the heartbeat: the marker gone, or a peer's token
        # in its place. A transient filesystem error (ESTALE / EIO on the NFS-style filesystems this lock targets)
        # keeps the heartbeat alive to retry next tick, the way the touch below already does, so one blip does not
        # silently drop a held lock.
        try:
            fd = _open_marker_fd(marker_name, dir_fd=dir_fd)
        except FileNotFoundError:
            return False
        except OSError:
            return True
        try:
            try:
                data = _read_marker_fd(fd)
            except OSError:  # a transient read error or EAGAIN from a hostile FIFO; retry rather than drop the lock
                return True
            info = _parse_marker_bytes(data)
            # Token mismatch means another process already evicted our marker and created its own; stop the
            # thread so it does not keep a stranger's file alive.
            if info is None or not hmac.compare_digest(info.token, token):
                return False
            # A transient touch failure (ESTALE / EIO on the NFS-style filesystems this lock targets) must not
            # kill the heartbeat thread: the read above just confirmed the marker is still ours, so swallow the
            # error and retry on the next tick rather than letting the lease lapse while we still hold the lock.
            with suppress(OSError):
                touch(marker_name, fd=fd)
            return True
        finally:
            os.close(fd)

    def _reset_after_fork_in_child(self) -> None:  # pragma: forked child
        self._locks = _Locks(
            internal=threading.Lock(),
            transaction=threading.Lock(),
            state=self._locks.state,
        )
        self._hold = None
        self._readers_dir_fd = None
        self._readers_dir_fd_token = None


class _HeartbeatThread(threading.Thread):
    def __init__(
        self,
        refresh: Callable[[], bool],
        interval: float,
        stop_event: threading.Event,
        name: str,
    ) -> None:
        super().__init__(name=name, daemon=True)
        self._refresh = refresh
        self._interval = interval
        self._stop_event = stop_event

    def run(self) -> None:
        while not self._stop_event.wait(self._interval):
            if not self._refresh():
                self._stop_event.set()
                return


def _read_marker(name: str, *, dir_fd: int | None = None) -> tuple[_MarkerInfo | None, float] | None:
    fd = _open_marker(name, dir_fd=dir_fd)
    if fd is None:
        return None
    try:
        st = os.fstat(fd)
        # A legitimate marker is a regular file, so anything else at the path (a FIFO, say) is reported as a
        # malformed marker (its mtime still drives stale eviction) without being read. Reading is where
        # platforms diverge: an empty non-blocking read yields 0 bytes on Linux/macOS but EAGAIN on FreeBSD,
        # and the EAGAIN used to abort the stale-break and wedge the acquire until timeout (#587).
        if not stat.S_ISREG(st.st_mode):  # pragma: needs fifo
            return None, st.st_mtime
        data = os.read(fd, _MAX_MARKER_SIZE + 1)
    except OSError:  # pragma: no cover - marker vanished or turned unreadable between open and read
        return None
    finally:
        os.close(fd)
    return _parse_marker_bytes(data), st.st_mtime


def _read_marker_fd(fd: int) -> bytes:
    return os.read(fd, _MAX_MARKER_SIZE + 1)


def _open_marker_fd(name: str, *, dir_fd: int | None = None) -> int:
    # The file is ours; these guard a hostile mid-flight swap. O_NOFOLLOW rejects a symlink; O_NONBLOCK keeps
    # a real FIFO from blocking the open forever, so it reads as a malformed marker instead of wedging a peer
    # that holds the state lock.
    flags = os.O_RDONLY | _O_NOFOLLOW | _O_NONBLOCK
    return os.open(name, flags, dir_fd=dir_fd) if _SUPPORTS_DIR_FD and dir_fd is not None else os.open(name, flags)


def _open_marker(name: str, *, dir_fd: int | None = None) -> int | None:
    try:
        return _open_marker_fd(name, dir_fd=dir_fd)
    except OSError:
        return None


def _parse_marker_bytes(data: bytes) -> _MarkerInfo | None:
    # Trust nothing about attacker-controlled markers; any deviation returns None so callers fall through
    # to stale cleanup. ``re.match`` caches compiled patterns internally, so the regex is built only once
    # despite being defined inline.
    if not data or len(data) > _MAX_MARKER_SIZE:
        return None
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError:
        return None
    match = re.match(
        r"""
        \A                                  # start of string
        (?P<token>    [0-9a-f]{32}     ) \n # 128-bit hex token
        (?P<pid>      [1-9][0-9]{0,9}  ) \n # decimal pid: no leading zero, ≤ 10 digits
        (?P<hostname> [\x21-\x7e]{1,253})   # printable non-whitespace ASCII (RFC 1123 hostname limit)
        \n*                                 # tolerate sloppy writers that append extra newlines
        \Z                                  # end of string
        """,
        text,
        re.VERBOSE,
    )
    if match is None:
        return None
    pid = int(match["pid"], 10)
    if pid > 2**31 - 1:
        return None
    return _MarkerInfo(token=match["token"], pid=pid, hostname=match["hostname"])


def _unlink(name: str, *, dir_fd: int | None = None) -> None:
    with suppress(FileNotFoundError):
        if _SUPPORTS_DIR_FD and dir_fd is not None:  # pragma: needs dir-fd
            # Path.unlink has no dir_fd support, so we stay on os.unlink for the dirfd path.
            os.unlink(name, dir_fd=dir_fd)
        else:
            Path(name).unlink()


def _break_stale_marker(  # ruff:ignore[too-many-return-statements]  # each return is a distinct abort/commit point in the break protocol
    name: str,
    *,
    stale_threshold: float,
    now: float,
    dir_fd: int | None = None,
) -> bool:
    # Atomic break pattern: read → rename to unique break-name → re-verify → unlink. The rename gives us a
    # private name nobody else can touch; if the re-verify sees a newer mtime or a different token, the
    # legitimate holder's heartbeat fired between read and rename and we must abort (leaving the .break.*
    # file behind rather than rollback-renaming, because rollback is itself racy).
    if (read_result := _read_marker(name, dir_fd=dir_fd)) is None:
        return False
    info_before, mtime_before = read_result
    if now - mtime_before <= stale_threshold:
        return False
    if info_before is None:
        _unlink(name, dir_fd=dir_fd)
        return True

    break_name = f"{name}{_BREAK_SUFFIX}.{os.getpid()}.{secrets.token_hex(16)}"
    try:
        if _SUPPORTS_DIR_FD and dir_fd is not None:  # pragma: needs dir-fd
            os.rename(name, break_name, src_dir_fd=dir_fd, dst_dir_fd=dir_fd)
        else:
            Path(name).rename(break_name)
    except OSError:  # pragma: no cover - race where the marker vanishes between read and rename
        return False

    read_after = _read_marker(break_name, dir_fd=dir_fd)
    if read_after is None:  # pragma: no cover - race where a peer unlinks the break-name file
        return False
    info_after, mtime_after = read_after
    if info_after is None:  # pragma: no cover - content replaced post-rename by a racing peer
        _unlink(break_name, dir_fd=dir_fd)
        return True
    if not hmac.compare_digest(info_before.token, info_after.token):  # pragma: no cover - race only
        return False
    if mtime_after > mtime_before:  # pragma: no cover - heartbeat raced our rename
        return False
    _unlink(break_name, dir_fd=dir_fd)
    return True


def _atomic_create_marker(name: str, token: str, *, dir_fd: int | None = None) -> None:
    # O_NOFOLLOW blocks the symlink-overwrite attack where an attacker pre-creates the marker path as a
    # symlink pointing at a victim file. Mode 0o600 keeps the token unreadable to other users.
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | _O_NOFOLLOW
    if _SUPPORTS_DIR_FD and dir_fd is not None:  # pragma: needs dir-fd
        fd = os.open(name, flags, 0o600, dir_fd=dir_fd)
    else:
        fd = os.open(name, flags, 0o600)
    # Write the whole record before the marker counts as created. On failure remove it only while the path still names
    # the file we opened, so a rollback never deletes a marker a concurrent reader recreated at this name.
    identity: tuple[int, int] | None = None
    try:
        st = os.fstat(fd)
        identity = st.st_dev, st.st_ino
        write_all(fd, f"{token}\n{os.getpid()}\n{socket.gethostname()}\n".encode("ascii"))
    except BaseException:
        os.close(fd)
        if identity is not None and _same_file(name, identity, dir_fd=dir_fd):
            _unlink(name, dir_fd=dir_fd)
        raise
    else:
        os.close(fd)


def _same_file(name: str, identity: tuple[int, int], *, dir_fd: int | None) -> bool:
    try:
        st = os.lstat(name, dir_fd=dir_fd) if _SUPPORTS_DIR_FD and dir_fd is not None else os.lstat(name)
    except OSError:
        return False
    return (st.st_dev, st.st_ino) == identity


def _file_exists(path: str) -> bool:
    try:
        st = os.lstat(path)
    except FileNotFoundError:
        return False
    return stat.S_ISREG(st.st_mode)


def _is_housekeeping_name(name: str) -> bool:
    return name.startswith(".") or _BREAK_SUFFIX in name


@dataclass(frozen=True)
class _Paths:
    state: str
    write: str
    readers: str


@dataclass
class _Locks:
    internal: threading.Lock
    transaction: threading.Lock
    state: SoftFileLock


@dataclass(frozen=True)
class _MarkerInfo:
    token: str
    pid: int
    hostname: str


@dataclass
class _Hold:
    """Everything that exists only while a lock is held; ``None`` when the instance has no lock."""

    level: int
    mode: _Mode
    write_thread_id: int | None
    marker_name: str
    is_reader: bool
    token: str
    heartbeat_thread: _HeartbeatThread
    heartbeat_stop: threading.Event


def _cleanup_all_instances() -> None:  # pragma: no cover - runs from atexit at interpreter shutdown
    for instance in list(_ALL_INSTANCES.values()):
        with suppress(Exception):
            instance.release(force=True)


atexit.register(_cleanup_all_instances)
_register_fork_class(SoftReadWriteLock)


__all__ = [
    "SoftReadWriteLock",
]
