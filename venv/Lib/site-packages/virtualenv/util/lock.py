"""holds locking functionality that works across processes."""

from __future__ import annotations

import logging
import os
from abc import ABC, abstractmethod
from contextlib import contextmanager, suppress
from pathlib import Path
from threading import Lock, RLock
from typing import TYPE_CHECKING

from filelock import FileLock, Timeout

if TYPE_CHECKING:
    from collections.abc import Iterator
    from types import TracebackType

LOGGER = logging.getLogger(__name__)


class _CountedFileLock(FileLock):
    def __init__(self, lock_file: str) -> None:
        parent = os.path.dirname(lock_file)
        with suppress(OSError):
            os.makedirs(parent, exist_ok=True)

        super().__init__(lock_file)
        self.count = 0
        self.thread_safe = RLock()

    def acquire(  # ty: ignore[invalid-method-override]
        self,
        timeout: float | None = None,
        poll_interval: float = 0.05,
    ) -> None:
        if not self.thread_safe.acquire(timeout=-1 if timeout is None else timeout):
            raise Timeout(self.lock_file)
        if self.count == 0:
            try:
                super().acquire(timeout, poll_interval)
            except BaseException:
                self.thread_safe.release()
                raise
        self.count += 1

    def release(self, force: bool = False) -> None:  # ruff:ignore[boolean-default-value-positional-argument]
        with self.thread_safe:
            if self.count > 0:
                if self.count == 1:
                    super().release(force=force)
                self.count -= 1
                if self.count == 0:
                    # if we have no more users of this lock, release the thread lock
                    self.thread_safe.release()


_lock_store = {}
_store_lock = Lock()


class PathLockBase(ABC):
    def __init__(self, folder: str | Path) -> None:
        path = Path(folder)
        self.path = path.resolve() if path.exists() else path

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.path})"

    def __truediv__(self, other: str) -> PathLockBase:
        return type(self)(self.path / other)

    @abstractmethod
    def __enter__(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None:
        raise NotImplementedError

    @abstractmethod
    @contextmanager
    def lock_for_key(self, name: str, no_block: bool = False) -> Iterator[None]:  # ruff:ignore[boolean-default-value-positional-argument]
        raise NotImplementedError

    @abstractmethod
    @contextmanager
    def non_reentrant_lock_for_key(self, name: str) -> Iterator[None]:
        raise NotImplementedError


class ReentrantFileLock(PathLockBase):
    def __init__(self, folder: str | Path) -> None:
        super().__init__(folder)
        self._lock = None

    def _create_lock(self, name: str = "") -> _CountedFileLock:
        lock_file = str(self.path / f"{name}.lock")
        with _store_lock:
            if lock_file not in _lock_store:
                _lock_store[lock_file] = _CountedFileLock(lock_file)
            return _lock_store[lock_file]

    @staticmethod
    def _del_lock(lock: _CountedFileLock | None) -> None:
        if lock is not None:
            with _store_lock, lock.thread_safe:
                if lock.count == 0:
                    _lock_store.pop(lock.lock_file, None)

    def __del__(self) -> None:
        self._del_lock(self._lock)

    def __enter__(self) -> None:
        self._lock = self._create_lock()
        self._lock_file(self._lock)

    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None:
        self._release(self._lock)  # ty: ignore[invalid-argument-type]
        self._del_lock(self._lock)
        self._lock = None

    def _lock_file(self, lock: _CountedFileLock, no_block: bool = False) -> None:  # ruff:ignore[boolean-default-value-positional-argument]
        # multiple processes might be trying to get a first lock... so we cannot check if this directory exist without
        # a lock, but that lock might then become expensive, and it's not clear where that lock should live.
        # Instead here we just ignore if we fail to create the directory.
        with suppress(OSError):
            os.makedirs(str(self.path), exist_ok=True)

        try:
            lock.acquire(0.0001)
        except Timeout:
            if no_block:
                raise
            LOGGER.debug("lock file %s present, will block until released", lock.lock_file)
            lock.release()  # release the acquire try from above
            lock.acquire()

    @staticmethod
    def _release(lock: _CountedFileLock) -> None:
        lock.release()

    @contextmanager
    def lock_for_key(self, name: str, no_block: bool = False) -> Iterator[None]:  # ruff:ignore[boolean-default-value-positional-argument]
        lock = self._create_lock(name)
        try:
            with self._lock_and_yield(lock, no_block):
                yield
        finally:
            self._del_lock(lock)
            lock = None

    @contextmanager
    def _lock_and_yield(self, lock: _CountedFileLock, no_block: bool) -> Iterator[None]:
        self._lock_file(lock, no_block)
        try:
            yield
        finally:
            self._release(lock)

    @contextmanager
    def non_reentrant_lock_for_key(self, name: str) -> Iterator[None]:
        with _CountedFileLock(str(self.path / f"{name}.lock")):
            yield


class NoOpFileLock(PathLockBase):
    def __enter__(self) -> None:
        raise NotImplementedError

    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None:
        raise NotImplementedError

    @contextmanager
    def lock_for_key(self, name: str, no_block: bool = False) -> Iterator[None]:  # ruff:ignore[unused-method-argument, boolean-default-value-positional-argument]
        yield

    @contextmanager
    def non_reentrant_lock_for_key(self, name: str) -> Iterator[None]:  # ruff:ignore[unused-method-argument]
        yield


__all__ = [
    "NoOpFileLock",
    "ReentrantFileLock",
    "Timeout",
]
