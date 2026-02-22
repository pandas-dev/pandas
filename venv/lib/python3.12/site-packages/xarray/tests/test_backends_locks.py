from __future__ import annotations

import threading

from xarray.backends import locks
from xarray.backends.locks import CombinedLock, SerializableLock


def test_threaded_lock() -> None:
    lock1 = locks._get_threaded_lock("foo")
    assert isinstance(lock1, type(threading.Lock()))
    lock2 = locks._get_threaded_lock("foo")
    assert lock1 is lock2

    lock3 = locks._get_threaded_lock("bar")
    assert lock1 is not lock3


def test_combined_lock_locked_returns_false_when_no_locks_acquired() -> None:
    """CombinedLock.locked() should return False when no locks are held."""
    lock1 = threading.Lock()
    lock2 = threading.Lock()
    combined = CombinedLock([lock1, lock2])

    assert combined.locked() is False
    assert lock1.locked() is False
    assert lock2.locked() is False


def test_combined_lock_locked_returns_true_when_one_lock_acquired() -> None:
    """CombinedLock.locked() should return True when any lock is held."""
    lock1 = threading.Lock()
    lock2 = threading.Lock()
    combined = CombinedLock([lock1, lock2])

    lock1.acquire()
    try:
        assert combined.locked() is True
    finally:
        lock1.release()

    assert combined.locked() is False


def test_combined_lock_locked_returns_true_when_all_locks_acquired() -> None:
    """CombinedLock.locked() should return True when all locks are held."""
    lock1 = threading.Lock()
    lock2 = threading.Lock()
    combined = CombinedLock([lock1, lock2])

    lock1.acquire()
    lock2.acquire()
    try:
        assert combined.locked() is True
    finally:
        lock1.release()
        lock2.release()

    assert combined.locked() is False


def test_combined_lock_locked_with_serializable_locks() -> None:
    """CombinedLock.locked() should work with SerializableLock instances."""
    lock1 = SerializableLock()
    lock2 = SerializableLock()
    combined = CombinedLock([lock1, lock2])

    assert combined.locked() is False

    lock1.acquire()
    try:
        assert combined.locked() is True
    finally:
        lock1.release()

    assert combined.locked() is False


def test_combined_lock_locked_with_context_manager() -> None:
    """CombinedLock.locked() should reflect state when using context manager."""
    lock1 = threading.Lock()
    lock2 = threading.Lock()
    combined = CombinedLock([lock1, lock2])

    assert combined.locked() is False

    with combined:
        assert combined.locked() is True

    assert combined.locked() is False
