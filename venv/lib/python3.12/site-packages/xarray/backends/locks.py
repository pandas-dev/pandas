from __future__ import annotations

import multiprocessing
import threading
import uuid
import weakref
from collections.abc import Hashable, MutableMapping
from typing import Any, ClassVar
from weakref import WeakValueDictionary


# SerializableLock is adapted from Dask:
# https://github.com/dask/dask/blob/74e898f0ec712e8317ba86cc3b9d18b6b9922be0/dask/utils.py#L1160-L1224
# Used under the terms of Dask's license, see licenses/DASK_LICENSE.
class SerializableLock:
    """A Serializable per-process Lock

    This wraps a normal ``threading.Lock`` object and satisfies the same
    interface.  However, this lock can also be serialized and sent to different
    processes.  It will not block concurrent operations between processes (for
    this you should look at ``dask.multiprocessing.Lock`` or ``locket.lock_file``
    but will consistently deserialize into the same lock.

    So if we make a lock in one process::

        lock = SerializableLock()

    And then send it over to another process multiple times::

        bytes = pickle.dumps(lock)
        a = pickle.loads(bytes)
        b = pickle.loads(bytes)

    Then the deserialized objects will operate as though they were the same
    lock, and collide as appropriate.

    This is useful for consistently protecting resources on a per-process
    level.

    The creation of locks is itself not threadsafe.
    """

    _locks: ClassVar[WeakValueDictionary[Hashable, threading.Lock]] = (
        WeakValueDictionary()
    )
    token: Hashable
    lock: threading.Lock

    def __init__(self, token: Hashable | None = None):
        self.token = token or str(uuid.uuid4())
        if self.token in SerializableLock._locks:
            self.lock = SerializableLock._locks[self.token]
        else:
            self.lock = threading.Lock()
            SerializableLock._locks[self.token] = self.lock

    def acquire(self, *args, **kwargs):
        return self.lock.acquire(*args, **kwargs)

    def release(self, *args, **kwargs):
        return self.lock.release(*args, **kwargs)

    def __enter__(self):
        self.lock.__enter__()

    def __exit__(self, *args):
        self.lock.__exit__(*args)

    def locked(self):
        return self.lock.locked()

    def __getstate__(self):
        return self.token

    def __setstate__(self, token):
        self.__init__(token)

    def __str__(self):
        return f"<{self.__class__.__name__}: {self.token}>"

    __repr__ = __str__


# Locks used by multiple backends.
# Neither HDF5 nor the netCDF-C library are thread-safe.
HDF5_LOCK = SerializableLock()
NETCDFC_LOCK = SerializableLock()


_FILE_LOCKS: MutableMapping[Any, threading.Lock] = weakref.WeakValueDictionary()


def _get_threaded_lock(key):
    try:
        lock = _FILE_LOCKS[key]
    except KeyError:
        lock = _FILE_LOCKS[key] = threading.Lock()
    return lock


def _get_multiprocessing_lock(key):
    # TODO: make use of the key -- maybe use locket.py?
    # https://github.com/mwilliamson/locket.py
    del key  # unused
    return multiprocessing.Lock()


def _get_lock_maker(scheduler=None):
    """Returns an appropriate function for creating resource locks.

    Parameters
    ----------
    scheduler : str or None
        Dask scheduler being used.

    See Also
    --------
    dask.utils.get_scheduler_lock
    """

    if scheduler is None or scheduler == "threaded":
        return _get_threaded_lock
    elif scheduler == "multiprocessing":
        return _get_multiprocessing_lock
    elif scheduler == "distributed":
        # Lazy import distributed since it is can add a significant
        # amount of time to import
        try:
            from dask.distributed import Lock as DistributedLock
        except ImportError:
            DistributedLock = None
        return DistributedLock
    else:
        raise KeyError(scheduler)


def get_dask_scheduler(get=None, collection=None) -> str | None:
    """Determine the dask scheduler that is being used.

    None is returned if no dask scheduler is active.

    See Also
    --------
    dask.base.get_scheduler
    """
    try:
        # Fix for bug caused by dask installation that doesn't involve the toolz library
        # Issue: 4164
        import dask
        from dask.base import get_scheduler

        actual_get = get_scheduler(get, collection)
    except ImportError:
        return None

    try:
        from dask.distributed import Client

        if isinstance(actual_get.__self__, Client):
            return "distributed"
    except (ImportError, AttributeError):
        pass

    try:
        # As of dask=2.6, dask.multiprocessing requires cloudpickle to be installed
        # Dependency removed in https://github.com/dask/dask/pull/5511
        if actual_get is dask.multiprocessing.get:
            return "multiprocessing"
    except AttributeError:
        pass

    return "threaded"


def get_write_lock(key):
    """Get a scheduler appropriate lock for writing to the given resource.

    Parameters
    ----------
    key : str
        Name of the resource for which to acquire a lock. Typically a filename.

    Returns
    -------
    Lock object that can be used like a threading.Lock object.
    """
    scheduler = get_dask_scheduler()
    lock_maker = _get_lock_maker(scheduler)
    return lock_maker(key)


def acquire(lock, blocking=True):
    """Acquire a lock, possibly in a non-blocking fashion.

    Includes backwards compatibility hacks for old versions of Python, dask
    and dask-distributed.
    """
    if blocking:
        # no arguments needed
        return lock.acquire()
    else:
        # "blocking" keyword argument not supported for:
        # - threading.Lock on Python 2.
        # - dask.SerializableLock with dask v1.0.0 or earlier.
        # - multiprocessing.Lock calls the argument "block" instead.
        # - dask.distributed.Lock uses the blocking argument as the first one
        return lock.acquire(blocking)


class CombinedLock:
    """A combination of multiple locks.

    Like a locked door, a CombinedLock is locked if any of its constituent
    locks are locked.
    """

    def __init__(self, locks):
        self.locks = tuple(set(locks))  # remove duplicates

    def acquire(self, blocking=True):
        return all(acquire(lock, blocking=blocking) for lock in self.locks)

    def release(self):
        for lock in self.locks:
            lock.release()

    def __enter__(self):
        for lock in self.locks:
            lock.__enter__()

    def __exit__(self, *args):
        for lock in self.locks:
            lock.__exit__(*args)

    def locked(self):
        return any(lock.locked for lock in self.locks)

    def __repr__(self):
        return f"CombinedLock({list(self.locks)!r})"


class DummyLock:
    """DummyLock provides the lock API without any actual locking."""

    def acquire(self, blocking=True):
        pass

    def release(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *args):
        pass

    def locked(self):
        return False


def combine_locks(locks):
    """Combine a sequence of locks into a single lock."""
    all_locks = []
    for lock in locks:
        if isinstance(lock, CombinedLock):
            all_locks.extend(lock.locks)
        elif lock is not None:
            all_locks.append(lock)

    num_locks = len(all_locks)
    if num_locks > 1:
        return CombinedLock(all_locks)
    elif num_locks == 1:
        return all_locks[0]
    else:
        return DummyLock()


def ensure_lock(lock):
    """Ensure that the given object is a lock."""
    if lock is None or lock is False:
        return DummyLock()
    return lock
