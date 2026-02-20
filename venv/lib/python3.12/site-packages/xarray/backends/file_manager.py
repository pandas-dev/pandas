from __future__ import annotations

import atexit
import threading
import uuid
import warnings
from collections.abc import Callable, Hashable, Iterator, Mapping, MutableMapping
from contextlib import AbstractContextManager, contextmanager
from typing import Any, Generic, Literal, TypeVar, cast

from xarray.backends.locks import acquire
from xarray.backends.lru_cache import LRUCache
from xarray.core import utils
from xarray.core.options import OPTIONS
from xarray.core.types import Closable, Lock

# Global cache for storing open files.
FILE_CACHE: LRUCache[Any, Closable] = LRUCache(
    maxsize=OPTIONS["file_cache_maxsize"], on_evict=lambda k, v: v.close()
)
assert FILE_CACHE.maxsize, "file cache must be at least size one"

T_File = TypeVar("T_File", bound=Closable)

REF_COUNTS: dict[Any, int] = {}

_OMIT_MODE = utils.ReprObject("<omitted>")


class FileManager(Generic[T_File]):
    """Manager for acquiring and closing a file object.

    Use FileManager subclasses (CachingFileManager in particular) on backend
    storage classes to automatically handle issues related to keeping track of
    many open files and transferring them between multiple processes.
    """

    def acquire(self, needs_lock: bool = True) -> T_File:
        """Acquire the file object from this manager."""
        raise NotImplementedError()

    def acquire_context(
        self, needs_lock: bool = True
    ) -> AbstractContextManager[T_File]:
        """Context manager for acquiring a file. Yields a file object.

        The context manager unwinds any actions taken as part of acquisition
        (i.e., removes it from any cache) if an exception is raised from the
        context. It *does not* automatically close the file.
        """
        raise NotImplementedError()

    def close(self, needs_lock: bool = True) -> None:
        """Close the file object associated with this manager, if needed."""
        raise NotImplementedError()


class CachingFileManager(FileManager[T_File]):
    """Wrapper for automatically opening and closing file objects.

    Unlike files, CachingFileManager objects can be safely pickled and passed
    between processes. They should be explicitly closed to release resources,
    but a per-process least-recently-used cache for open files ensures that you
    can safely create arbitrarily large numbers of FileManager objects.

    Don't directly close files acquired from a FileManager. Instead, call
    FileManager.close(), which ensures that closed files are removed from the
    cache as well.

    Example usage::

        manager = FileManager(open, "example.txt", mode="w")
        f = manager.acquire()
        f.write(...)
        manager.close()  # ensures file is closed

    Note that as long as previous files are still cached, acquiring a file
    multiple times from the same FileManager is essentially free::

        f1 = manager.acquire()
        f2 = manager.acquire()
        assert f1 is f2

    """

    def __init__(
        self,
        opener: Callable[..., T_File],
        *args: Any,
        mode: Any = _OMIT_MODE,
        kwargs: Mapping[str, Any] | None = None,
        lock: Lock | None | Literal[False] = None,
        cache: MutableMapping[Any, T_File] | None = None,
        manager_id: Hashable | None = None,
        ref_counts: dict[Any, int] | None = None,
    ):
        """Initialize a CachingFileManager.

        The cache, manager_id and ref_counts arguments exist solely to
        facilitate dependency injection, and should only be set for tests.

        Parameters
        ----------
        opener : callable
            Function that when called like ``opener(*args, **kwargs)`` returns
            an open file object. The file object must implement a ``close()``
            method.
        *args
            Positional arguments for opener. A ``mode`` argument should be
            provided as a keyword argument (see below). All arguments must be
            hashable.
        mode : optional
            If provided, passed as a keyword argument to ``opener`` along with
            ``**kwargs``. ``mode='w' `` has special treatment: after the first
            call it is replaced by ``mode='a'`` in all subsequent function to
            avoid overriding the newly created file.
        kwargs : dict, optional
            Keyword arguments for opener, excluding ``mode``. All values must
            be hashable.
        lock : duck-compatible threading.Lock, optional
            Lock to use when modifying the cache inside acquire() and close().
            By default, uses a new threading.Lock() object. If set, this object
            should be pickleable.
        cache : MutableMapping, optional
            Mapping to use as a cache for open files. By default, uses xarray's
            global LRU file cache. Because ``cache`` typically points to a
            global variable and contains non-picklable file objects, an
            unpickled FileManager objects will be restored with the default
            cache.
        manager_id : hashable, optional
            Identifier for this CachingFileManager.
        ref_counts : dict, optional
            Optional dict to use for keeping track the number of references to
            the same file.
        """
        self._opener = opener
        self._args = args
        self._mode = mode
        self._kwargs = {} if kwargs is None else dict(kwargs)

        if lock is None or lock is False:
            self._use_default_lock = True
            self._lock: Lock = threading.Lock()
        else:
            self._use_default_lock = False
            self._lock = lock

        # cache[self._key] stores the file associated with this object.
        if cache is None:
            cache = cast(MutableMapping[Any, T_File], FILE_CACHE)
        self._cache: MutableMapping[Any, T_File] = cache
        if manager_id is None:
            # Each call to CachingFileManager should separately open files.
            manager_id = str(uuid.uuid4())
        self._manager_id = manager_id
        self._key = self._make_key()

        # ref_counts[self._key] stores the number of CachingFileManager objects
        # in memory referencing this same file. We use this to know if we can
        # close a file when the manager is deallocated.
        if ref_counts is None:
            ref_counts = REF_COUNTS
        self._ref_counter = _RefCounter(ref_counts)
        self._ref_counter.increment(self._key)

    def _make_key(self) -> _HashedSequence:
        """Make a key for caching files in the LRU cache."""
        value = (
            self._opener,
            self._args,
            "a" if self._mode == "w" else self._mode,
            tuple(sorted(self._kwargs.items())),
            self._manager_id,
        )
        return _HashedSequence(value)

    @contextmanager
    def _optional_lock(self, needs_lock: bool):
        """Context manager for optionally acquiring a lock."""
        if needs_lock:
            with self._lock:
                yield
        else:
            yield

    def acquire(self, needs_lock: bool = True) -> T_File:
        """Acquire a file object from the manager.

        A new file is only opened if it has expired from the
        least-recently-used cache.

        This method uses a lock, which ensures that it is thread-safe. You can
        safely acquire a file in multiple threads at the same time, as long as
        the underlying file object is thread-safe.

        Returns
        -------
        file-like
            An open file object, as returned by ``opener(*args, **kwargs)``.
        """
        file, _ = self._acquire_with_cache_info(needs_lock)
        return file

    @contextmanager
    def acquire_context(self, needs_lock: bool = True) -> Iterator[T_File]:
        """Context manager for acquiring a file."""
        file, cached = self._acquire_with_cache_info(needs_lock)
        try:
            yield file
        except Exception:
            if not cached:
                self.close(needs_lock)
            raise

    def _acquire_with_cache_info(self, needs_lock: bool = True) -> tuple[T_File, bool]:
        """Acquire a file, returning the file and whether it was cached."""
        with self._optional_lock(needs_lock):
            try:
                file = self._cache[self._key]
            except KeyError:
                kwargs = self._kwargs
                if self._mode is not _OMIT_MODE:
                    kwargs = kwargs.copy()
                    kwargs["mode"] = self._mode
                file = self._opener(*self._args, **kwargs)
                if self._mode == "w":
                    # ensure file doesn't get overridden when opened again
                    self._mode = "a"
                self._cache[self._key] = file
                return file, False
            else:
                return file, True

    def close(self, needs_lock: bool = True) -> None:
        """Explicitly close any associated file object (if necessary)."""
        # TODO: remove needs_lock if/when we have a reentrant lock in
        # dask.distributed: https://github.com/dask/dask/issues/3832
        with self._optional_lock(needs_lock):
            default = None
            file = self._cache.pop(self._key, default)
            if file is not None:
                file.close()

    def __del__(self) -> None:
        # If we're the only CachingFileManger referencing an unclosed file,
        # remove it from the cache upon garbage collection.
        #
        # We keep track of our own reference count because we don't want to
        # close files if another identical file manager needs it. This can
        # happen if a CachingFileManager is pickled and unpickled without
        # closing the original file.
        ref_count = self._ref_counter.decrement(self._key)

        if not ref_count and self._key in self._cache:
            if acquire(self._lock, blocking=False):
                # Only close files if we can do so immediately.
                try:
                    self.close(needs_lock=False)
                finally:
                    self._lock.release()

            if OPTIONS["warn_for_unclosed_files"]:
                warnings.warn(
                    f"deallocating {self}, but file is not already closed. "
                    "This may indicate a bug.",
                    RuntimeWarning,
                    stacklevel=2,
                )

    def __getstate__(self):
        """State for pickling."""
        # cache is intentionally omitted: we don't want to try to serialize
        # these global objects.
        lock = None if self._use_default_lock else self._lock
        return (
            self._opener,
            self._args,
            self._mode,
            self._kwargs,
            lock,
            self._manager_id,
        )

    def __setstate__(self, state) -> None:
        """Restore from a pickle."""
        opener, args, mode, kwargs, lock, manager_id = state
        self.__init__(  # type: ignore[misc]
            opener, *args, mode=mode, kwargs=kwargs, lock=lock, manager_id=manager_id
        )

    def __repr__(self) -> str:
        args_string = ", ".join(map(repr, self._args))
        if self._mode is not _OMIT_MODE:
            args_string += f", mode={self._mode!r}"
        return (
            f"{type(self).__name__}({self._opener!r}, {args_string}, "
            f"kwargs={self._kwargs}, manager_id={self._manager_id!r})"
        )


class _RefCounter:
    """Class for keeping track of reference counts."""

    def __init__(self, counts):
        self._counts = counts
        self._lock = threading.Lock()

    def increment(self, name):
        with self._lock:
            count = self._counts[name] = self._counts.get(name, 0) + 1
        return count

    def decrement(self, name):
        with self._lock:
            count = self._counts[name] - 1
            if count:
                self._counts[name] = count
            else:
                del self._counts[name]
        return count


class _HashedSequence(list):
    """Speedup repeated look-ups by caching hash values.

    Based on what Python uses internally in functools.lru_cache.

    Python doesn't perform this optimization automatically:
    https://bugs.python.org/issue1462796
    """

    def __init__(self, tuple_value):
        self[:] = tuple_value
        self.hashvalue = hash(tuple_value)

    def __hash__(self) -> int:  # type: ignore[override]
        return self.hashvalue


def _get_none() -> None:
    return None


class PickleableFileManager(FileManager[T_File]):
    """File manager that supports pickling by reopening a file object.

    Use PickleableFileManager for wrapping file-like objects that do not natively
    support pickling (e.g., netCDF4.Dataset and h5netcdf.File) in cases where a
    global cache is not desirable (e.g., for netCDF files opened from bytes in
    memory, or from existing file objects).
    """

    def __init__(
        self,
        opener: Callable[..., T_File],
        *args: Any,
        mode: Any = _OMIT_MODE,
        lock: Lock | None | Literal[False] = None,
        kwargs: Mapping[str, Any] | None = None,
    ):
        kwargs = {} if kwargs is None else dict(kwargs)
        self._opener = opener
        self._args = args
        self._mode = "a" if mode == "w" else mode
        self._kwargs = kwargs
        self._lock = lock

        # Note: No need for locking with PickleableFileManager, because all
        # opening of files happens in the constructor.
        if mode != _OMIT_MODE:
            kwargs = kwargs | {"mode": mode}
        self._file: T_File | None = opener(*args, **kwargs)

    @property
    def _closed(self) -> bool:
        # If opener() raised an error in the constructor, _file may not be set
        return getattr(self, "_file", None) is None

    def _get_unclosed_file(self) -> T_File:
        if self._closed:
            raise RuntimeError("file is closed")
        file = self._file
        assert file is not None
        return file

    def acquire(self, needs_lock: bool = True) -> T_File:
        del needs_lock  # unused
        return self._get_unclosed_file()

    @contextmanager
    def acquire_context(self, needs_lock: bool = True) -> Iterator[T_File]:
        del needs_lock  # unused
        yield self._get_unclosed_file()

    def close(self, needs_lock: bool = True) -> None:
        if not self._closed:
            file = self._get_unclosed_file()
            if needs_lock and self._lock:
                with self._lock:
                    file.close()
            else:
                file.close()
            self._file = None
            # Remove all references to opener arguments, so they can be garbage
            # collected.
            self._args = ()
            self._mode = _OMIT_MODE
            self._kwargs = {}

    def __del__(self) -> None:
        if not self._closed:
            self.close()

            if OPTIONS["warn_for_unclosed_files"]:
                warnings.warn(
                    f"deallocating {self}, but file is not already closed. "
                    "This may indicate a bug.",
                    RuntimeWarning,
                    stacklevel=2,
                )

    def __getstate__(self):
        # file is intentionally omitted: we want to open it again
        opener = _get_none if self._closed else self._opener
        return (opener, self._args, self._mode, self._lock, self._kwargs)

    def __setstate__(self, state) -> None:
        opener, args, mode, lock, kwargs = state
        self.__init__(opener, *args, mode=mode, lock=lock, kwargs=kwargs)  # type: ignore[misc]

    def __repr__(self) -> str:
        if self._closed:
            return f"<closed {type(self).__name__}>"
        args_string = ", ".join(map(repr, self._args))
        if self._mode is not _OMIT_MODE:
            args_string += f", mode={self._mode!r}"
        kwargs = (
            self._kwargs | {"memory": utils.ReprObject("...")}
            if "memory" in self._kwargs
            else self._kwargs
        )
        return f"{type(self).__name__}({self._opener!r}, {args_string}, {kwargs=})"


@atexit.register
def _remove_del_methods():
    # We don't need to close unclosed files at program exit, and may not be able
    # to, because Python is cleaning up imports / globals.
    del CachingFileManager.__del__
    del PickleableFileManager.__del__


class DummyFileManager(FileManager[T_File]):
    """FileManager that simply wraps an open file in the FileManager interface."""

    def __init__(
        self,
        value: T_File,
        *,
        close: Callable[[], None] | None = None,
        lock: Lock | None | Literal[False] = None,
    ):
        if close is None:
            close = value.close
        self._lock = lock
        self._value = value
        self._close = close

    def acquire(self, needs_lock: bool = True) -> T_File:
        del needs_lock  # unused
        return self._value

    @contextmanager
    def acquire_context(self, needs_lock: bool = True) -> Iterator[T_File]:
        del needs_lock  # unused
        yield self._value

    def close(self, needs_lock: bool = True) -> None:
        if needs_lock and self._lock:
            with self._lock:
                self._close()
        else:
            self._close()
