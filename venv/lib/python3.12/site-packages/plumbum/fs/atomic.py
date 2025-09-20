"""
Atomic file operations
"""

from __future__ import annotations

import atexit
import contextlib
import os
import threading

from plumbum.machines.local import local

try:
    import fcntl
except ImportError:
    import msvcrt

    try:
        from pywintypes import error as WinError
        from win32con import LOCKFILE_EXCLUSIVE_LOCK, LOCKFILE_FAIL_IMMEDIATELY
        from win32file import OVERLAPPED, LockFileEx, UnlockFile
    except ImportError:
        print(  # noqa: T201
            "On Windows, Plumbum requires Python for Windows Extensions (pywin32)"
        )
        raise

    @contextlib.contextmanager
    def locked_file(fileno, blocking=True):
        hndl = msvcrt.get_osfhandle(fileno)
        try:
            LockFileEx(
                hndl,
                LOCKFILE_EXCLUSIVE_LOCK
                | (0 if blocking else LOCKFILE_FAIL_IMMEDIATELY),
                0xFFFFFFFF,
                0xFFFFFFFF,
                OVERLAPPED(),
            )
        except WinError as ex:
            raise OSError(*ex.args) from None
        try:
            yield
        finally:
            UnlockFile(hndl, 0, 0, 0xFFFFFFFF, 0xFFFFFFFF)

else:
    if hasattr(fcntl, "lockf"):

        @contextlib.contextmanager
        def locked_file(fileno, blocking=True):
            fcntl.lockf(fileno, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
            try:
                yield
            finally:
                fcntl.lockf(fileno, fcntl.LOCK_UN)

    else:

        @contextlib.contextmanager
        def locked_file(fileno, blocking=True):
            fcntl.flock(fileno, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
            try:
                yield
            finally:
                fcntl.flock(fileno, fcntl.LOCK_UN)


class AtomicFile:
    """
    Atomic file operations implemented using file-system advisory locks (``flock`` on POSIX,
    ``LockFile`` on Windows).

    .. note::
        On Linux, the manpage says ``flock`` might have issues with NFS mounts. You should
        take this into account.

    .. versionadded:: 1.3
    """

    CHUNK_SIZE = 32 * 1024

    def __init__(self, filename, ignore_deletion=False):
        self.path = local.path(filename)
        self._ignore_deletion = ignore_deletion
        self._thdlock = threading.Lock()
        self._owned_by = None
        self._fileobj = None
        self.reopen()

    def __repr__(self):
        return f"<AtomicFile: {self.path}>" if self._fileobj else "<AtomicFile: closed>"

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, t, v, tb):
        self.close()

    def close(self):
        if self._fileobj is not None:
            self._fileobj.close()
            self._fileobj = None

    def reopen(self):
        """
        Close and reopen the file; useful when the file was deleted from the file system
        by a different process
        """
        self.close()
        self._fileobj = os.fdopen(
            os.open(str(self.path), os.O_CREAT | os.O_RDWR, 384), "r+b", 0
        )

    @contextlib.contextmanager
    def locked(self, blocking=True):
        """
        A context manager that locks the file; this function is reentrant by the thread currently
        holding the lock.

        :param blocking: if ``True``, the call will block until we can grab the file system lock.
                         if ``False``, the call may fail immediately with the underlying exception
                         (``IOError`` or ``WindowsError``)
        """
        if self._owned_by == threading.get_ident():
            yield
            return
        with self._thdlock, locked_file(self._fileobj.fileno(), blocking):
            if not self.path.exists() and not self._ignore_deletion:
                raise ValueError("Atomic file removed from filesystem")
            self._owned_by = threading.get_ident()
            try:
                yield
            finally:
                self._owned_by = None

    def delete(self):
        """
        Atomically delete the file (holds the lock while doing it)
        """
        with self.locked():
            self.path.delete()

    def _read_all(self):
        self._fileobj.seek(0)
        data = []
        while True:
            buf = self._fileobj.read(self.CHUNK_SIZE)
            data.append(buf)
            if len(buf) < self.CHUNK_SIZE:
                break
        return b"".join(data)

    def read_atomic(self):
        """Atomically read the entire file"""
        with self.locked():
            return self._read_all()

    def read_shared(self):
        """Read the file **without** holding the lock"""
        return self._read_all()

    def write_atomic(self, data):
        """Writes the given data atomically to the file. Note that it overwrites the entire file;
        ``write_atomic("foo")`` followed by ``write_atomic("bar")`` will result in only ``"bar"``.
        """
        with self.locked():
            self._fileobj.seek(0)
            while data:
                chunk = data[: self.CHUNK_SIZE]
                self._fileobj.write(chunk)
                data = data[len(chunk) :]
            self._fileobj.flush()
            self._fileobj.truncate()


class AtomicCounterFile:
    """
    An atomic counter based on AtomicFile. Each time you call ``next()``, it will
    atomically read and increment the counter's value, returning its previous value

    Example::

        acf = AtomicCounterFile.open("/some/file")
        print(acf.next())  # e.g., 7
        print(acf.next())  # 8
        print(acf.next())  # 9

    .. versionadded:: 1.3
    """

    def __init__(self, atomicfile, initial=0):
        """
        :param atomicfile: an :class:`AtomicFile <plumbum.atomic.AtomicFile>` instance
        :param initial: the initial value (used when the first time the file is created)
        """
        self.atomicfile = atomicfile
        self.initial = initial

    def __enter__(self):
        return self

    def __exit__(self, t, v, tb):
        self.close()

    def close(self):
        self.atomicfile.close()

    @classmethod
    def open(cls, filename):
        """
        Shortcut for ``AtomicCounterFile(AtomicFile(filename))``
        """
        return cls(AtomicFile(filename))

    def reset(self, value=None):
        """
        Reset the counter's value to the one given. If ``None``, it will default to the
        initial value provided to the constructor
        """
        if value is None:
            value = self.initial
        if not isinstance(value, int):
            raise TypeError(f"value must be an integer, not {type(value)!r}")
        self.atomicfile.write_atomic(str(value).encode("utf8"))

    def next(self):
        """
        Read and increment the counter, returning its previous value
        """
        with self.atomicfile.locked():
            curr = self.atomicfile.read_atomic().decode("utf8")
            curr = self.initial if not curr else int(curr)
            self.atomicfile.write_atomic(str(curr + 1).encode("utf8"))
            return curr


class PidFileTaken(SystemExit):
    """
    This exception is raised when PidFile.acquire fails to lock the pid file. Note that it
    derives from ``SystemExit``, so unless explicitly handled, it will terminate the process
    cleanly
    """

    def __init__(self, msg, pid):
        SystemExit.__init__(self, msg)
        self.pid = pid


class PidFile:
    """
    A PID file is a file that's locked by some process from the moment it starts until it dies
    (the OS will clear the lock when the process exits). It is used to prevent two instances
    of the same process (normally a daemon) from running concurrently. The PID file holds its
    process' PID, so you know who's holding it.

    .. versionadded:: 1.3
    """

    def __init__(self, filename):
        self.atomicfile = AtomicFile(filename)
        self._ctx = None

    def __enter__(self):
        self.acquire()

    def __exit__(self, t, v, tb):
        self.release()

    def __del__(self):
        with contextlib.suppress(Exception):
            self.release()

    def close(self):
        self.atomicfile.close()

    def acquire(self):
        """
        Attempt to acquire the PID file. If it's already locked, raises
        :class:`PidFileTaken <plumbum.atomic.PidFileTaken>`. You should normally acquire
        the file as early as possible when the program starts
        """
        if self._ctx is not None:
            return
        self._ctx = self.atomicfile.locked(blocking=False)
        try:
            self._ctx.__enter__()  # pylint: disable=unnecessary-dunder-call
        except OSError:
            self._ctx = None
            try:
                pid = self.atomicfile.read_shared().strip().decode("utf8")
            except OSError:
                pid = "Unknown"
            raise PidFileTaken(
                f"PID file {self.atomicfile.path!r} taken by process {pid}",
                pid,
            ) from None
        self.atomicfile.write_atomic(str(os.getpid()).encode("utf8"))
        atexit.register(self.release)

    def release(self):
        """
        Release the PID file (should only happen when the program terminates)
        """
        if self._ctx is None:
            return
        self.atomicfile.delete()
        try:
            self._ctx.__exit__(None, None, None)
        finally:
            self._ctx = None
