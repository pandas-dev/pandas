"""
Atomic file operations
"""

from __future__ import annotations

__lazy_modules__ = {"atexit", "msvcrt", "plumbum.machines.local", "threading"}

import atexit
import contextlib
import errno
import os
import sys
import threading
import typing

from plumbum.machines.local import local

if typing.TYPE_CHECKING:
    from collections.abc import Generator
    from io import FileIO

    from plumbum._compat.typing import Self

if not sys.platform.startswith("win32"):
    import fcntl

    if hasattr(fcntl, "lockf"):

        @contextlib.contextmanager
        def locked_file(
            fileno: int, blocking: bool = True
        ) -> Generator[None, None, None]:
            fcntl.lockf(fileno, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
            try:
                yield
            finally:
                fcntl.lockf(fileno, fcntl.LOCK_UN)

    else:

        @contextlib.contextmanager
        def locked_file(
            fileno: int, blocking: bool = True
        ) -> Generator[None, None, None]:
            fcntl.flock(fileno, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
            try:
                yield
            finally:
                fcntl.flock(fileno, fcntl.LOCK_UN)
else:
    import ctypes
    import msvcrt
    from ctypes.wintypes import BOOL, DWORD, HANDLE
    from ctypes.wintypes import LPVOID as PVOID
    from ctypes.wintypes import WPARAM as ULONG_PTR

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)

    # Refer: https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-lockfileex
    LOCKFILE_FAIL_IMMEDIATELY = 0x01
    LOCKFILE_EXCLUSIVE_LOCK = 0x02

    # Refer: https://learn.microsoft.com/en-us/windows/win32/api/minwinbase/ns-minwinbase-overlapped
    class OVERLAPPED(ctypes.Structure):
        class DUMMYUNIONNAME(ctypes.Union):
            class DUMMYSTRUCTNAME(ctypes.Structure):
                _fields_ = (
                    ("Offset", DWORD),
                    ("OffsetHigh", DWORD),
                )

            _fields_ = (
                ("_offsets", DUMMYSTRUCTNAME),
                ("Pointer", PVOID),
            )

        _fields_ = (
            ("Internal", ULONG_PTR),
            ("InternalHigh", ULONG_PTR),
            ("_offsets_or_ptr", DUMMYUNIONNAME),
            ("hEvent", HANDLE),
        )

    LPOVERLAPPED = ctypes.POINTER(OVERLAPPED)

    # Refer: https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-lockfile
    LockFile = kernel32.LockFile
    LockFile.restype = BOOL
    LockFile.argtypes = [
        HANDLE,  # hFile
        DWORD,  # dwFileOffsetLow
        DWORD,  # dwFileOffsetHigh
        DWORD,  # nNumberOfBytesToLockLow
        DWORD,  # nNumberOfBytesToLockHigh
    ]

    # Refer: https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-unlockfile
    UnlockFile = kernel32.UnlockFile
    UnlockFile.restype = BOOL
    UnlockFile.argtypes = [
        HANDLE,  # hFile,
        DWORD,  # dwFileOffsetLow,
        DWORD,  # dwFileOffsetHigh,
        DWORD,  # nNumberOfBytesToUnlockLow,
        DWORD,  # nNumberOfBytesToUnlockHigh
    ]

    # Refer: https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-lockfileex
    LockFileEx = kernel32.LockFileEx
    LockFileEx.restype = BOOL
    LockFileEx.argtypes = [
        HANDLE,  # hFile
        DWORD,  # dwFlags
        DWORD,  # dwReserved - must be set to zero
        DWORD,  # nNumberOfBytesToLockLow
        DWORD,  # nNumberOfBytesToLockHigh
        LPOVERLAPPED,  # lpOverlapped
    ]

    # Refer: https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-unlockfileex
    UnlockFileEx = kernel32.UnlockFileEx
    UnlockFileEx.restype = BOOL
    UnlockFileEx.argtypes = [
        HANDLE,  # hFile
        DWORD,  # dwReserved - must be set to zero
        DWORD,  # nNumberOfBytesToUnlockLow
        DWORD,  # nNumberOfBytesToUnlockHigh
        LPOVERLAPPED,  # lpOverlapped
    ]

    @contextlib.contextmanager
    def locked_file(fileno: int, blocking: bool = True) -> Generator[None, None, None]:
        hndl = msvcrt.get_osfhandle(fileno)
        overlapped = OVERLAPPED()
        ok = LockFileEx(
            hndl,
            LOCKFILE_EXCLUSIVE_LOCK | (0 if blocking else LOCKFILE_FAIL_IMMEDIATELY),
            0,
            0xFFFFFFFF,
            0xFFFFFFFF,
            ctypes.byref(overlapped),
        )
        if not ok:
            raise ctypes.WinError(ctypes.get_last_error())
        try:
            yield
        finally:
            exc = sys.exc_info()[1]
            ok = UnlockFile(hndl, 0, 0, 0xFFFFFFFF, 0xFFFFFFFF)
            if not ok:
                next_exc = ctypes.WinError(ctypes.get_last_error())
                if exc is None:
                    exc, next_exc = next_exc, None
                raise exc from next_exc


class AtomicFile:
    """
    Atomic file operations implemented using file-system advisory locks (``flock`` on POSIX,
    ``LockFile`` on Windows).

    .. note::
        On Linux, the manpage says ``flock`` might have issues with NFS mounts. You should
        take this into account.

    .. versionadded:: 1.3
    """

    __slots__ = ("_fileobj", "_ignore_deletion", "_owned_by", "_thdlock", "path")

    CHUNK_SIZE = 32 * 1024

    def __init__(self, filename: str, ignore_deletion: bool = False):
        self.path = local.path(filename)
        self._ignore_deletion = ignore_deletion
        self._thdlock = threading.Lock()
        self._owned_by: int | None = None
        self._fileobj: FileIO | None = None

        self.reopen()

    def __repr__(self) -> str:
        return f"<AtomicFile: {self.path}>" if self._fileobj else "<AtomicFile: closed>"

    def __del__(self) -> None:
        self.close()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, t: object, v: object, tb: object) -> None:
        self.close()

    def close(self) -> None:
        if self._fileobj is not None:
            self._fileobj.close()
            self._fileobj = None

    def reopen(self) -> None:
        """
        Close and reopen the file; useful when the file was deleted from the file system
        by a different process
        """
        self.close()
        self._fileobj = os.fdopen(
            os.open(str(self.path), os.O_CREAT | os.O_RDWR, 0o600), "r+b", 0
        )

    @contextlib.contextmanager
    def locked(self, blocking: bool = True) -> Generator[None, None, None]:
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

        assert self._fileobj is not None
        # Honor ``blocking`` for the in-process thread lock too: a plain
        # ``with self._thdlock`` would block a non-blocking caller when another
        # thread of this process holds the lock. Mirror the filesystem-lock
        # failure path by raising an ``OSError`` that "would block".
        if not self._thdlock.acquire(blocking=blocking):
            raise OSError(errno.EAGAIN, "Atomic file lock is held by another thread")
        try:
            with locked_file(self._fileobj.fileno(), blocking):
                if not self.path.exists() and not self._ignore_deletion:
                    raise ValueError("Atomic file removed from filesystem")
                self._owned_by = threading.get_ident()
                try:
                    yield
                finally:
                    self._owned_by = None
        finally:
            self._thdlock.release()

    def delete(self) -> None:
        """
        Atomically delete the file (holds the lock while doing it)
        """
        with self.locked():
            self.path.delete()

    def _read_all(self) -> bytes:
        assert self._fileobj is not None
        self._fileobj.seek(0)
        data = []
        while True:
            buf = self._fileobj.read(self.CHUNK_SIZE)
            data.append(buf)
            if len(buf) < self.CHUNK_SIZE:
                break
        return b"".join(data)

    def read_atomic(self) -> bytes:
        """Atomically read the entire file"""
        with self.locked():
            return self._read_all()

    def read_shared(self) -> bytes:
        """Read the file **without** holding the lock"""
        return self._read_all()

    def write_atomic(self, data: bytes) -> None:
        """Writes the given data atomically to the file. Note that it overwrites the entire file;
        ``write_atomic("foo")`` followed by ``write_atomic("bar")`` will result in only ``"bar"``.
        """
        with self.locked():
            assert self._fileobj is not None
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

    __slots__ = ("atomicfile", "initial")

    def __init__(self, atomicfile: AtomicFile, initial: int = 0):
        """
        :param atomicfile: an :class:`AtomicFile <plumbum.fs.atomic.AtomicFile>` instance
        :param initial: the initial value (used when the first time the file is created)
        """
        self.atomicfile = atomicfile
        self.initial = initial

    def __enter__(self) -> Self:
        return self

    def __exit__(self, t: object, v: object, tb: object) -> None:
        self.close()

    def close(self) -> None:
        self.atomicfile.close()

    @classmethod
    def open(cls, filename: str) -> Self:
        """
        Shortcut for ``AtomicCounterFile(AtomicFile(filename))``
        """
        return cls(AtomicFile(filename))

    def reset(self, value: int | None = None) -> None:
        """
        Reset the counter's value to the one given. If ``None``, it will default to the
        initial value provided to the constructor
        """
        if value is None:
            value = self.initial
        if not isinstance(value, int):
            raise TypeError(f"value must be an integer, not {type(value)!r}")
        self.atomicfile.write_atomic(str(value).encode("utf8"))

    def next(self) -> int:
        """
        Read and increment the counter, returning its previous value
        """
        with self.atomicfile.locked():
            curr_str = self.atomicfile.read_atomic().decode("utf8")
            curr = self.initial if not curr_str else int(curr_str)
            self.atomicfile.write_atomic(str(curr + 1).encode("utf8"))
            return curr


class PidFileTaken(SystemExit):
    """
    This exception is raised when PidFile.acquire fails to lock the pid file. Note that it
    derives from ``SystemExit``, so unless explicitly handled, it will terminate the process
    cleanly.
    """

    def __init__(self, msg: str, pid: str):
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

    __slots__ = ("_ctx", "atomicfile")

    def __init__(self, filename: str):
        self.atomicfile = AtomicFile(filename)
        self._ctx: contextlib.AbstractContextManager[None] | None = None

    def __enter__(self) -> None:
        self.acquire()

    def __exit__(self, t: object, v: object, tb: object) -> None:
        self.release()

    def __del__(self) -> None:
        with contextlib.suppress(Exception):
            self.release()

    def close(self) -> None:
        self.atomicfile.close()

    def acquire(self) -> None:
        """
        Attempt to acquire the PID file. If it's already locked, raises
        :class:`PidFileTaken <plumbum.fs.atomic.PidFileTaken>`. You should normally acquire
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

    def release(self) -> None:
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


__all__ = [
    "AtomicCounterFile",
    "AtomicFile",
    "PidFile",
    "PidFileTaken",
]


def __dir__() -> list[str]:
    return list(__all__)
