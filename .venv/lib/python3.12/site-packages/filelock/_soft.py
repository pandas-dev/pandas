from __future__ import annotations

import os
import socket
import sys
import time
from contextlib import suppress
from errno import EACCES, EEXIST, EPERM, ESRCH
from pathlib import Path

from ._api import BaseFileLock
from ._util import ensure_directory_exists, raise_on_not_writable_file

_WIN_SYNCHRONIZE = 0x100000
_WIN_ERROR_INVALID_PARAMETER = 87


class SoftFileLock(BaseFileLock):
    """
    Portable file lock based on file existence.

    Unlike :class:`UnixFileLock <filelock.UnixFileLock>` and :class:`WindowsFileLock <filelock.WindowsFileLock>`, this
    lock does not use OS-level locking primitives. Instead, it creates the lock file with ``O_CREAT | O_EXCL`` and
    treats its existence as the lock indicator. This makes it work on any filesystem but leaves stale lock files behind
    if the process crashes without releasing the lock.

    To mitigate stale locks, the lock file contains the PID and hostname of the holding process. On contention, if the
    holder is on the same host and its PID no longer exists, the stale lock is broken automatically.

    """

    def _acquire(self) -> None:
        raise_on_not_writable_file(self.lock_file)
        ensure_directory_exists(self.lock_file)
        flags = (
            os.O_WRONLY  # open for writing only
            | os.O_CREAT
            | os.O_EXCL  # together with above raise EEXIST if the file specified by filename exists
            | os.O_TRUNC  # truncate the file to zero byte
        )
        if (o_nofollow := getattr(os, "O_NOFOLLOW", None)) is not None:
            flags |= o_nofollow
        try:
            file_handler = os.open(self.lock_file, flags, self._open_mode())
        except OSError as exception:
            if not (
                exception.errno == EEXIST or (exception.errno == EACCES and sys.platform == "win32")
            ):  # pragma: win32 no cover
                raise
            if exception.errno == EEXIST and sys.platform != "win32":  # pragma: win32 no cover
                self._try_break_stale_lock()
        else:
            self._write_lock_info(file_handler)
            self._context.lock_file_fd = file_handler

    def _try_break_stale_lock(self) -> None:
        with suppress(OSError):
            content = Path(self.lock_file).read_text(encoding="utf-8")
            lines = content.strip().splitlines()
            if len(lines) != 2:  # noqa: PLR2004
                return
            pid_str, hostname = lines
            if hostname != socket.gethostname():
                return
            pid = int(pid_str)
            if self._is_process_alive(pid):
                return
            break_path = f"{self.lock_file}.break.{os.getpid()}"
            Path(self.lock_file).rename(break_path)
            Path(break_path).unlink()

    @staticmethod
    def _is_process_alive(pid: int) -> bool:
        if sys.platform == "win32":  # pragma: win32 cover
            import ctypes  # noqa: PLC0415

            kernel32 = ctypes.windll.kernel32
            handle = kernel32.OpenProcess(_WIN_SYNCHRONIZE, 0, pid)
            if handle:
                kernel32.CloseHandle(handle)
                return True
            return kernel32.GetLastError() != _WIN_ERROR_INVALID_PARAMETER
        try:
            os.kill(pid, 0)
        except OSError as exc:
            if exc.errno == ESRCH:
                return False
            if exc.errno == EPERM:
                return True
            raise
        return True

    @staticmethod
    def _write_lock_info(fd: int) -> None:
        with suppress(OSError):
            os.write(fd, f"{os.getpid()}\n{socket.gethostname()}\n".encode())

    def _release(self) -> None:
        assert self._context.lock_file_fd is not None  # noqa: S101
        os.close(self._context.lock_file_fd)
        self._context.lock_file_fd = None
        if sys.platform == "win32":
            self._windows_unlink_with_retry()
        else:
            with suppress(OSError):
                Path(self.lock_file).unlink()

    def _windows_unlink_with_retry(self) -> None:
        max_retries = 10
        retry_delay = 0.001
        for attempt in range(max_retries):
            # Windows doesn't immediately release file handles after close, causing EACCES/EPERM on unlink
            try:
                Path(self.lock_file).unlink()
            except OSError as exc:  # noqa: PERF203
                if exc.errno not in {EACCES, EPERM}:
                    return
                if attempt < max_retries - 1:
                    time.sleep(retry_delay)
                    retry_delay *= 2
            else:
                return


__all__ = [
    "SoftFileLock",
]
