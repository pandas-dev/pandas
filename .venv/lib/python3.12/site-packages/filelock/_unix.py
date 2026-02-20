from __future__ import annotations

import os
import sys
import warnings
from contextlib import suppress
from errno import EAGAIN, ENOSYS, EWOULDBLOCK
from pathlib import Path
from typing import cast

from ._api import BaseFileLock
from ._util import ensure_directory_exists

#: a flag to indicate if the fcntl API is available
has_fcntl = False
if sys.platform == "win32":  # pragma: win32 cover

    class UnixFileLock(BaseFileLock):
        """Uses the :func:`fcntl.flock` to hard lock the lock file on unix systems."""

        def _acquire(self) -> None:
            raise NotImplementedError

        def _release(self) -> None:
            raise NotImplementedError

else:  # pragma: win32 no cover
    try:
        import fcntl

        _ = (fcntl.flock, fcntl.LOCK_EX, fcntl.LOCK_NB, fcntl.LOCK_UN)
    except (ImportError, AttributeError):
        pass
    else:
        has_fcntl = True

    class UnixFileLock(BaseFileLock):
        """Uses the :func:`fcntl.flock` to hard lock the lock file on unix systems."""

        def _acquire(self) -> None:  # noqa: C901, PLR0912
            ensure_directory_exists(self.lock_file)
            open_flags = os.O_RDWR | os.O_TRUNC
            o_nofollow = getattr(os, "O_NOFOLLOW", None)
            if o_nofollow is not None:
                open_flags |= o_nofollow
            open_flags |= os.O_CREAT
            open_mode = self._open_mode()
            try:
                fd = os.open(self.lock_file, open_flags, open_mode)
            except FileNotFoundError:
                # On FUSE/NFS, os.open(O_CREAT) is not atomic: LOOKUP + CREATE can be split, allowing a concurrent
                # unlink() to delete the file between them. For valid paths, treat ENOENT as transient contention.
                # For invalid paths (e.g., empty string), re-raise to avoid infinite retry loops.
                if self.lock_file and Path(self.lock_file).parent.exists():
                    return
                raise
            except PermissionError:
                # Sticky-bit dirs (e.g. /tmp): O_CREAT fails if the file is owned by another user (#317).
                # Fall back to opening the existing file without O_CREAT.
                if not Path(self.lock_file).exists():
                    raise
                try:
                    fd = os.open(self.lock_file, open_flags & ~os.O_CREAT, open_mode)
                except FileNotFoundError:
                    return
            if self.has_explicit_mode:
                with suppress(PermissionError):
                    os.fchmod(fd, self._context.mode)
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as exception:
                os.close(fd)
                if exception.errno == ENOSYS:
                    with suppress(OSError):
                        Path(self.lock_file).unlink()
                    self._fallback_to_soft_lock()
                    self._acquire()
                    return
                if exception.errno not in {EAGAIN, EWOULDBLOCK}:
                    raise
            else:
                # The file may have been unlinked by a concurrent _release() between our open() and flock().
                # A lock on an unlinked inode is useless — discard and let the retry loop start fresh.
                if os.fstat(fd).st_nlink == 0:
                    os.close(fd)
                else:
                    self._context.lock_file_fd = fd

        def _fallback_to_soft_lock(self) -> None:
            from ._soft import SoftFileLock  # noqa: PLC0415

            warnings.warn("flock not supported on this filesystem, falling back to SoftFileLock", stacklevel=2)
            from .asyncio import AsyncSoftFileLock, BaseAsyncFileLock  # noqa: PLC0415

            self.__class__ = AsyncSoftFileLock if isinstance(self, BaseAsyncFileLock) else SoftFileLock

        def _release(self) -> None:
            fd = cast("int", self._context.lock_file_fd)
            self._context.lock_file_fd = None
            with suppress(OSError):
                Path(self.lock_file).unlink()
            fcntl.flock(fd, fcntl.LOCK_UN)
            os.close(fd)


__all__ = [
    "UnixFileLock",
    "has_fcntl",
]
