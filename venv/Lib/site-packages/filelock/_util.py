from __future__ import annotations

import os
import secrets
import stat
import sys
from errno import EACCES, EIO, EISDIR
from pathlib import Path
from typing import Final


def write_all(fd: int, data: bytes) -> None:
    """
    Write the whole buffer to *fd*, looping over the short writes ``os.write`` is allowed to make.

    A marker written with a bare ``os.write`` can land partially: a peer reading it mid-write parses a truncated record
    as malformed or as a foreign holder. Looping until the buffer drains keeps the record atomic in the process and
    kernel view. No ``fsync``: filelock needs a complete record, not crash-durable storage.

    :param fd: file descriptor open for writing.
    :param data: bytes to write in full.

    :raises OSError: if a write reports zero progress before the record is complete.

    """
    remaining = memoryview(data)
    while remaining:
        if (written := os.write(fd, remaining)) == 0:
            raise OSError(EIO, "os.write wrote 0 bytes before the record was complete")
        remaining = remaining[written:]


def raise_on_not_writable_file(filename: str) -> None:
    """
    Raise an exception if attempting to open the file for writing would fail.

    Separates files that can never be written from files that are writable but currently locked.

    :param filename: file to check

    :raises OSError: as if the file was opened for writing.

    """
    try:
        # lstat, not stat: settles exists-and-writable in one syscall, and a hostile symlink at the lock path would
        # make stat inspect the link target, letting an attacker turn a contended acquire into a misleading
        # PermissionError / IsADirectoryError and probe that target's attributes. The real open passes O_NOFOLLOW and
        # refuses the symlink anyway.
        file_stat = os.lstat(filename)
    except OSError:
        return  # does not exist, or an error the caller cannot act on

    # No mtime guard: the old `if st_mtime != 0` skip covered NFS/Linux quirks where os.lstat returned an all-zero
    # struct, which it no longer does. Skipping on mtime 0 let a read-only file or a directory at the lock path pass
    # as missing, so acquire() blocked forever on an open that cannot succeed.
    if not (file_stat.st_mode & stat.S_IWUSR):
        raise PermissionError(EACCES, "Permission denied", filename)

    if stat.S_ISDIR(file_stat.st_mode):
        if sys.platform == "win32":  # pragma: win32 cover
            raise PermissionError(EACCES, "Permission denied", filename)
        raise IsADirectoryError(EISDIR, "Is a directory", filename)  # pragma: win32 no cover


def ensure_directory_exists(filename: Path | str) -> None:
    """
    Ensure the directory containing the file exists (create it if necessary).

    :param filename: file.

    """
    Path(filename).parent.mkdir(parents=True, exist_ok=True)


def break_lock_file(lock_file: str, mtime_before: float, ino_before: int) -> None:
    """
    Atomically break a stale lock file judged stale at modification time *mtime_before*.

    Rename the file to a process-private name before unlinking it, so two processes breaking the same lock cannot
    delete each other's work: only one rename of a given inode wins, the loser gets ``OSError``. After the rename,
    re-check the file. A newer modification time, or a different inode than *ino_before*, means a peer recreated the
    lock between the stale decision and the rename, so we grabbed a live file and abort, leaving the renamed file in
    place. A rollback rename is itself racy, the same trade-off as the soft read/write marker break. The inode check
    matters because filesystems with coarse modification-time granularity (NFS, FAT) can give a same-second recreation
    the old mtime, so mtime alone would miss it and unlink a live lock; the inode is the reliable identity, mirroring
    the token re-check in the soft read/write marker break. ``lstat`` avoids following a hostile symlink swapped in
    after the decision.

    The break name carries a random token so it is unguessable and unique per attempt. Without it two breakers in the
    same process share ``<lock>.break.<pid>``, and a second break can rename a recreated live lock onto that path in
    the window between the re-verify ``lstat`` above and the ``unlink`` below, deleting a live lock the inode check
    just approved. A private name keeps anyone else from targeting our break path, matching the soft read/write marker
    break.

    :param lock_file: path to the lock file to break.
    :param mtime_before: modification time observed when the lock was judged stale.
    :param ino_before: inode number observed when the lock was judged stale.

    :raises OSError: if the rename fails (e.g. the file vanished or is not owned in a sticky directory).

    """
    break_path = f"{lock_file}.break.{os.getpid()}.{secrets.token_hex(16)}"
    Path(lock_file).rename(break_path)
    try:
        st_after = os.lstat(break_path)
    except OSError:
        return
    if st_after.st_mtime > mtime_before or st_after.st_ino != ino_before:
        return
    Path(break_path).unlink()


def touch(name: str, *, fd: int | None = None) -> None:
    # Prefer the already-open, already-verified fd so a peer that swaps a symlink or a different file in at the
    # path after our O_NOFOLLOW read cannot redirect the touch: utime then targets the inode behind the fd.
    # Where the platform cannot utime an fd, fall back to a path-based touch that still refuses to follow a
    # symlink where supported, matching the O_NOFOLLOW reads used elsewhere here.
    if fd is not None and _SUPPORTS_UTIME_FD:  # pragma: needs utime-fd
        os.utime(fd, None)
        return
    os.utime(name, None, follow_symlinks=not _SUPPORTS_UTIME_NOFOLLOW)


# Retargeting os.utime to an open fd lets a heartbeat refresh the exact inode it verified instead of whatever the
# pathname now names.
_SUPPORTS_UTIME_FD: Final[bool] = sys.platform != "win32" and os.utime in os.supports_fd
# os.utime follows symlinks unless told not to; not every platform can refuse the follow, so probe support.
_SUPPORTS_UTIME_NOFOLLOW: Final[bool] = os.utime in os.supports_follow_symlinks


__all__ = [
    "break_lock_file",
    "ensure_directory_exists",
    "raise_on_not_writable_file",
    "touch",
    "write_all",
]
