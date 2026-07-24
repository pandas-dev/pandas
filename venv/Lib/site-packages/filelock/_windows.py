from __future__ import annotations

import os
import sys
from contextlib import suppress
from pathlib import Path
from typing import Final, cast

from ._api import BaseFileLock
from ._util import ensure_directory_exists, raise_on_not_writable_file

if sys.platform == "win32":  # pragma: win32 cover
    import ctypes
    import msvcrt
    from ctypes import wintypes

    _GENERIC_READ: Final[int] = 0x80000000
    _GENERIC_WRITE: Final[int] = 0x40000000
    _SYNCHRONIZE: Final[int] = 0x00100000
    _DESIRED_ACCESS: Final[int] = _GENERIC_READ | _GENERIC_WRITE | _SYNCHRONIZE
    _FILE_SHARE_READ_WRITE: Final[int] = (
        0x00000001 | 0x00000002
    )  # read | write; matches os.open (_SH_DENYNO), no delete
    _FILE_OPEN_IF: Final[int] = 3  # open the file if it exists, create it otherwise; the NtCreateFile OPEN_ALWAYS
    _FILE_ATTRIBUTE_READONLY: Final[int] = 0x00000001
    _FILE_ATTRIBUTE_NORMAL: Final[int] = 0x00000080
    _FILE_ATTRIBUTE_REPARSE_POINT: Final[int] = 0x00000400
    # CreateOptions: keep the handle synchronous (the CRT and msvcrt.locking rely on a maintained file position),
    # refuse a directory, and open a reparse point rather than following it so the check below acts on the link itself.
    _FILE_SYNCHRONOUS_IO_NONALERT: Final[int] = 0x00000020
    _FILE_NON_DIRECTORY_FILE: Final[int] = 0x00000040
    _FILE_OPEN_REPARSE_POINT: Final[int] = 0x00200000
    _CREATE_OPTIONS: Final[int] = _FILE_SYNCHRONOUS_IO_NONALERT | _FILE_NON_DIRECTORY_FILE | _FILE_OPEN_REPARSE_POINT
    _OBJ_CASE_INSENSITIVE: Final[int] = 0x00000040  # Win32 name lookups are case-insensitive
    _OWNER_WRITE: Final[int] = 0o200

    # LockFileEx locks a byte range at an offset carried in OVERLAPPED, independent of the descriptor's file position.
    # msvcrt.locking starts at the current position instead, so a metadata write between lock and unlock could shift
    # the byte a later unlock targets; the explicit offset removes that hazard for both the path lock and #608's
    # descriptor lock.
    _LOCKFILE_FAIL_IMMEDIATELY: Final[int] = 0x00000001
    _LOCKFILE_EXCLUSIVE_LOCK: Final[int] = 0x00000002
    _ERROR_LOCK_VIOLATION: Final[int] = 33  # another handle holds the byte range

    # NtCreateFile returns the raw NTSTATUS as its value, where CreateFileW collapses several of these into one
    # ERROR_ACCESS_DENIED. Telling them apart is the point (#604): a name pending deletion or a share conflict is
    # transient and worth a retry, a real access denial is not.
    _STATUS_SUCCESS: Final[int] = 0x00000000
    _STATUS_ACCESS_DENIED: Final[int] = 0xC0000022
    _STATUS_SHARING_VIOLATION: Final[int] = 0xC0000043
    _STATUS_DELETE_PENDING: Final[int] = 0xC0000056

    _ntdll: Final[ctypes.WinDLL] = ctypes.WinDLL("ntdll")
    _kernel32: Final[ctypes.WinDLL] = ctypes.WinDLL("kernel32", use_last_error=True)

    class _UNICODE_STRING(ctypes.Structure):  # ruff:ignore[invalid-class-name]  # mirrors the Win32 struct name
        _fields_ = (
            ("Length", wintypes.USHORT),  # byte length, not character count
            ("MaximumLength", wintypes.USHORT),
            ("Buffer", wintypes.LPWSTR),
        )

    class _OBJECT_ATTRIBUTES(ctypes.Structure):  # ruff:ignore[invalid-class-name]  # mirrors the Win32 struct name
        _fields_ = (
            ("Length", wintypes.ULONG),
            ("RootDirectory", wintypes.HANDLE),
            ("ObjectName", ctypes.POINTER(_UNICODE_STRING)),
            ("Attributes", wintypes.ULONG),
            ("SecurityDescriptor", ctypes.c_void_p),
            ("SecurityQualityOfService", ctypes.c_void_p),
        )

    class _IO_STATUS_BLOCK(ctypes.Structure):  # ruff:ignore[invalid-class-name]  # mirrors the Win32 struct name
        _fields_ = (
            ("Status", ctypes.c_void_p),  # a union of NTSTATUS and PVOID, so it is pointer-sized
            ("Information", ctypes.c_void_p),
        )

    class _OVERLAPPED(ctypes.Structure):  # mirrors the Win32 struct name
        _fields_ = (
            ("Internal", ctypes.c_void_p),  # ULONG_PTR: pointer-sized, not DWORD, or the x64 layout corrupts Offset
            ("InternalHigh", ctypes.c_void_p),
            ("Offset", wintypes.DWORD),  # the DUMMYUNIONNAME struct, flattened: low 32 bits of the byte offset
            ("OffsetHigh", wintypes.DWORD),
            ("hEvent", wintypes.HANDLE),
        )

    class _BY_HANDLE_FILE_INFORMATION(ctypes.Structure):  # ruff:ignore[invalid-class-name]  # mirrors the Win32 struct name
        _fields_ = (
            ("dwFileAttributes", wintypes.DWORD),
            ("ftCreationTime", wintypes.FILETIME),
            ("ftLastAccessTime", wintypes.FILETIME),
            ("ftLastWriteTime", wintypes.FILETIME),
            ("dwVolumeSerialNumber", wintypes.DWORD),
            ("nFileSizeHigh", wintypes.DWORD),
            ("nFileSizeLow", wintypes.DWORD),
            ("nNumberOfLinks", wintypes.DWORD),
            ("nFileIndexHigh", wintypes.DWORD),
            ("nFileIndexLow", wintypes.DWORD),
        )

    _ntdll.NtCreateFile.restype = wintypes.LONG  # NTSTATUS
    _ntdll.NtCreateFile.argtypes = [
        ctypes.POINTER(wintypes.HANDLE),
        wintypes.DWORD,
        ctypes.POINTER(_OBJECT_ATTRIBUTES),
        ctypes.POINTER(_IO_STATUS_BLOCK),
        ctypes.POINTER(ctypes.c_longlong),  # PLARGE_INTEGER AllocationSize
        wintypes.ULONG,
        wintypes.ULONG,
        wintypes.ULONG,
        wintypes.ULONG,
        ctypes.c_void_p,
        wintypes.ULONG,
    ]
    _ntdll.RtlDosPathNameToNtPathName_U_WithStatus.restype = wintypes.LONG  # NTSTATUS
    _ntdll.RtlDosPathNameToNtPathName_U_WithStatus.argtypes = [
        wintypes.LPCWSTR,
        ctypes.POINTER(_UNICODE_STRING),
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    _ntdll.RtlFreeUnicodeString.restype = None
    _ntdll.RtlFreeUnicodeString.argtypes = [ctypes.POINTER(_UNICODE_STRING)]
    _ntdll.RtlNtStatusToDosError.restype = wintypes.ULONG
    _ntdll.RtlNtStatusToDosError.argtypes = [wintypes.LONG]

    _kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    _kernel32.CloseHandle.restype = wintypes.BOOL
    _kernel32.GetFileInformationByHandle.argtypes = [wintypes.HANDLE, ctypes.POINTER(_BY_HANDLE_FILE_INFORMATION)]
    _kernel32.GetFileInformationByHandle.restype = wintypes.BOOL
    _kernel32.LockFileEx.argtypes = [
        wintypes.HANDLE,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.DWORD,
        ctypes.POINTER(_OVERLAPPED),
    ]
    _kernel32.LockFileEx.restype = wintypes.BOOL
    _kernel32.UnlockFileEx.argtypes = [
        wintypes.HANDLE,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.DWORD,
        ctypes.POINTER(_OVERLAPPED),
    ]
    _kernel32.UnlockFileEx.restype = wintypes.BOOL

    def _lock_fd_nonblocking(fd: int) -> bool:
        # One nonblocking exclusive LockFileEx attempt shared by WindowsFileLock and lock_descriptor, over the one-byte
        # range at offset 0. True on acquisition, False on contention, raise otherwise. The caller owns fd; the handle
        # from get_osfhandle belongs to the CRT descriptor and must not be closed here.
        overlapped = _OVERLAPPED()  # zero-initialized, so Offset/OffsetHigh/hEvent are 0
        flags = _LOCKFILE_EXCLUSIVE_LOCK | _LOCKFILE_FAIL_IMMEDIATELY
        if _kernel32.LockFileEx(msvcrt.get_osfhandle(fd), flags, 0, 1, 0, ctypes.byref(overlapped)):
            return True
        err = ctypes.get_last_error()
        if err == _ERROR_LOCK_VIOLATION:
            return False
        # A non-contention LockFileEx failure is not reproducible in-process.
        raise ctypes.WinError(err)  # pragma: no cover

    def _unlock_fd(fd: int) -> None:
        overlapped = _OVERLAPPED()  # the same offset 0 and one-byte length the lock used
        # Unlocking the exact range we hold does not fail.
        if not _kernel32.UnlockFileEx(msvcrt.get_osfhandle(fd), 0, 1, 0, ctypes.byref(overlapped)):  # pragma: no cover
            raise ctypes.WinError(ctypes.get_last_error())

    class WindowsFileLock(BaseFileLock):
        """
        Uses ``LockFileEx`` to hard lock a byte range of the lock file on Windows systems.

        Lock file cleanup: Windows attempts to delete the lock file after release, but deletion is
        not guaranteed in multi-threaded scenarios where another thread holds an open handle. The lock
        file may persist on disk, which does not affect lock correctness.
        """

        def _acquire(self) -> None:
            raise_on_not_writable_file(self.lock_file)
            ensure_directory_exists(self.lock_file)

            # The reparse test is bound to the opened handle, so a symlink or junction swapped in cannot defeat it
            # through a check-then-open TOCTOU race.
            fd = _open_non_reparse_fd(self.lock_file, self._open_mode())
            if fd is None:
                return  # open contention (share conflict or a name pending deletion); let the retry loop try again
            try:
                locked = _lock_fd_nonblocking(fd)
                if locked:
                    self._mark_descriptor_owned(fd)
            except BaseException:  # pragma: no cover  # cleanup only if the lock attempt itself raises
                os.close(fd)
                raise
            if not locked:
                os.close(fd)  # another holder owns the byte-range lock; let the retry loop try again

        def _release(self) -> None:
            fd = cast("int", self._context.lock_file_fd)
            # Retain the descriptor until the OS unlock succeeds: if UnlockFileEx raises, the byte-range lock is still
            # held, so is_locked must keep reporting held rather than losing the fd. Only after the unlock commits do
            # close and unlink run as post-unlock cleanup; their failure cannot make the lock held again.
            _unlock_fd(fd)
            self._mark_descriptor_released()
            self._close_released_fd(fd, default_suppresses=False)
            if not self._preserve_lock_file:  # preserve_lock_file keeps a stable file identity for the caller (#605)
                with suppress(OSError):
                    Path(self.lock_file).unlink()

    def _open_non_reparse_fd(path: str, mode: int) -> int | None:
        """
        Open *path* for locking while refusing reparse points, bound to the handle actually locked.

        The file is opened through ``NtCreateFile`` with ``FILE_OPEN_REPARSE_POINT`` so a symlink or junction planted
        at the path is not followed, and the reparse decision is read from *that* handle via
        ``GetFileInformationByHandle`` rather than from a prior pathname query. Reading the held handle closes the
        check-then-open race: an attacker cannot swap the path between validation and use because both act on the same
        handle. Share mode omits delete so a peer cannot unlink or rename the file out from under a live holder,
        matching ``os.open``'s ``_SH_DENYNO``.

        ``NtCreateFile`` is used instead of ``CreateFileW`` because its return value carries the raw ``NTSTATUS``.
        Windows collapses a transient delete-pending name and a permanent access denial into the same Win32
        ``ERROR_ACCESS_DENIED``; the status keeps them apart, so a real denial fails fast instead of spinning until the
        caller's timeout (#604).

        The reparse option only guards the final path component; Windows still follows reparse points in intermediate
        directories. This assumes the lock file sits in a lock directory untrusted users cannot modify. A path with
        attacker-controlled parent directories would need component-by-component handle validation.

        :param path: the lock file path.
        :param mode: the permission mode; as ``os.open`` does on Windows, a cleared owner-write bit creates the file
            read-only. The attribute only takes effect when the file is created, not when an existing one is opened.

        :returns: a file descriptor owning the opened handle, or ``None`` on a sharing violation or a delete-pending
            name the caller should treat as contention and retry.

        :raises OSError: if the path resolves to a reparse point, or the open fails for any other reason, raised with
            the Win32 error the status maps to.

        """
        # Emit the audit event os.open would, so consumers watching "open" still see the path-level open and can veto.
        sys.audit("open", path, None, os.O_RDWR | os.O_CREAT)
        handle, status = _nt_open(path, read_only=not mode & _OWNER_WRITE)
        if status != _STATUS_SUCCESS:
            if status in {_STATUS_SHARING_VIOLATION, _STATUS_DELETE_PENDING}:
                return None
            winerror = _ntdll.RtlNtStatusToDosError(status)
            raise OSError(None, ctypes.FormatError(winerror).strip(), path, winerror)

        info = _BY_HANDLE_FILE_INFORMATION()
        # Querying an open handle we just created does not fail.
        if not _kernel32.GetFileInformationByHandle(handle, ctypes.byref(info)):  # pragma: no cover
            err = ctypes.get_last_error()
            _kernel32.CloseHandle(handle)
            raise ctypes.WinError(err)
        if info.dwFileAttributes & _FILE_ATTRIBUTE_REPARSE_POINT:
            _kernel32.CloseHandle(handle)
            msg = f"Lock file is a reparse point (symlink/junction): {path}"
            raise OSError(msg)

        try:
            # O_NOINHERIT mirrors os.open on Windows: the lock fd must not leak into child processes.
            return msvcrt.open_osfhandle(handle, os.O_RDWR | os.O_NOINHERIT)
        except BaseException:  # pragma: no cover  # open_osfhandle audits too; a hook raising must not leak the handle
            _kernel32.CloseHandle(handle)
            raise

    def _nt_open(path: str, *, read_only: bool) -> tuple[int, int]:
        """
        Open *path* through ``NtCreateFile`` and return ``(handle, status)``.

        ``RtlDosPathNameToNtPathName_U_WithStatus`` translates the Win32 path to the NT namespace, handling relative,
        drive, UNC and extended-length path forms as Win32 itself would, and allocates a buffer that
        ``RtlFreeUnicodeString`` releases. The handle is ``0`` unless the status is ``STATUS_SUCCESS``.
        """
        nt_name = _UNICODE_STRING()
        status = _ntdll.RtlDosPathNameToNtPathName_U_WithStatus(path, ctypes.byref(nt_name), None, None) & 0xFFFFFFFF
        if status != _STATUS_SUCCESS:
            return 0, status
        try:
            attributes = _OBJECT_ATTRIBUTES()
            attributes.Length = ctypes.sizeof(_OBJECT_ATTRIBUTES)
            attributes.ObjectName = ctypes.pointer(nt_name)
            attributes.Attributes = _OBJ_CASE_INSENSITIVE
            handle = wintypes.HANDLE()
            io_status = _IO_STATUS_BLOCK()
            status = (
                _ntdll.NtCreateFile(
                    ctypes.byref(handle),
                    _DESIRED_ACCESS,
                    ctypes.byref(attributes),
                    ctypes.byref(io_status),
                    None,
                    _FILE_ATTRIBUTE_READONLY if read_only else _FILE_ATTRIBUTE_NORMAL,
                    _FILE_SHARE_READ_WRITE,
                    _FILE_OPEN_IF,
                    _CREATE_OPTIONS,
                    None,
                    0,
                )
                & 0xFFFFFFFF
            )
        finally:
            _ntdll.RtlFreeUnicodeString(ctypes.byref(nt_name))
        if status != _STATUS_SUCCESS:
            return 0, status
        return handle.value or 0, status

else:  # pragma: win32 no cover

    class WindowsFileLock(BaseFileLock):
        """Uses ``LockFileEx`` to hard lock a byte range of the lock file on Windows systems."""

        def _acquire(self) -> None:
            raise NotImplementedError

        def _release(self) -> None:
            raise NotImplementedError


__all__ = [
    "WindowsFileLock",
]
