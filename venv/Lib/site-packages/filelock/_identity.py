from __future__ import annotations

import os
import socket
import sys
from errno import EPERM, ESRCH
from pathlib import Path
from typing import Final


def host_name() -> str:
    """The hostname recorded alongside an owner, so a marker written on another machine is never probed here."""
    return socket.gethostname()


def owner_is_stale(pid: int, hostname: str, start_token: int | None) -> bool:
    """
    Whether the recorded owner is provably gone, so reclaiming its marker cannot detach a live holder.

    Fail closed: return ``True`` only when this process can prove the exact recorded owner is dead. A marker from
    another host cannot be probed; a live PID whose start token still matches, or whose token cannot be read, is the
    holder or is indistinguishable from it; a live PID whose start token differs is a recycled PID, so the process that
    wrote the marker is gone. PostgreSQL, Qt ``QLockFile`` and Mercurial all break a stale lock only on proof of death
    and treat an unreadable or foreign owner as still holding.
    """
    if hostname != host_name():
        return False
    if not process_alive(pid):
        return True
    if start_token is None:
        return False
    current = process_start_token(pid)
    return current is not None and current != start_token


if sys.platform == "win32":  # pragma: win32 cover
    import ctypes
    from ctypes import wintypes

    _KERNEL32: Final[ctypes.WinDLL] = ctypes.WinDLL("kernel32", use_last_error=True)
    _KERNEL32.CloseHandle.argtypes = [wintypes.HANDLE]
    _KERNEL32.CloseHandle.restype = wintypes.BOOL
    _KERNEL32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    _KERNEL32.OpenProcess.restype = wintypes.HANDLE
    _KERNEL32.GetProcessTimes.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(wintypes.FILETIME),
        ctypes.POINTER(wintypes.FILETIME),
        ctypes.POINTER(wintypes.FILETIME),
        ctypes.POINTER(wintypes.FILETIME),
    ]
    _KERNEL32.GetProcessTimes.restype = wintypes.BOOL

    _WIN_SYNCHRONIZE: Final[int] = 0x100000
    _WIN_PROCESS_QUERY_LIMITED_INFORMATION: Final[int] = 0x1000
    _WIN_ERROR_INVALID_PARAMETER: Final[int] = 87
    _WIN_INHERIT_HANDLE: Final[bool] = False

    def process_alive(pid: int) -> bool:
        """Whether a process with this PID exists, treating an access denial as proof it does."""
        handle = _KERNEL32.OpenProcess(_WIN_SYNCHRONIZE, _WIN_INHERIT_HANDLE, pid)
        if handle:
            _KERNEL32.CloseHandle(handle)
            return True
        return ctypes.get_last_error() != _WIN_ERROR_INVALID_PARAMETER

    def process_start_token(pid: int) -> int | None:
        """The process creation FILETIME as a 100ns tick count, or ``None`` when it cannot be read."""
        handle = _KERNEL32.OpenProcess(_WIN_PROCESS_QUERY_LIMITED_INFORMATION, _WIN_INHERIT_HANDLE, pid)
        if not handle:
            return None
        creation, exit_time, kernel_time, user_time = (wintypes.FILETIME() for _ in range(4))
        try:
            if not _KERNEL32.GetProcessTimes(
                handle,
                ctypes.byref(creation),
                ctypes.byref(exit_time),
                ctypes.byref(kernel_time),
                ctypes.byref(user_time),
            ):
                return None  # pragma: no cover  # win32 GetProcessTimes failure path; not reproducible on a live handle
        finally:
            _KERNEL32.CloseHandle(handle)
        return (creation.dwHighDateTime << 32) | creation.dwLowDateTime

else:  # pragma: win32 no cover

    def process_alive(pid: int) -> bool:
        """Whether a process with this PID exists, treating an access denial (``EPERM``) as proof it does."""
        try:
            os.kill(pid, 0)
        except OSError as error:
            if error.errno == ESRCH:
                return False
            if error.errno == EPERM:
                return True
            raise
        return True

    if sys.platform in {"linux", "android"}:  # pragma: linux cover
        # Termux/Android reports sys.platform == "android" but runs the Linux kernel, so /proc/<pid>/stat and the boot
        # id are the same reliable start-time source; treat it exactly like Linux rather than the tokenless fallback.
        # comm (field 2) is wrapped in parentheses and may itself contain spaces or a ')', so the fixed fields start
        # after the final ')'. starttime is field 22 overall, the twentieth of those trailing fields (index 19).
        _STARTTIME_INDEX: Final[int] = 19

        def _read_boot_id() -> int:
            # starttime is measured in clock ticks since boot, so on its own it repeats across a reboot. Folding the
            # boot id into the high bits makes the Linux token reboot-safe like the absolute clocks macOS and Windows
            # expose, while staying a single integer so a 3.29 reader still parses the third marker line. 0 when the
            # kernel does not expose a boot id degrades to bare starttime, which stays safe (a reboot collision fails
            # closed rather than reclaiming a live marker).
            try:
                boot_id = Path("/proc/sys/kernel/random/boot_id").read_text(encoding="ascii")
                return int(boot_id.strip().replace("-", ""), 16)
            except (OSError, ValueError):  # pragma: no cover  # the kernel always exposes boot_id as a UUID on Linux
                return 0

        _BOOT_ID: Final[int] = _read_boot_id()

        def process_start_token(pid: int) -> int | None:
            """The ``/proc/<pid>/stat`` ``starttime`` folded with the boot id, or ``None`` when the process is gone."""
            try:
                data = Path(f"/proc/{pid}/stat").read_bytes()
            except OSError:
                return None
            # psutil identifies a process by the same (pid, starttime) pair; the boot id extends that across reboots.
            fields = data[data.rfind(b")") + 1 :].split()
            if len(fields) <= _STARTTIME_INDEX:  # pragma: no cover  # a truncated /proc read, never seen in practice
                return None
            try:
                starttime = int(fields[_STARTTIME_INDEX])
            except ValueError:  # pragma: no cover  # /proc always renders starttime as an integer
                return None
            return (_BOOT_ID << 64) | starttime

    elif sys.platform == "darwin":  # pragma: darwin cover
        import ctypes
        import struct

        _LIBC: Final[ctypes.CDLL] = ctypes.CDLL(None, use_errno=True)
        _CTL_KERN: Final[int] = 1
        _KERN_PROC: Final[int] = 14
        _KERN_PROC_PID: Final[int] = 1
        # kinfo_proc opens with kp_proc.p_starttime (a struct timeval) at offset 0: int64 seconds, int32 microseconds.
        # The offset is fixed by the struct chain kinfo_proc -> extern_proc -> p_un, so a read at 0 is not a guess.
        _TIMEVAL_AT_ZERO: Final[str] = "<qi"
        _TIMEVAL_SIZE: Final[int] = struct.calcsize(_TIMEVAL_AT_ZERO)

        def process_start_token(pid: int) -> int | None:
            """The process start time in microseconds from ``sysctl(KERN_PROC_PID)``, or ``None`` when it is gone."""
            mib = (ctypes.c_int * 4)(_CTL_KERN, _KERN_PROC, _KERN_PROC_PID, pid)
            length = ctypes.c_size_t(0)
            # The size probe reports the kinfo_proc size for any PID, so the read below, not the probe, tells a live
            # process from a gone one: a gone PID leaves the fetch a zero-length success, so re-check the length after.
            if _LIBC.sysctl(mib, 4, None, ctypes.byref(length), None, 0) != 0:
                return None
            buffer = (ctypes.c_char * length.value)()
            if _LIBC.sysctl(mib, 4, buffer, ctypes.byref(length), None, 0) != 0 or length.value < _TIMEVAL_SIZE:
                return None
            seconds, microseconds = struct.unpack_from(_TIMEVAL_AT_ZERO, buffer.raw, 0)
            return seconds * 1_000_000 + microseconds

    else:  # pragma: no cover  # a POSIX platform without a proven start-time source falls back to fail-closed liveness

        def process_start_token(pid: int) -> int | None:
            """No proven start-time source, so the owner carries no token and liveness rests on the PID alone."""
            del pid
            return None


__all__ = [
    "host_name",
    "owner_is_stale",
    "process_alive",
    "process_start_token",
]
