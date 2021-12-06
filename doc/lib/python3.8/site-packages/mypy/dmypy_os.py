import sys

from typing import Any, Callable

if sys.platform == 'win32':
    import ctypes
    from ctypes.wintypes import DWORD, HANDLE
    import subprocess

    PROCESS_QUERY_LIMITED_INFORMATION = ctypes.c_ulong(0x1000)

    kernel32 = ctypes.windll.kernel32
    OpenProcess = kernel32.OpenProcess  # type: Callable[[DWORD, int, int], HANDLE]
    GetExitCodeProcess = kernel32.GetExitCodeProcess  # type: Callable[[HANDLE, Any], int]
else:
    import os
    import signal


def alive(pid: int) -> bool:
    """Is the process alive?"""
    if sys.platform == 'win32':
        # why can't anything be easy...
        status = DWORD()
        handle = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION,
                             0,
                             pid)
        GetExitCodeProcess(handle, ctypes.byref(status))
        return status.value == 259  # STILL_ACTIVE
    else:
        try:
            os.kill(pid, 0)
        except OSError:
            return False
        return True


def kill(pid: int) -> None:
    """Kill the process."""
    if sys.platform == 'win32':
        subprocess.check_output("taskkill /pid {pid} /f /t".format(pid=pid))
    else:
        os.kill(pid, signal.SIGKILL)
