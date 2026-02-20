"""
Terminal size utility
---------------------
"""

from __future__ import annotations

import contextlib
import os
import platform
import warnings
from struct import Struct

from plumbum import local


def get_terminal_size(default: tuple[int, int] = (80, 25)) -> tuple[int, int]:
    """
    Get width and height of console; works on linux, os x, windows and cygwin

    Adapted from https://gist.github.com/jtriley/1108174
    Originally from: http://stackoverflow.com/questions/566746/how-to-get-console-window-width-in-python
    """
    current_os = platform.system()
    if current_os == "Windows":  # pragma: no cover
        size = _get_terminal_size_windows()
        if not size:
            # needed for window's python in cygwin's xterm!
            size = _get_terminal_size_tput()
    elif current_os in ("Linux", "Darwin", "FreeBSD", "SunOS") or current_os.startswith(
        "CYGWIN"
    ):
        size = _get_terminal_size_linux()

    else:  # pragma: no cover
        warnings.warn(
            "Plumbum does not know the type of the current OS for term size, defaulting to UNIX",
            stacklevel=2,
        )
        size = _get_terminal_size_linux()

    # we'll assume the standard 80x25 if for any reason we don't know the terminal size
    if size is None:
        return default
    return size


def _get_terminal_size_windows():  # pragma: no cover
    try:
        from ctypes import create_string_buffer, windll

        STDERR_HANDLE = -12
        h = windll.kernel32.GetStdHandle(STDERR_HANDLE)
        csbi_struct = Struct("hhhhHhhhhhh")
        csbi = create_string_buffer(csbi_struct.size)
        res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)
        if res:
            _, _, _, _, _, left, top, right, bottom, _, _ = csbi_struct.unpack(csbi.raw)
            return right - left + 1, bottom - top + 1
        return None
    except Exception:
        return None


def _get_terminal_size_tput():  # pragma: no cover
    # get terminal width
    # src: http://stackoverflow.com/questions/263890/how-do-i-find-the-width-height-of-a-terminal-window
    try:
        tput = local["tput"]
        cols = int(tput("cols"))
        rows = int(tput("lines"))
        return (cols, rows)
    except Exception:
        return None


def _ioctl_GWINSZ(fd: int) -> tuple[int, int] | None:
    # struct winsize has 4 unsigned shorts: ws_row, ws_col, ws_xpixel, ws_ypixel
    # We only need the first two (rows and cols), but must provide full struct size
    winsize = Struct("HHHH")  # 4 unsigned shorts = 8 bytes
    try:
        import fcntl
        import termios

        # TODO: Clean this up. Problems could be hidden by the broad except.
        result = fcntl.ioctl(fd, termios.TIOCGWINSZ, b"\x00" * winsize.size)
        rows, cols, _, _ = winsize.unpack(result)
        return rows, cols
    except Exception:
        return None


def _get_terminal_size_linux() -> tuple[int, int] | None:
    cr = _ioctl_GWINSZ(0) or _ioctl_GWINSZ(1) or _ioctl_GWINSZ(2)
    if not cr:
        with contextlib.suppress(Exception):
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = _ioctl_GWINSZ(fd)
            os.close(fd)
    if not cr:
        try:
            cr = (int(os.environ["LINES"]), int(os.environ["COLUMNS"]))
        except Exception:
            return None
    return cr[1], cr[0]
