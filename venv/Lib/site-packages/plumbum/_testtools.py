from __future__ import annotations

__lazy_modules__ = {"tempfile"}

import os
import platform
import sys
import tempfile

import pytest


def _new_files_inherit_process_gid() -> bool:
    """Whether a newly created file's group is the process gid.

    On BSD/macOS a new file inherits its parent directory's group rather than
    the creating process's gid, so chown tests that assert ``p.gid ==
    os.getgid()`` do not hold (e.g. when the temp dir is group-owned by wheel).
    """
    if not hasattr(os, "getgid"):
        return False
    try:
        with tempfile.NamedTemporaryFile() as f:
            return os.fstat(f.fileno()).st_gid == os.getgid()
    except OSError:
        return False


skip_without_chown = pytest.mark.skipif(
    not hasattr(os, "chown") or not _new_files_inherit_process_gid(),
    reason="os.chown not supported, or new files do not inherit the process gid (e.g. BSD/macOS)",
)

skip_without_tty = pytest.mark.skipif(not sys.stdin.isatty(), reason="Not a TTY")

skip_on_windows = pytest.mark.skipif(
    sys.platform == "win32", reason="Windows not supported for this test (yet)"
)

xfail_on_windows = pytest.mark.xfail(
    sys.platform == "win32", reason="Windows not supported for this test (yet)"
)

xfail_on_pypy = pytest.mark.xfail(
    platform.python_implementation() == "PyPy",
    reason="PyPy is currently not working on this test!",
)

__all__ = [
    "skip_on_windows",
    "skip_without_chown",
    "skip_without_tty",
    "xfail_on_pypy",
    "xfail_on_windows",
]


def __dir__() -> list[str]:
    return list(__all__)
