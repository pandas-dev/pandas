from __future__ import annotations

import os
import platform
import sys

import pytest

skip_without_chown = pytest.mark.skipif(
    not hasattr(os, "chown"), reason="os.chown not supported"
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
