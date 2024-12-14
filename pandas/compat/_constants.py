"""
_constants
======

Constants relevant for the Python implementation.
"""

from __future__ import annotations

import platform
import sys
import sysconfig

IS64 = sys.maxsize > 2**32

PY311 = sys.version_info >= (3, 11)
PY312 = sys.version_info >= (3, 12)
PY313 = sys.version_info >= (3, 13)
PYPY = platform.python_implementation() == "PyPy"
WASM = (sys.platform == "emscripten") or (platform.machine() in ["wasm32", "wasm64"])
IS_FREE_THREADING = False if not PY313 else sys._is_gil_enabled()  # type: ignore[attr-defined]
ISMUSL = "musl" in (sysconfig.get_config_var("HOST_GNU_TYPE") or "")
REF_COUNT = 2 if PY311 else 3

__all__ = [
    "IS64",
    "ISMUSL",
    "IS_FREE_THREADING",
    "PY311",
    "PY312",
    "PY313",
    "PYPY",
    "WASM",
]
