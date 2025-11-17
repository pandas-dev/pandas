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

PY312 = sys.version_info >= (3, 12)
PY314 = sys.version_info >= (3, 14)
PYPY = platform.python_implementation() == "PyPy"
WASM = (sys.platform == "emscripten") or (platform.machine() in ["wasm32", "wasm64"])
ISMUSL = "musl" in (sysconfig.get_config_var("HOST_GNU_TYPE") or "")
REF_COUNT = 2
CHAINED_WARNING_DISABLED = PYPY or (PY314 and not sys._is_gil_enabled())  # type: ignore[attr-defined]
CHAINED_WARNING_DISABLED_INPLACE_METHOD = PYPY or PY314


__all__ = [
    "IS64",
    "ISMUSL",
    "PY312",
    "PY314",
    "PYPY",
    "WASM",
]
