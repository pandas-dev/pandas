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
# the refcount for self in a chained __setitem__/.(i)loc indexing/method call
REF_COUNT = 2 if PY314 else 3
REF_COUNT_IDX = 2
REF_COUNT_METHOD = 1 if PY314 else 2
CHAINED_WARNING_DISABLED = PYPY


__all__ = [
    "IS64",
    "ISMUSL",
    "PY312",
    "PY314",
    "PYPY",
    "WASM",
]
