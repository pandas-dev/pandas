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
PY314 = sys.version_info >= (3, 14)
PYPY = platform.python_implementation() == "PyPy"
WASM = (sys.platform == "emscripten") or (platform.machine() in ["wasm32", "wasm64"])
ISMUSL = "musl" in (sysconfig.get_config_var("HOST_GNU_TYPE") or "")
if PY311:
    REF_COUNT = 2
elif PY314:
    REF_COUNT = 1
else:
    # Python 3.10 and older
    REF_COUNT = 3

__all__ = [
    "IS64",
    "ISMUSL",
    "PY311",
    "PY312",
    "PYPY",
    "WASM",
]
