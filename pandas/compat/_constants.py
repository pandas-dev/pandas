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

PY310 = sys.version_info >= (3, 10)
PY311 = sys.version_info >= (3, 11)
PYPY = platform.python_implementation() == "PyPy"
ISMUSL = "musl" in (sysconfig.get_config_var("HOST_GNU_TYPE") or "")


__all__ = [
    "IS64",
    "ISMUSL",
    "PY310",
    "PY311",
    "PYPY",
]
