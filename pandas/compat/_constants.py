"""
_constants
======

Constants relevant for the Python implementation.
"""

from __future__ import annotations

import platform
import sys

IS64 = sys.maxsize > 2**32

PY39 = sys.version_info >= (3, 9)
PY310 = sys.version_info >= (3, 10)
PY311 = sys.version_info >= (3, 11)
PYPY = platform.python_implementation() == "PyPy"


__all__ = [
    "IS64",
    "PY39",
    "PY310",
    "PY311",
    "PYPY",
]
