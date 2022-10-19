"""
_constants
======

Constants relevant for the Python implementation.
"""

from __future__ import annotations

import sys
import platform

PY39 = sys.version_info >= (3, 9)
PY310 = sys.version_info >= (3, 10)
PY311 = sys.version_info >= (3, 11)
PYPY = platform.python_implementation() == "PyPy"
IS64 = sys.maxsize > 2**32
