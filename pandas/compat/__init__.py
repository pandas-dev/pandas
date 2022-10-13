"""
compat
======

Cross-compatible functions for different versions of Python.

Other items:
* platform checker
"""
from __future__ import annotations

import os
from pickle import PickleBuffer
import platform
import sys

try:
    import lzma

    has_lzma = True
except ImportError:
    has_lzma = False

from pandas._typing import F
from pandas.compat.numpy import (
    is_numpy_dev,
    np_version_under1p21,
)
from pandas.compat.pyarrow import (
    pa_version_under1p01,
    pa_version_under2p0,
    pa_version_under3p0,
    pa_version_under4p0,
    pa_version_under5p0,
    pa_version_under6p0,
    pa_version_under7p0,
    pa_version_under8p0,
    pa_version_under9p0,
)

PY39 = sys.version_info >= (3, 9)
PY310 = sys.version_info >= (3, 10)
PY311 = sys.version_info >= (3, 11)
PYPY = platform.python_implementation() == "PyPy"
IS64 = sys.maxsize > 2**32


if has_lzma:

    class _LZMAFile(lzma.LZMAFile):
        def write(self, b) -> int:
            if isinstance(b, PickleBuffer):
                # Workaround issue where `lzma.LZMAFile` expects `len`
                # to return the number of bytes in `b` by converting
                # `b` into something that meets that constraint with
                # minimal copying.
                try:
                    # coerce to 1-D `uint8` C-contiguous `memoryview` zero-copy
                    b = b.raw()
                except BufferError:
                    # perform in-memory copy if buffer is not contiguous
                    b = memoryview(b).tobytes()
            return super().write(b)


def set_function_name(f: F, name: str, cls) -> F:
    """
    Bind the name/qualname attributes of the function.
    """
    f.__name__ = name
    f.__qualname__ = f"{cls.__name__}.{name}"
    f.__module__ = cls.__module__
    return f


def is_platform_little_endian() -> bool:
    """
    Checking if the running platform is little endian.

    Returns
    -------
    bool
        True if the running platform is little endian.
    """
    return sys.byteorder == "little"


def is_platform_windows() -> bool:
    """
    Checking if the running platform is windows.

    Returns
    -------
    bool
        True if the running platform is windows.
    """
    return sys.platform in ["win32", "cygwin"]


def is_platform_linux() -> bool:
    """
    Checking if the running platform is linux.

    Returns
    -------
    bool
        True if the running platform is linux.
    """
    return sys.platform == "linux"


def is_platform_mac() -> bool:
    """
    Checking if the running platform is mac.

    Returns
    -------
    bool
        True if the running platform is mac.
    """
    return sys.platform == "darwin"


def is_platform_arm() -> bool:
    """
    Checking if the running platform use ARM architecture.

    Returns
    -------
    bool
        True if the running platform uses ARM architecture.
    """
    return platform.machine() in ("arm64", "aarch64") or platform.machine().startswith(
        "armv"
    )


def is_ci_environment() -> bool:
    """
    Checking if running in a continuous integration environment by checking
    the PANDAS_CI environment variable.

    Returns
    -------
    bool
        True if the running in a continuous integration environment.
    """
    return os.environ.get("PANDAS_CI", "0") == "1"


def get_lzma_file() -> type[_LZMAFile]:
    """
    Importing the `LZMAFile` class from the `lzma` module.

    Returns
    -------
    class
        The `LZMAFile` class from the `lzma` module.

    Raises
    ------
    RuntimeError
        If the `lzma` module was not imported correctly, or didn't exist.
    """
    if not has_lzma:
        raise RuntimeError(
            "lzma module not available. "
            "A Python re-install with the proper dependencies, "
            "might be required to solve this issue."
        )
    return _LZMAFile


__all__ = [
    "is_numpy_dev",
    "np_version_under1p21",
    "pa_version_under1p01",
    "pa_version_under2p0",
    "pa_version_under3p0",
    "pa_version_under4p0",
    "pa_version_under5p0",
    "pa_version_under6p0",
    "pa_version_under7p0",
    "pa_version_under8p0",
    "pa_version_under9p0",
]
