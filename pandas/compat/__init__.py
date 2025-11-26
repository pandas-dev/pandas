"""
compat
======

Cross-compatible functions for different versions of Python.

Other items:
* platform checker
"""

from __future__ import annotations

import os
import platform
import sys
from typing import TYPE_CHECKING

from pandas.compat._constants import (
    CHAINED_WARNING_DISABLED,
    IS64,
    ISMUSL,
    PY312,
    PY314,
    PYPY,
    WASM,
)
from pandas.compat.numpy import is_numpy_dev
from pandas.compat.pyarrow import (
    HAS_PYARROW,
    PYARROW_MIN_VERSION,
    pa_version_under14p0,
    pa_version_under14p1,
    pa_version_under16p0,
    pa_version_under17p0,
    pa_version_under18p0,
    pa_version_under19p0,
    pa_version_under20p0,
    pa_version_under21p0,
)

if TYPE_CHECKING:
    from pandas._typing import F


def set_function_name(f: F, name: str, cls: type) -> F:
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


def is_platform_power() -> bool:
    """
    Checking if the running platform use Power architecture.

    Returns
    -------
    bool
        True if the running platform uses ARM architecture.
    """
    return platform.machine() in ("ppc64", "ppc64le")


def is_platform_riscv64() -> bool:
    """
    Checking if the running platform use riscv64 architecture.

    Returns
    -------
    bool
        True if the running platform uses riscv64 architecture.
    """
    return platform.machine() == "riscv64"


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


__all__ = [
    "CHAINED_WARNING_DISABLED",
    "HAS_PYARROW",
    "IS64",
    "ISMUSL",
    "PY312",
    "PY314",
    "PYARROW_MIN_VERSION",
    "PYPY",
    "WASM",
    "is_numpy_dev",
    "pa_version_under14p0",
    "pa_version_under14p1",
    "pa_version_under16p0",
    "pa_version_under17p0",
    "pa_version_under18p0",
    "pa_version_under19p0",
    "pa_version_under20p0",
    "pa_version_under21p0",
]
