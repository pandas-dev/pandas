"""
compat
======

Cross-compatible functions for different versions of Python.

Other items:
* platform checker
"""

from __future__ import annotations

import builtins
import platform
import sys
from typing import (
    TYPE_CHECKING,
    Any,
    TypeGuard,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

from pandas.compat._constants import (
    CHAINED_WARNING_DISABLED,
    IS64,
    ISMUSL,
    PY312,
    PY313,
    PY314,
    PY315,
    PYPY,
    WASM,
)

# built-in frozendict in Python 3.15+ (GH#71853)
DICT_TYPES: tuple[type, ...]
frozendict_type = getattr(builtins, "frozendict", None)
if frozendict_type is not None:
    DICT_TYPES = (dict, frozendict_type)
else:
    DICT_TYPES = (dict,)


def is_dict_or_frozendict(val: object) -> TypeGuard[Mapping[Any, Any]]:
    return isinstance(val, DICT_TYPES)


from pandas.compat.numpy import is_numpy_dev
from pandas.compat.pyarrow import (
    HAS_PYARROW,
    PYARROW_INSTALLED,
    PYARROW_MIN_VERSION,
    pa_version_under17p0,
    pa_version_under18p0,
    pa_version_under19p0,
    pa_version_under20p0,
    pa_version_under21p0,
    pa_version_under22p0,
    pa_version_under23p0,
    pa_version_under25p0,
    pa_version_under26p0,
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


__all__ = [
    "CHAINED_WARNING_DISABLED",
    "DICT_TYPES",
    "HAS_PYARROW",
    "IS64",
    "ISMUSL",
    "PY312",
    "PY313",
    "PY314",
    "PY315",
    "PYARROW_INSTALLED",
    "PYARROW_MIN_VERSION",
    "PYPY",
    "WASM",
    "is_dict_or_frozendict",
    "is_numpy_dev",
    "pa_version_under17p0",
    "pa_version_under18p0",
    "pa_version_under19p0",
    "pa_version_under20p0",
    "pa_version_under21p0",
    "pa_version_under22p0",
    "pa_version_under23p0",
    "pa_version_under25p0",
    "pa_version_under26p0",
]
