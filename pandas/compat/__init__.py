"""
compat
======

Cross-compatible functions for different versions of Python.

Other items:
* platform checker
"""
import platform
import struct
import sys
import warnings

PY35 = sys.version_info[:2] == (3, 5)
PY36 = sys.version_info >= (3, 6)
PY37 = sys.version_info >= (3, 7)
PY38 = sys.version_info >= (3, 8)
PYPY = platform.python_implementation() == "PyPy"


# ----------------------------------------------------------------------------
# functions largely based / taken from the six module

# Much of the code in this module comes from Benjamin Peterson's six library.
# The license for this library can be found in LICENSES/SIX and the code can be
# found at https://bitbucket.org/gutworth/six


def set_function_name(f, name, cls):
    """
    Bind the name/qualname attributes of the function.
    """
    f.__name__ = name
    f.__qualname__ = "{klass}.{name}".format(klass=cls.__name__, name=name)
    f.__module__ = cls.__module__
    return f


# https://github.com/pandas-dev/pandas/pull/9123
def is_platform_little_endian() -> bool:
    """
    Checks if the running platform is little endian.

    Returns
    -------
    bool
        True if the running platform is little endian.
    """
    return sys.byteorder == "little"


def is_platform_windows() -> bool:
    """
    Checks if the running platform is windows.

    Returns
    -------
    bool
        True if the running platform is windows.
    """
    return sys.platform == "win32" or sys.platform == "cygwin"


def is_platform_linux() -> bool:
    """
    Checks if the running platform is linux.

    Returns
    -------
    bool
        True if the running platform is linux.
    """
    return sys.platform == "linux2"


def is_platform_mac() -> bool:
    """
    Checks if the running platform is mac.

    Returns
    -------
    bool
        True if the running platform is mac.
    """
    return sys.platform == "darwin"


def is_platform_32bit() -> bool:
    """
    Checks if the running platform is 32-bit.

    Returns
    -------
    bool
        True if the running platform is 32-bit.
    """
    return struct.calcsize("P") * 8 < 64


def _import_lzma():
    """
    Attempts to import the lzma module.

    Warns
    -----
    When the lzma module is not available.
    """
    try:
        import lzma

        return lzma
    except ImportError:
        msg = (
            "Could not import the lzma module. "
            "Your installed Python is incomplete. "
            "Attempting to use lzma compression will result in a RuntimeError."
        )
        warnings.warn(msg)


def _get_lzma_file(lzma):
    """
    Returns
    -------
    method
        The lzma.LZMAFile method.

    Raises
    ------
    RuntimeError
        If the module lzma was not imported correctly, or didn't exist.
    """
    if lzma is None:
        raise RuntimeError(
            "lzma module not available. "
            "A Python re-install with the proper "
            "dependencies might be required to solve this issue."
        )
    return lzma.LZMAFile
