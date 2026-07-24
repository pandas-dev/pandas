# coding: utf-8
"""Compatibility tricks for Python 3. Mainly to do with unicode.

This file is deprecated and will be removed in a future version.
"""

import platform
import builtins as builtin_mod

from .encoding import DEFAULT_ENCODING
from typing import Optional


def decode(s: bytes, encoding: str | None = None) -> str:
    encoding = encoding or DEFAULT_ENCODING
    return s.decode(encoding, "replace")


def encode(u: str, encoding: Optional[str]=None) -> bytes:
    encoding = encoding or DEFAULT_ENCODING
    return u.encode(encoding, "replace")


def cast_unicode(s: str | bytes, encoding: Optional[str]=None) -> str:
    if isinstance(s, bytes):
        return decode(s, encoding)
    return s


def safe_unicode(e):
    """unicode(e) with various fallbacks. Used for exceptions, which may not be
    safe to call unicode() on.
    """
    try:
        return str(e)
    except UnicodeError:
        pass

    try:
        return repr(e)
    except UnicodeError:
        pass

    return "Unrecoverably corrupt evalue"


# keep reference to builtin_mod because the kernel overrides that value
# to forward requests to a frontend.
def input(prompt=""):
    return builtin_mod.input(prompt)


def execfile(fname, glob, loc=None, compiler=None):
    __tracebackhide__ = "__ipython_bottom__"
    loc = loc if (loc is not None) else glob
    with open(fname, "rb") as f:
        compiler = compiler or compile
        exec(compiler(f.read(), fname, "exec"), glob, loc)


PYPY = platform.python_implementation() == "PyPy"
