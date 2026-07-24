from __future__ import annotations

__lazy_modules__ = {"inspect", "io"}

import inspect
import os
import sys
from contextlib import contextmanager
from io import StringIO
from typing import IO, TYPE_CHECKING, Any, TextIO

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

IS_WIN32 = sys.platform == "win32"


class ProcInfo:
    __slots__ = ("args", "pid", "stat", "uid")

    def __init__(self, pid: int, uid: int | str, stat: str, args: str):
        self.pid = pid
        self.uid = uid
        self.stat = stat
        self.args = args

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.pid!r}, {self.uid!r}, {self.stat!r}, {self.args!r})"


@contextmanager
def captured_stdout(stdin: str = "") -> Generator[TextIO, None, None]:
    """
    Captures stdout (similar to the redirect_stdout in Python 3.4+, but with slightly different arguments)
    """
    prevstdin = sys.stdin
    prevstdout = sys.stdout
    sys.stdin = StringIO(stdin)
    sys.stdout = StringIO()
    try:
        yield sys.stdout
    finally:
        sys.stdin = prevstdin
        sys.stdout = prevstdout


class StaticProperty:
    # This acts like a static property, allowing access via class or object.
    # This is a non-data descriptor.

    __slots__ = ("__doc__", "_function")

    def __init__(self, function: Callable[[], Any]):
        self._function = function
        self.__doc__ = function.__doc__

    def __get__(self, obj: object, klass: object = None) -> Any:
        return self._function()


def getdoc(obj: object) -> str | None:
    """
    This gets a docstring if available, and cleans it, but does not look up docs in
    inheritance tree (Pre Python 3.5 behavior of ``inspect.getdoc``).
    """
    try:
        doc = obj.__doc__
    except AttributeError:
        return None
    if not isinstance(doc, str):
        return None
    return inspect.cleandoc(doc)


def read_fd_decode_safely(fd: IO[str], size: int = 4096) -> tuple[bytes, str]:
    """
    This reads a utf-8 file descriptor and returns a chunk, growing up to
    three bytes if needed to decode the character at the end.

    Returns the data and the decoded text.
    """
    data = os.read(fd.fileno(), size)
    for _ in range(3):
        try:
            return data, data.decode("utf-8")
        except UnicodeDecodeError as e:  # noqa: PERF203
            if e.reason != "unexpected end of data":
                raise
            data += os.read(fd.fileno(), 1)

    return data, data.decode("utf-8")


__all__ = [
    "IS_WIN32",
    "ProcInfo",
    "StaticProperty",
    "captured_stdout",
    "getdoc",
    "read_fd_decode_safely",
]


def __dir__() -> list[str]:
    return list(__all__)
