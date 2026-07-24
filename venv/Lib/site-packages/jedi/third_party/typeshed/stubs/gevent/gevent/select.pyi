import sys
from _typeshed import FileDescriptorLike
from collections.abc import Iterable
from select import error as error
from typing import Any

def select(
    rlist: Iterable[Any], wlist: Iterable[Any], xlist: Iterable[Any], timeout: float | None = None
) -> tuple[list[Any], list[Any], list[Any]]: ...

if sys.platform != "win32":
    __all__ = ["error", "poll", "select"]
else:
    __all__ = ["error", "select"]

class poll:
    def register(self, fd: FileDescriptorLike, eventmask: int = ...) -> None: ...
    def modify(self, fd: FileDescriptorLike, eventmask: int) -> None: ...
    def poll(self, timeout: float | None = None) -> list[tuple[int, int]]: ...
    def unregister(self, fd: FileDescriptorLike) -> None: ...
