from _typeshed import SupportsRead
from collections.abc import Callable
from typing import Any

__all__ = ("loads", "load", "TOMLDecodeError")

class TOMLDecodeError(ValueError): ...

def load(__fp: SupportsRead[bytes], *, parse_float: Callable[[str], Any] = ...) -> dict[str, Any]: ...
def loads(__s: str, *, parse_float: Callable[[str], Any] = ...) -> dict[str, Any]: ...
