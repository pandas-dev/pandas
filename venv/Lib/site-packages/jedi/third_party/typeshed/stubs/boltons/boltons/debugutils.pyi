from collections.abc import Callable
from typing import Any

def pdb_on_signal(signalnum: int | None = None) -> None: ...
def pdb_on_exception(limit: int = 100) -> None: ...
def wrap_trace(
    obj, hook: Callable[..., Any] = ..., which: str | None = None, events: str | None = None, label: str | None = None
): ...

__all__ = ["pdb_on_signal", "pdb_on_exception", "wrap_trace"]
