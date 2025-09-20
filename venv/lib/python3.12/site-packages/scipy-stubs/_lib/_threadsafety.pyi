from collections.abc import Callable
from typing_extensions import TypeVar

from scipy._typing import ExitMixin

__all__ = ["ReentrancyError", "ReentrancyLock", "non_reentrant"]

_FT = TypeVar("_FT", bound=Callable[..., object])

class ReentrancyError(RuntimeError): ...

class ReentrancyLock(ExitMixin):
    def __init__(self, /, err_msg: str) -> None: ...
    def __enter__(self, /) -> None: ...
    def decorate(self, /, func: _FT) -> _FT: ...

def non_reentrant(err_msg: str | None = None) -> Callable[[_FT], _FT]: ...
