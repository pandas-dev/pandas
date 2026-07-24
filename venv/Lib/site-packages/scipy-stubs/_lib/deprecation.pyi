from collections.abc import Callable
from types import ModuleType
from typing import TypeVar

__all__ = ["_deprecated"]

_F = TypeVar("_F", bound=Callable[..., object])

class _DeprecationHelperStr:
    def __init__(self, /, content: str, message: str) -> None: ...

def _deprecated(msg: str, stacklevel: int = 2) -> Callable[[_F], _F]: ...
def deprecate_cython_api(
    module: ModuleType, routine_name: str, new_name: str | None = None, message: str | None = None
) -> None: ...
