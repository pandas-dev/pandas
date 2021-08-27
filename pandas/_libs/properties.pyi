from typing import (
    Callable,
    Protocol,
    TypeVar,
)

_F = TypeVar("_F", contravariant=True)
_G = TypeVar("_G")

class cache_readonly(Protocol[_F, _G]):
    def __init__(self, func: Callable[[_F], _G]) -> None: ...
    def __get__(self, obj: _F, typ) -> _G: ...
    def __set__(self, obj: _F, value: _G) -> None: ...
