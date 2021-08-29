from typing import (
    Callable,
    Protocol,
    TypeVar,
    overload,
)

_F = TypeVar("_F")
_G = TypeVar("_G")

class cache_readonly(Protocol[_F, _G]):
    def __init__(self, func: Callable[[_F], _G]) -> None: ...
    @overload
    def __get__(self, obj: _F, typ) -> _G: ...
    @overload
    def __get__(self, obj: None, typ) -> cache_readonly[_F, _G]: ...
    def __set__(self, obj: _F, value: _G) -> None: ...
