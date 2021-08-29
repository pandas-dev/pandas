from __future__ import annotations

from typing import (
    Callable,
    Protocol,
    Type,
    TypeVar,
    overload,
)

_F = TypeVar("_F")
_G = TypeVar("_G")

class cache_readonly(Protocol[_F, _G]):
    def __init__(self, func: Callable[[_F], _G]) -> None: ...
    @overload
    def __get__(self, obj: _F, typ: Type[_F] | None) -> _G: ...
    @overload
    def __get__(self, obj: None, typ: Type[_F] | None) -> cache_readonly[_F, _G]: ...
    def __set__(self, obj: _F, value: _G) -> None: ...
