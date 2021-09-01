from __future__ import annotations

from typing import (
    Any,
    Callable,
    Generic,
    Type,
    TypeVar,
    overload,
)

T = TypeVar("T")
F = TypeVar("F", contravariant=True)
G = TypeVar("G", covariant=True)

class CachedProperty(Generic[F, G]):
    def __init__(self, func: Callable[[F], G]) -> None: ...
    @overload
    def __get__(self, obj: F, typ: Type[F]) -> G: ...
    @overload
    def __get__(self: T, obj: None, typ: Type[F]) -> T: ...
    def __set__(self, obj: F, value: Any) -> None: ...

cache_readonly = CachedProperty

class AxisProperty:
    def __init__(self, axis: int = ..., doc: str = ...) -> None: ...
    def __get__(self, obj: Any, typ: Type) -> Any: ...
    def __set__(self, obj: Any, value: Any) -> None: ...
