from typing import (
    Any,
    Type,
)

CachedProperty: Type[property]
cache_readonly = CachedProperty

class AxisProperty:
    def __init__(self, axis: int = ..., doc: str = ...) -> None: ...
    def __get__(self, obj: Any, typ: Type) -> Any: ...
    def __set__(self, obj: Any, value: Any) -> None: ...
