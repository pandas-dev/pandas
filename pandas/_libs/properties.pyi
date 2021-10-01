from typing import (
    Any,
    Type,
)

cache_readonly = property

class AxisProperty:
    def __init__(self, axis: int = ..., doc: str = ...) -> None: ...
    def __get__(self, obj: Any, typ: Type) -> Any: ...  # TODO
    def __set__(self, obj: Any, value: Any) -> None: ...  # TODO
