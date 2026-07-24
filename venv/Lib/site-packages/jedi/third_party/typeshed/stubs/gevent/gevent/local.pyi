from typing import Any
from typing_extensions import Self

class local:
    __slots__ = (
        "_local__impl",
        "_local_type_set_descriptors",
        "_local_type_get_descriptors",
        "_local_type_vars",
        "_local_type_del_descriptors",
        "_local_type",
        "_local_type_set_or_del_descriptors",
    )
    def __new__(cls, *args: object, **kwargs: object) -> Self: ...
    def __copy__(self) -> Self: ...
    def __getattribute__(self, name: str, /) -> Any: ...
    def __delattr__(self, name: str, /) -> None: ...
    def __setattr__(self, name: str, value: Any, /) -> None: ...

__all__ = ["local"]
