from typing import Any, final
from weakref import ref

@final
class ValuedWeakRef(ref):
    __slots__ = ("value",)
    value: Any

@final
class IdentRegistry:
    def __init__(self) -> None: ...
    def get_ident(self, obj: object) -> int: ...
    def __len__(self) -> int: ...

__all__ = ["IdentRegistry"]
