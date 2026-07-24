from _typeshed import Incomplete
from collections.abc import Callable
from typing import ClassVar

class XRayMiddleware:
    name: ClassVar[str]
    api: ClassVar[int]
    def __init__(self, recorder) -> None: ...
    def apply(self, callback: Callable[..., Incomplete], route) -> Callable[..., Incomplete]: ...
