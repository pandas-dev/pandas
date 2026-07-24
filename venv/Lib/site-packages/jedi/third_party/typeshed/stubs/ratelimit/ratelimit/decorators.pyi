from collections.abc import Callable
from typing import TypeVar
from typing_extensions import ParamSpec

_P = ParamSpec("_P")
_T = TypeVar("_T")

class RateLimitDecorator:
    def __init__(
        self, calls: int = 15, period: float = 900, clock: Callable[[], float] = ..., raise_on_limit: bool = True
    ) -> None: ...
    def __call__(self, func: Callable[_P, _T]) -> Callable[_P, _T]: ...

def sleep_and_retry(func: Callable[_P, _T]) -> Callable[_P, _T]: ...
