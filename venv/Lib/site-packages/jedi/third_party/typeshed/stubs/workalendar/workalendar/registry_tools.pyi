from collections.abc import Callable
from typing import TypeVar

_T = TypeVar("_T", bound=type)

def iso_register(iso_code: str) -> Callable[[_T], _T]: ...
