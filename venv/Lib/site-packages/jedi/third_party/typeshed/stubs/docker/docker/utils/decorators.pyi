from collections.abc import Callable
from typing import TypeVar
from typing_extensions import ParamSpec

_P = ParamSpec("_P")
_R = TypeVar("_R")

def check_resource(resource_name: str) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]: ...
def minimum_version(version: str) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]: ...
def update_headers(f: Callable[_P, _R]) -> Callable[_P, _R]: ...
