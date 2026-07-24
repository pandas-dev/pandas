from collections.abc import Iterable, Iterator
from typing import Any, Generic, TypeVar

_T = TypeVar("_T", default=Any)

class Filter(Generic[_T]):
    source: Iterable[_T]
    def __init__(self, source: Iterable[_T]) -> None: ...
    def __iter__(self) -> Iterator[_T]: ...
    def __getattr__(self, name: str) -> Any: ...  # Depends on `source`
