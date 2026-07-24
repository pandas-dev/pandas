from collections.abc import Iterator
from typing import Generic, TypeVar
from typing_extensions import Self

from requests import Response

_T_co = TypeVar("_T_co", covariant=True)

class CancellableStream(Generic[_T_co]):
    def __init__(self, stream: Iterator[_T_co], response: Response) -> None: ...
    def __iter__(self) -> Self: ...
    def __next__(self) -> _T_co: ...
    next = __next__
    def close(self) -> None: ...
