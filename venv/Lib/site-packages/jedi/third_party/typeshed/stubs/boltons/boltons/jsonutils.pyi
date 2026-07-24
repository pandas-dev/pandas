from collections.abc import Generator
from typing import IO, Any, overload
from typing_extensions import Self

@overload
def reverse_iter_lines(
    file_obj: IO[bytes], blocksize: int = 4096, preseek: bool = True, encoding: None = None
) -> Generator[bytes]: ...
@overload
def reverse_iter_lines(file_obj: IO[str], blocksize: int = 4096, preseek: bool = True, *, encoding: str) -> Generator[str]: ...
@overload
def reverse_iter_lines(file_obj: IO[str], blocksize: int, preseek: bool, encoding: str) -> Generator[str]: ...

class JSONLIterator:
    ignore_errors: bool
    def __init__(
        self, file_obj: IO[str], ignore_errors: bool = False, reverse: bool = False, rel_seek: float | None = None
    ) -> None: ...
    @property
    def cur_byte_pos(self) -> int: ...
    def __iter__(self) -> Self: ...
    def next(self) -> Any: ...
    __next__ = next

__all__ = ["JSONLIterator", "reverse_iter_lines"]
