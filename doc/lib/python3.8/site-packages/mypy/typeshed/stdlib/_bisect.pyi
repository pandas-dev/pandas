import sys
from _typeshed import SupportsLessThan
from typing import Callable, MutableSequence, Optional, Sequence, TypeVar

_T = TypeVar("_T")

if sys.version_info >= (3, 10):
    def bisect_left(
        a: Sequence[_T], x: _T, lo: int = ..., hi: Optional[int] = ..., *, key: Optional[Callable[[_T], SupportsLessThan]] = ...
    ) -> int: ...
    def bisect_right(
        a: Sequence[_T], x: _T, lo: int = ..., hi: Optional[int] = ..., *, key: Optional[Callable[[_T], SupportsLessThan]] = ...
    ) -> int: ...
    def insort_left(
        a: MutableSequence[_T],
        x: _T,
        lo: int = ...,
        hi: Optional[int] = ...,
        *,
        key: Optional[Callable[[_T], SupportsLessThan]] = ...,
    ) -> None: ...
    def insort_right(
        a: MutableSequence[_T],
        x: _T,
        lo: int = ...,
        hi: Optional[int] = ...,
        *,
        key: Optional[Callable[[_T], SupportsLessThan]] = ...,
    ) -> None: ...

else:
    def bisect_left(a: Sequence[_T], x: _T, lo: int = ..., hi: Optional[int] = ...) -> int: ...
    def bisect_right(a: Sequence[_T], x: _T, lo: int = ..., hi: Optional[int] = ...) -> int: ...
    def insort_left(a: MutableSequence[_T], x: _T, lo: int = ..., hi: Optional[int] = ...) -> None: ...
    def insort_right(a: MutableSequence[_T], x: _T, lo: int = ..., hi: Optional[int] = ...) -> None: ...
