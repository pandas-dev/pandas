from collections.abc import Callable, Iterable
from typing import Any, TypeVar, overload

_T = TypeVar("_T")
_S = TypeVar("_S")

__license__: str
__title__: str

@overload
def first(iterable: Iterable[_T]) -> _T | None: ...
@overload
def first(iterable: Iterable[_T], default: _S) -> _T | _S: ...
@overload
def first(iterable: Iterable[_T], default: _S, key: Callable[[_T], Any] | None) -> _T | _S: ...
@overload
def first(iterable: Iterable[_T], default: None, key: Callable[[_T], Any] | None) -> _T | None: ...
@overload
def first(iterable: Iterable[_T], *, key: Callable[[_T], Any] | None) -> _T | None: ...
