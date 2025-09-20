"""
Runtime-protocols for the `pickle` standard library.
https://docs.python.org/3/library/pickle.html
"""

import sys
from collections.abc import Callable, Iterable
from typing import Protocol, Self, SupportsIndex, TypeAlias

if sys.version_info >= (3, 13):
    from typing import ParamSpec, TypeVar, override, runtime_checkable
else:
    from typing_extensions import ParamSpec, TypeVar, override, runtime_checkable


__all__ = (
    "CanGetnewargs",
    "CanGetnewargsEx",
    "CanGetstate",
    "CanReduce",
    "CanReduceEx",
    "CanSetstate",
)


def __dir__() -> tuple[str, ...]:
    return __all__


###


_Tss = ParamSpec("_Tss", default=...)
_ArgT = TypeVar("_ArgT", default=object)
_ArgT_co = TypeVar("_ArgT_co", covariant=True, default=object)
_StateT_co = TypeVar("_StateT_co", covariant=True)
_StateT_contra = TypeVar("_StateT_contra", contravariant=True)

_Tuple: TypeAlias = tuple[object, ...]
_Iter1: TypeAlias = Iterable[object]
_Iter2: TypeAlias = Iterable[tuple[object, object]]
_Callable: TypeAlias = Callable[_Tss, object]

_ReduceValue: TypeAlias = (
    str
    | tuple[_Callable, _Tuple]
    | tuple[_Callable, _Tuple, object]
    | tuple[_Callable, _Tuple, object, _Iter1 | None]
    | tuple[_Callable, _Tuple, object, _Iter1 | None, _Iter2 | None]
    | tuple[
        _Callable,
        _Tuple,
        object,
        _Iter1 | None,
        _Iter2 | None,
        _Callable[[object, object]] | None,
    ]
)
_RT_co = TypeVar("_RT_co", bound=_ReduceValue, default=_ReduceValue, covariant=True)


###


@runtime_checkable
class CanReduce(Protocol[_RT_co]):
    """https://docs.python.org/3/library/pickle.html#object.__reduce__"""

    @override
    def __reduce__(self, /) -> _RT_co: ...


@runtime_checkable
class CanReduceEx(Protocol[_RT_co]):
    """https://docs.python.org/3/library/pickle.html#object.__reduce_ex__"""

    @override
    def __reduce_ex__(self, v: SupportsIndex, /) -> _RT_co: ...


@runtime_checkable
class CanGetstate(Protocol[_StateT_co]):
    """
    https://docs.python.org/3/library/pickle.html#object.__getstate__
    """

    @override
    def __getstate__(self, /) -> _StateT_co: ...


@runtime_checkable
class CanSetstate(Protocol[_StateT_contra]):
    """https://docs.python.org/3/library/pickle.html#object.__setstate__"""

    def __setstate__(self, state: _StateT_contra, /) -> None: ...


@runtime_checkable
class CanGetnewargs(Protocol[_ArgT_co]):
    """https://docs.python.org/3/library/pickle.html#object.__getnewargs__"""

    def __new__(cls, /, *args: _ArgT_co) -> Self: ...
    def __getnewargs__(self, /) -> tuple[_ArgT_co, ...]: ...


@runtime_checkable
class CanGetnewargsEx(Protocol[_ArgT_co, _ArgT]):
    """https://docs.python.org/3/library/pickle.html#object.__getnewargs_ex__"""

    def __new__(cls, /, *args: _ArgT_co, **kwargs: _ArgT) -> Self: ...
    def __getnewargs_ex__(self, /) -> tuple[tuple[_ArgT_co, ...], dict[str, _ArgT]]: ...
