"""
Runtime-protocols for the `copy` standard library.
https://docs.python.org/3/library/copy.html
"""

import sys
from typing import Any, Protocol, Self, TypeVar

if sys.version_info >= (3, 13):
    from typing import override, runtime_checkable
else:
    from typing_extensions import override, runtime_checkable

__all__ = (
    "CanCopy",
    "CanCopySelf",
    "CanDeepcopy",
    "CanDeepcopySelf",
    "CanReplace",
    "CanReplaceSelf",
)


def __dir__() -> tuple[str, ...]:
    return __all__


_T_co = TypeVar("_T_co", covariant=True)


@runtime_checkable
class CanCopy(Protocol[_T_co]):
    """Anything that can be used as `copy.copy(_: CanCopy[T]) -> T`."""

    def __copy__(self, /) -> _T_co: ...


@runtime_checkable  # https://github.com/astral-sh/ty/issues/1800
class CanCopySelf(CanCopy["CanCopySelf"], Protocol):  # ty:ignore[unsupported-base]
    """Runtime-checkable alias `CanCopySelf = CanCopy[Self]`."""

    @override
    def __copy__(self, /) -> Self: ...


@runtime_checkable
class CanDeepcopy(Protocol[_T_co]):
    """Anything that can be used as `copy.deepcopy(_: CanDeepcopy[T]) -> T`."""

    def __deepcopy__(self, memo: dict[int, object], /) -> _T_co: ...


@runtime_checkable  # https://github.com/astral-sh/ty/issues/1800
class CanDeepcopySelf(CanDeepcopy["CanDeepcopySelf"], Protocol):  # ty:ignore[unsupported-base]
    """Runtime-checkable alias `CanDeepcopySelf = CanDeepcopy[Self]`."""

    @override
    def __deepcopy__(self, memo: dict[int, object], /) -> Self: ...


@runtime_checkable
class CanReplace(Protocol[_T_co]):
    """
    Represents anything that will be accepted by the Python 3.13+ `copy.replace()`
    stdlib function. If `x: CanReplace[T, *Ts]` and `y = copy.replace(x, **kwargs)`,
    then `y: T` and `**kwargs: *_Ts`.

    See https://docs.python.org/3/library/copy.html#copy.replace for details.

    Note the documented `**kwargs` cannot be annotated, because typeshed does not
    always include them in the `__replace__` signatures of e.g. `datetime.date`.
    This LSP violation forces us to be definsive and use what in `Callable` is `...`,
    but which looks like `*_: Any, **changes: Any` here.
    """

    def __replace__(self, /, *_: Any, **changes: Any) -> _T_co: ...


@runtime_checkable  # https://github.com/astral-sh/ty/issues/1800
class CanReplaceSelf(CanReplace["CanReplaceSelf"], Protocol):  # ty:ignore[unsupported-base]
    """
    Runtime-checkable alias `CanReplaceSelf = CanReplace[Self]`, i.e.
    `copy.replace: def[S <: CanReplaceSelf](_: S) -> S`.
    """

    @override
    def __replace__(self, /, *_: Any, **changes: Any) -> Self: ...
