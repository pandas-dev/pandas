"""
Runtime-protocols for the `copy` standard library.
https://docs.python.org/3/library/copy.html
"""

import sys
from typing import Any, Protocol, Self, override

if sys.version_info >= (3, 13):
    from typing import runtime_checkable
else:
    from typing_extensions import runtime_checkable

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


@runtime_checkable
class CanCopy[T_co](Protocol):
    """Anything that can be used as `copy.copy(_: CanCopy[T]) -> T`."""

    def __copy__(self, /) -> T_co: ...


@runtime_checkable  # https://github.com/astral-sh/ty/issues/1800
class CanCopySelf(CanCopy["CanCopySelf"], Protocol):
    """Runtime-checkable alias `CanCopySelf = CanCopy[Self]`."""

    @override
    def __copy__(self, /) -> Self: ...


@runtime_checkable
class CanDeepcopy[T_co](Protocol):
    """Anything that can be used as `copy.deepcopy(_: CanDeepcopy[T]) -> T`."""

    def __deepcopy__(self, memo: dict[int, object], /) -> T_co: ...


@runtime_checkable  # https://github.com/astral-sh/ty/issues/1800
class CanDeepcopySelf(CanDeepcopy["CanDeepcopySelf"], Protocol):
    """Runtime-checkable alias `CanDeepcopySelf = CanDeepcopy[Self]`."""

    @override
    def __deepcopy__(self, memo: dict[int, object], /) -> Self: ...


@runtime_checkable
class CanReplace[T_co](Protocol):
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

    def __replace__(self, /, *_: Any, **changes: Any) -> T_co: ...


@runtime_checkable  # https://github.com/astral-sh/ty/issues/1800
class CanReplaceSelf(CanReplace["CanReplaceSelf"], Protocol):
    """
    Runtime-checkable alias `CanReplaceSelf = CanReplace[Self]`, i.e.
    `copy.replace: def[S <: CanReplaceSelf](_: S) -> S`.
    """

    @override
    def __replace__(self, /, *_: Any, **changes: Any) -> Self: ...
