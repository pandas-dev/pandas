import datetime as dt
import sys
from typing import (
    Any,
    Generic,
    Protocol,
    Self,
    TypeAlias,
    _ProtocolMeta,  # noqa: PLC2701
    final,
)

if sys.version_info >= (3, 13):
    from typing import TypeIs, TypeVar, override, runtime_checkable
else:
    from typing_extensions import TypeIs, TypeVar, override, runtime_checkable

from ._can import CanFloat, CanIndex

__all__ = [
    "Just",
    "JustBytes",
    "JustComplex",
    "JustDate",
    "JustFloat",
    "JustInt",
    "JustObject",
]


def __dir__() -> list[str]:
    return __all__


###


_T = TypeVar("_T")
_TypeT = TypeVar("_TypeT", bound=type)
_ObjectT = TypeVar("_ObjectT", default=object)

_CanFloatOrIndex: TypeAlias = CanFloat | CanIndex


###


# NOTE: Both mypy and pyright incorrectly report LSP violations in `@final` protocols,
# even though these are purely structural, and therefore the LSP does not apply.

# mypy: disable-error-code="override, explicit-override"
# pyright: reportIncompatibleMethodOverride=false


@final  # https://github.com/python/mypy/issues/17288
class Just(Protocol[_T]):  # type: ignore[misc]
    """
    An runtime-checkable invariant type "wrapper", where `Just[T]` only accepts
    instances of `T`, and but rejects instances of any strict subtypes of `T`.

    Note that e.g. `Literal[""]` and `LiteralString` are not a strict `str` subtypes,
    and are therefore assignable to `Just[str]`, but instances of `class S(str): ...`
    are **not** assignable to `Just[str]`.

    Note:
        This protocol is not marked with `@runtime_checkable`, because using
        `isinstance(value, Just)` would be `True` for any type or object, and has no
        semantic meaning.
    """

    @property
    @override
    def __class__(self, /) -> type[_T]: ...  # pyrefly: ignore[bad-override]
    @__class__.setter
    def __class__(self, t: type[_T], /) -> None: ...


@final
class _JustMeta(_ProtocolMeta, Generic[_ObjectT]):
    __just_class__: type[_ObjectT]  # pyright: ignore[reportUninitializedInstanceVariable]

    def __new__(  # noqa: PYI019
        mcls: type[_TypeT],
        /,
        *args: *tuple[str, tuple[type, ...], dict[str, Any]],
        just: type[_ObjectT],
    ) -> _TypeT:
        self = super().__new__(mcls, *args)  # type: ignore[misc]
        self.__just_class__ = just  # pyrefly: ignore[missing-attribute]
        return self

    @override
    def __instancecheck__(self, instance: object) -> TypeIs[_ObjectT]:
        return self.__subclasscheck__(type(instance))

    @override
    def __subclasscheck__(self, subclass: Any) -> TypeIs[type[_ObjectT]]:
        tp = self.__just_class__

        if isinstance(subclass, int):
            # "bare" bool and int literals
            subclass = type(subclass)

        if not isinstance(subclass, type):
            # unwrap subscripted generics, with special-casing for `Just[...]`
            from types import GenericAlias  # noqa: I001, PLC0415
            from optype.types._typeforms import GenericType  # noqa: PLC0415

            while type(subclass) in {GenericType, GenericAlias}:
                origin = subclass.__origin__
                if origin is Just:
                    (subclass,) = subclass.__args__
                else:
                    subclass = origin

        if not isinstance(subclass, type):
            raise TypeError("issubclass() arg 1 must be a class") from None

        return subclass is tp


@runtime_checkable
@final  # https://github.com/python/mypy/issues/17288
class JustBytes(Protocol, metaclass=_JustMeta, just=bytes):  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]
    """
    A runtime checkable `Just[bytes]`, that also works on `pyright<1.1.390`.

    Useful as workaround for `mypy`'s `bytes` promotion (which can be disabled with the
    undocumented `--disable-bytearray-promotion` and  `--disable-memoryview-promotion`
    flags). See https://github.com/python/mypy/issues/15313s for more info.
    Note that this workaround requires `mypy >=1.15` or the `--disable-*-promotion`
    flags to work.
    """

    @property
    @override
    def __class__(self, /) -> type[bytes]: ...  # pyrefly: ignore[bad-override]
    @__class__.setter
    def __class__(self, t: type[bytes], /) -> None: ...


@runtime_checkable
@final  # https://github.com/python/mypy/issues/17288
class JustInt(Protocol, metaclass=_JustMeta, just=int):  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]
    """
    A runtime-checkable `Just[int]` that's compatible with `pyright<1.1.390`.

    Example:
        This example shows a situation you want to accept instances of just `int`,
        and reject and other instances, even they're subtypes of `int` such as `bool`,

        ```python
        def f(x: int, /) -> int:
            assert type(x) is int
            return x


        f(1337)  # accepted
        f(True)  # accepted :(
        ```

        But because `bool <: int`, booleans are also assignable to `int`, which we
        don't want to allow.

        Previously, it was considered impossible to achieve this using just python
        typing, and it was told that a PEP would be necessary to make it work.

        But this is not the case at all, and `optype` provides a clean counter-example
        that works with pyright and (even) mypy:

        ```python
        import optype as op


        def f_fixed(x: op.JustInt, /) -> int:
            assert type(x) is int
            return x  # accepted


        f_fixed(1337)  # accepted
        f_fixed(True)  # rejected :)
        ```

    Note:
        On `mypy` this behaves in the same way as `Just[int]`, which accidentally
        already works because of a mypy bug.
    """

    @property
    @override
    def __class__(self, /) -> type[int]: ...  # pyrefly: ignore[bad-override]
    @__class__.setter
    def __class__(self, t: type[int], /) -> None: ...

    # workaround for `pyright<1.1.390` and `basedpyright<1.22.1`
    def __new__(cls, x: str | bytes | bytearray, /, base: CanIndex) -> Self: ...


@runtime_checkable
@final  # https://github.com/python/mypy/issues/17288
class JustFloat(Protocol, metaclass=_JustMeta, just=float):  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]
    """
    A runtime checkable `Just[float]` that also works on `pyright<1.1.390`.

    Warning:
        Unlike `JustInt`, this does not work on `mypy<1.41.2`.
    """

    @property
    @override
    def __class__(self, /) -> type[float]: ...  # pyrefly: ignore[bad-override]
    @__class__.setter
    def __class__(self, t: type[float], /) -> None: ...

    # workaround for `pyright<1.1.390` and `basedpyright<1.22.1`
    def hex(self, /) -> str: ...


@runtime_checkable
@final  # https://github.com/python/mypy/issues/17288
class JustComplex(Protocol, metaclass=_JustMeta, just=complex):  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]
    """
    A runtime checkable `Just[complex]`, that also works on `pyright<1.1.390`.

    Warning:
        Unlike `JustInt`, this does not work on `mypy<1.41.2`.
    """

    @property
    @override
    def __class__(self, /) -> type[complex]: ...  # pyrefly: ignore[bad-override]
    @__class__.setter
    def __class__(self, t: type[complex], /) -> None: ...

    # workaround for `pyright<1.1.390` and `basedpyright<1.22.1`
    def __new__(cls, /, real: _CanFloatOrIndex, imag: _CanFloatOrIndex) -> Self: ...


@runtime_checkable
@final  # https://github.com/python/mypy/issues/17288
class JustDate(Protocol, metaclass=_JustMeta, just=dt.date):  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]
    """
    A runtime checkable `Just[datetime.Date]`, that rejects `datetime.datetime`.
    """

    @property
    @override
    def __class__(self, /) -> type[dt.date]: ...  # pyrefly: ignore[bad-override]
    @__class__.setter
    def __class__(self, t: type[dt.date], /) -> None: ...


@runtime_checkable
@final  # https://github.com/python/mypy/issues/17288
class JustObject(Protocol, metaclass=_JustMeta, just=object):  # type: ignore[misc]  # pyright: ignore[reportGeneralTypeIssues]
    """
    A runtime checkable `Just[object]`, that also works on `pyright<1.1.390`.

    Useful for typing `object()` sentinels.
    """

    @property
    @override
    def __class__(self, /) -> type[object]: ...
    @__class__.setter
    def __class__(self, t: type[object], /) -> None: ...
