import sys
import types
from collections.abc import Generator, Iterable
from typing import Any, Never, Protocol, Self, SupportsIndex, TypeAlias, overload

if sys.version_info >= (3, 13):
    from typing import ParamSpec, TypeVar, override, runtime_checkable
else:
    from typing_extensions import ParamSpec, TypeVar, override, runtime_checkable

__all__ = [
    "CanAEnter",
    "CanAEnterSelf",
    "CanAExit",
    "CanAIter",
    "CanAIterSelf",
    "CanANext",
    "CanAbs",
    "CanAbsSelf",
    "CanAdd",
    "CanAddSame",
    "CanAddSelf",
    "CanAnd",
    "CanAndSame",
    "CanAndSelf",
    "CanAsyncWith",
    "CanAsyncWithSelf",
    "CanAwait",
    "CanBool",
    "CanBuffer",
    "CanBytes",
    "CanCall",
    "CanCeil",
    "CanComplex",
    "CanContains",
    "CanDelattr",
    "CanDelete",
    "CanDelitem",
    "CanDir",
    "CanDivmod",
    "CanEnter",
    "CanEnterSelf",
    "CanEq",
    "CanExit",
    "CanFloat",
    "CanFloor",
    "CanFloordiv",
    "CanFloordivSame",
    "CanFloordivSelf",
    "CanFormat",
    "CanGe",
    "CanGet",
    "CanGetMissing",
    "CanGetattr",
    "CanGetattribute",
    "CanGetitem",
    "CanGt",
    "CanHash",
    "CanIAdd",
    "CanIAddSame",
    "CanIAddSelf",
    "CanIAnd",
    "CanIAndSame",
    "CanIAndSelf",
    "CanIFloordiv",
    "CanIFloordivSame",
    "CanIFloordivSelf",
    "CanILshift",
    "CanILshiftSame",
    "CanILshiftSelf",
    "CanIMatmul",
    "CanIMatmulSame",
    "CanIMatmulSelf",
    "CanIMod",
    "CanIModSame",
    "CanIModSelf",
    "CanIMul",
    "CanIMulSame",
    "CanIMulSelf",
    "CanIOr",
    "CanIOrSame",
    "CanIOrSelf",
    "CanIPow",
    "CanIPowSame",
    "CanIPowSelf",
    "CanIRshift",
    "CanIRshiftSame",
    "CanIRshiftSelf",
    "CanISub",
    "CanISubSame",
    "CanISubSelf",
    "CanITruediv",
    "CanITruedivSame",
    "CanITruedivSelf",
    "CanIXor",
    "CanIXorSame",
    "CanIXorSelf",
    "CanIndex",
    "CanInt",
    "CanInvert",
    "CanInvertSelf",
    "CanIter",
    "CanIterSelf",
    "CanLe",
    "CanLen",
    "CanLengthHint",
    "CanLshift",
    "CanLshiftSame",
    "CanLshiftSelf",
    "CanLt",
    "CanMatmul",
    "CanMatmulSame",
    "CanMatmulSelf",
    "CanMissing",
    "CanMod",
    "CanModSame",
    "CanModSelf",
    "CanMul",
    "CanMulSame",
    "CanMulSelf",
    "CanNe",
    "CanNeg",
    "CanNegSelf",
    "CanNext",
    "CanOr",
    "CanOrSame",
    "CanOrSelf",
    "CanPos",
    "CanPosSelf",
    "CanPow",
    "CanPow2",
    "CanPow3",
    "CanPowSame",
    "CanPowSelf",
    "CanRAdd",
    "CanRAddSelf",
    "CanRAnd",
    "CanRAndSelf",
    "CanRDivmod",
    "CanRFloordiv",
    "CanRFloordivSelf",
    "CanRLshift",
    "CanRLshiftSelf",
    "CanRMatmul",
    "CanRMatmulSelf",
    "CanRMod",
    "CanRModSelf",
    "CanRMul",
    "CanRMulSelf",
    "CanROr",
    "CanROrSelf",
    "CanRPow",
    "CanRPowSelf",
    "CanRRshift",
    "CanRRshiftSelf",
    "CanRSub",
    "CanRSubSelf",
    "CanRTruediv",
    "CanRTruedivSelf",
    "CanRXor",
    "CanRXorSelf",
    "CanReleaseBuffer",
    "CanRepr",
    "CanReversed",
    "CanRound",
    "CanRound1",
    "CanRound2",
    "CanRshift",
    "CanRshiftSame",
    "CanRshiftSelf",
    "CanSequence",
    "CanSet",
    "CanSetName",
    "CanSetattr",
    "CanSetitem",
    "CanStr",
    "CanSub",
    "CanSubSame",
    "CanSubSelf",
    "CanTruediv",
    "CanTruedivSame",
    "CanTruedivSelf",
    "CanTrunc",
    "CanWith",
    "CanWithSelf",
    "CanXor",
    "CanXorSame",
    "CanXorSelf",
]


def __dir__() -> list[str]:
    return __all__


###


# this should (but can't) be contravariant
# https://github.com/astral-sh/ty/issues/1798
_Tss = ParamSpec("_Tss", default=...)  # ty:ignore[invalid-paramspec]

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)
_T_contra = TypeVar("_T_contra", contravariant=True)
_TT_co = TypeVar("_TT_co", covariant=True, default=_T_contra)

_K_contra = TypeVar("_K_contra", contravariant=True)
_V_contra = TypeVar("_V_contra", contravariant=True)
_V_co = TypeVar("_V_co", covariant=True)
_VV_co = TypeVar("_VV_co", default=_V_co, covariant=True)

_BoolT_co = TypeVar("_BoolT_co", default=bool, covariant=True)
_ExcT = TypeVar("_ExcT", bound=BaseException)

_T_object_contra = TypeVar("_T_object_contra", contravariant=True, default=object)
_T_object_co = TypeVar("_T_object_co", covariant=True, default=object)
_T_bool_co = TypeVar("_T_bool_co", default=bool, covariant=True)
_T_int_contra = TypeVar("_T_int_contra", default=int, contravariant=True)
_T_int_co = TypeVar("_T_int_co", default=int, covariant=True)
_T_float_co = TypeVar("_T_float_co", default=float, covariant=True)
_T_None_co = TypeVar("_T_None_co", default=None, covariant=True)
_T_Never_contra = TypeVar("_T_Never_contra", default=Never, contravariant=True)
_T_Never_co = TypeVar("_T_Never_co", default=Never, covariant=True)

_IterT_str_co = TypeVar(
    "_IterT_str_co",
    # we use `Any` instead of `object` to side-step LSP errors in `CanDir`
    bound=Iterable[Any],
    default=Iterable[str],
    covariant=True,
)

# we can't use `CanIndex` here, because of a recent regression in pyright 1.1.392
_IndexT_contra = TypeVar(
    "_IndexT_contra",
    bound=SupportsIndex | slice,  # pyrefly: ignore[invalid-annotation]
    contravariant=True,
)

# return type that is usually `None`, but can be anything, as it is ignored at runtime
_Ignored: TypeAlias = object
# This should be `asyncio.Future[typing.Any] | None`. But that would make this
# incompatible with `Awaitable` -- it (annoyingly) uses `Any`:
# https://github.com/python/typeshed/blob/587ad6/stdlib/asyncio/futures.pyi#L51
_FutureOrNone: TypeAlias = object
_AsyncGen: TypeAlias = Generator[_FutureOrNone, None, _T]


###

# Type conversion


@runtime_checkable
class CanBool(Protocol[_BoolT_co]):
    def __bool__(self, /) -> _BoolT_co: ...


@runtime_checkable
class CanInt(Protocol):
    def __int__(self, /) -> int: ...


@runtime_checkable
class CanFloat(Protocol):
    def __float__(self, /) -> float: ...


@runtime_checkable
class CanComplex(Protocol):
    def __complex__(self, /) -> complex: ...


@runtime_checkable
class CanBytes(Protocol):
    def __bytes__(self, /) -> bytes: ...


@runtime_checkable
class CanStr(Protocol):
    @override
    def __str__(self, /) -> str: ...


# Object representation


@runtime_checkable
class CanHash(Protocol):
    @override
    def __hash__(self, /) -> int: ...


@runtime_checkable
class CanIndex(Protocol):
    def __index__(self, /) -> int: ...


@runtime_checkable
class CanRepr(Protocol):
    @override
    def __repr__(self, /) -> str: ...


@runtime_checkable
class CanFormat(Protocol):
    @override
    def __format__(self, format_spec: str, /) -> str: ...


# Iteration


@runtime_checkable
class CanNext(Protocol[_V_co]):
    """
    Similar to `Iterator[V]`, but without the requirement to
    also have a `__iter__` method, which isn't needed in most cases (at least
    not in cpython).
    """

    def __next__(self, /) -> _V_co: ...


_CanNextT_co = TypeVar("_CanNextT_co", bound=CanNext[object], covariant=True)


@runtime_checkable
class CanIter(Protocol[_CanNextT_co]):
    """Like `Iterable[V]`, but with a flexible return type."""

    def __iter__(self, /) -> _CanNextT_co: ...


@runtime_checkable
class CanIterSelf(CanNext[_V_co], CanIter[CanNext[_V_co]], Protocol[_V_co]):
    """Like `Iterator[V]`, but without the `abc` nonsense."""

    @override
    def __iter__(self, /) -> Self: ...


# Async Iteration


@runtime_checkable
class CanANext(Protocol[_V_co]):
    def __anext__(self, /) -> _V_co: ...


_CanANextT_co = TypeVar("_CanANextT_co", bound=CanANext[object], covariant=True)


@runtime_checkable
class CanAIter(Protocol[_CanANextT_co]):
    def __aiter__(self, /) -> _CanANextT_co: ...


@runtime_checkable
class CanAIterSelf(CanAIter["CanAIterSelf[_V_co]"], CanANext[_V_co], Protocol[_V_co]):  # ty:ignore[unsupported-base]
    """Like `AsyncIterator[T]`, but without the `abc` nonsense."""

    @override
    def __aiter__(self, /) -> Self: ...


# Rich comparison operands


@runtime_checkable
class CanEq(Protocol[_T_object_contra, _T_bool_co]):  # noqa: PLW1641
    """
    Unfortunately, `typeshed` (incorrectly) annotates `object.__eq__` as
    `(Self, object) -> bool`.
    As a counter-example, consider `numpy.ndarray`. It's `__eq__` method
    returns a boolean (mask) array of the same shape as the input array.
    Moreover, `numpy.ndarray` doesn't even implement `CanBool` (`bool()`
    raises a `TypeError` for shapes of size > 1).
    There is nothing wrong with this implementation, even though `typeshed`
    (incorrectly) won't allow it (because `numpy.ndarray <: object`).

    So in reality, it should be `__eq__: (Self, X, /) -> Y`, with `X` unbounded
    and *contra*variant, and `+Y` unbounded and *co*variant.
    """

    @override
    def __eq__(self, rhs: _T_object_contra, /) -> _T_bool_co: ...  # type: ignore[override]  # pyright:ignore[reportIncompatibleMethodOverride]  # pyrefly: ignore[bad-override]


@runtime_checkable
class CanNe(Protocol[_T_object_contra, _T_bool_co]):
    """
    Just like `__eq__`, the `__ne__` method is incorrectly annotated in
    typeshed. Refer to `CanEq` for why this is.
    """

    @override
    def __ne__(self, rhs: _T_object_contra, /) -> _T_bool_co: ...  # type: ignore[override]  # pyright:ignore[reportIncompatibleMethodOverride]  # pyrefly: ignore[bad-override]


@runtime_checkable
class CanLt(Protocol[_T_object_contra, _T_bool_co]):
    def __lt__(self, rhs: _T_object_contra, /) -> _T_bool_co: ...


@runtime_checkable
class CanLe(Protocol[_T_object_contra, _T_bool_co]):
    def __le__(self, rhs: _T_object_contra, /) -> _T_bool_co: ...


@runtime_checkable
class CanGt(Protocol[_T_object_contra, _T_bool_co]):
    def __gt__(self, rhs: _T_object_contra, /) -> _T_bool_co: ...


@runtime_checkable
class CanGe(Protocol[_T_object_contra, _T_bool_co]):
    def __ge__(self, rhs: _T_object_contra, /) -> _T_bool_co: ...


# Callables


@runtime_checkable
class CanCall(Protocol[_Tss, _T_object_co]):
    """Equivalent to `typing.Callable` in theory, but not always in practice."""

    def __call__(self, /, *args: _Tss.args, **kwargs: _Tss.kwargs) -> _T_object_co: ...


# Dynamic attribute access


@runtime_checkable
class CanGetattr(Protocol[_T_object_co]):
    def __getattr__(self, name: str, /) -> _T_object_co: ...


@runtime_checkable
class CanGetattribute(Protocol[_T_object_co]):
    """Note that `isinstance(x, CanGetattribute)` is always `True`."""

    @override
    def __getattribute__(self, name: str, /) -> _T_object_co: ...


@runtime_checkable
class CanSetattr(Protocol[_T_object_contra]):
    """Note that `isinstance(x, CanSetattr)` is always true."""

    @override
    def __setattr__(self, name: str, value: _T_object_contra, /) -> _Ignored: ...  # type: ignore[misc, override]  # pyright: ignore[reportIncompatibleMethodOverride]


@runtime_checkable
class CanDelattr(Protocol):
    @override
    def __delattr__(self, name: str, /) -> _Ignored: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # pyrefly: ignore[bad-override]


@runtime_checkable
class CanDir(Protocol[_IterT_str_co]):
    @override
    def __dir__(self, /) -> _IterT_str_co: ...


# Descriptors


@runtime_checkable
class CanGet(Protocol[_T_contra, _V_co, _VV_co]):
    @overload
    def __get__(self, obj: _T_contra, cls: type | None = ..., /) -> _V_co: ...
    @overload
    def __get__(self, obj: None, cls: type[_T_contra], /) -> _VV_co: ...


@runtime_checkable
class CanGetSelf(Protocol[_T_contra, _V_co]):
    """CanGetSelf[-T, +V] = CanGet[T, V, Self]"""

    @overload
    def __get__(self, obj: _T_contra, cls: type | None = ..., /) -> _V_co: ...
    @overload
    def __get__(self, obj: None, cls: type[_T_contra], /) -> Self: ...


@runtime_checkable
class CanSet(Protocol[_T_contra, _V_contra]):
    def __set__(self, owner: _T_contra, value: _V_contra, /) -> _Ignored: ...


@runtime_checkable
class CanDelete(Protocol[_T_contra]):
    def __delete__(self, owner: _T_contra, /) -> _Ignored: ...


@runtime_checkable
class CanSetName(Protocol[_T_contra]):
    def __set_name__(self, cls: type[_T_contra], name: str, /) -> _Ignored: ...


# Collection type operands.


@runtime_checkable
class CanLen(Protocol):
    def __len__(self, /) -> int: ...


@runtime_checkable
class CanLengthHint(Protocol):
    def __length_hint__(self, /) -> int: ...


@runtime_checkable
class CanGetitem(Protocol[_K_contra, _V_co]):
    def __getitem__(self, key: _K_contra, /) -> _V_co: ...


@runtime_checkable
class CanSetitem(Protocol[_K_contra, _V_contra]):
    def __setitem__(self, key: _K_contra, value: _V_contra, /) -> _Ignored: ...


@runtime_checkable
class CanDelitem(Protocol[_K_contra]):
    def __delitem__(self, key: _K_contra, /) -> None: ...


@runtime_checkable
class CanReversed(Protocol[_T_co]):
    # `builtin.reversed` can return anything, but in practice it's always
    # something that can be iterated over (e.g. iterable or sequence-like)
    def __reversed__(self, /) -> _T_co: ...


@runtime_checkable
class CanContains(Protocol[_T_object_contra]):
    # usually the key is required to also be a hashable object, but this
    # isn't strictly required
    def __contains__(self, key: _T_object_contra, /) -> bool: ...


@runtime_checkable
class CanMissing(Protocol[_K_contra, _V_co]):
    def __missing__(self, key: _K_contra, /) -> _V_co: ...


@runtime_checkable
class CanGetMissing(
    CanGetitem[_K_contra, _V_co],
    CanMissing[_K_contra, _V_co],
    Protocol[_K_contra, _V_co, _VV_co],
): ...


@runtime_checkable
class CanSequence(
    CanGetitem[_IndexT_contra, _V_co],
    CanLen,
    Protocol[_IndexT_contra, _V_co],
):
    """
    A sequence is an object with a __len__ method and a
    __getitem__ method that takes int(-like) argument as key (the index).
    Additionally, it is expected to be 0-indexed (the first element is at
    index 0) and "dense" (i.e. the indices are consecutive integers, and are
    obtainable with e.g. `range(len(_))`).
    """


###
# Arithmetic operands

# __add__


@runtime_checkable
class CanAdd(Protocol[_T_contra, _TT_co]):
    def __add__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanAddSelf(Protocol[_T_contra]):
    """CanAddSelf[-T] = CanAdd[T, Self]"""

    def __add__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanAddSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanAddSame[-T = Never, +R = Never] = CanAdd[Self | T, Self | R]"""

    def __add__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __sub__


@runtime_checkable
class CanSub(Protocol[_T_contra, _TT_co]):
    def __sub__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanSubSelf(Protocol[_T_contra]):
    """CanSubSelf[-T] = CanSub[T, Self]"""

    def __sub__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanSubSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanSubSame[-T = Never, +R = Never] = CanSub[Self | T, Self | R]"""

    def __sub__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __mul__


@runtime_checkable
class CanMul(Protocol[_T_contra, _TT_co]):
    def __mul__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanMulSelf(Protocol[_T_contra]):
    """CanMulSelf[-T] = CanMul[T, Self]"""

    def __mul__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanMulSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanMulSame[-T = Never, +R = Never] = CanMul[Self | T, Self | R]"""

    def __mul__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __matmul__


@runtime_checkable
class CanMatmul(Protocol[_T_contra, _TT_co]):
    def __matmul__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanMatmulSelf(Protocol[_T_contra]):
    """CanMatmul[-T, Self]"""

    def __matmul__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanMatmulSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanMatmulSame[-T = Never, +R = Never] = CanMatmul[Self | T, Self | R]"""

    def __matmul__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __truediv__


@runtime_checkable
class CanTruediv(Protocol[_T_contra, _TT_co]):
    def __truediv__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanTruedivSelf(Protocol[_T_contra]):
    """CanTruedivSelf[-T] = CanTruediv[T, Self]"""

    def __truediv__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanTruedivSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanTruedivSame[-T = Never, +R = Never] = CanTruediv[Self | T, Self | R]"""

    def __truediv__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __floordiv__


@runtime_checkable
class CanFloordiv(Protocol[_T_contra, _TT_co]):
    def __floordiv__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanFloordivSelf(Protocol[_T_contra]):
    """CanFloordivSelf[-T] = CanFloordiv[T, Self]"""

    def __floordiv__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanFloordivSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanFloordivSame[-T = Never, +R = Never] = CanFloordiv[Self | T, Self | R]"""

    def __floordiv__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __mod__


@runtime_checkable
class CanMod(Protocol[_T_contra, _TT_co]):
    def __mod__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanModSelf(Protocol[_T_contra]):
    """CanModSelf[-T] = CanMod[T, Self]"""

    def __mod__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanModSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanModSame[-T = Never, +R = Never] = CanMod[Self | T, Self | R]"""

    def __mod__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __divmod__


@runtime_checkable
class CanDivmod(Protocol[_T_contra, _T_co]):
    def __divmod__(self, rhs: _T_contra, /) -> _T_co: ...


# __pow__


@runtime_checkable
class CanPow2(Protocol[_T_contra, _TT_co]):
    def __pow__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanPow3(Protocol[_T_contra, _V_contra, _T_int_co]):
    def __pow__(self, exp: _T_contra, mod: _V_contra, /) -> _T_int_co: ...


@runtime_checkable
class CanPow(
    CanPow2[_T_contra, _TT_co],
    CanPow3[_T_contra, _V_contra, _T_int_co],
    Protocol[_T_contra, _V_contra, _TT_co, _T_int_co],
):
    @overload
    @override
    def __pow__(self, exp: _T_contra, /) -> _TT_co: ...
    @overload
    def __pow__(self, exp: _T_contra, mod: _V_contra, /) -> _T_int_co: ...


@runtime_checkable
class CanPowSelf(Protocol[_T_contra]):
    """CanPowSelf[-T] = CanPow2[T, Self]"""

    def __pow__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanPowSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanPowSame[-T = Never, +R = Never] = CanPow2[Self | T, Self | R]"""

    def __pow__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __lshift__


@runtime_checkable
class CanLshift(Protocol[_T_contra, _TT_co]):
    def __lshift__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanLshiftSelf(Protocol[_T_contra]):
    """CanLshiftSelf[-T] = CanLshift[T, Self]"""

    def __lshift__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanLshiftSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanLshiftSame[-T = Never, +R = Never] = CanLshift[Self | T, Self | R]"""

    def __lshift__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __rshift__


@runtime_checkable
class CanRshift(Protocol[_T_contra, _TT_co]):
    def __rshift__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRshiftSelf(Protocol[_T_contra]):
    """CanRshiftSelf[-T] = CanRshift[T, Self]"""

    def __rshift__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanRshiftSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanRshiftSame[-T = Never, +R = Never] = CanRshift[Self | T, Self | R]"""

    def __rshift__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __and__


@runtime_checkable
class CanAnd(Protocol[_T_contra, _TT_co]):
    def __and__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanAndSelf(Protocol[_T_contra]):
    """CanAndSelf[-T] = CanAnd[T, Self]"""

    def __and__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanAndSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanAndSame[-T = Never, +R = Never] = CanAnd[Self | T, Self | R]"""

    def __and__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __xor__


@runtime_checkable
class CanXor(Protocol[_T_contra, _TT_co]):
    def __xor__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanXorSelf(Protocol[_T_contra]):
    """CanXorSelf[-T] = CanXor[T, Self]"""

    def __xor__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanXorSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanXorSame[-T = Never, +R = Never] = CanXor[Self | T, Self | R]"""

    def __xor__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


# __or__


@runtime_checkable
class CanOr(Protocol[_T_contra, _TT_co]):
    def __or__(self, rhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanOrSelf(Protocol[_T_contra]):
    """CanOrSelf[-T] = CanOr[T, Self]"""

    def __or__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanOrSame(Protocol[_T_Never_contra, _T_Never_co]):
    """CanOrSame[-T = Never, +R = Never] = CanOr[Self | T, Self | R]"""

    def __or__(self, rhs: Self | _T_Never_contra, /) -> Self | _T_Never_co: ...


###
# Reflected arithmetic operands

# __radd__


@runtime_checkable
class CanRAdd(Protocol[_T_contra, _TT_co]):
    def __radd__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRAddSelf(Protocol[_T_contra]):
    """CanRAddSelf[-T] = CanRAdd[T, Self]"""

    def __radd__(self, lhs: _T_contra, /) -> Self: ...


# __rsub__


@runtime_checkable
class CanRSub(Protocol[_T_contra, _TT_co]):
    def __rsub__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRSubSelf(Protocol[_T_contra]):
    """CanRSubSelf[-T] = CanRSub[T, Self]"""

    def __rsub__(self, lhs: _T_contra, /) -> Self: ...


# __rmul__


@runtime_checkable
class CanRMul(Protocol[_T_contra, _TT_co]):
    def __rmul__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRMulSelf(Protocol[_T_contra]):
    """CanRMulSelf[-T] = CanRMul[T, Self]"""

    def __rmul__(self, lhs: _T_contra, /) -> Self: ...


# __rmatmul__


@runtime_checkable
class CanRMatmul(Protocol[_T_contra, _TT_co]):
    def __rmatmul__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRMatmulSelf(Protocol[_T_contra]):
    """CanRMatmulSelf[-T] = CanRMatmul[T, Self]"""

    def __rmatmul__(self, lhs: _T_contra, /) -> Self: ...


# __rtruediv__


@runtime_checkable
class CanRTruediv(Protocol[_T_contra, _TT_co]):
    def __rtruediv__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRTruedivSelf(Protocol[_T_contra]):
    """CanRTruedivSelf[-T] = CanRTruediv[T, Self]"""

    def __rtruediv__(self, lhs: _T_contra, /) -> Self: ...


# __rfloordiv__


@runtime_checkable
class CanRFloordiv(Protocol[_T_contra, _TT_co]):
    def __rfloordiv__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRFloordivSelf(Protocol[_T_contra]):
    """CanRFloordivSelf[-T] = CanRFloordiv[T, Self]"""

    def __rfloordiv__(self, lhs: _T_contra, /) -> Self: ...


# __rmod__


@runtime_checkable
class CanRMod(Protocol[_T_contra, _TT_co]):
    def __rmod__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRModSelf(Protocol[_T_contra]):
    """CanRModSelf[-T] = CanRMod[T, Self]"""

    def __rmod__(self, lhs: _T_contra, /) -> Self: ...


# __rdivmod__


@runtime_checkable
class CanRDivmod(Protocol[_T_contra, _T_co]):
    # can return anything, but is almost always a 2-tuple
    def __rdivmod__(self, lhs: _T_contra, /) -> _T_co: ...


# __rpow__


@runtime_checkable
class CanRPow(Protocol[_T_contra, _TT_co]):
    def __rpow__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRPowSelf(Protocol[_T_contra]):
    """CanRPowSelf[-T] = CanRPow[T, Self]"""

    def __rpow__(self, lhs: _T_contra, /) -> Self: ...


# __rlshift__


@runtime_checkable
class CanRLshift(Protocol[_T_contra, _TT_co]):
    def __rlshift__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRLshiftSelf(Protocol[_T_contra]):
    """CanRLshiftSelf[-T] = CanRLshift[T, Self]"""

    def __rlshift__(self, lhs: _T_contra, /) -> Self: ...


# __rrshift__


@runtime_checkable
class CanRRshift(Protocol[_T_contra, _TT_co]):
    def __rrshift__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRRshiftSelf(Protocol[_T_contra]):
    """CanRRshiftSelf[-T] = CanRRshift[T, Self]"""

    def __rrshift__(self, lhs: _T_contra, /) -> Self: ...


# __rand__


@runtime_checkable
class CanRAnd(Protocol[_T_contra, _TT_co]):
    def __rand__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRAndSelf(Protocol[_T_contra]):
    """CanRAndSelf[-T] = CanRAnd[T, Self]"""

    def __rand__(self, lhs: _T_contra, /) -> Self: ...


# __rxor__


@runtime_checkable
class CanRXor(Protocol[_T_contra, _TT_co]):
    def __rxor__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanRXorSelf(Protocol[_T_contra]):
    """CanRXorSelf[-T] = CanRXor[T, Self]"""

    def __rxor__(self, lhs: _T_contra, /) -> Self: ...


# __ror__


@runtime_checkable
class CanROr(Protocol[_T_contra, _TT_co]):
    def __ror__(self, lhs: _T_contra, /) -> _TT_co: ...


@runtime_checkable
class CanROrSelf(Protocol[_T_contra]):
    """CanROrSelf[-T] = CanROr[T, Self]"""

    def __ror__(self, lhs: _T_contra, /) -> Self: ...


###
# Augmented arithmetic operands

# ruff: noqa: PYI034

# __iadd__


@runtime_checkable
class CanIAdd(Protocol[_T_contra, _T_co]):
    def __iadd__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIAddSelf(Protocol[_T_contra]):
    """CanIAddSelf[-T] = CanIAdd[T, Self]"""

    def __iadd__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIAddSame(Protocol[_T_Never_contra]):
    """CanIAddSame[-T = Never] = CanIAdd[Self | T, Self]"""

    def __iadd__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __isub__


@runtime_checkable
class CanISub(Protocol[_T_contra, _T_co]):
    def __isub__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanISubSelf(Protocol[_T_contra]):
    """CanISubSelf[-T] = CanISub[T, Self]"""

    def __isub__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanISubSame(Protocol[_T_Never_contra]):
    """CanISubSame[-T = Never] = CanISub[Self | T, Self]"""

    def __isub__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __imul__


@runtime_checkable
class CanIMul(Protocol[_T_contra, _T_co]):
    def __imul__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIMulSelf(Protocol[_T_contra]):
    """CanIMulSelf[-T] = CanIMul[T, Self]"""

    def __imul__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIMulSame(Protocol[_T_Never_contra]):
    """CanIMulSame[-T = Never] = CanIMul[Self | T, Self]"""

    def __imul__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __imatmul__


@runtime_checkable
class CanIMatmul(Protocol[_T_contra, _T_co]):
    def __imatmul__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIMatmulSelf(Protocol[_T_contra]):
    """CanIMatmulSelf[-T] = CanIMatmul[T, Self]"""

    def __imatmul__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIMatmulSame(Protocol[_T_Never_contra]):
    """CanIMatmulSame[-T = Never] = CanIMatmul[Self | T, Self]"""

    def __imatmul__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __itruediv__


@runtime_checkable
class CanITruediv(Protocol[_T_contra, _T_co]):
    def __itruediv__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanITruedivSelf(Protocol[_T_contra]):
    """CanITruedivSelf[-T] = CanITruediv[T, Self]"""

    def __itruediv__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanITruedivSame(Protocol[_T_Never_contra]):
    """CanITruedivSame[-T = Never] = CanITruediv[Self | T, Self]"""

    def __itruediv__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __ifloordiv__


@runtime_checkable
class CanIFloordiv(Protocol[_T_contra, _T_co]):
    def __ifloordiv__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIFloordivSelf(Protocol[_T_contra]):
    """CanIFloordivSelf[-T] = CanIFloordiv[T, Self]"""

    def __ifloordiv__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIFloordivSame(Protocol[_T_Never_contra]):
    """CanIFloordivSame[-T = Never] = CanIFloordiv[Self | T, Self]"""

    def __ifloordiv__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __imod__


@runtime_checkable
class CanIMod(Protocol[_T_contra, _T_co]):
    def __imod__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIModSelf(Protocol[_T_contra]):
    """CanIModSelf[-T] = CanIMod[T, Self]"""

    def __imod__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIModSame(Protocol[_T_Never_contra]):
    """CanIModSame[-T = Never] = CanIMod[Self | T, Self]"""

    def __imod__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __ipow__


@runtime_checkable
class CanIPow(Protocol[_T_contra, _T_co]):
    # no augmented pow/3 exists
    def __ipow__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIPowSelf(Protocol[_T_contra]):
    """CanIPowSelf[-T] = CanIPow[T, Self]"""

    def __ipow__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIPowSame(Protocol[_T_Never_contra]):
    """CanIPowSame[-T = Never] = CanIPow[Self | T, Self]"""

    def __ipow__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __ilshift__


@runtime_checkable
class CanILshift(Protocol[_T_contra, _T_co]):
    def __ilshift__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanILshiftSelf(Protocol[_T_contra]):
    """CanILshiftSelf[-T] = CanILshift[T, Self]"""

    def __ilshift__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanILshiftSame(Protocol[_T_Never_contra]):
    """CanILshiftSame[-T = Never] = CanILshift[Self | T, Self]"""

    def __ilshift__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __irshift__


@runtime_checkable
class CanIRshift(Protocol[_T_contra, _T_co]):
    def __irshift__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIRshiftSelf(Protocol[_T_contra]):
    """CanIRshiftSelf[-T] = CanIRshift[T, Self]"""

    def __irshift__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIRshiftSame(Protocol[_T_Never_contra]):
    """CanIRshiftSame[-T = Never] = CanIRshift[Self | T, Self]"""

    def __irshift__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __iand__


@runtime_checkable
class CanIAnd(Protocol[_T_contra, _T_co]):
    def __iand__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIAndSelf(Protocol[_T_contra]):
    """CanIAndSelf[-T] = CanIAnd[T, Self]"""

    def __iand__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIAndSame(Protocol[_T_Never_contra]):
    """CanIAndSame[-T = Never] = CanIAnd[Self | T, Self]"""

    def __iand__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __ixor__


@runtime_checkable
class CanIXor(Protocol[_T_contra, _T_co]):
    def __ixor__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIXorSelf(Protocol[_T_contra]):
    """CanIXorSelf[-T] = CanIXor[T, Self]"""

    def __ixor__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIXorSame(Protocol[_T_Never_contra]):
    """CanIXorSame[-T = Never] = CanIXor[Self | T, Self]"""

    def __ixor__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


# __ior__


@runtime_checkable
class CanIOr(Protocol[_T_contra, _T_co]):
    def __ior__(self, rhs: _T_contra, /) -> _T_co: ...


@runtime_checkable
class CanIOrSelf(Protocol[_T_contra]):
    """CanIOrSelf[-T] = CanIOr[T, Self]"""

    def __ior__(self, rhs: _T_contra, /) -> Self: ...


@runtime_checkable
class CanIOrSame(Protocol[_T_Never_contra]):
    """CanIOrSame[-T = Never] = CanIOr[Self | T, Self]"""

    def __ior__(self, rhs: Self | _T_Never_contra, /) -> Self: ...


###
# Unary arithmetic ops

# __neg__


@runtime_checkable
class CanNeg(Protocol[_T_co]):
    def __neg__(self, /) -> _T_co: ...


@runtime_checkable
class CanNegSelf(Protocol[_T_Never_co]):
    """CanNegSelf[+T = Never] = CanNeg[Self | T]"""

    def __neg__(self, /) -> Self | _T_Never_co: ...


# __pos__


@runtime_checkable
class CanPos(Protocol[_T_co]):
    def __pos__(self, /) -> _T_co: ...


@runtime_checkable
class CanPosSelf(Protocol[_T_Never_co]):
    """CanPosSelf[+T = Never] = CanPos[Self | T]"""

    def __pos__(self, /) -> Self | _T_Never_co: ...


# __abs__


@runtime_checkable
class CanAbs(Protocol[_T_co]):
    def __abs__(self, /) -> _T_co: ...


@runtime_checkable
class CanAbsSelf(Protocol[_T_Never_co]):
    """CanAbsSelf[+T = Never] = CanAbs[Self | T]"""

    def __abs__(self, /) -> Self | _T_Never_co: ...


# __invert__


@runtime_checkable
class CanInvert(Protocol[_T_co]):
    def __invert__(self, /) -> _T_co: ...


@runtime_checkable
class CanInvertSelf(Protocol[_T_Never_co]):
    """CanInvertSelf[+T = Never] = CanInvert[Self | T]"""

    def __invert__(self, /) -> Self | _T_Never_co: ...


###
# Rounding


@runtime_checkable
class CanRound1(Protocol[_T_int_co]):
    def __round__(self, /) -> _T_int_co: ...


@runtime_checkable
class CanRound2(Protocol[_T_int_contra, _T_float_co]):
    def __round__(self, /, ndigits: _T_int_contra) -> _T_float_co: ...


@runtime_checkable
class CanRound(
    CanRound1[_T_int_co],
    CanRound2[_T_int_contra, _T_float_co],
    Protocol[_T_int_contra, _T_int_co, _T_float_co],
):
    @overload
    @override
    def __round__(self, /) -> _T_int_co: ...
    @overload
    def __round__(self, /, ndigits: _T_int_contra) -> _T_float_co: ...


@runtime_checkable
class CanTrunc(Protocol[_T_int_co]):
    def __trunc__(self, /) -> _T_int_co: ...


@runtime_checkable
class CanFloor(Protocol[_T_int_co]):
    def __floor__(self, /) -> _T_int_co: ...


@runtime_checkable
class CanCeil(Protocol[_T_int_co]):
    def __ceil__(self, /) -> _T_int_co: ...


# Awaitables


@runtime_checkable
class CanAwait(Protocol[_T_co]):
    # The "return" value of a `yield from _` is attached to the `StopIteration`
    # exception when raised in `__next__()`. However, there is currently no way to
    # express this in Python's type system. So because `__await__()` could return any
    # iterator of `None | asyncio.Future`, it's theoretically impossible to annotate
    # awaitables using the Python type system.
    # In practice, the `collections.abc.Generator` type can be used for this. But note
    # that this type only works because type-checkers special-cased it. It's not
    # possible to write a custom type that behaves like `Generator` in this regard.

    @overload
    def __await__(self: "CanAwait[_T_co]", /) -> _AsyncGen[_T_co]: ...
    @overload
    def __await__(self: "CanAwait[None]", /) -> CanNext[_FutureOrNone]: ...


# Context managers


@runtime_checkable
class CanEnter(Protocol[_T_co]):
    def __enter__(self, /) -> _T_co: ...


@runtime_checkable
class CanEnterSelf(Protocol):
    """CanEnterSelf = CanEnter[Self]"""

    def __enter__(self, /) -> Self: ...


@runtime_checkable
class CanExit(Protocol[_T_None_co]):
    @overload
    def __exit__(self, exc_type: None, exc: None, tb: None, /) -> None: ...
    @overload
    def __exit__(  # noqa: PYI036
        self,
        exc_type: type[_ExcT],
        exc: _ExcT,
        tb: types.TracebackType,
        /,
    ) -> _T_None_co: ...


@runtime_checkable
class CanWith(CanEnter[_T_co], CanExit[_T_None_co], Protocol[_T_co, _T_None_co]): ...


@runtime_checkable
class CanWithSelf(CanEnterSelf, CanExit[_T_None_co], Protocol[_T_None_co]):
    """CanWithSelf[+R = None] = CanWith[Self, R]"""


# Async context managers


@runtime_checkable
class CanAEnter(Protocol[_T_co]):
    def __aenter__(self, /) -> CanAwait[_T_co]: ...


@runtime_checkable
class CanAEnterSelf(Protocol):
    """CanAEnterSelf = CanAEnter[Self]"""

    def __aenter__(self, /) -> CanAwait[Self]: ...


@runtime_checkable
class CanAExit(Protocol[_T_None_co]):
    @overload
    def __aexit__(self, exc_type: None, exc: None, tb: None, /) -> CanAwait[None]: ...
    @overload
    def __aexit__(
        self,
        exc_type: type[_ExcT],
        exc_: _ExcT,
        tb: types.TracebackType,
        /,
    ) -> CanAwait[_T_None_co]: ...


@runtime_checkable
class CanAsyncWith(
    CanAEnter[_T_co],
    CanAExit[_T_None_co],
    Protocol[_T_co, _T_None_co],
): ...


@runtime_checkable
class CanAsyncWithSelf(CanAEnterSelf, CanAExit[_T_None_co], Protocol[_T_None_co]):
    """CanAsyncWithSelf[+R = None] = CanAsyncWith[Self, R]"""


# Buffer protocol


@runtime_checkable
class CanBuffer(Protocol):
    def __buffer__(self, buffer: int, /) -> memoryview: ...


@runtime_checkable
class CanReleaseBuffer(Protocol):
    def __release_buffer__(self, buffer: memoryview, /) -> None: ...
