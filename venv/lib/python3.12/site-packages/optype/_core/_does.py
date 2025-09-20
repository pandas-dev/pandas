import sys
from collections.abc import AsyncIterator, Callable, Iterable, Iterator
from typing import Literal, Protocol, TypeAlias, overload

if sys.version_info >= (3, 13):
    from typing import ParamSpec, TypeVar
else:
    from typing_extensions import ParamSpec, TypeVar

import optype._core._can as _c

__all__ = [
    "DoesAIter",
    "DoesANext",
    "DoesAbs",
    "DoesAdd",
    "DoesAnd",
    "DoesBool",
    "DoesBytes",
    "DoesCall",
    "DoesCeil",
    "DoesComplex",
    "DoesContains",
    "DoesDelattr",
    "DoesDelitem",
    "DoesDir",
    "DoesDivmod",
    "DoesEq",
    "DoesFloat",
    "DoesFloor",
    "DoesFloordiv",
    "DoesFormat",
    "DoesGe",
    "DoesGetattr",
    "DoesGetitem",
    "DoesGt",
    "DoesHash",
    "DoesIAdd",
    "DoesIAnd",
    "DoesIFloordiv",
    "DoesILshift",
    "DoesIMatmul",
    "DoesIMod",
    "DoesIMul",
    "DoesIOr",
    "DoesIPow",
    "DoesIRshift",
    "DoesISub",
    "DoesITruediv",
    "DoesIXor",
    "DoesIndex",
    "DoesInt",
    "DoesInvert",
    "DoesIter",
    "DoesLe",
    "DoesLen",
    "DoesLengthHint",
    "DoesLshift",
    "DoesLt",
    "DoesMatmul",
    "DoesMissing",
    "DoesMod",
    "DoesMul",
    "DoesNe",
    "DoesNeg",
    "DoesNext",
    "DoesOr",
    "DoesPos",
    "DoesPow",
    "DoesRAdd",
    "DoesRAnd",
    "DoesRDivmod",
    "DoesRFloordiv",
    "DoesRLshift",
    "DoesRMatmul",
    "DoesRMod",
    "DoesRMul",
    "DoesROr",
    "DoesRPow",
    "DoesRRshift",
    "DoesRSub",
    "DoesRTruediv",
    "DoesRXor",
    "DoesRepr",
    "DoesReversed",
    "DoesRound",
    "DoesRshift",
    "DoesSetattr",
    "DoesSetitem",
    "DoesStr",
    "DoesSub",
    "DoesTruediv",
    "DoesTrunc",
    "DoesXor",
]


def __dir__() -> list[str]:
    return __all__


###

_JustFalse: TypeAlias = Literal[False]
_JustTrue: TypeAlias = Literal[True]
_Just0: TypeAlias = Literal[0]
# cannot use `optype.typing.LiteralByte` here, as it starts at 0
_PosInt: TypeAlias = Literal[
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
]  # fmt: skip


_Tss = ParamSpec("_Tss")
_KeyT = TypeVar("_KeyT")
_ValT = TypeVar("_ValT")
_AttrT = TypeVar("_AttrT")
_LeftT = TypeVar("_LeftT")
_RightT = TypeVar("_RightT")
_ModT = TypeVar("_ModT")
_NDigitsT = TypeVar("_NDigitsT")
_OutT = TypeVar("_OutT")
_DefaultT = TypeVar("_DefaultT")
_SentinelT = TypeVar("_SentinelT")
_IteratorT = TypeVar("_IteratorT", bound=Iterator[object] | _c.CanNext[object])
_AIteratorT = TypeVar("_AIteratorT", bound=AsyncIterator[object] | _c.CanANext[object])
_IterT = TypeVar("_IterT", bound=Iterable[object])
_BoolT = TypeVar("_BoolT", bound=bool)
_IntT = TypeVar("_IntT", bound=int)
_StrT = TypeVar("_StrT", bound=str)
_FormatT = TypeVar("_FormatT", bound=str)
_BytesT = TypeVar("_BytesT", bound=bytes)


###


# iteration


class DoesNext(Protocol):
    @overload
    def __call__(self, iterator: _c.CanNext[_ValT], /) -> _ValT: ...
    @overload
    def __call__(
        self,
        iterator: _c.CanNext[_ValT],
        default: _DefaultT,
        /,
    ) -> _ValT | _DefaultT: ...


class DoesANext(Protocol):
    @overload
    def __call__(self, aiterator: _c.CanANext[_ValT], /) -> _ValT: ...
    @overload
    async def __call__(
        self,
        aiterator: _c.CanANext[_c.CanAwait[_ValT]],
        default: _DefaultT,
        /,
    ) -> _ValT | _DefaultT: ...


class DoesIter(Protocol):
    @overload
    def __call__(self, iterable: _c.CanIter[_IteratorT], /) -> _IteratorT: ...
    @overload
    def __call__(
        self,
        sequence: _c.CanGetitem[_c.CanIndex, _ValT],
        /,
    ) -> _c.CanIterSelf[_ValT]: ...
    @overload
    def __call__(
        self,
        callable_: _c.CanCall[[], _ValT | None],
        sentinel: None,
        /,
    ) -> _c.CanIterSelf[_ValT]: ...
    @overload
    def __call__(
        self,
        callable_: _c.CanCall[[], _ValT | _SentinelT],
        sentinel: _SentinelT,
        /,
    ) -> _c.CanIterSelf[_ValT]: ...


class DoesAIter(Protocol):
    def __call__(self, aiterable: _c.CanAIter[_AIteratorT], /) -> _AIteratorT: ...


# type conversion


class DoesComplex(Protocol):
    def __call__(self, obj: _c.CanComplex, /) -> complex: ...


class DoesFloat(Protocol):
    def __call__(self, obj: _c.CanFloat, /) -> float: ...


class DoesInt(Protocol):
    def __call__(self, obj: _c.CanInt[_IntT], /) -> _IntT: ...


class DoesBool(Protocol):
    @overload
    def __call__(self, obj: _c.CanBool[_BoolT], /) -> _BoolT: ...
    @overload
    def __call__(self, obj: _c.CanLen[_Just0], /) -> _JustFalse: ...
    @overload
    def __call__(self, obj: _c.CanLen[_PosInt], /) -> _JustTrue: ...
    @overload
    def __call__(self, obj: object, /) -> bool: ...


class DoesStr(Protocol):
    def __call__(self, obj: _c.CanStr[_StrT], /) -> _StrT: ...


class DoesBytes(Protocol):
    def __call__(self, obj: _c.CanBytes[_BytesT], /) -> _BytesT: ...


# formatting


class DoesRepr(Protocol):
    def __call__(self, obj: _c.CanRepr[_StrT], /) -> _StrT: ...


class DoesFormat(Protocol):
    def __call__(
        self,
        obj: _c.CanFormat[_FormatT, _StrT],
        format_spec: _FormatT = ...,
        /,
    ) -> _StrT: ...


# rich comparison


class DoesLt(Protocol):
    @overload
    def __call__(self, lhs: _c.CanLt[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanGt[_LeftT, _OutT], /) -> _OutT: ...


class DoesLe(Protocol):
    @overload
    def __call__(self, lhs: _c.CanLe[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanGe[_LeftT, _OutT], /) -> _OutT: ...


class DoesEq(Protocol):
    @overload
    def __call__(self, lhs: _c.CanEq[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanEq[_LeftT, _OutT], /) -> _OutT: ...  # type: ignore[overload-cannot-match]  # pyright: ignore[reportOverlappingOverload]


class DoesNe(Protocol):
    @overload
    def __call__(self, lhs: _c.CanNe[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanNe[_LeftT, _OutT], /) -> _OutT: ...  # type: ignore[overload-cannot-match]  # pyright: ignore[reportOverlappingOverload]


class DoesGt(Protocol):
    @overload
    def __call__(self, lhs: _c.CanGt[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanLt[_LeftT, _OutT], /) -> _OutT: ...


class DoesGe(Protocol):
    @overload
    def __call__(self, lhs: _c.CanGe[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanLe[_LeftT, _OutT], /) -> _OutT: ...


# dynamic attribute access


class DoesGetattr(Protocol):
    @overload
    def __call__(self, obj: _c.CanGetattr[_StrT, _AttrT], name: _StrT, /) -> _AttrT: ...
    @overload
    def __call__(
        self,
        obj: _c.CanGetattribute[_StrT, _AttrT],
        name: _StrT,
        /,
    ) -> _AttrT: ...
    @overload
    def __call__(
        self,
        obj: _c.CanGetattr[_StrT, _AttrT],
        name: _StrT,
        default: _DefaultT,
        /,
    ) -> _AttrT | _DefaultT: ...
    @overload
    def __call__(
        self,
        obj: _c.CanGetattribute[_StrT, _AttrT],
        name: _StrT,
        default: _DefaultT,
        /,
    ) -> _AttrT | _DefaultT: ...


class DoesSetattr(Protocol):
    def __call__(
        self,
        obj: _c.CanSetattr[_StrT, _AttrT],
        name: _StrT,
        value: _AttrT,
        /,
    ) -> None: ...


class DoesDelattr(Protocol):
    def __call__(self, obj: _c.CanDelattr[_StrT], name: _StrT, /) -> None: ...


class DoesDir(Protocol):
    @overload
    def __call__(self, /) -> list[str]: ...
    @overload
    def __call__(self, obj: _c.CanDir[_IterT], /) -> _IterT: ...


# callables


class DoesCall(Protocol):
    def __call__(
        self,
        callable_: Callable[_Tss, _OutT],
        /,
        *args: _Tss.args,
        **kwargs: _Tss.kwargs,
    ) -> _OutT: ...


# containers and subscriptable types


class DoesLen(Protocol):
    def __call__(self, obj: _c.CanLen[_IntT], /) -> _IntT: ...


class DoesLengthHint(Protocol):
    def __call__(self, obj: _c.CanLengthHint[_IntT], /) -> _IntT: ...


class DoesGetitem(Protocol):
    def __call__(
        self,
        obj: _c.CanGetitem[_KeyT, _ValT] | _c.CanGetMissing[_KeyT, _ValT, _DefaultT],
        key: _KeyT,
        /,
    ) -> _ValT | _DefaultT: ...


class DoesSetitem(Protocol):
    def __call__(
        self,
        obj: _c.CanSetitem[_KeyT, _ValT],
        key: _KeyT,
        value: _ValT,
        /,
    ) -> None: ...


class DoesDelitem(Protocol):
    def __call__(self, obj: _c.CanDelitem[_KeyT], key: _KeyT, /) -> None: ...


class DoesMissing(Protocol):
    def __call__(
        self,
        obj: _c.CanMissing[_KeyT, _DefaultT],
        key: _KeyT,
        /,
    ) -> _DefaultT: ...


class DoesContains(Protocol):
    def __call__(self, obj: _c.CanContains[_KeyT, _BoolT], key: _KeyT, /) -> _BoolT: ...


class DoesReversed(Protocol):
    """
    This is correct type of `builtins.reversed`.

    Note that typeshed's annotations for `reversed` are completely wrong:
    https://github.com/python/typeshed/issues/11645
    """

    @overload
    def __call__(self, reversible: _c.CanReversed[_OutT], /) -> _OutT: ...
    @overload
    def __call__(
        self,
        sequence: _c.CanSequence[_c.CanIndex, _ValT],
        /,
    ) -> "reversed[_ValT]": ...


# binary infix operators


class DoesAdd(Protocol):
    @overload
    def __call__(self, lhs: _c.CanAdd[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRAdd[_LeftT, _OutT], /) -> _OutT: ...


class DoesSub(Protocol):
    @overload
    def __call__(self, lhs: _c.CanSub[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRSub[_LeftT, _OutT], /) -> _OutT: ...


class DoesMul(Protocol):
    @overload
    def __call__(self, lhs: _c.CanMul[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRMul[_LeftT, _OutT], /) -> _OutT: ...


class DoesMatmul(Protocol):
    @overload
    def __call__(self, lhs: _c.CanMatmul[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRMatmul[_LeftT, _OutT], /) -> _OutT: ...


class DoesTruediv(Protocol):
    @overload
    def __call__(
        self,
        lhs: _c.CanTruediv[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRTruediv[_LeftT, _OutT], /) -> _OutT: ...


class DoesFloordiv(Protocol):
    @overload
    def __call__(
        self,
        lhs: _c.CanFloordiv[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(
        self,
        lhs: _LeftT,
        rhs: _c.CanRFloordiv[_LeftT, _OutT],
        /,
    ) -> _OutT: ...


class DoesMod(Protocol):
    @overload
    def __call__(self, lhs: _c.CanMod[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRMod[_LeftT, _OutT], /) -> _OutT: ...


class DoesDivmod(Protocol):
    @overload
    def __call__(self, lhs: _c.CanDivmod[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRDivmod[_LeftT, _OutT], /) -> _OutT: ...


class DoesPow(Protocol):
    @overload
    def __call__(self, base: _c.CanPow2[_RightT, _OutT], exp: _RightT, /) -> _OutT: ...
    @overload
    def __call__(
        self,
        base: _c.CanPow3[_RightT, _ModT, _OutT],
        exp: _RightT,
        mod: _ModT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(self, base: _LeftT, exp: _c.CanRPow[_LeftT, _OutT], /) -> _OutT: ...


class DoesLshift(Protocol):
    @overload
    def __call__(self, lhs: _c.CanLshift[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRLshift[_LeftT, _OutT], /) -> _OutT: ...


class DoesRshift(Protocol):
    @overload
    def __call__(self, lhs: _c.CanRshift[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRRshift[_LeftT, _OutT], /) -> _OutT: ...


class DoesAnd(Protocol):
    @overload
    def __call__(self, lhs: _c.CanAnd[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRAnd[_LeftT, _OutT], /) -> _OutT: ...


class DoesXor(Protocol):
    @overload
    def __call__(self, lhs: _c.CanXor[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRXor[_LeftT, _OutT], /) -> _OutT: ...


class DoesOr(Protocol):
    @overload
    def __call__(self, lhs: _c.CanOr[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanROr[_LeftT, _OutT], /) -> _OutT: ...


# binary reflected operators


class DoesRAdd(Protocol):
    def __call__(self, rhs: _c.CanRAdd[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRSub(Protocol):
    def __call__(self, rhs: _c.CanRSub[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRMul(Protocol):
    def __call__(self, rhs: _c.CanRMul[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRMatmul(Protocol):
    def __call__(self, rhs: _c.CanRMatmul[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRTruediv(Protocol):
    def __call__(self, rhs: _c.CanRTruediv[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRFloordiv(Protocol):
    def __call__(
        self,
        rhs: _c.CanRFloordiv[_LeftT, _OutT],
        lhs: _LeftT,
        /,
    ) -> _OutT: ...


class DoesRMod(Protocol):
    def __call__(self, rhs: _c.CanRMod[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRDivmod(Protocol):
    def __call__(self, rhs: _c.CanRDivmod[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRPow(Protocol):
    def __call__(self, rhs: _c.CanRPow[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRLshift(Protocol):
    def __call__(self, rhs: _c.CanRLshift[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRRshift(Protocol):
    def __call__(self, rhs: _c.CanRRshift[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRAnd(Protocol):
    def __call__(self, rhs: _c.CanRAnd[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesRXor(Protocol):
    def __call__(self, rhs: _c.CanRXor[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


class DoesROr(Protocol):
    def __call__(self, rhs: _c.CanROr[_LeftT, _OutT], lhs: _LeftT, /) -> _OutT: ...


# augmented / in-place operators


class DoesIAdd(Protocol):
    @overload
    def __call__(self, lhs: _c.CanIAdd[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanAdd[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRAdd[_LeftT, _OutT], /) -> _OutT: ...


class DoesISub(Protocol):
    @overload
    def __call__(self, lhs: _c.CanISub[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanSub[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRSub[_LeftT, _OutT], /) -> _OutT: ...


class DoesIMul(Protocol):
    @overload
    def __call__(self, lhs: _c.CanIMul[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanMul[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRMul[_LeftT, _OutT], /) -> _OutT: ...


class DoesIMatmul(Protocol):
    @overload
    def __call__(
        self,
        lhs: _c.CanIMatmul[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanMatmul[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRMatmul[_LeftT, _OutT], /) -> _OutT: ...


class DoesITruediv(Protocol):
    @overload
    def __call__(
        self,
        lhs: _c.CanITruediv[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(
        self,
        lhs: _c.CanTruediv[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRTruediv[_LeftT, _OutT], /) -> _OutT: ...


class DoesIFloordiv(Protocol):
    @overload
    def __call__(
        self,
        lhs: _c.CanIFloordiv[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(
        self,
        lhs: _c.CanFloordiv[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(
        self,
        lhs: _LeftT,
        rhs: _c.CanRFloordiv[_LeftT, _OutT],
        /,
    ) -> _OutT: ...


class DoesIMod(Protocol):
    @overload
    def __call__(self, lhs: _c.CanIMod[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanMod[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRMod[_LeftT, _OutT], /) -> _OutT: ...


class DoesIPow(Protocol):
    @overload
    def __call__(self, lhs: _c.CanIPow[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanPow2[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRPow[_LeftT, _OutT], /) -> _OutT: ...


class DoesILshift(Protocol):
    @overload
    def __call__(
        self,
        lhs: _c.CanILshift[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanLshift[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRLshift[_LeftT, _OutT], /) -> _OutT: ...


class DoesIRshift(Protocol):
    @overload
    def __call__(
        self,
        lhs: _c.CanIRshift[_RightT, _OutT],
        rhs: _RightT,
        /,
    ) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanRshift[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRRshift[_LeftT, _OutT], /) -> _OutT: ...


class DoesIAnd(Protocol):
    @overload
    def __call__(self, lhs: _c.CanIAnd[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanAnd[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRAnd[_LeftT, _OutT], /) -> _OutT: ...


class DoesIXor(Protocol):
    @overload
    def __call__(self, lhs: _c.CanIXor[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanXor[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanRXor[_LeftT, _OutT], /) -> _OutT: ...


class DoesIOr(Protocol):
    @overload
    def __call__(self, lhs: _c.CanIOr[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _c.CanOr[_RightT, _OutT], rhs: _RightT, /) -> _OutT: ...
    @overload
    def __call__(self, lhs: _LeftT, rhs: _c.CanROr[_LeftT, _OutT], /) -> _OutT: ...


# unary arithmetic


class DoesNeg(Protocol):
    def __call__(self, obj: _c.CanNeg[_OutT], /) -> _OutT: ...


class DoesPos(Protocol):
    def __call__(self, obj: _c.CanPos[_OutT], /) -> _OutT: ...


class DoesAbs(Protocol):
    def __call__(self, obj: _c.CanAbs[_OutT], /) -> _OutT: ...


class DoesInvert(Protocol):
    def __call__(self, obj: _c.CanInvert[_OutT], /) -> _OutT: ...


# object identification


class DoesIndex(Protocol):
    def __call__(self, obj: _c.CanIndex[_IntT], /) -> _IntT: ...


class DoesHash(Protocol):
    def __call__(self, obj: _c.CanHash, /) -> int: ...


# rounding


class DoesRound(Protocol):
    @overload
    def __call__(self, obj: _c.CanRound1[_OutT], /) -> _OutT: ...
    @overload
    def __call__(self, obj: _c.CanRound1[_OutT], /, ndigits: None = None) -> _OutT: ...
    @overload
    def __call__(
        self,
        obj: _c.CanRound2[_NDigitsT, _OutT],
        /,
        ndigits: _NDigitsT,
    ) -> _OutT: ...


class DoesTrunc(Protocol):
    def __call__(self, obj: _c.CanTrunc[_OutT], /) -> _OutT: ...


class DoesFloor(Protocol):
    def __call__(self, obj: _c.CanFloor[_OutT], /) -> _OutT: ...


class DoesCeil(Protocol):
    def __call__(self, obj: _c.CanCeil[_OutT], /) -> _OutT: ...
