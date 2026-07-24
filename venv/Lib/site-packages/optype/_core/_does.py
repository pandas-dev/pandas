from collections.abc import AsyncIterator, Callable, Iterable, Iterator
from typing import Protocol, SupportsIndex, overload

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


# iteration


class DoesNext(Protocol):
    @overload
    def __call__[ValT](self, iterator: _c.CanNext[ValT], /) -> ValT: ...
    @overload
    def __call__[ValT, DefaultT](
        self,
        iterator: _c.CanNext[ValT],
        default: DefaultT,
        /,
    ) -> ValT | DefaultT: ...


class DoesANext(Protocol):
    @overload
    def __call__[ValT](self, aiterator: _c.CanANext[ValT], /) -> ValT: ...
    @overload
    async def __call__[ValT, DefaultT](
        self,
        aiterator: _c.CanANext[_c.CanAwait[ValT]],
        default: DefaultT,
        /,
    ) -> ValT | DefaultT: ...


class DoesIter(Protocol):
    @overload
    def __call__[IteratorT: Iterator[object] | _c.CanNext[object]](
        self,
        iterable: _c.CanIter[IteratorT],
        /,
    ) -> IteratorT: ...
    @overload
    def __call__[ValT](
        self,
        sequence: _c.CanGetitem[_c.CanIndex, ValT],
        /,
    ) -> _c.CanIterSelf[ValT]: ...
    @overload
    def __call__[ValT](
        self,
        callable_: _c.CanCall[[], ValT | None],
        sentinel: None,
        /,
    ) -> _c.CanIterSelf[ValT]: ...
    @overload
    def __call__[ValT, SentinelT](
        self,
        callable_: _c.CanCall[[], ValT | SentinelT],
        sentinel: SentinelT,
        /,
    ) -> _c.CanIterSelf[ValT]: ...


class DoesAIter(Protocol):
    def __call__[AIteratorT: AsyncIterator[object] | _c.CanANext[object]](
        self,
        aiterable: _c.CanAIter[AIteratorT],
        /,
    ) -> AIteratorT: ...


# type conversion


class DoesComplex(Protocol):
    def __call__(self, obj: _c.CanComplex, /) -> complex: ...


class DoesFloat(Protocol):
    def __call__(self, obj: _c.CanFloat, /) -> float: ...


class DoesInt(Protocol):
    def __call__(self, obj: _c.CanInt, /) -> int: ...


class DoesBool(Protocol):
    @overload
    def __call__[BoolT: bool](self, obj: _c.CanBool[BoolT], /) -> BoolT: ...
    @overload
    def __call__(self, obj: _c.CanLen, /) -> bool: ...
    @overload
    def __call__(self, obj: object, /) -> bool: ...


class DoesStr(Protocol):
    def __call__(self, obj: _c.CanStr, /) -> str: ...


class DoesBytes(Protocol):
    def __call__(self, obj: _c.CanBytes, /) -> bytes: ...


# formatting


class DoesRepr(Protocol):
    def __call__(self, obj: _c.CanRepr, /) -> str: ...


class DoesFormat(Protocol):
    def __call__(self, obj: _c.CanFormat, format_spec: str = ..., /) -> str: ...


# rich comparison


class DoesLt(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanLt[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanGt[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesLe(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanLe[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanGe[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesEq(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanEq[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](  # type: ignore[overload-cannot-match] # pyright: ignore[reportOverlappingOverload]
        self,
        lhs: LeftT,
        rhs: _c.CanEq[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesNe(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanNe[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](  # type: ignore[overload-cannot-match] # pyright: ignore[reportOverlappingOverload]
        self,
        lhs: LeftT,
        rhs: _c.CanNe[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesGt(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanGt[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanLt[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesGe(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanGe[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanLe[LeftT, OutT],
        /,
    ) -> OutT: ...


# dynamic attribute access


class DoesGetattr(Protocol):
    @overload
    def __call__[AttrT](self, obj: _c.CanGetattr[AttrT], name: str, /) -> AttrT: ...
    @overload
    def __call__[AttrT](
        self,
        obj: _c.CanGetattribute[AttrT],
        name: str,
        /,
    ) -> AttrT: ...
    @overload
    def __call__[AttrT, DefaultT](
        self,
        obj: _c.CanGetattr[AttrT],
        name: str,
        default: DefaultT,
        /,
    ) -> AttrT | DefaultT: ...
    @overload
    def __call__[AttrT, DefaultT](
        self,
        obj: _c.CanGetattribute[AttrT],
        name: str,
        default: DefaultT,
        /,
    ) -> AttrT | DefaultT: ...


class DoesSetattr(Protocol):
    def __call__[AttrT](
        self,
        obj: _c.CanSetattr[AttrT],
        name: str,
        value: AttrT,
        /,
    ) -> None: ...


class DoesDelattr(Protocol):
    def __call__(self, obj: _c.CanDelattr, name: str, /) -> None: ...


class DoesDir(Protocol):
    @overload
    def __call__(self, /) -> list[str]: ...
    @overload
    def __call__[IterT: Iterable[object]](self, obj: _c.CanDir[IterT], /) -> IterT: ...


# callables


class DoesCall(Protocol):
    def __call__[**Tss, OutT](
        self,
        callable_: Callable[Tss, OutT],
        /,
        *args: Tss.args,
        **kwargs: Tss.kwargs,
    ) -> OutT: ...


# containers and subscriptable types


class DoesLen(Protocol):
    def __call__(self, obj: _c.CanLen, /) -> int: ...


class DoesLengthHint(Protocol):
    def __call__(self, obj: _c.CanLengthHint, /) -> int: ...


class DoesGetitem(Protocol):
    def __call__[KeyT, ValT, DefaultT](
        self,
        obj: _c.CanGetitem[KeyT, ValT] | _c.CanGetMissing[KeyT, ValT, DefaultT],
        key: KeyT,
        /,
    ) -> ValT | DefaultT: ...


class DoesSetitem(Protocol):
    def __call__[KeyT, ValT](
        self,
        obj: _c.CanSetitem[KeyT, ValT],
        key: KeyT,
        value: ValT,
        /,
    ) -> None: ...


class DoesDelitem(Protocol):
    def __call__[KeyT](self, obj: _c.CanDelitem[KeyT], key: KeyT, /) -> None: ...


class DoesMissing(Protocol):
    def __call__[KeyT, DefaultT](
        self,
        obj: _c.CanMissing[KeyT, DefaultT],
        key: KeyT,
        /,
    ) -> DefaultT: ...


class DoesContains(Protocol):
    def __call__[KeyT](self, obj: _c.CanContains[KeyT], key: KeyT, /) -> bool: ...


class DoesReversed(Protocol):
    """
    This is correct type of `builtins.reversed`.

    Note that typeshed's annotations for `reversed` are completely wrong:
    https://github.com/python/typeshed/issues/11645
    """

    @overload
    def __call__[OutT](self, reversible: _c.CanReversed[OutT], /) -> OutT: ...
    @overload
    def __call__[ValT](
        self,
        sequence: _c.CanSequence[_c.CanIndex, ValT],
        /,
    ) -> "reversed[ValT]": ...


# binary infix operators


class DoesAdd(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanAdd[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRAdd[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesSub(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanSub[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRSub[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesMul(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanMul[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRMul[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesMatmul(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanMatmul[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRMatmul[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesTruediv(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanTruediv[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRTruediv[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesFloordiv(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanFloordiv[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRFloordiv[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesMod(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanMod[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRMod[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesDivmod(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanDivmod[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRDivmod[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesPow(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        base: _c.CanPow2[RightT, OutT],
        exp: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, ModT, OutT](
        self,
        base: _c.CanPow3[RightT, ModT, OutT],
        exp: RightT,
        mod: ModT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        base: LeftT,
        exp: _c.CanRPow[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesLshift(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanLshift[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRLshift[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesRshift(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanRshift[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRRshift[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesAnd(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanAnd[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRAnd[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesXor(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanXor[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRXor[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesOr(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanOr[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanROr[LeftT, OutT],
        /,
    ) -> OutT: ...


# binary reflected operators


class DoesRAdd(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRAdd[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRSub(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRSub[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRMul(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRMul[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRMatmul(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRMatmul[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRTruediv(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRTruediv[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRFloordiv(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRFloordiv[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRMod(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRMod[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRDivmod(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRDivmod[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRPow(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRPow[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRLshift(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRLshift[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRRshift(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRRshift[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRAnd(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRAnd[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesRXor(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanRXor[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


class DoesROr(Protocol):
    def __call__[LeftT, OutT](
        self,
        rhs: _c.CanROr[LeftT, OutT],
        lhs: LeftT,
        /,
    ) -> OutT: ...


# augmented / in-place operators


class DoesIAdd(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIAdd[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanAdd[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRAdd[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesISub(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanISub[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanSub[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRSub[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIMul(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIMul[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanMul[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRMul[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIMatmul(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIMatmul[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanMatmul[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRMatmul[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesITruediv(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanITruediv[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanTruediv[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRTruediv[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIFloordiv(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIFloordiv[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanFloordiv[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRFloordiv[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIMod(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIMod[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanMod[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRMod[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIPow(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIPow[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanPow2[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRPow[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesILshift(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanILshift[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanLshift[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRLshift[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIRshift(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIRshift[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanRshift[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRRshift[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIAnd(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIAnd[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanAnd[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRAnd[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIXor(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIXor[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanXor[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanRXor[LeftT, OutT],
        /,
    ) -> OutT: ...


class DoesIOr(Protocol):
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanIOr[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[RightT, OutT](
        self,
        lhs: _c.CanOr[RightT, OutT],
        rhs: RightT,
        /,
    ) -> OutT: ...
    @overload
    def __call__[LeftT, OutT](
        self,
        lhs: LeftT,
        rhs: _c.CanROr[LeftT, OutT],
        /,
    ) -> OutT: ...


# unary arithmetic


class DoesNeg(Protocol):
    def __call__[OutT](self, obj: _c.CanNeg[OutT], /) -> OutT: ...


class DoesPos(Protocol):
    def __call__[OutT](self, obj: _c.CanPos[OutT], /) -> OutT: ...


class DoesAbs(Protocol):
    def __call__[OutT](self, obj: _c.CanAbs[OutT], /) -> OutT: ...


class DoesInvert(Protocol):
    def __call__[OutT](self, obj: _c.CanInvert[OutT], /) -> OutT: ...


# object identification


class DoesIndex(Protocol):
    def __call__(self, obj: _c.CanIndex, /) -> int: ...


class DoesHash(Protocol):
    def __call__(self, obj: _c.CanHash, /) -> int: ...


# rounding


class DoesRound(Protocol):
    @overload
    def __call__[OutT](
        self,
        number: _c.CanRound1[OutT],
        ndigits: None = None,
        /,
    ) -> OutT: ...
    @overload  # this unnecessary overload works around the issue described below
    def __call__[OutT](
        self,
        number: _c.CanRound2[SupportsIndex, OutT],
        ndigits: SupportsIndex,
        /,
    ) -> OutT: ...
    @overload  # all type-checkers except ty fail to correctly resolve this overload
    def __call__[OutT, NDigitsT](
        self,
        number: _c.CanRound2[NDigitsT, OutT],
        ndigits: NDigitsT,
        /,
    ) -> OutT: ...


class DoesTrunc(Protocol):
    def __call__[OutT](self, obj: _c.CanTrunc[OutT], /) -> OutT: ...


class DoesFloor(Protocol):
    def __call__[OutT](self, obj: _c.CanFloor[OutT], /) -> OutT: ...


class DoesCeil(Protocol):
    def __call__[OutT](self, obj: _c.CanCeil[OutT], /) -> OutT: ...
