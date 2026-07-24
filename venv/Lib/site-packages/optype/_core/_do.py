# pyright: reportInvalidCast=false

import math
import operator as _o
from typing import Final, cast, overload

import optype._core._can as _c
import optype._core._does as _d

__all__ = [
    "do_abs",
    "do_add",
    "do_aiter",
    "do_and",
    "do_anext",
    "do_bool",
    "do_bytes",
    "do_call",
    "do_ceil",
    "do_complex",
    "do_contains",
    "do_delattr",
    "do_delitem",
    "do_dir",
    "do_divmod",
    "do_eq",
    "do_float",
    "do_floor",
    "do_floordiv",
    "do_format",
    "do_ge",
    "do_getattr",
    "do_getitem",
    "do_gt",
    "do_hash",
    "do_iadd",
    "do_iand",
    "do_ifloordiv",
    "do_ilshift",
    "do_imatmul",
    "do_imod",
    "do_imul",
    "do_index",
    "do_int",
    "do_invert",
    "do_ior",
    "do_ipow",
    "do_irshift",
    "do_isub",
    "do_iter",
    "do_itruediv",
    "do_ixor",
    "do_le",
    "do_len",
    "do_length_hint",
    "do_lshift",
    "do_lt",
    "do_matmul",
    "do_missing",
    "do_mod",
    "do_mul",
    "do_ne",
    "do_neg",
    "do_next",
    "do_or",
    "do_pos",
    "do_pow",
    "do_radd",
    "do_rand",
    "do_rdivmod",
    "do_repr",
    "do_reversed",
    "do_rfloordiv",
    "do_rlshift",
    "do_rmatmul",
    "do_rmod",
    "do_rmul",
    "do_ror",
    "do_round",
    "do_rpow",
    "do_rrshift",
    "do_rshift",
    "do_rsub",
    "do_rtruediv",
    "do_rxor",
    "do_setattr",
    "do_setitem",
    "do_str",
    "do_sub",
    "do_truediv",
    "do_trunc",
    "do_xor",
]


def __dir__() -> list[str]:
    return __all__


###


# type conversion
do_bool: Final[_d.DoesBool] = cast("_d.DoesBool", bool)
do_int: Final[_d.DoesInt] = cast("_d.DoesInt", int)
do_float: Final[_d.DoesFloat] = cast("_d.DoesFloat", float)
do_complex: Final[_d.DoesComplex] = cast("_d.DoesComplex", complex)
do_bytes: Final[_d.DoesBytes] = cast("_d.DoesBytes", bytes)
do_str: Final[_d.DoesStr] = cast("_d.DoesStr", str)

# formatting
do_repr: Final[_d.DoesRepr] = repr
do_format: Final[_d.DoesFormat] = format

# iteration
do_next: Final[_d.DoesNext] = next
do_iter: Final[_d.DoesIter] = cast("_d.DoesIter", iter)

# async iteration
do_anext: Final[_d.DoesANext] = cast("_d.DoesANext", anext)
do_aiter: Final[_d.DoesAIter] = cast("_d.DoesAIter", aiter)

# rich comparison
do_lt: Final[_d.DoesLt] = cast("_d.DoesLt", _o.lt)
do_le: Final[_d.DoesLe] = cast("_d.DoesLe", _o.le)
do_eq: Final[_d.DoesEq] = cast("_d.DoesEq", _o.eq)
do_ne: Final[_d.DoesNe] = cast("_d.DoesNe", _o.ne)
do_gt: Final[_d.DoesGt] = cast("_d.DoesGt", _o.gt)
do_ge: Final[_d.DoesGe] = cast("_d.DoesGe", _o.ge)

# attributes
do_getattr: Final[_d.DoesGetattr] = cast("_d.DoesGetattr", getattr)
do_setattr: Final[_d.DoesSetattr] = cast("_d.DoesSetattr", setattr)
do_delattr: Final[_d.DoesDelattr] = delattr
do_dir: Final[_d.DoesDir] = cast("_d.DoesDir", dir)

# callables

do_call: Final[_d.DoesCall] = _o.call

# containers and sequences

do_len: Final[_d.DoesLen] = len
do_length_hint: Final[_d.DoesLengthHint] = cast("_d.DoesLengthHint", _o.length_hint)


# `operator.getitem` isn't used, because it has an (unreasonably loose, and
# redundant) overload for `(Sequence[T], slice) -> Sequence[T]`
# https://github.com/python/typeshed/blob/587ad6b/stdlib/_operator.pyi#L84-L86


@overload
def do_getitem[KT, VT, DT](
    obj: _c.CanGetMissing[KT, VT, DT],
    key: KT,
    /,
) -> VT | DT: ...
@overload
def do_getitem[KT, VT](obj: _c.CanGetitem[KT, VT], key: KT, /) -> VT: ...
def do_getitem[KT, VT, DT](
    obj: _c.CanGetitem[KT, VT] | _c.CanGetMissing[KT, VT, DT],
    key: KT,
    /,
) -> VT | DT:
    """Same as `value = obj[key]`."""
    return obj[key]


def do_setitem[KT, VT](obj: _c.CanSetitem[KT, VT], key: KT, value: VT, /) -> None:
    """Same as `obj[key] = value`."""
    obj[key] = value


def do_delitem[KT](obj: _c.CanDelitem[KT], key: KT, /) -> None:
    """Same as `del obj[key]`."""
    # https://github.com/astral-sh/ty/issues/1799
    del obj[key]


def do_missing[KT, DT](obj: _c.CanMissing[KT, DT], key: KT, /) -> DT:
    return obj.__missing__(key)


# `operator.contains` cannot be used, as it incorrectly requires `key`
# to be exactly of type `object`, so that it only accepts `object()`...
def do_contains[KT](obj: _c.CanContains[KT], key: KT, /) -> bool:
    """Same as `key in obj`."""
    return key in obj


# `builtins.reversed` is annotated incorrectly within typeshed:
# https://github.com/python/typeshed/issues/11645
do_reversed: Final[_d.DoesReversed] = cast("_d.DoesReversed", reversed)


# infix ops
do_add: Final[_d.DoesAdd] = cast("_d.DoesAdd", _o.add)
do_sub: Final[_d.DoesSub] = cast("_d.DoesSub", _o.sub)
do_mul: Final[_d.DoesMul] = cast("_d.DoesMul", _o.mul)
do_matmul: Final[_d.DoesMatmul] = cast("_d.DoesMatmul", _o.matmul)
do_truediv: Final[_d.DoesTruediv] = cast("_d.DoesTruediv", _o.truediv)
do_floordiv: Final[_d.DoesFloordiv] = cast("_d.DoesFloordiv", _o.floordiv)
do_mod: Final[_d.DoesMod] = cast("_d.DoesMod", _o.mod)
do_divmod: Final[_d.DoesDivmod] = divmod
do_pow: Final[_d.DoesPow] = cast("_d.DoesPow", pow)
do_lshift: Final[_d.DoesLshift] = cast("_d.DoesLshift", _o.lshift)
do_rshift: Final[_d.DoesRshift] = cast("_d.DoesRshift", _o.rshift)
do_and: Final[_d.DoesAnd] = cast("_d.DoesAnd", _o.and_)
do_xor: Final[_d.DoesXor] = cast("_d.DoesXor", _o.xor)
do_or: Final[_d.DoesOr] = cast("_d.DoesOr", _o.or_)


# reflected ops
# (a DRY `do_r*` decorator won't work; the overloads get lost during casting to
# `CanCall` or `Callable`, within the decorator function signature).


def do_radd[LeftT, OutT](a: _c.CanRAdd[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b + a`."""
    return b + a


def do_rsub[LeftT, OutT](a: _c.CanRSub[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b - a`."""
    return b - a


def do_rmul[LeftT, OutT](a: _c.CanRMul[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b * a`."""
    return b * a


def do_rmatmul[LeftT, OutT](a: _c.CanRMatmul[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b @ a`."""
    return b @ a


def do_rtruediv[LeftT, OutT](a: _c.CanRTruediv[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b / a`."""
    return b / a


def do_rfloordiv[LeftT, OutT](a: _c.CanRFloordiv[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b // a`."""
    return b // a


def do_rmod[LeftT, OutT](a: _c.CanRMod[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b % a`."""
    return b % a


def do_rdivmod[LeftT, OutT](a: _c.CanRDivmod[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `divmod(b, a)`."""
    return divmod(b, a)


def do_rpow[LeftT, OutT](a: _c.CanRPow[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b ** a`."""
    return b**a


def do_rlshift[LeftT, OutT](a: _c.CanRLshift[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b << a`."""
    return b << a


def do_rrshift[LeftT, OutT](a: _c.CanRRshift[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b >> a`."""
    return b >> a


def do_rand[LeftT, OutT](a: _c.CanRAnd[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b & a`."""
    return b & a


def do_rxor[LeftT, OutT](a: _c.CanRXor[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b ^ a`."""
    return b ^ a


def do_ror[LeftT, OutT](a: _c.CanROr[LeftT, OutT], b: LeftT, /) -> OutT:
    """Same as `b | a`."""
    return b | a


# augmented ops
do_iadd: Final[_d.DoesIAdd] = cast("_d.DoesIAdd", _o.iadd)
do_isub: Final[_d.DoesISub] = cast("_d.DoesISub", _o.isub)
do_imul: Final[_d.DoesIMul] = cast("_d.DoesIMul", _o.imul)
do_imatmul: Final[_d.DoesIMatmul] = cast("_d.DoesIMatmul", _o.imatmul)
do_itruediv: Final[_d.DoesITruediv] = cast("_d.DoesITruediv", _o.itruediv)
do_ifloordiv: Final[_d.DoesIFloordiv] = cast("_d.DoesIFloordiv", _o.ifloordiv)
do_imod: Final[_d.DoesIMod] = cast("_d.DoesIMod", _o.imod)
do_ipow: Final[_d.DoesIPow] = cast("_d.DoesIPow", _o.ipow)
do_ilshift: Final[_d.DoesILshift] = cast("_d.DoesILshift", _o.ilshift)
do_irshift: Final[_d.DoesIRshift] = cast("_d.DoesIRshift", _o.irshift)
do_iand: Final[_d.DoesIAnd] = cast("_d.DoesIAnd", _o.iand)
do_ixor: Final[_d.DoesIXor] = cast("_d.DoesIXor", _o.ixor)
do_ior: Final[_d.DoesIOr] = cast("_d.DoesIOr", _o.ior)

# unary ops
do_neg: Final[_d.DoesNeg] = _o.neg
do_pos: Final[_d.DoesPos] = _o.pos
do_abs: Final[_d.DoesAbs] = abs
do_invert: Final[_d.DoesInvert] = _o.invert

# fingerprinting
do_hash: Final[_d.DoesHash] = hash
do_index: Final[_d.DoesIndex] = _o.index

# rounding
# (the typeshed stubs for `round` are unnecessarily strict)
do_round: Final[_d.DoesRound] = cast("_d.DoesRound", round)
do_trunc: Final[_d.DoesTrunc] = math.trunc
do_floor: Final[_d.DoesFloor] = cast("_d.DoesFloor", math.floor)
do_ceil: Final[_d.DoesCeil] = cast("_d.DoesCeil", math.ceil)
