from collections.abc import Callable
from typing import Literal, overload

import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ordqz", "qz"]

###

type _FloatND = onp.ArrayND[npc.floating]
type _ComplexND = onp.ArrayND[npc.complexfloating]
type _InexactND = onp.ArrayND[npc.inexact]

type _Tuple4[T] = tuple[T, T, T, T]
type _Tuple2C12[T2, T1] = tuple[T2, T2, _ComplexND, T1, T2, T2]

type _OutputReal = Literal["real", "r"]
type _OutputComplex = Literal["complex", "c"]

type _Sort = Literal["lhp", "rhp", "iuc", "ouc"] | Callable[[float, float], bool]

###

# NOTE: `sort` will raise `ValueError` if not `None`.
@overload  # float, {"real"}
def qz(
    A: onp.ToFloatND,
    B: onp.ToFloatND,
    output: _OutputReal = "real",
    lwork: onp.ToJustInt | None = None,
    sort: None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> _Tuple4[_FloatND]: ...
@overload  # complex, {"real"}
def qz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    output: _OutputReal = "real",
    lwork: onp.ToJustInt | None = None,
    sort: None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> _Tuple4[_InexactND]: ...
@overload  # complex, {"complex"}
def qz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    output: _OutputComplex,
    lwork: onp.ToJustInt | None = None,
    sort: None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> _Tuple4[_ComplexND]: ...

#
@overload  # float, {"real"}
def ordqz(
    A: onp.ToFloatND,
    B: onp.ToFloatND,
    sort: _Sort = "lhp",
    output: _OutputReal = "real",
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> _Tuple2C12[_FloatND, _FloatND]: ...
@overload  # complex, {"real"}
def ordqz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    sort: _Sort = "lhp",
    output: _OutputReal = "real",
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> _Tuple2C12[_InexactND, _InexactND]: ...
@overload  # complex, {"complex"} (positional)
def ordqz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    sort: _Sort,
    output: _OutputComplex,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> _Tuple2C12[_ComplexND, _ComplexND]: ...
@overload  # complex, {"complex"} (keyword)
def ordqz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    sort: _Sort = "lhp",
    *,
    output: _OutputComplex,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> _Tuple2C12[_ComplexND, _ComplexND]: ...
