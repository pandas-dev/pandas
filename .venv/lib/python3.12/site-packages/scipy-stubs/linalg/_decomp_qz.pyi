from collections.abc import Callable
from typing import Literal, TypeAlias, TypeVar, overload

import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ordqz", "qz"]

_T1 = TypeVar("_T1")
_T2 = TypeVar("_T2")
_Tuple4: TypeAlias = tuple[_T2, _T2, _T2, _T2]
_Tuple222: TypeAlias = tuple[_T2, _T2, _T1, _T1, _T2, _T2]

_FloatND: TypeAlias = onp.ArrayND[npc.floating]
_ComplexND: TypeAlias = onp.ArrayND[npc.complexfloating]
_InexactND: TypeAlias = onp.ArrayND[npc.inexact]

_OutputReal: TypeAlias = Literal["real", "r"]
_OutputComplex: TypeAlias = Literal["complex", "c"]

_Sort: TypeAlias = Literal["lhp", "rhp", "iuc", "ouc"] | Callable[[float, float], onp.ToBool]

###

# NOTE: `sort` will raise `ValueError` if not `None`.
@overload  # float, {"real"}
def qz(
    A: onp.ToFloatND,
    B: onp.ToFloatND,
    output: _OutputReal = "real",
    lwork: onp.ToJustInt | None = None,
    sort: None = None,
    overwrite_a: onp.ToBool = False,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Tuple4[_FloatND]: ...
@overload  # complex, {"real"}
def qz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    output: _OutputReal = "real",
    lwork: onp.ToJustInt | None = None,
    sort: None = None,
    overwrite_a: onp.ToBool = False,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Tuple4[_InexactND]: ...
@overload  # complex, {"complex"}
def qz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    output: _OutputComplex,
    lwork: onp.ToJustInt | None = None,
    sort: None = None,
    overwrite_a: onp.ToBool = False,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Tuple4[_ComplexND]: ...

#
@overload  # float, {"real"}
def ordqz(
    A: onp.ToFloatND,
    B: onp.ToFloatND,
    sort: _Sort = "lhp",
    output: _OutputReal = "real",
    overwrite_a: onp.ToBool = False,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Tuple222[_FloatND, _FloatND]: ...
@overload  # complex, {"real"}
def ordqz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    sort: _Sort = "lhp",
    output: _OutputReal = "real",
    overwrite_a: onp.ToBool = False,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Tuple222[_InexactND, _InexactND]: ...
@overload  # complex, {"complex"} (positional)
def ordqz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    sort: _Sort,
    output: _OutputComplex,
    overwrite_a: onp.ToBool = False,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Tuple222[_ComplexND, _ComplexND]: ...
@overload  # complex, {"complex"} (keyword)
def ordqz(
    A: onp.ToComplexND,
    B: onp.ToComplexND,
    sort: _Sort = "lhp",
    *,
    output: _OutputComplex,
    overwrite_a: onp.ToBool = False,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Tuple222[_ComplexND, _ComplexND]: ...
