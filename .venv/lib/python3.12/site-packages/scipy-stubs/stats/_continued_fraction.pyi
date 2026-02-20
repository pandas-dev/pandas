from _typeshed import Incomplete
from collections.abc import Callable
from types import ModuleType
from typing import Any, Concatenate, Generic, TypeAlias, TypeVar, TypedDict, overload, type_check_only

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._lib._util import _RichResult

###

_T = TypeVar("_T", bound=npt.ArrayLike)
_AT = TypeVar("_AT", bound=npt.ArrayLike)
_BT = TypeVar("_BT", bound=npt.ArrayLike)

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

_FnCoef: TypeAlias = Callable[Concatenate[int, ...], _T]

@type_check_only
class _Tolerances(TypedDict, total=False):
    eps: float | None
    tiny: float | None

@type_check_only
class _ContinuedFraction0D(_RichResult):
    success: np.bool_
    status: np.int32
    f: np.float64
    nit: np.int32
    nfev: np.int32

@type_check_only
class _ContinuedFractionND(_RichResult, Generic[_ShapeT_co]):
    success: onp.ArrayND[np.bool_, _ShapeT_co]
    status: onp.ArrayND[np.int32, _ShapeT_co]
    f: onp.ArrayND[np.float64, _ShapeT_co]
    nit: onp.ArrayND[np.int32, _ShapeT_co]
    nfev: onp.ArrayND[np.int32, _ShapeT_co]

###

# undocumented
@overload
def _logaddexp(x: onp.ToFloat, y: onp.ToFloat, xp: None = None) -> np.float64 | Any: ...
@overload
def _logaddexp(x: onp.ToComplex, y: onp.ToComplex, xp: None = None) -> np.complex128 | Any: ...
@overload
def _logaddexp(x: onp.ToFloatND, y: onp.ToFloat | onp.ToFloatND, xp: None = None) -> onp.ArrayND[np.float64 | Any]: ...
@overload
def _logaddexp(x: onp.ToComplexND, y: onp.ToComplex | onp.ToComplexND, xp: None = None) -> onp.ArrayND[np.complex128 | Any]: ...
@overload
def _logaddexp(x: Incomplete, y: Incomplete, xp: ModuleType) -> Incomplete: ...

# undocumented
def _continued_fraction_iv(
    a: _FnCoef[_AT], b: _FnCoef[_BT], args: tuple[object, ...], tolerances: _Tolerances | None, maxiter: int, log: bool
) -> tuple[
    _FnCoef[_AT],  # a
    _FnCoef[_BT],  # b
    tuple[Any, ...],  # args
    float | None,  # eps
    float | None,  # tiny
    int,  # maxiter
    bool,  # log
    _AT,  # a0
    _BT,  # b0
    tuple[Any, ...],  # shape
    np.dtype[Any],  # dtype
    ModuleType,  # xp
]: ...

# undocumented
@overload
def _continued_fraction(
    a: _FnCoef[onp.ToFloat],
    b: _FnCoef[onp.ToFloat],
    *,
    args: tuple[object, ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: int = 100,
    log: bool = False,
) -> _ContinuedFraction0D: ...
@overload
def _continued_fraction(
    a: _FnCoef[onp.ArrayND[npc.floating | npc.integer | np.bool_, _ShapeT]],
    b: _FnCoef[onp.ArrayND[npc.floating | npc.integer | np.bool_, _ShapeT] | onp.ToFloat],
    *,
    args: tuple[object, ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: int = 100,
    log: bool = False,
) -> _ContinuedFractionND[_ShapeT]: ...
@overload
def _continued_fraction(
    a: _FnCoef[onp.ArrayND[npc.floating | npc.integer | np.bool_, _ShapeT] | onp.ToFloat],
    b: _FnCoef[onp.ArrayND[npc.floating | npc.integer | np.bool_, _ShapeT]],
    *,
    args: tuple[object, ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: int = 100,
    log: bool = False,
) -> _ContinuedFractionND[_ShapeT]: ...
