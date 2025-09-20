from collections.abc import Callable
from typing import Any, Concatenate, Generic, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc
import optype.typing as opt

from scipy._lib._util import _RichResult

_FloatT = TypeVar("_FloatT", bound=npc.floating, default=np.float64)
_FloatT_co = TypeVar("_FloatT_co", bound=npc.floating, default=np.float64, covariant=True)
_ShapeT = TypeVar("_ShapeT", bound=onp.AtLeast1D, default=onp.AtLeast0D[Any])
_ShapeT_co = TypeVar("_ShapeT_co", bound=onp.AtLeast1D, default=onp.AtLeast0D[Any], covariant=True)
_ShapeT2_co = TypeVar("_ShapeT2_co", bound=onp.AtLeast2D, default=onp.AtLeast0D[Any], covariant=True)

_Ignored: TypeAlias = object

_Function00: TypeAlias = Callable[Concatenate[_FloatT, ...], onp.ToFloat]
_Function11: TypeAlias = Callable[Concatenate[onp.Array1D[_FloatT], ...], onp.ToFloat1D]
_FunctionNN: TypeAlias = Callable[Concatenate[onp.ArrayND[_FloatT, Any], ...], onp.ToFloatND]

@type_check_only
class _Tolerances(TypedDict, total=False):
    rtol: onp.ToFloat
    atol: onp.ToFloat

@type_check_only
class _DerivativeResult0D(_RichResult, Generic[_FloatT_co]):
    success: np.bool_
    status: np.int32
    nfev: np.int32
    nit: np.int32
    x: _FloatT_co
    df: _FloatT_co
    error: _FloatT_co

@type_check_only
class _DerivativeResultND(_RichResult, Generic[_FloatT_co, _ShapeT_co]):
    success: onp.Array[_ShapeT_co, np.bool_]
    status: onp.Array[_ShapeT_co, np.int32]
    nfev: onp.Array[_ShapeT_co, np.int32]
    nit: onp.Array[_ShapeT_co, np.int32]
    x: onp.Array[_ShapeT_co, _FloatT_co]
    df: onp.Array[_ShapeT_co, _FloatT_co]
    error: onp.Array[_ShapeT_co, _FloatT_co]

@type_check_only
class _JacobianResult(_RichResult, Generic[_FloatT_co, _ShapeT_co]):
    status: onp.Array[_ShapeT_co, np.int32]
    df: onp.Array[_ShapeT_co, _FloatT_co]
    error: onp.Array[_ShapeT_co, _FloatT_co]
    nit: onp.Array[_ShapeT_co, np.int32]
    nfev: onp.Array[_ShapeT_co, np.int32]
    success: onp.Array[_ShapeT_co, np.bool_]

@type_check_only
class _HessianResult(_RichResult, Generic[_FloatT_co, _ShapeT2_co]):
    status: onp.Array[_ShapeT2_co, np.int32]
    error: onp.Array[_ShapeT2_co, _FloatT_co]
    nfev: onp.Array[_ShapeT2_co, np.int32]
    success: onp.Array[_ShapeT2_co, np.bool_]
    ddf: onp.Array[_ShapeT2_co, _FloatT_co]

###

@overload  # 0-d float64
def derivative(
    f: _Function00[np.float64],
    x: float | np.float64 | npc.integer | onp.CanArray0D[np.float64 | npc.integer],
    *,
    args: tuple[onp.ToScalar, ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: opt.AnyInt = 10,
    order: opt.AnyInt = 8,
    initial_step: onp.ToFloat = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt = 0,
    preserve_shape: bool = False,
    callback: Callable[[_DerivativeResult0D[np.float64]], _Ignored] | None = None,
) -> _DerivativeResult0D[np.float64]: ...
@overload  # 0-d <known>
def derivative(
    f: _Function00[_FloatT],
    x: _FloatT | onp.CanArray0D[_FloatT],
    *,
    args: tuple[onp.ToScalar, ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: opt.AnyInt = 10,
    order: opt.AnyInt = 8,
    initial_step: onp.ToFloat = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt = 0,
    preserve_shape: bool = False,
    callback: Callable[[_DerivativeResult0D[_FloatT]], _Ignored] | None = None,
) -> _DerivativeResult0D[_FloatT]: ...
@overload  # n-d <known>
def derivative(
    f: _FunctionNN[_FloatT],
    x: onp.CanArrayND[_FloatT, _ShapeT],
    *,
    args: tuple[onp.ToScalar | onp.ToArray1D | onp.CanArray[_ShapeT], ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: opt.AnyInt = 10,
    order: opt.AnyInt = 8,
    initial_step: onp.ToFloat | onp.ToFloat1D | onp.CanArrayND[npc.floating, _ShapeT] = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustInt1D | onp.CanArrayND[npc.integer, _ShapeT] = 0,
    preserve_shape: bool = False,
    callback: Callable[[_DerivativeResultND[_FloatT, _ShapeT]], _Ignored] | None = None,
) -> _DerivativeResultND[_FloatT, _ShapeT]: ...
@overload  # 1-d <unknown>
def derivative(
    f: _Function11[_FloatT],
    x: onp.ToFloatStrict1D,
    *,
    args: tuple[onp.ToScalar | onp.ToArrayStrict1D, ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: opt.AnyInt = 10,
    order: opt.AnyInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatStrict1D = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustIntStrict1D = 0,
    preserve_shape: bool = False,
    callback: Callable[[_DerivativeResultND[_FloatT, tuple[int]]], _Ignored] | None = None,
) -> _DerivativeResultND[_FloatT, tuple[int]]: ...
@overload  # n-d <unknown>
def derivative(
    f: _FunctionNN[_FloatT],
    x: onp.ToFloatND,
    *,
    args: tuple[onp.ToScalar | onp.ToArrayND, ...] = (),
    tolerances: _Tolerances | None = None,
    maxiter: opt.AnyInt = 10,
    order: opt.AnyInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatND = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustIntND = 0,
    preserve_shape: bool = False,
    callback: Callable[[_DerivativeResultND[_FloatT, _ShapeT]], _Ignored] | None = None,
) -> _DerivativeResultND[_FloatT, _ShapeT]: ...

#
def jacobian(
    f: Callable[[onp.Array[Any, _FloatT]], onp.ToFloat | onp.ToFloatND],
    x: onp.ToFloatND,
    *,
    tolerances: _Tolerances | None = None,
    maxiter: opt.AnyInt = 10,
    order: opt.AnyInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatND = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustIntND = 0,
) -> _JacobianResult[_FloatT, onp.AtLeast1D]: ...

#
def hessian(
    f: Callable[[onp.Array[Any, _FloatT]], onp.ToFloat | onp.ToFloatND],
    x: onp.ToFloatND,
    *,
    tolerances: _Tolerances | None = None,
    maxiter: opt.AnyInt = 10,
    order: opt.AnyInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatND = 0.5,
    step_factor: onp.ToFloat = 2.0,
) -> _HessianResult[_FloatT, onp.AtLeast2D]: ...
