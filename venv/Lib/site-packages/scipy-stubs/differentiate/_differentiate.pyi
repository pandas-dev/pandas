from _typeshed import ConvertibleToInt, Unused
from collections.abc import Callable, Mapping
from typing import Any, Concatenate, Generic, Literal, TypedDict, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._lib._util import _RichResult

_FloatT_co = TypeVar("_FloatT_co", bound=npc.floating, default=np.float64, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, *tuple[int, ...]], default=tuple[Any, ...], covariant=True)
_ShapeT2_co = TypeVar("_ShapeT2_co", bound=tuple[int, int, *tuple[int, ...]], default=tuple[Any, ...], covariant=True)

type _Function00[FloatT: npc.floating] = Callable[Concatenate[FloatT, ...], onp.ToFloat]
type _Function11[FloatT: npc.floating] = Callable[Concatenate[onp.Array1D[FloatT], ...], onp.ToFloat1D]
type _FunctionNN[FloatT: npc.floating] = Callable[Concatenate[onp.ArrayND[FloatT, Any], ...], onp.ToFloatND]

@type_check_only
class _Tolerances(TypedDict, total=False):
    rtol: onp.ToFloat
    atol: onp.ToFloat

@type_check_only
class _DerivativeResult0D(_RichResult, Generic[_FloatT_co]):
    success: np.bool
    status: np.int32
    nfev: np.int32
    nit: np.int32
    x: _FloatT_co
    df: _FloatT_co
    error: _FloatT_co

@type_check_only
class _DerivativeResultND(_RichResult, Generic[_FloatT_co, _ShapeT_co]):
    success: onp.Array[_ShapeT_co, np.bool]
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
    success: onp.Array[_ShapeT_co, np.bool]

@type_check_only
class _HessianResult(_RichResult, Generic[_FloatT_co, _ShapeT2_co]):
    status: onp.Array[_ShapeT2_co, np.int32]
    error: onp.Array[_ShapeT2_co, _FloatT_co]
    nfev: onp.Array[_ShapeT2_co, np.int64]
    success: onp.Array[_ShapeT2_co, np.bool]
    ddf: onp.Array[_ShapeT2_co, _FloatT_co]

###

@overload  # 0-d float64
def derivative(
    f: _Function00[np.float64],
    x: float | np.float64 | npc.integer | onp.CanArray0[np.float64 | npc.integer],
    *,
    args: tuple[onp.ToScalar, ...] = (),
    kwargs: Mapping[str, onp.ToScalar] | None = None,
    tolerances: _Tolerances | None = None,
    maxiter: ConvertibleToInt = 10,
    order: ConvertibleToInt = 8,
    initial_step: onp.ToFloat = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt = 0,
    preserve_shape: Literal[False] = False,
    callback: Callable[[_DerivativeResult0D[np.float64]], Unused] | None = None,
) -> _DerivativeResult0D[np.float64]: ...
@overload  # 0-d <known>
def derivative[FloatT: npc.floating](
    f: _Function00[FloatT],
    x: FloatT | onp.CanArray0[FloatT | npc.integer],
    *,
    args: tuple[onp.ToScalar, ...] = (),
    kwargs: Mapping[str, onp.ToScalar] | None = None,
    tolerances: _Tolerances | None = None,
    maxiter: ConvertibleToInt = 10,
    order: ConvertibleToInt = 8,
    initial_step: onp.ToFloat = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt = 0,
    preserve_shape: Literal[False] = False,
    callback: Callable[[_DerivativeResult0D[FloatT]], Unused] | None = None,
) -> _DerivativeResult0D[FloatT]: ...
@overload  # 1-d <unknown>
def derivative[FloatT: npc.floating](
    f: _Function11[FloatT],
    x: onp.ToFloatStrict1D,
    *,
    args: tuple[onp.ToScalar | onp.ToArrayStrict1D, ...] = (),
    kwargs: Mapping[str, onp.ToScalar | onp.ToArrayStrict1D] | None = None,
    tolerances: _Tolerances | None = None,
    maxiter: ConvertibleToInt = 10,
    order: ConvertibleToInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatStrict1D = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustIntStrict1D = 0,
    preserve_shape: Literal[False] = False,
    callback: Callable[[_DerivativeResultND[FloatT, tuple[int]]], Unused] | None = None,
) -> _DerivativeResultND[FloatT, tuple[int]]: ...
@overload  # n-d <known>
def derivative[FloatT: npc.floating](
    f: _FunctionNN[FloatT],
    x: FloatT | onp.ToArrayND[FloatT],
    *,
    args: tuple[onp.ToScalar | onp.ToArrayND, ...] = (),
    kwargs: Mapping[str, onp.ToScalar | onp.ToArrayND] | None = None,
    tolerances: _Tolerances | None = None,
    maxiter: ConvertibleToInt = 10,
    order: ConvertibleToInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatND = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustIntND = 0,
    preserve_shape: bool = False,
    callback: Callable[[_DerivativeResultND[FloatT]], Unused] | None = None,
) -> _DerivativeResultND[FloatT]: ...
@overload  # n-d <unknown>
def derivative[FloatT: npc.floating](
    f: _FunctionNN[FloatT],
    x: onp.ToFloat | onp.ToFloatND,
    *,
    args: tuple[onp.ToScalar | onp.ToArrayND, ...] = (),
    kwargs: Mapping[str, onp.ToScalar | onp.ToArrayND] | None = None,
    tolerances: _Tolerances | None = None,
    maxiter: ConvertibleToInt = 10,
    order: ConvertibleToInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatND = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustIntND = 0,
    preserve_shape: bool = False,
    callback: Callable[[_DerivativeResultND[FloatT]], Unused] | None = None,
) -> _DerivativeResultND[FloatT]: ...

#
def jacobian[FloatT: npc.floating](
    f: Callable[[onp.Array[Any, FloatT]], onp.ToFloat | onp.ToFloatND],
    x: onp.ToFloatND,
    *,
    tolerances: _Tolerances | None = None,
    maxiter: ConvertibleToInt = 10,
    order: ConvertibleToInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatND = 0.5,
    step_factor: onp.ToFloat = 2.0,
    step_direction: onp.ToJustInt | onp.ToJustIntND = 0,
) -> _JacobianResult[FloatT, onp.AtLeast1D]: ...

#
def hessian[FloatT: npc.floating](
    f: Callable[[onp.Array[Any, FloatT]], onp.ToFloat | onp.ToFloatND],
    x: onp.ToFloatND,
    *,
    tolerances: _Tolerances | None = None,
    maxiter: ConvertibleToInt = 10,
    order: ConvertibleToInt = 8,
    initial_step: onp.ToFloat | onp.ToFloatND = 0.5,
    step_factor: onp.ToFloat = 2.0,
) -> _HessianResult[FloatT, onp.AtLeast2D]: ...
