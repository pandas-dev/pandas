from collections.abc import Callable, Sequence
from typing import Final, Literal, TypeAlias, TypedDict, Unpack, overload, type_check_only
from typing_extensions import TypeAliasType, TypeVar, TypeVarTuple, deprecated

import numpy as np
import optype as op
import optype.numpy as onp

from scipy.sparse.linalg import LinearOperator

__all__ = ["LbfgsInvHessProduct", "fmin_l_bfgs_b"]

_T = TypeVar("_T")
_Ts = TypeVarTuple("_Ts", default=Unpack[tuple[()]])

_Fn = TypeAliasType("_Fn", Callable[[onp.Array1D[np.float64], *_Ts], _T], type_params=(_T, _Ts))

_Bounds: TypeAlias = Sequence[tuple[onp.ToFloat | None, onp.ToFloat | None]]
_FMinResult: TypeAlias = tuple[onp.Array1D[np.float64], float, _InfoDict]

_ToFloatOr1D: TypeAlias = onp.ToFloat | onp.ToFloat1D
_ToFloatAnd1D: TypeAlias = tuple[onp.ToFloat, onp.ToFloat1D]

_Ignored: TypeAlias = object
_NoValueType: TypeAlias = op.JustObject

@type_check_only
class _InfoDict(TypedDict):
    grad: onp.Array1D[np.float64]
    task: str  # undocumented
    funcalls: int
    nit: int
    warnflag: Literal[0, 1, 2]

###

status_messages: Final[dict[int, str]] = ...
task_messages: Final[dict[int, str]] = ...

class LbfgsInvHessProduct(LinearOperator[np.float64]):
    sk: Final[onp.Array2D[np.float64]]
    yk: Final[onp.Array2D[np.float64]]
    n_corrs: Final[int]
    rho: Final[float]

    def __init__(self, /, sk: onp.ToFloat2D, yk: onp.ToFloat2D) -> None: ...
    def todense(self, /) -> onp.Array2D[np.float64]: ...

@overload  # no args, no fprime, no approx_grad
def fmin_l_bfgs_b(
    func: _Fn[_ToFloatAnd1D],
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[()] = (),
    approx_grad: onp.ToFalse = 0,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    iprint: _NoValueType = ...,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    disp: _NoValueType = ...,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
@overload  # args, no fprime, no approx_grad
def fmin_l_bfgs_b(
    func: _Fn[_ToFloatAnd1D, *_Ts],
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    approx_grad: onp.ToFalse = 0,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    iprint: _NoValueType = ...,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    disp: _NoValueType = ...,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
@overload  # fprime, no approx_grad
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts],
    x0: _ToFloatOr1D,
    fprime: _Fn[onp.ToFloat1D, *_Ts],
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    approx_grad: onp.ToFalse = 0,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    iprint: _NoValueType = ...,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    disp: _NoValueType = ...,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
@overload  # no fprime, approx_grad (keyword)
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts],
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    *,
    approx_grad: onp.ToTrue,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    iprint: _NoValueType = ...,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    disp: _NoValueType = ...,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
@overload  # no fprime, unknown approx_grad  (keyword)
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts] | _Fn[_ToFloatAnd1D, *_Ts],
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    *,
    approx_grad: onp.ToBool,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    iprint: _NoValueType = ...,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    disp: _NoValueType = ...,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
@overload  # iprint
@deprecated("The `iprint` keyword is deprecated and will be removed from SciPy 1.18.0.")
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts] | _Fn[_ToFloatAnd1D, *_Ts],
    x0: _ToFloatOr1D,
    fprime: _Fn[onp.ToFloat1D, *_Ts] | None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    approx_grad: onp.ToBool = 0,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    *,
    iprint: int,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    disp: _NoValueType = ...,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
@overload  # disp
@deprecated("The `disp` keyword is deprecated and will be removed from SciPy 1.18.0.")
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts] | _Fn[_ToFloatAnd1D, *_Ts],
    x0: _ToFloatOr1D,
    fprime: _Fn[onp.ToFloat1D, *_Ts] | None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    approx_grad: onp.ToBool = 0,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    iprint: _NoValueType = ...,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    *,
    disp: int,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
@overload  # iprint and disp
@deprecated("The `iprint` and `disp` keywords are deprecated and will be removed from SciPy 1.18.0.")
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts] | _Fn[_ToFloatAnd1D, *_Ts],
    x0: _ToFloatOr1D,
    fprime: _Fn[onp.ToFloat1D, *_Ts] | None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    approx_grad: onp.ToBool = 0,
    bounds: _Bounds | None = None,
    m: onp.ToJustInt = 10,
    factr: onp.ToFloat = 1e7,
    pgtol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat = 1e-8,
    *,
    iprint: int,
    maxfun: onp.ToJustInt = 15_000,
    maxiter: onp.ToJustInt = 15_000,
    disp: int,
    callback: _Fn[_Ignored] | None = None,
    maxls: onp.ToJustInt = 20,
) -> _FMinResult: ...
