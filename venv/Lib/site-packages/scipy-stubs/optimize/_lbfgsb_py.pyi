from _typeshed import Unused
from collections.abc import Callable, Sequence
from typing import Final, Literal, Self, TypedDict, Unpack, overload, override, type_check_only
from typing_extensions import TypeVarTuple

import numpy as np
import optype.numpy as onp

from scipy.sparse.linalg import LinearOperator

__all__ = ["LbfgsInvHessProduct", "fmin_l_bfgs_b"]

_Ts = TypeVarTuple("_Ts", default=Unpack[tuple[()]])

type _Fn[T, *Ts] = Callable[[onp.Array1D[np.float64], *Ts], T]

type _Bounds = Sequence[tuple[onp.ToFloat | None, onp.ToFloat | None]]
type _FMinResult = tuple[onp.Array1D[np.float64], float, _InfoDict]

type _ToFloatOr1D = onp.ToFloat | onp.ToFloat1D
type _ToFloatAnd1D = tuple[onp.ToFloat, onp.ToFloat1D]

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

    @override
    def __new__(cls, /, sk: onp.ToFloat2D, yk: onp.ToFloat2D) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, sk: onp.ToFloat2D, yk: onp.ToFloat2D) -> None: ...  # pyrefly:ignore[bad-override]

    #
    def todense(self, /) -> onp.Array2D[np.float64]: ...

@overload  # no args, no fprime, no approx_grad
def fmin_l_bfgs_b(
    func: _Fn[_ToFloatAnd1D],
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[()] = (),
    approx_grad: onp.ToFalse = 0,
    bounds: _Bounds | None = None,
    m: int = 10,
    factr: float = 1e7,
    pgtol: float = 1e-5,
    epsilon: float = 1e-8,
    # *,
    maxfun: int = 15_000,
    maxiter: int = 15_000,
    callback: _Fn[Unused] | None = None,
    maxls: int = 20,
) -> _FMinResult: ...
@overload  # args, no fprime, no approx_grad
def fmin_l_bfgs_b(
    func: _Fn[_ToFloatAnd1D, *_Ts],  # ty:ignore[invalid-type-arguments]
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    approx_grad: onp.ToFalse = 0,
    bounds: _Bounds | None = None,
    m: int = 10,
    factr: float = 1e7,
    pgtol: float = 1e-5,
    epsilon: float = 1e-8,
    # *,
    maxfun: int = 15_000,
    maxiter: int = 15_000,
    callback: _Fn[Unused] | None = None,
    maxls: int = 20,
) -> _FMinResult: ...
@overload  # fprime, no approx_grad
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts],  # ty:ignore[invalid-type-arguments]
    x0: _ToFloatOr1D,
    fprime: _Fn[onp.ToFloat1D, *_Ts],  # ty:ignore[invalid-type-arguments]
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    approx_grad: onp.ToFalse = 0,
    bounds: _Bounds | None = None,
    m: int = 10,
    factr: float = 1e7,
    pgtol: float = 1e-5,
    epsilon: float = 1e-8,
    # *,
    maxfun: int = 15_000,
    maxiter: int = 15_000,
    callback: _Fn[Unused] | None = None,
    maxls: int = 20,
) -> _FMinResult: ...
@overload  # no fprime, approx_grad (keyword)
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts],  # ty:ignore[invalid-type-arguments]
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    *,
    approx_grad: onp.ToTrue,
    bounds: _Bounds | None = None,
    m: int = 10,
    factr: float = 1e7,
    pgtol: float = 1e-5,
    epsilon: float = 1e-8,
    maxfun: int = 15_000,
    maxiter: int = 15_000,
    callback: _Fn[Unused] | None = None,
    maxls: int = 20,
) -> _FMinResult: ...
@overload  # no fprime, unknown approx_grad  (keyword)
def fmin_l_bfgs_b(
    func: _Fn[onp.ToFloat, *_Ts] | _Fn[_ToFloatAnd1D, *_Ts],  # ty:ignore[invalid-type-arguments]
    x0: _ToFloatOr1D,
    fprime: None = None,
    args: tuple[*_Ts] = ...,  # stubdefaulter: ignore[missing-default]
    *,
    approx_grad: onp.ToBool,
    bounds: _Bounds | None = None,
    m: int = 10,
    factr: float = 1e7,
    pgtol: float = 1e-5,
    epsilon: float = 1e-8,
    maxfun: int = 15_000,
    maxiter: int = 15_000,
    callback: _Fn[Unused] | None = None,
    maxls: int = 20,
) -> _FMinResult: ...
