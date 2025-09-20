from collections.abc import Callable, Mapping
from typing import Concatenate, Generic, Literal, TypeAlias, TypeVar, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._nonlin import InverseJacobian
from ._optimize import OptimizeResult as _OptimizeResult
from scipy.sparse.linalg import LinearOperator

__all__ = ["root"]

_RootMethod: TypeAlias = Literal[
    "hybr",
    "lm",
    "broyden1",
    "broyden2",
    "anderson",
    "linearmixing",
    "diagbroyden",
    "excitingmixing",
    "krylov",
    "df-sane",
]  # fmt: skip

_JacOptionsT = TypeVar("_JacOptionsT", bound=Mapping[str, object])

@type_check_only
class _RootOptionsHybr(TypedDict, total=False):
    col_deriv: onp.ToBool
    xtol: onp.ToFloat
    maxffev: onp.ToJustInt
    band: tuple[onp.ToJustInt, onp.ToJustInt]
    eps: onp.ToFloat
    factor: onp.ToFloat
    diag: onp.ToFloat1D

@type_check_only
class _RootOptionsLM(TypedDict, total=False):
    col_deriv: onp.ToBool
    ftol: onp.ToFloat
    xtol: onp.ToFloat
    gtol: onp.ToFloat
    maxiter: onp.ToJustInt
    eps: onp.ToFloat
    factor: onp.ToFloat
    diag: onp.ToFloat1D

@type_check_only
class _RootOptionsJacBase(TypedDict, Generic[_JacOptionsT], total=False):
    nit: onp.ToJustInt
    disp: onp.ToBool
    maxiter: onp.ToJustInt
    ftol: onp.ToFloat
    fatol: onp.ToFloat
    xtol: onp.ToFloat
    xatol: onp.ToFloat
    tol_norm: Callable[[onp.Array1D[np.float64]], onp.ToFloat]
    line_search: Literal["armijo", "wolfe"] | None
    jac_options: _JacOptionsT

@type_check_only
class _JacOptionsBase(TypedDict, total=False):
    alpha: onp.ToFloat

@type_check_only
class _JacOptionsBroyden(_JacOptionsBase, total=False):
    reduction_method: Literal["restart", "simple", "svd"] | tuple[Literal["svd"], onp.ToJustInt]
    max_rank: onp.ToJustInt

@type_check_only
class _JacOptionsAnderson(_JacOptionsBase, total=False):
    M: onp.ToFloat
    w0: onp.ToFloat

@type_check_only
class _JacOptionsExcitingMixing(_JacOptionsBase, total=False):
    alphamax: onp.ToFloat

@type_check_only
class _JacOptionsKrylov(TypedDict, total=False):
    rdiff: onp.ToFloat
    method: Literal["lgmres", "gmres", "bicgstab", "cgs", "minres", "tfqmr"] | Callable[..., onp.ToFloatND]
    inner_M: LinearOperator | InverseJacobian
    inner_rtol: onp.ToFloat
    inner_atol: onp.ToFloat
    inner_maxiter: onp.ToJustInt
    outer_k: onp.ToJustInt

_RootOptionsBroyden: TypeAlias = _RootOptionsJacBase[_JacOptionsBroyden]
_RootOptionsAnderson: TypeAlias = _RootOptionsJacBase[_JacOptionsAnderson]
_RootOptionsLinearMixing: TypeAlias = _RootOptionsJacBase[_JacOptionsBase]
_RootOptionsDiagBroyden: TypeAlias = _RootOptionsJacBase[_JacOptionsBase]
_RootOptionsExcitingMixing: TypeAlias = _RootOptionsJacBase[_JacOptionsExcitingMixing]
_RootOptionsKrylov: TypeAlias = _RootOptionsJacBase[_JacOptionsKrylov]

@type_check_only
class _RootOptionsDFSane(TypedDict, total=False):
    ftol: onp.ToFloat
    fatol: onp.ToFloat
    fnorm: Callable[[onp.ArrayND[np.float64]], onp.ToFloat]
    maxfev: onp.ToJustInt
    disp: onp.ToBool
    eta_strategy: Callable[[int, onp.ArrayND[np.float64], onp.ArrayND[np.float64]], onp.ToFloat1D]
    sigma_eps: onp.ToFloat
    sigma_0: onp.ToFloat
    M: onp.ToJustInt
    line_search: Literal["cruz", "cheng"]

_RootOptions: TypeAlias = (
    _RootOptionsHybr
    | _RootOptionsLM
    | _RootOptionsBroyden
    | _RootOptionsAnderson
    | _RootOptionsLinearMixing
    | _RootOptionsDiagBroyden
    | _RootOptionsExcitingMixing
    | _RootOptionsKrylov
    | _RootOptionsDFSane
)

###

class OptimizeResult(_OptimizeResult):
    x: onp.ArrayND[npc.floating]
    success: bool
    message: str
    nfev: int

@overload  # `jac` not given (None), False, or a callable
def root(
    fun: Callable[Concatenate[onp.ArrayND[np.float64], ...], onp.ToFloatND],
    x0: onp.ToFloatND,
    args: tuple[object, ...] = (),
    method: _RootMethod = "hybr",
    jac: Callable[Concatenate[onp.ArrayND[np.float64], ...], onp.ToFloatND] | onp.ToFalse | None = None,
    tol: onp.ToFloat | None = None,
    callback: Callable[[onp.ArrayND[np.float64], onp.ArrayND[np.float64]], None] | None = None,
    options: _RootOptions | None = None,
) -> _OptimizeResult: ...
@overload  # `jac` truthy (positional)
def root(
    fun: Callable[Concatenate[onp.ArrayND[np.float64], ...], tuple[onp.ToFloatND, onp.ToFloatND]],
    x0: onp.ToFloatND,
    args: tuple[object, ...],
    method: _RootMethod,
    jac: onp.ToTrue,
    tol: onp.ToFloat | None = None,
    callback: Callable[[onp.ArrayND[np.float64], onp.ArrayND[np.float64]], None] | None = None,
    options: _RootOptions | None = None,
) -> _OptimizeResult: ...
@overload  # `jac` truthy (keyword)
def root(
    fun: Callable[Concatenate[onp.ArrayND[np.float64], ...], tuple[onp.ToFloatND, onp.ToFloatND]],
    x0: onp.ToFloatND,
    args: tuple[object, ...] = (),
    method: _RootMethod = "hybr",
    *,
    jac: onp.ToTrue,
    tol: onp.ToFloat | None = None,
    callback: Callable[[onp.ArrayND[np.float64], onp.ArrayND[np.float64]], None] | None = None,
    options: _RootOptions | None = None,
) -> _OptimizeResult: ...
