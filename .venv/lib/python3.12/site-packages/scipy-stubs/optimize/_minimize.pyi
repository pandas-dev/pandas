from collections.abc import Callable, Mapping, Sequence
from typing import Concatenate, Final, Literal, LiteralString, Protocol, TypeAlias, TypeVar, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc
from numpy_typing_compat import ABCPolyBase

from ._hessian_update_strategy import HessianUpdateStrategy
from ._optimize import OptimizeResult as _OptimizeResult
from ._typing import Bound, Bounds, Constraint, Constraints, MethodMimimize, MethodMinimizeScalar
from scipy.sparse.linalg import LinearOperator

__all__ = ["minimize", "minimize_scalar"]

###

_T = TypeVar("_T")
_Float1DT = TypeVar("_Float1DT", bound=_Float1D)

_Tuple2: TypeAlias = tuple[_T, _T]
_Tuple3: TypeAlias = tuple[_T, _T, _T]
_Args: TypeAlias = tuple[object, ...]

_Floating: TypeAlias = float | npc.floating
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]

_RT = TypeVar("_RT")
# NOTE: `ABCPolyBase` is required to work around https://github.com/scipy/scipy-stubs/issues/465
_Fun0D: TypeAlias = Callable[Concatenate[float, ...], _RT] | Callable[Concatenate[np.float64, ...], _RT] | ABCPolyBase
_Fun1D: TypeAlias = Callable[Concatenate[_Float1D, ...], _RT]
_Fun1Dp: TypeAlias = Callable[Concatenate[_Float1D, _Float1D, ...], _RT]

_FDMethod: TypeAlias = Literal["2-point", "3-point", "cs"]

_ToBracket: TypeAlias = Sequence[_Tuple2[onp.ToFloat] | _Tuple3[onp.ToFloat]]
_ToBound: TypeAlias = _Tuple2[onp.ToFloat]
_Ignored: TypeAlias = object

_MinimizeScalarResultT = TypeVar("_MinimizeScalarResultT", bound=_MinimizeScalarResultBase)
_MinimizeScalarResultT_co = TypeVar("_MinimizeScalarResultT_co", bound=_MinimizeScalarResultBase, covariant=True)

@type_check_only
class _CallbackResult(Protocol):
    def __call__(self, /, intermediate_result: OptimizeResult) -> None: ...

@type_check_only
class _CallbackVector(Protocol):
    def __call__(self, /, xk: _Float1D) -> None: ...

@type_check_only
class _MinimizeMethodFun(Protocol):
    def __call__(self, fun: _Fun1D[onp.ToFloat], x0: onp.ToFloat1D, /, args: _Args) -> OptimizeResult: ...

@type_check_only
class _MinimizeScalarMethodFun(Protocol[_MinimizeScalarResultT_co]):
    def __call__(
        self, fun: _Fun0D[onp.ToFloat], /, *, args: _Args, bracket: _ToBracket, bound: _ToBound
    ) -> _MinimizeScalarResultT_co: ...

@type_check_only
class _MinimizeOptions(TypedDict, total=False):
    # Nelder-Mead, Powell, CG, BFGS, Newton-CG
    return_all: onp.ToBool
    # Nelder-Mead, Powell, CG, BFGS, L-BFGS-B, Newton-CG, TNC, COBYLA, COBYQA, SLSQP, trust-constr
    disp: onp.ToBool
    # Nelder-Mead, Powell, CG, BFGS, L-BFGS-B, Newton-CG, COBYLA, SLSQP, trust-constr
    maxiter: int
    # Nelder-Mead, Powell, COBYQA
    maxfev: int
    # TNC
    maxCGit: int
    offset: _Floating
    stepmx: _Floating
    accuracy: _Floating
    minfev: _Floating
    rescale: _Floating
    # L-BFGS-B, TNC
    maxfun: int
    # L-BFGS-B
    maxcor: int
    iprint: int
    maxls: int
    # Nelder-Mead
    initial_simplex: onp.ToFloatND
    adaptive: onp.ToBool
    xatol: _Floating
    fatol: _Floating
    # CG, BFGS, L-BFGS-B, dogleg, trust-ncg, trust-exact, TNC, trust-constr
    gtol: _Floating
    # Powell, Newton-CG, TNC, trust-constr
    xtol: _Floating
    # Powell, L-BFGS-B, TNC, SLSQP
    ftol: _Floating
    # BFGS
    xrtol: _Floating
    hess_inv0: onp.ArrayND[npc.floating]
    # COBYLA
    tol: _Floating
    catool: _Floating
    rhobeg: _Floating
    f_target: _Floating
    # COBYQA
    feasibility_tol: _Floating
    final_tr_radius: _Floating
    # Powell
    direc: onp.ArrayND[npc.floating]
    # CG, BFGS, Newton-CG, L-BFGS-B, TNC, SLSQP
    eps: _Floating | onp.ArrayND[npc.floating]
    # CG, BFGS, Newton-CG
    c1: _Floating
    c2: _Floating
    # CG, BFGS
    norm: _Floating
    # CG, BFGS, L-BFGS-B, TNC, SLSQP, trust-constr
    finite_diff_rel_step: onp.ToFloat | onp.ToFloatND
    # dogleg, trust-ncg, trust-exact
    initial_trust_radius: _Floating
    max_trust_radius: _Floating
    # COBYQA, trust-constr
    initial_tr_radius: _Floating
    # trust-constr
    barrier_tol: _Floating
    sparse_jacobian: onp.ToBool
    initial_constr_penalty: _Floating
    initial_barrier_parameter: _Floating
    initial_barrier_tolerance: _Floating
    factorization_method: Literal["NormalEquation", "AugmentedSystem", "QRFactorization", "SVDFactorization"]
    verbose: Literal[0, 1, 2, 3]
    # dogleg, trust-ncg, trust-exact, TNC
    eta: _Floating
    # trust-krylov
    inexact: onp.ToBool
    # TNC (list of floats), COBYQA (bool)
    scale: Sequence[_Floating] | onp.ToBool
    # trust-exact
    subproblem_maxiter: float

@type_check_only
class _MinimizeScalarOptionsCommon(TypedDict, total=False):
    maxiter: int
    disp: Literal[0, 1, 2, 3]

@type_check_only
class _MinimizeScalarOptionsBracketed(_MinimizeScalarOptionsCommon, TypedDict, total=False):
    xtol: _Floating

@type_check_only
class _MinimizeScalarOptionsBounded(_MinimizeScalarOptionsCommon, TypedDict, total=False):
    xatol: _Floating

@type_check_only
class _MinimizeScalarResultBase(_OptimizeResult):
    x: np.float64
    fun: np.float64

@type_check_only
class _MinimizeScalarResult(_MinimizeScalarResultBase):
    success: bool
    message: LiteralString
    nit: int
    nfev: int

###

MINIMIZE_METHODS: Final[list[MethodMimimize]] = ...
MINIMIZE_METHODS_NEW_CB: Final[list[MethodMimimize]] = ...
MINIMIZE_SCALAR_METHODS: Final[list[MethodMinimizeScalar]] = ...

# NOTE: This `OptimizeResult` "flavor" is specific to `minimize`
class OptimizeResult(_OptimizeResult):
    success: bool
    status: int
    message: LiteralString
    x: _Float1D
    nit: int
    maxcv: float  # requires `bounds`
    fun: float
    nfev: int
    jac: _Float1D  # requires `jac`
    njev: int  # requires `jac`
    hess: _Float2D  # requires `hess` or `hessp`
    hess_inv: _Float2D | LinearOperator  # requires `hess` or `hessp`, depends on solver
    nhev: int  # requires `hess` or `hessp`

@overload  # identity function with and one parameter, `jac` not truthy
def minimize(
    fun: Callable[Concatenate[_Float1DT, ...], _Float1DT],
    x0: onp.ToFloat,
    args: _Args = (),
    method: MethodMimimize | _MinimizeMethodFun | None = None,
    jac: _Fun1D[onp.ToFloat1D] | _FDMethod | onp.ToFalse | None = None,
    hess: _Fun1D[onp.ToFloat2D] | _FDMethod | HessianUpdateStrategy | None = None,
    hessp: _Fun1Dp[onp.ToFloat1D] | None = None,
    bounds: Bounds | None = None,
    constraints: Constraints = (),
    tol: onp.ToFloat | None = None,
    callback: _CallbackResult | _CallbackVector | None = None,
    options: _MinimizeOptions | None = None,
) -> OptimizeResult: ...
@overload  # `fun` return scalar, `jac` not truthy
def minimize(
    fun: _Fun1D[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    method: MethodMimimize | _MinimizeMethodFun | None = None,
    jac: _Fun1D[onp.ToFloat1D] | _FDMethod | onp.ToFalse | None = None,
    hess: _Fun1D[onp.ToFloat2D] | _FDMethod | HessianUpdateStrategy | None = None,
    hessp: _Fun1Dp[onp.ToFloat1D] | None = None,
    bounds: Bounds | None = None,
    constraints: Constraints = (),
    tol: onp.ToFloat | None = None,
    callback: _CallbackResult | _CallbackVector | None = None,
    options: _MinimizeOptions | None = None,
) -> OptimizeResult: ...
@overload  # fun` return (scalar, vector), `jac` truthy  (positional)
def minimize(
    fun: _Fun1D[tuple[onp.ToFloat, onp.ToFloat1D]],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args,
    method: MethodMimimize | _MinimizeMethodFun | None,
    jac: onp.ToTrue,
    hess: _Fun1D[onp.ToFloat2D] | _FDMethod | HessianUpdateStrategy | None = None,
    hessp: _Fun1Dp[onp.ToFloat1D] | None = None,
    bounds: Bounds | None = None,
    constraints: Constraints = (),
    tol: onp.ToFloat | None = None,
    callback: _CallbackResult | _CallbackVector | None = None,
    options: _MinimizeOptions | None = None,
) -> OptimizeResult: ...
@overload  # fun` return (scalar, vector), `jac` truthy  (keyword)
def minimize(
    fun: _Fun1D[tuple[onp.ToFloat, onp.ToFloat1D]],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    method: MethodMimimize | _MinimizeMethodFun | None = None,
    *,
    jac: onp.ToTrue,
    hess: _Fun1D[onp.ToFloat2D] | _FDMethod | HessianUpdateStrategy | None = None,
    hessp: _Fun1Dp[onp.ToFloat1D] | None = None,
    bounds: Bounds | None = None,
    constraints: Constraints = (),
    tol: onp.ToFloat | None = None,
    callback: _CallbackResult | _CallbackVector | None = None,
    options: _MinimizeOptions | None = None,
) -> OptimizeResult: ...

#
@overload  # method="brent" or method="golden"
def minimize_scalar(
    fun: _Fun0D[onp.ToFloat],
    bracket: _ToBracket | None = None,
    bounds: None = None,
    args: _Args = (),
    method: Literal["brent", "golden"] | None = None,  # default: "brent"
    tol: onp.ToFloat | None = None,
    options: _MinimizeScalarOptionsBracketed | None = None,
) -> _MinimizeScalarResult: ...
@overload  # bound=<given>  (positional)
def minimize_scalar(
    fun: _Fun0D[onp.ToFloat],
    bracket: _Ignored | None,
    bounds: _ToBound,
    args: _Args = (),
    method: Literal["bounded"] | None = None,
    tol: onp.ToFloat | None = None,
    options: _MinimizeScalarOptionsBounded | None = None,
) -> _MinimizeScalarResult: ...
@overload  # bound=<given>  (keyword)
def minimize_scalar(
    fun: _Fun0D[onp.ToFloat],
    bracket: _Ignored | None = None,
    *,
    bounds: _ToBound,
    args: _Args = (),
    method: Literal["bounded"] | None = None,
    tol: onp.ToFloat | None = None,
    options: _MinimizeScalarOptionsBounded | None = None,
) -> _MinimizeScalarResult: ...
@overload  # method=<custom>  (positional)
def minimize_scalar(
    fun: _Fun0D[onp.ToFloat],
    bracket: _ToBracket | None,
    bounds: _ToBound | None,
    args: _Args,
    method: _MinimizeScalarMethodFun[_MinimizeScalarResultT],
    tol: onp.ToFloat | None = None,
    options: Mapping[str, object] | None = None,
) -> _MinimizeScalarResultT: ...
@overload  # method=<custom>  (keyword)
def minimize_scalar(
    fun: _Fun0D[onp.ToFloat],
    bracket: _ToBracket | None = None,
    bounds: _ToBound | None = None,
    args: _Args = (),
    *,
    method: _MinimizeScalarMethodFun[_MinimizeScalarResultT],
    tol: onp.ToFloat | None = None,
    options: Mapping[str, object] | None = None,
) -> _MinimizeScalarResultT: ...

# undocumented
def standardize_bounds(bounds: Bounds, x0: onp.ToFloat1D, meth: MethodMimimize) -> Bounds | list[Bound]: ...
def standardize_constraints(constraints: Constraints, x0: onp.ToFloat1D, meth: MethodMimimize) -> list[Constraint]: ...
