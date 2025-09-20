from collections.abc import Callable, Iterable
from typing import Concatenate, Final, Literal, Protocol, TypeAlias, TypedDict, type_check_only

import numpy as np
import optype.numpy as onp

from .canonical_constraint import CanonicalConstraint
from scipy.optimize._constraints import Bounds, LinearConstraint, NonlinearConstraint, PreparedConstraint
from scipy.optimize._differentiable_functions import _DoesMap
from scipy.optimize._optimize import OptimizeResult as _OptimizeResult
from scipy.sparse import sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

_StopCond: TypeAlias = Literal[0, 1, 2, 3, 4]

_Float: TypeAlias = float | np.float64
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]
_Sparse: TypeAlias = spmatrix | sparray
_Matrix: TypeAlias = _Float2D | _Sparse

_HessPFunc: TypeAlias = Callable[Concatenate[_Float1D, _Float1D, ...], _Float1D]
_ObjectiveHessFunc: TypeAlias = Callable[[_Float1D], _Matrix]
_ConstraintsHessFunc: TypeAlias = Callable[[_Float1D, _Float1D, _Float1D], _Matrix]

@type_check_only
class _CGInfo(TypedDict):
    niter: int
    stop_cond: _StopCond
    hits_boundary: bool

@type_check_only
class _Objective(Protocol):
    f: _Float
    g: _Float1D
    nfev: int
    njev: int
    nhev: int

###

TERMINATION_MESSAGES: dict[Literal[0, 1, 2, 3], str]  # undocumented

class OptimizeResult(_OptimizeResult):
    x: _Float1D
    optimality: _Float
    const_violation: _Float
    fun: _Float
    grad: _Float1D
    lagrangian_grad: _Float1D
    nit: int
    nfev: int
    njev: int
    nhev: int
    cg_niter: int
    method: Literal["equality_constrained_sqp", "tr_interior_point"]
    constr: list[_Float]
    jac: list[_Matrix]
    v: list[_FloatND]
    constr_nfev: list[int]
    constr_njev: list[int]
    constr_nhev: list[int]
    tr_radius: _Float
    constr_penalty: _Float
    barrier_tolerance: _Float
    barrier_parameter: _Float
    execution_time: _Float
    message: str
    status: Literal[0, 1, 2, 3]
    cg_stop_cond: _StopCond

# undocumented
class HessianLinearOperator:
    n: Final[int]
    hessp: Final[_HessPFunc]

    def __init__(self, /, hessp: _HessPFunc, n: int) -> None: ...
    def __call__(self, /, x: _Float1D, *args: object) -> LinearOperator: ...

# undocumented
class LagrangianHessian:
    n: Final[int]
    objective_hess: Final[_ObjectiveHessFunc]
    constraints_hess: Final[_ConstraintsHessFunc]

    def __init__(self, /, n: int, objective_hess: _ObjectiveHessFunc, constraints_hess: _ConstraintsHessFunc) -> None: ...
    def __call__(self, /, x: _Float1D, v_eq: _Float1D, v_ineq: _Float1D | None = None) -> LinearOperator: ...

# undocumented
def update_state_sqp(
    state: OptimizeResult,
    x: _Float1D,
    last_iteration_failed: bool,
    objective: _Objective,
    prepared_constraints: Iterable[CanonicalConstraint | PreparedConstraint],
    start_time: float,
    tr_radius: _Float,
    constr_penalty: _Float,
    cg_info: _CGInfo,
) -> OptimizeResult: ...

# undocumented
def update_state_ip(
    state: OptimizeResult,
    x: _Float1D,
    last_iteration_failed: bool,
    objective: _Objective,
    prepared_constraints: Iterable[CanonicalConstraint | PreparedConstraint],
    start_time: float,
    tr_radius: _Float,
    constr_penalty: _Float,
    cg_info: _CGInfo,
    barrier_parameter: _Float,
    barrier_tolerance: _Float,
) -> OptimizeResult: ...

#
def _minimize_trustregion_constr(
    fun: Callable[Concatenate[_Float1D, ...], _Float],
    x0: onp.ToFloat1D,
    args: tuple[object, ...],
    grad: Callable[Concatenate[_Float1D, ...], onp.ToFloat1D] | None,
    hess: Callable[Concatenate[_Float1D, ...], onp.ToFloat2D | _Sparse | LinearOperator] | None,
    hessp: _HessPFunc | None,
    bounds: Bounds | None,
    constraints: LinearConstraint | NonlinearConstraint | None,
    xtol: float = 1e-8,
    gtol: float = 1e-8,
    barrier_tol: float = 1e-8,
    sparse_jacobian: bool | None = None,
    callback: Callable[[OptimizeResult], None] | None = None,
    maxiter: int = 1000,
    verbose: Literal[0, 1, 2] = 0,
    finite_diff_rel_step: onp.ToFloat1D | None = None,
    initial_constr_penalty: float = 1.0,
    initial_tr_radius: float = 1.0,
    initial_barrier_parameter: float = 0.1,
    initial_barrier_tolerance: float = 0.1,
    factorization_method: str | None = None,
    disp: bool = False,
    workers: int | _DoesMap | None = None,
) -> OptimizeResult: ...
