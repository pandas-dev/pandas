from typing import Final, Literal, TypeAlias, TypeVar

import numpy as np
import optype as op
import optype.numpy as onp

from scipy.sparse import sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

_T = TypeVar("_T")

_Int1D: TypeAlias = onp.Array1D[np.int_]

_Float: TypeAlias = float | np.float64
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]

_Sparse: TypeAlias = sparray | spmatrix
_ToMatrix: TypeAlias = onp.ToFloat2D | _Sparse | LinearOperator
_Matrix: TypeAlias = _Float2D | _Sparse | LinearOperator

###

EPS: Final[float] = ...

def intersect_trust_region(x: onp.ToFloat1D, s: onp.ToFloatND, Delta: onp.ToFloat) -> tuple[_Float, _Float]: ...
def solve_trust_region_2d(B: onp.ToFloat2D, g: onp.ToFloat1D, Delta: onp.ToFloat) -> tuple[_Float1D, bool]: ...
def solve_lsq_trust_region(
    n: onp.ToInt,
    m: onp.ToInt,
    uf: onp.ToFloat1D,
    s: onp.ToFloatND,
    V: onp.ToFloat2D,
    Delta: onp.ToFloat,
    initial_alpha: onp.ToFloat | None = None,
    rtol: onp.ToFloat = 0.01,
    max_iter: onp.ToInt = 10,
) -> tuple[_Float1D, _Float, int]: ...
def update_tr_radius(
    Delta: onp.ToFloat,
    actual_reduction: onp.ToFloat,
    predicted_reduction: onp.ToFloat,
    step_norm: onp.ToFloat,
    bound_hit: op.CanBool,
) -> tuple[_Float, _Float]: ...

#
def build_quadratic_1d(
    J: _ToMatrix, g: onp.ToFloat1D, s: onp.ToFloat1D, diag: onp.ToFloat1D | None = None, s0: onp.ToFloat1D | None = None
) -> tuple[_Float, _Float, _Float]: ...
def minimize_quadratic_1d(
    a: onp.ToFloat, b: onp.ToFloat, lb: onp.ToFloat1D, ub: onp.ToFloat1D, c: onp.ToFloat = 0
) -> tuple[_Float, _Float]: ...
def evaluate_quadratic(
    J: _ToMatrix, g: onp.ToFloat1D, s: onp.ToFloat1D | onp.ToFloat2D, diag: onp.ToFloat1D | None = None
) -> _Float | _Float1D: ...

#
def in_bounds(x: onp.ToFloatND, lb: onp.ToFloatND, ub: onp.ToFloatND) -> np.bool_: ...
def step_size_to_bound(x: onp.ToFloat1D, s: onp.ToFloat1D, lb: onp.ToFloat1D, ub: onp.ToFloat1D) -> tuple[_Float, _Int1D]: ...
def find_active_constraints(x: onp.ToFloat1D, lb: onp.ToFloat1D, ub: onp.ToFloat1D, rtol: onp.ToFloat = 1e-10) -> _Int1D: ...
def make_strictly_feasible(x: onp.ToFloat1D, lb: onp.ToFloat1D, ub: onp.ToFloat1D, rstep: onp.ToFloat = 1e-10) -> _Float1D: ...
def CL_scaling_vector(x: onp.ToFloat1D, g: onp.ToFloat1D, lb: onp.ToFloat1D, ub: onp.ToFloat1D) -> tuple[_Float1D, _Float1D]: ...
def reflective_transformation(y: onp.ToFloat1D, lb: onp.ToFloat1D, ub: onp.ToFloat1D) -> tuple[_Float1D, _Float1D]: ...

#
def print_header_nonlinear() -> None: ...
def print_iteration_nonlinear(
    iteration: int, nfev: int, cost: float, cost_reduction: float, step_norm: float, optimality: float
) -> None: ...
def print_header_linear() -> None: ...
def print_iteration_linear(iteration: int, cost: float, cost_reduction: float, step_norm: float, optimality: float) -> None: ...

#
def compute_grad(J: _ToMatrix, f: onp.ToFloat1D) -> _Float1D | _Sparse: ...
def compute_jac_scale(J: _ToMatrix, scale_inv_old: onp.ToFloat | onp.ToFloat2D | None = None) -> tuple[_Float2D, _Float2D]: ...

#
def left_multiplied_operator(J: _ToMatrix, d: onp.ToFloat1D) -> LinearOperator: ...
def right_multiplied_operator(J: _ToMatrix, d: onp.ToFloat1D) -> LinearOperator: ...
def regularized_lsq_operator(J: _ToMatrix, diag: onp.ToFloat1D) -> LinearOperator: ...

#
def right_multiply(J: _ToMatrix, d: onp.ToFloat1D, copy: bool = True) -> _Matrix: ...
def left_multiply(J: _ToMatrix, d: onp.ToFloat1D, copy: bool = True) -> _Matrix: ...

#
def check_termination(
    dF: onp.ToFloat,
    F: onp.ToFloat,
    dx_norm: onp.ToFloat,
    x_norm: onp.ToFloat,
    ratio: onp.ToFloat,
    ftol: onp.ToFloat,
    xtol: onp.ToFloat,
) -> Literal[2, 3, 4] | None: ...

#
def scale_for_robust_loss_function(J: _ToMatrix, f: _T, rho: onp.ToFloat1D) -> tuple[_ToMatrix, _T]: ...
