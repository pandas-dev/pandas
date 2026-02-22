from typing import Literal, TypeAlias

import numpy as np
import optype as op
import optype.numpy as onp

from scipy.optimize import OptimizeResult as _OptimizeResult

_Float: TypeAlias = float | np.float64
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Int1D: TypeAlias = onp.Array1D[np.int_]

###

class OptimizeResult(_OptimizeResult):
    x: _Float1D
    fun: _Float
    cost: _Float
    initial_cost: _Float
    optimality: _Float
    active_mask: _Int1D
    nit: int
    status: int

def regularized_lsq_with_qr(
    m: onp.ToJustInt,
    n: onp.ToJustInt,
    R: onp.ToFloat2D,
    QTb: onp.ToFloat1D,
    perm: onp.ToFloat1D,
    diag: onp.ToFloat1D,
    copy_R: op.CanBool = True,
) -> _Float1D: ...

# undocumented
def backtracking(
    A: onp.ToFloat2D,
    g: onp.ToFloat1D,
    x: onp.ToFloat1D,
    p: onp.ToFloat | onp.ToFloat1D,
    theta: onp.ToFloat | onp.ToFloat1D,
    p_dot_g: onp.ToFloat,
    lb: onp.ToFloat1D,
    ub: onp.ToFloat1D,
) -> tuple[_Float1D, _Float1D, _Float]: ...

# undocumented
def select_step(
    x: onp.ToFloat1D,
    A_h: onp.ToFloat1D,
    g_h: onp.ToFloat1D,
    c_h: onp.ToFloat1D,
    p: onp.ToFloat | onp.ToFloat1D,
    p_h: onp.ToFloat1D,
    d: onp.ToFloat | onp.ToFloat1D,
    lb: onp.ToFloat1D,
    ub: onp.ToFloat1D,
    theta: onp.ToFloat | onp.ToFloat1D,
) -> _Float | _Float1D: ...

# undocumented
def trf_linear(
    A: onp.ToFloat2D,
    b: onp.ToFloat2D,
    x_lsq: onp.ToFloat1D,
    lb: onp.ToFloat1D,
    ub: onp.ToFloat1D,
    tol: onp.ToFloat,
    lsq_solver: Literal["exact", "lsmr"],
    lsmr_tol: onp.ToFloat | None,
    max_iter: onp.ToJustInt,
    verbose: Literal[0, 1, 2],
    *,
    lsmr_maxiter: onp.ToJustInt | None = None,
) -> OptimizeResult: ...
