from collections.abc import Callable, Mapping
from typing import Literal

import numpy as np
import optype.numpy as onp

from scipy.optimize import OptimizeResult as _OptimizeResult
from scipy.optimize._typing import TRSolver
from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

###

type _Vector_b1 = onp.Array1D[np.bool]
type _Vector_i0 = onp.Array1D[np.int_]

type _Scalar_f8 = float | np.float64
type _Vector_f8 = onp.Array1D[np.float64]
type _Matrix_f8 = onp.Array2D[np.float64] | _spbase | LinearOperator

type _FunResid = Callable[[_Vector_f8], _Vector_f8]
type _FunJac = Callable[[_Vector_f8, _Vector_f8], _Matrix_f8]
type _FunLoss = Callable[[_Vector_f8], _Scalar_f8]

###

class OptimizeResult(_OptimizeResult):
    x: _Vector_f8
    cost: _Scalar_f8
    fun: _Scalar_f8
    jac: _Vector_f8
    grad: _Vector_f8
    optimality: _Scalar_f8
    active_mask: _Vector_b1
    nfev: int
    njev: int
    status: Literal[0, 1, 2, 3, 4]

def lsmr_operator(Jop: LinearOperator, d: _Vector_f8, active_set: _Vector_b1) -> LinearOperator: ...  # undocumented

#
def find_intersection(
    x: _Vector_f8, tr_bounds: _Vector_f8, lb: _Vector_f8, ub: _Vector_f8
) -> tuple[_Vector_f8, _Vector_f8, _Vector_b1, _Vector_b1, _Vector_b1, _Vector_b1]: ...

#
def dogleg_step(
    x: _Vector_f8,
    newton_step: _Vector_f8,
    g: _Vector_f8,
    a: _Scalar_f8,
    b: _Scalar_f8,
    tr_bounds: _Vector_f8,
    lb: _Vector_f8,
    ub: _Vector_f8,
) -> tuple[_Vector_f8, _Vector_i0, np.bool]: ...

#
def dogbox(
    fun: _FunResid,
    jac: _FunJac,
    x0: _Vector_f8,
    f0: _Vector_f8,
    J0: _Matrix_f8,
    lb: _Vector_f8,
    ub: _Vector_f8,
    ftol: _Scalar_f8,
    xtol: _Scalar_f8,
    gtol: _Scalar_f8,
    max_nfev: int,
    x_scale: Literal["jac"] | _Scalar_f8 | _Vector_f8,
    loss_function: _FunLoss,
    tr_solver: TRSolver | None,
    tr_options: Mapping[str, object],
    verbose: bool,
    callback: Callable[[OptimizeResult], object] | None = None,
) -> OptimizeResult: ...  # undocumented
