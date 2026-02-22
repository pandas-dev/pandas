from collections.abc import Callable
from typing import Literal, TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["tr_interior_point"]

_StateT = TypeVar("_StateT")

###

def tr_interior_point(
    fun: Callable[[onp.Array1D[np.float64]], onp.ToFloat],
    grad: Callable[[onp.Array1D[np.float64]], onp.ToFloat1D],
    lagr_hess: Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat2D],
    n_vars: int,
    n_ineq: int,
    n_eq: int,
    constr: Callable[[onp.Array1D[np.float64]], tuple[onp.ToFloat1D, onp.ToFloat1D]],
    jac: Callable[[onp.Array1D[np.float64]], tuple[onp.ToFloat2D, onp.ToFloat2D]],
    x0: onp.ToFloat1D,
    fun0: onp.ToFloat,
    grad0: onp.ToFloat1D,
    constr_ineq0: onp.ArrayND[npc.floating],
    jac_ineq0: onp.ToArray2D,
    constr_eq0: onp.ArrayND[npc.floating],
    jac_eq0: onp.ToArray2D,
    stop_criteria: Callable[
        [
            _StateT,  # state
            onp.Array1D[np.float64],  # x
            bool,  # last_iteration_failed
            np.float64,  # optimality
            np.float64,  # const_violation
            np.float64,  # trust_radius
            np.float64,  # penalty
            dict[str, int],  # cg_info
        ],
        onp.ToBool,
    ],
    enforce_feasibility: onp.ToArray1D,
    xtol: float,
    state: _StateT,
    initial_barrier_parameter: float,
    initial_tolerance: float,
    initial_penalty: onp.ToFloat,
    initial_trust_radius: onp.ToFloat,
    factorization_method: Literal["NormalEquation", "AugmentedSystem", "QRFactorization", "SVDFactorization"],
    finite_diff_bounds: tuple[onp.ToFloat1D, onp.ToFloat1D],
) -> tuple[onp.Array1D[npc.floating], _StateT]: ...
