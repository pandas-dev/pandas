from collections.abc import Callable
from typing import Literal, TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import dia_matrix

__all__ = ["equality_constrained_sqp"]

_StateT = TypeVar("_StateT")

def default_scaling(x: onp.ToArray1D) -> dia_matrix: ...
def equality_constrained_sqp(
    fun_and_constr: Callable[[onp.Array1D[np.float64]], tuple[float | npc.floating, onp.ToFloat1D]],
    grad_and_jac: Callable[[onp.Array1D[np.float64]], tuple[onp.ToFloat1D, onp.ToFloat2D]],
    lagr_hess: Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat2D],
    x0: onp.ToFloat1D,
    fun0: onp.ToFloat,
    grad0: onp.ToFloat1D,
    constr0: onp.ToFloat1D,
    jac0: onp.ToFloat2D,
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
    state: _StateT,
    initial_penalty: onp.ToFloat,
    initial_trust_radius: onp.ToFloat,
    factorization_method: Literal["NormalEquation", "AugmentedSystem", "QRFactorization", "SVDFactorization"],
    trust_lb: onp.Array1D[npc.floating] | None = None,
    trust_ub: onp.Array1D[npc.floating] | None = None,
    scaling: Callable[[onp.Array1D[np.float64]], dia_matrix] = ...,
) -> tuple[onp.Array1D[npc.floating], _StateT]: ...
