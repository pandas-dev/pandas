from _typeshed import Unused
from collections.abc import Callable
from typing import Concatenate

import numpy as np
import optype.numpy as onp

from ._minimize import OptimizeResult
from ._typing import Bounds, Constraints

def _minimize_cobyqa(
    fun: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat],
    x0: onp.ToFloat1D,
    args: tuple[object, ...] = (),
    bounds: Bounds | None = None,
    constraints: Constraints = (),
    callback: Callable[[onp.Array1D[np.float64]], Unused] | None = None,
    disp: bool = False,
    maxfev: int | None = None,
    maxiter: int | None = None,
    f_target: float = ...,  # = -np.inf
    feasibility_tol: float = 1e-8,
    initial_tr_radius: float = 1.0,
    final_tr_radius: float = 1e-6,
    scale: bool = False,
) -> OptimizeResult: ...
