from collections.abc import Callable
from typing import Protocol, type_check_only

import numpy as np
import optype.numpy as onp

from ._trlib import TRLIBQuadraticSubproblem

__all__ = ["TRLIBQuadraticSubproblem", "get_trlib_quadratic_subproblem"]

@type_check_only
class _SubproblemFactory(Protocol):
    def __call__(
        self,
        /,
        x: onp.ToFloat1D,
        fun: Callable[[onp.Array1D[np.float64]], onp.ToFloat],
        jac: Callable[[onp.Array1D[np.float64]], onp.ToFloat1D],
        hess: Callable[[onp.Array1D[np.float64]], onp.ToFloat2D] | None = None,
        hessp: Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat1D] | None = None,
    ) -> TRLIBQuadraticSubproblem: ...

###

def get_trlib_quadratic_subproblem(
    tol_rel_i: onp.ToFloat = -2.0, tol_rel_b: onp.ToFloat = -3.0, disp: onp.ToBool = False
) -> _SubproblemFactory: ...
