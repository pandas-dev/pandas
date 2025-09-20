# https://github.com/scipy/scipy/blob/v1.14.1/scipy/optimize/_trlib/_trlib.pyx

from collections.abc import Callable, Mapping
from typing import Final, Never

import numpy as np
import optype.numpy as onp

from scipy.optimize._trustregion import BaseQuadraticSubproblem

###

__test__: Final[Mapping[Never, Never]]  # undocumented

class TRLIBQuadraticSubproblem(BaseQuadraticSubproblem):  # undocumented
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        fun: Callable[[onp.Array1D[np.float64]], onp.ToFloat],
        jac: Callable[[onp.Array1D[np.float64]], onp.ToFloat1D],
        hess: Callable[[onp.Array1D[np.float64]], onp.ToFloat2D] | None,
        hessp: Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat1D] | None,
        tol_rel_i: onp.ToFloat = -2.0,
        tol_rel_b: onp.ToFloat = -3.0,
        disp: onp.ToBool = False,
    ) -> None: ...
