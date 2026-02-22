from collections.abc import Callable
from typing import Never

import numpy as np
import optype.numpy as onp

from .base import DenseOutput, OdeSolver

class LSODA(OdeSolver[np.float64]):
    def __init__(
        self,
        /,
        fun: Callable[[float, onp.Array1D[np.float64]], onp.Array1D[np.float64]],
        t0: float,
        y0: onp.Array1D[np.float64],
        t_bound: float,
        first_step: float | None = None,
        min_step: float = 0.0,
        max_step: float = ...,
        rtol: onp.ToFloat | onp.ToFloat1D = 0.001,
        atol: onp.ToFloat | onp.ToFloat1D = 1e-06,
        jac: Callable[[float, onp.Array1D[np.float64]], onp.Array2D[np.float64]] | None = None,
        lband: int | None = None,
        uband: int | None = None,
        vectorized: bool = False,
        **extraneous: Never,
    ) -> None: ...

class LsodaDenseOutput(DenseOutput[np.float64]):
    h: float
    yh: onp.Array1D[np.float64]
    p: onp.Array1D[np.intp]

    def __init__(self, /, t_old: float, t: float, h: float, order: int, yh: onp.Array1D[np.float64]) -> None: ...
