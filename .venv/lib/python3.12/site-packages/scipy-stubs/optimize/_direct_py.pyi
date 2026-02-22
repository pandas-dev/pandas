from collections.abc import Callable
from typing import Concatenate

import numpy as np
import optype.numpy as onp

from ._constraints import Bounds
from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["direct"]

class OptimizeResult(_OptimizeResult):
    x: onp.Array1D[np.float64]
    fun: float | np.float64
    status: int
    success: bool
    message: str
    nfev: int
    nit: int

###

ERROR_MESSAGES: tuple[str, str, str, str, str, str, str, str, str, str, str] = ...
SUCCESS_MESSAGES: tuple[str, str, str] = ...

def direct(
    func: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat],
    bounds: tuple[onp.ToFloat1D, onp.ToFloat1D] | Bounds,
    *,
    args: tuple[object, ...] = (),
    eps: onp.ToFloat = 1e-4,
    maxfun: onp.ToJustInt | None = None,
    maxiter: onp.ToJustInt = 1000,
    locally_biased: onp.ToBool = True,
    f_min: onp.ToFloat = ...,
    f_min_rtol: onp.ToFloat = 0.0001,
    vol_tol: onp.ToFloat = 1e-16,
    len_tol: onp.ToFloat = 1e-06,
    callback: Callable[[onp.ToArray1D], None] | None = None,
) -> OptimizeResult: ...
