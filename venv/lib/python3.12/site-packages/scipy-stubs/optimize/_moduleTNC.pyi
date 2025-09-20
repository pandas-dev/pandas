# defined in scipy/optimize/tnc/_moduleTNC.pyx

from collections.abc import Callable

import numpy as np
import optype.numpy as onp

def tnc_minimize(
    func_and_grad: Callable[..., tuple[float, onp.Array1D[np.float64]]],
    x0: onp.Array1D[np.float64],
    low: onp.Array1D[np.float64],
    up: onp.Array1D[np.float64],
    scale: onp.Array1D[np.float64],
    offset: onp.Array1D[np.float64],
    messages: int,
    maxCGit: int,
    maxfun: int,
    eta: float,
    stepmx: float,
    accuracy: float,
    fmin: float,
    ftol: float,
    xtol: float,
    pgtol: float,
    rescale: float,
    callback: Callable[..., object] | None = None,
) -> tuple[int, int, int, onp.Array1D[np.float64], float, onp.Array1D[np.float64]]: ...  # -> (rc, nfeval, niter, x, f, g)
