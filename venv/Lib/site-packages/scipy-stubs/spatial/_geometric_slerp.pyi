from typing import overload

import numpy as np
import optype.numpy as onp

__all__ = ["geometric_slerp"]

@overload
def geometric_slerp(start: onp.ToFloat1D, end: onp.ToFloat1D, t: onp.ToFloat, tol: float = 1e-7) -> onp.Array1D[np.float64]: ...
@overload
def geometric_slerp(start: onp.ToFloat1D, end: onp.ToFloat1D, t: onp.ToFloat1D, tol: float = 1e-7) -> onp.Array2D[np.float64]: ...
