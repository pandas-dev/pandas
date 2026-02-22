from typing import overload
from typing_extensions import deprecated

import numpy as np
import optype.numpy as onp

__all__ = ["nnls"]

@overload
def nnls(
    A: onp.ToFloat2D, b: onp.ToFloat1D, *, maxiter: onp.ToInt | None = None
) -> tuple[onp.ArrayND[np.float64], np.float64]: ...
@overload
@deprecated("The atol parameter is deprecated and will be removed in SciPy 1.18.0. It is not used in the implementation.")
def nnls(
    A: onp.ToFloat2D, b: onp.ToFloat1D, *, maxiter: onp.ToInt | None = None, atol: onp.ToFloat
) -> tuple[onp.ArrayND[np.float64], np.float64]: ...
