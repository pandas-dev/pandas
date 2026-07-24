import numpy as np
import optype.numpy as onp

__all__ = ["nnls"]

def nnls(A: onp.ToFloat2D, b: onp.ToFloat1D, *, maxiter: int | None = None) -> tuple[onp.ArrayND[np.float64], float]: ...
