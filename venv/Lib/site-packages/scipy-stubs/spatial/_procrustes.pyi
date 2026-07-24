import numpy as np
import optype.numpy as onp

__all__ = ["procrustes"]

def procrustes(
    data1: onp.ToFloat2D, data2: onp.ToFloat2D
) -> tuple[onp.Array2D[np.float64], onp.Array2D[np.float64], np.float64]: ...
