import numpy as np
import optype.numpy as onp

__all__ = ["lebedev_rule"]

def lebedev_rule(n: int) -> tuple[onp.Array2D[np.float64], onp.Array1D[np.float64]]: ...
