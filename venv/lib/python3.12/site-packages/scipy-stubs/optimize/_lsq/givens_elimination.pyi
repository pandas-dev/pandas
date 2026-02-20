# scipy/optimize/_lsq/givens_elimination.pyx

import numpy as np
import optype.numpy as onp

def givens_elimination(S: onp.Array2D[np.float64], v: onp.Array1D[np.float64], diag: onp.Array1D[np.float64]) -> None: ...
