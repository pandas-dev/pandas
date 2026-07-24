# defined in scipy/interpolate/_rgi_cython.pyx

import numpy as np
import optype.numpy as onp

###

def evaluate_linear_2d[Inexact64T: (np.float64, np.complex128)](
    values: onp.Array2D[Inexact64T],
    indices: onp.Array2D[np.intp],
    norm_distances: onp.Array2D[np.float64],
    grid: tuple[onp.Array1D[np.float64], ...],
    out: onp.Array2D[Inexact64T],
) -> onp.Array2D[Inexact64T]: ...
def find_indices(
    grid: tuple[onp.Array1D[np.float64], ...], xi: onp.Array2D[np.float64]
) -> tuple[onp.Array2D[np.intp], onp.Array2D[np.float64]]: ...
