# defined in scipy/interpolate/_rgi_cython.pyx

from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

_Inexact64T = TypeVar("_Inexact64T", np.float64, np.complex128)

###

def evaluate_linear_2d(
    values: onp.Array2D[_Inexact64T],
    indices: onp.Array2D[np.intp],
    norm_distances: onp.Array2D[np.float64],
    grid: tuple[onp.Array1D[np.float64], ...],
    out: onp.Array2D[_Inexact64T],
) -> onp.Array2D[_Inexact64T]: ...
def find_indices(
    grid: tuple[onp.Array1D[np.float64], ...], xi: onp.Array2D[np.float64]
) -> tuple[onp.Array2D[np.intp], onp.Array2D[np.float64]]: ...
