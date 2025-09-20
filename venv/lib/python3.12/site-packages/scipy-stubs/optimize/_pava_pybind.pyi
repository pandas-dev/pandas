# defined in scipy/optimize/_pava/pava_pybind.cpp

import numpy as np
import optype.numpy as onp

def pava(
    xa: onp.Array1D[np.float64], wa: onp.Array1D[np.float64], ra: onp.Array1D[np.intp], /
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.float64], onp.Array1D[np.intp], np.intp]: ...  # -> (x, w, r, b)
