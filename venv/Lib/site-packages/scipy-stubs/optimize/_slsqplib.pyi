# scipy/optimize/__slsqp.h

from typing import Any

import numpy as np
import optype.numpy as onp

class error(Exception): ...  # undocumented

def nnls(A: onp.Array2D[np.float64], b: onp.Array1D[np.float64], maxiter: int, /) -> onp.Array1D[np.float64]: ...  # undocumented
def slsqp(
    state_dict: dict[str, Any],
    funx: float,
    gradx: onp.Array1D[np.float64],
    C: onp.Array2D[np.float64],
    d: onp.Array1D[np.float64],
    sol: onp.Array1D[np.float64],
    mult: onp.Array1D[np.float64],
    xl: onp.Array1D[np.float64],
    xu: onp.Array1D[np.float64],
    buffer: onp.Array1D[np.float64],
    indices: onp.Array1D[np.int32],
) -> None: ...  # undocumented
