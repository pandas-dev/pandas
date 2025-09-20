# see scipy/linalg/_solve_toeplitz.pyx

from typing import TypeVar

import numpy as np
import optype.numpy as onp

_dz = TypeVar("_dz", np.float64, np.complex128)  # static-typing analogue of the `cdef fused dz: ...`

def levinson(
    a: onp.ArrayND[_dz],  # shape: (2n - 1,)
    b: onp.ArrayND[_dz],  # shape: (n,)
) -> tuple[onp.Array1D[_dz], onp.Array1D[_dz]]: ...  # undocumented
