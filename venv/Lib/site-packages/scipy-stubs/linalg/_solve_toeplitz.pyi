# see scipy/linalg/_solve_toeplitz.pyx

import numpy as np
import optype.numpy as onp

def levinson[
    dz: (np.float64, np.complex128)  # static-typing analogue of the `cdef fused dz: ...`
](
    a: onp.ArrayND[dz],  # shape: (2n - 1,)
    b: onp.ArrayND[dz],  # shape: (n,)
) -> tuple[onp.Array1D[dz], onp.Array1D[dz]]: ...  # undocumented
