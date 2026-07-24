# defined in scipy/spatial/_hausdorff.pyx

import numpy as np
import optype.numpy as onp

__all__ = ["directed_hausdorff"]

def directed_hausdorff(
    ar1: onp.Array2D[np.float64], ar2: onp.Array2D[np.float64], seed: onp.random.ToRNG = 0
) -> tuple[float | complex, np.int64, np.int64]: ...
