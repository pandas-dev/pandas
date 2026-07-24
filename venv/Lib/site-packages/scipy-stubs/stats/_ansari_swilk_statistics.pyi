# defined in scipy/stats/_ansari_swilk_statistics.pyx

import numpy as np
import optype.numpy as onp

def gscale(test: int, other: int) -> tuple[int, onp.Array1D[np.float32], int]: ...  # undocumented
def swilk(
    x: onp.Array1D[np.float64], a: onp.Array1D[np.float64], init: bool = False, n1: int = -1
) -> tuple[float, float, int]: ...  # undocumented
