import numpy as np
import optype.numpy as onp

def rvs_rcont1(
    row: onp.Array1D[np.int64], col: onp.Array1D[np.int64], ntot: int | np.int64, size: int, random_state: onp.random.ToRNG
) -> onp.Array3D[np.int64]: ...
def rvs_rcont2(
    row: onp.Array1D[np.int64], col: onp.Array1D[np.int64], ntot: int | np.int64, size: int, random_state: onp.random.ToRNG
) -> onp.Array3D[np.int64]: ...
