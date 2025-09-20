from _typeshed import Incomplete

import numpy as np
import optype.numpy as onp

def setulb(
    m: int,
    x: onp.Array1D[np.float64],
    l: onp.Array1D[np.float64],
    u: onp.Array1D[np.float64],
    nbd: onp.Array1D[np.int32 | np.int64],
    f: onp.Array0D[np.int32 | np.int64],
    g: onp.Array1D[np.int32 | np.int64],
    factr: float,
    pgtol: float,
    wa: onp.Array1D[np.float64],
    iwa: onp.Array1D[np.int32 | np.int64],
    task: onp.Array1D[np.int32 | np.int64],
    lsave: onp.Array1D[np.int32 | np.int64],
    isave: onp.Array1D[np.int32 | np.int64],
    dsave: onp.Array1D[np.float64],
    maxls: int | None,
    ln_task: onp.Array1D[np.int32 | np.int64],
) -> tuple[Incomplete, ...]: ...
