# defined in scipy/stats/_qmvnt_cy.pyx

from typing import Final

import numpy as np
import optype.numpy as onp

###

gammaincinv: Final[np.ufunc] = ...  # implicit re-export from `scipy.special`, required by stubtest

def _qmvn_inner(
    q: onp.Array1D[np.float64],
    rndm: onp.Array2D[np.float64],
    n_qmc_samples: int,
    n_batches: int,
    cho: onp.Array2D[np.float64],
    lo: onp.Array1D[np.float64],
    hi: onp.Array1D[np.float64],
) -> tuple[float, np.float64, int]: ...  # undocumented
def _qmvt_inner(
    q: onp.Array1D[np.float64],
    rndm: onp.Array2D[np.float64],
    n_qmc_samples: int,
    n_batches: int,
    cho: onp.Array2D[np.float64],
    lo: onp.Array1D[np.float64],
    hi: onp.Array1D[np.float64],
    nu: float,
) -> tuple[float, np.float64, int]: ...  # undocumented
