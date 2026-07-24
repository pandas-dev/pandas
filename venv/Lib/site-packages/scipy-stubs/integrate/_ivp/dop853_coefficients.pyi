from typing import Final

import numpy as np
import optype.numpy as onp

N_STAGES: Final = 12
N_STAGES_EXTENDED: Final = 16
INTERPOLATOR_POWER: Final = 7
C: Final[onp.Array1D[np.float64]]
A: Final[onp.Array2D[np.float64]]
B: Final[onp.Array1D[np.float64]]
E3: Final[onp.Array1D[np.float64]]
E5: Final[onp.Array1D[np.float64]]
D: Final[onp.Array2D[np.float64]]
