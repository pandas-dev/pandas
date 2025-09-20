from typing import Literal

import numpy as np
import optype.numpy as onp

def get_poly_vinit(kind: Literal["poly", "vinit"], dtype: type[np.uint32 | np.uint64]) -> onp.Array2D[np.uint32 | np.uint64]: ...
