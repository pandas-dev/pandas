from typing import Literal

import numpy as np
import optype.numpy as onp

def ellip_harm(
    h2: onp.ToFloat,
    k2: onp.ToFloat,
    n: onp.ToInt,
    p: onp.ToFloat,
    s: onp.ToFloat,
    signm: Literal[-1, 1] = 1,
    signn: Literal[-1, 1] = 1,
) -> np.float64: ...
def ellip_harm_2(h2: onp.ToFloat, k2: onp.ToFloat, n: onp.ToInt, p: onp.ToInt, s: onp.ToFloat) -> onp.Array0D[np.float64]: ...
def ellip_normal(h2: onp.ToFloat, k2: onp.ToFloat, n: onp.ToFloat, p: onp.ToFloat) -> onp.Array0D[np.float64]: ...
