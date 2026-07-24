from typing import Literal, final, type_check_only

import numpy as np
import optype.numpy as onp

from scipy._lib._util import _RichResult

###

@type_check_only
@final
class _WhittakerHendersonResult(_RichResult):
    x: onp.Array1D[np.float64]
    lamb: float

###

def whittaker_henderson(
    signal: onp.ToFloat1D, *, lamb: Literal["reml"] | float = "reml", order: int = 2, weights: onp.ToFloat1D | None = None
) -> _WhittakerHendersonResult: ...
