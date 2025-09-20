from typing import Literal, TypeAlias

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

_Mode: TypeAlias = Literal["mirror", "constant", "nearest", "wrap", "interp"]

def savgol_coeffs(
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    pos: int | None = None,
    use: Literal["conv", "dot"] = "conv",
) -> onp.Array1D[npc.floating]: ...

#
def savgol_filter(
    x: onp.ToFloatND,
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: op.CanIndex = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.ArrayND[np.float32 | np.float64]: ...
