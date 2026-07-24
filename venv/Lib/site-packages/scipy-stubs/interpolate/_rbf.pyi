from collections.abc import Callable
from typing import Literal

import numpy as np
import optype.numpy as onp

__all__ = ["Rbf"]

###

type _Mode = Literal["1-D", "N-D"]
type _Function = (
    Literal["multiquadric", "inverse", "gaussian", "linear", "cubic", "quintic", "thin_plate"]
    | Callable[[Rbf, float], onp.ToFloat]
)

###

# legacy
class Rbf:
    N: int
    di: onp.Array1D[np.float64]
    xi: onp.Array2D[np.float64]
    function: _Function
    epsilon: float
    smooth: float
    norm: str | Callable[..., onp.ToFloat2D]
    mode: _Mode
    nodes: onp.Array1D[np.float64]

    @property
    def A(self, /) -> onp.Array2D[np.float64]: ...  # undocumented

    #
    def __init__(
        self,
        /,
        *args: onp.ToFloat1D,
        function: _Function = ...,
        epsilon: onp.ToFloat = ...,
        smooth: onp.ToFloat = ...,
        norm: str | Callable[..., onp.ToFloat2D] = ...,
        mode: _Mode = ...,
    ) -> None: ...
    def __call__(self, /, *args: onp.ToFloat | onp.ToFloatND) -> onp.ArrayND[np.float64]: ...
