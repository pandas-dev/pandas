from collections.abc import Callable
from typing import overload
from typing_extensions import TypeVarTuple

import numpy as np
import optype.numpy as onp

_Ts = TypeVarTuple("_Ts")

###

def _central_diff_weights(Np: int, ndiv: int = 1) -> onp.Array1D[np.float64]: ...  # undocumented

#
@overload
def _derivative(
    func: Callable[[float], onp.ToFloat], x0: float, dx: float = 1.0, n: int = 1, args: tuple[()] = (), order: int = 3
) -> np.float64: ...
@overload
def _derivative(
    func: Callable[[float, *_Ts], onp.ToFloat], x0: float, dx: float = 1.0, n: int = 1, *, args: tuple[*_Ts], order: int = 3
) -> np.float64: ...  # undocumented
