from collections.abc import Callable, Sequence
from typing import Concatenate, Literal, TypeAlias

import numpy as np
import optype.numpy as onp

__all__ = ["fmin_cobyla"]

_Ignored: TypeAlias = object

###

def fmin_cobyla(
    func: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat],
    x0: onp.ToArray1D,
    cons: Sequence[Callable[[onp.Array1D[np.float64]], onp.ToFloat | onp.ToFloat1D]],
    args: tuple[object, ...] = (),
    consargs: tuple[object, ...] | None = None,
    rhobeg: onp.ToFloat = 1.0,
    rhoend: onp.ToFloat = 0.0001,
    maxfun: onp.ToInt = 1000,
    disp: Literal[0, 1, 2, 3] | None = None,
    catol: onp.ToFloat = 0.0002,
    *,
    callback: Callable[[onp.Array1D[np.float64]], _Ignored] | None = None,
) -> onp.Array1D[np.float64]: ...
