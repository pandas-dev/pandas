# defined in scipy/optimize/_directmodule.c

from _typeshed import Incomplete
from collections.abc import Callable

import numpy as np
import optype.numpy as onp

from scipy._lib._ccallback import LowLevelCallable

def direct(
    f: Callable[..., Incomplete] | LowLevelCallable,
    lv: onp.Array1D[np.float64],
    ub: onp.Array1D[np.float64],
    f_args: tuple[object, ...],
    disp: bool,
    magic_eps: float,
    max_feval: int,
    max_iter: int,
    algorithm: str,
    fglobal: onp.Array1D[np.float64],
    fglobal_reltol: float,
    volume_reltol: float,
    sigma_reltol: float,
    callback: Callable[..., object],
    /,
) -> tuple[list[Incomplete], float, int, int, int]: ...  # -> (x_seq, minf, ret_code, numfunc, numiter)
