# defined in scipy/integrate/_odepackmodule.c

from collections.abc import Callable
from typing import Literal, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp

class error(Exception): ...

@type_check_only
class _InfoDict(TypedDict):
    hu: onp.Array1D[np.float64]
    tcur: onp.Array1D[np.float64]
    tolsf: onp.Array1D[np.float64]
    tsw: onp.Array1D[np.float64]
    nst: onp.Array1D[np.int32]
    nfe: onp.Array1D[np.int32]
    nje: onp.Array1D[np.int32]
    nqu: onp.Array1D[np.int32]
    imxer: int
    lenrw: int
    leniw: int
    mused: onp.Array1D[np.int32]

###

# NOTE: This (private) `odeint` differs from the one in `_odepack_py` because it
# doesn't have a `printmessg` kwarg, and has a different return type.
@overload
def odeint(
    fun: Callable[..., onp.ToFloat1D | float],
    y0: onp.ToFloatND | float,
    t: onp.ToFloat1D,
    args: tuple[object, ...] = (),
    Dfun: Callable[..., onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    ml: int = -1,
    mu: int = -1,
    full_output: Literal[0, False] = 0,
    rtol: float | None = None,
    atol: float | None = None,
    tcrit: onp.ToFloat1D | None = None,
    h0: float = 0.0,
    hmax: float = 0.0,
    hmin: float = 0.0,
    ixpr: onp.ToBool = 0,
    mxstep: int = 0,
    mxhnil: int = 0,
    mxordn: int = 12,
    mxords: int = 5,
    tfirst: int = 0,
) -> tuple[onp.Array2D[np.float64], int]: ...
@overload
def odeint(
    fun: Callable[..., onp.ToFloat1D | float],
    y0: onp.ToFloatND | float,
    t: onp.ToFloat1D,
    args: tuple[object, ...],
    Dfun: Callable[..., onp.ToFloat2D] | None,
    col_deriv: onp.ToBool,
    ml: int,
    mu: int,
    full_output: Literal[True, 1],
    rtol: float | None = None,
    atol: float | None = None,
    tcrit: onp.ToFloat1D | None = None,
    h0: float = 0.0,
    hmax: float = 0.0,
    hmin: float = 0.0,
    ixpr: onp.ToBool = 0,
    mxstep: int = 0,
    mxhnil: int = 0,
    mxordn: int = 12,
    mxords: int = 5,
    tfirst: int = 0,
) -> tuple[onp.Array2D[np.float64], _InfoDict, int]: ...
