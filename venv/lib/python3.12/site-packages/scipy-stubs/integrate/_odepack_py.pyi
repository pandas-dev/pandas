import threading
from collections.abc import Callable
from typing import Final, Literal, TypeAlias, TypeVar, TypeVarTuple, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp

__all__ = ["ODEintWarning", "odeint"]

_Ts = TypeVarTuple("_Ts")
_R = TypeVar("_R")

_FuncYT: TypeAlias = Callable[[onp.Array1D[np.float64], float, *_Ts], _R]
_FuncTY: TypeAlias = Callable[[float, onp.Array1D[np.float64], *_Ts], _R]

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
    message: str

###

ODE_LOCK: Final[threading.Lock] = ...

class ODEintWarning(Warning): ...

@overload  # args=() (default), full_output=False (default), tfirst=False (default)
def odeint(
    func: _FuncYT[*tuple[()], onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[()] = (),
    Dfun: _FuncYT[*tuple[()], onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    full_output: Literal[False, 0] = 0,
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    tfirst: Literal[False, 0] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # args=() (default), full_output=False (default), *, tfirst=True
def odeint(
    func: _FuncTY[*tuple[()], onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[()] = (),
    Dfun: _FuncTY[*tuple[()], onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    full_output: Literal[False, 0] = 0,
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    *,
    tfirst: Literal[True, 1],
) -> onp.Array2D[np.float64]: ...
@overload  # args=() (default), *, full_output=True, tfirst=False (default)
def odeint(
    func: _FuncYT[*tuple[()], onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[()] = (),
    Dfun: _FuncYT[*tuple[()], onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    *,
    full_output: Literal[True, 1],
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    tfirst: Literal[False, 0] = False,
) -> tuple[onp.Array2D[np.float64], _InfoDict]: ...
@overload  # args=() (default), full_output=True, *, tfirst=True
def odeint(
    func: _FuncTY[*tuple[()], onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[()] = (),
    Dfun: _FuncTY[*tuple[()], onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    *,
    full_output: Literal[True, 1],
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    tfirst: Literal[True, 1],
) -> tuple[onp.Array2D[np.float64], _InfoDict]: ...
@overload  # args=<given>, full_output=False (default), tfirst=False (default)
def odeint(
    func: _FuncYT[*_Ts, onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[*_Ts],
    Dfun: _FuncYT[*_Ts, onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    full_output: Literal[False, 0] = 0,
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    tfirst: Literal[False, 0] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # args=<given>, full_output=False (default), *, tfirst=True
def odeint(
    func: _FuncTY[*_Ts, onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[*_Ts],
    Dfun: _FuncTY[*_Ts, onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    full_output: Literal[False, 0] = 0,
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    *,
    tfirst: Literal[True, 1],
) -> onp.Array2D[np.float64]: ...
@overload  # args=<given>, *, full_output=True, tfirst=False (default)
def odeint(
    func: _FuncYT[*_Ts, onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[*_Ts],
    Dfun: _FuncYT[*_Ts, onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    *,
    full_output: Literal[True, 1],
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    tfirst: Literal[False, 0] = False,
) -> tuple[onp.Array2D[np.float64], _InfoDict]: ...
@overload  # args=<given>, *, full_output=True, tfirst=True
def odeint(
    func: _FuncTY[*_Ts, onp.ToFloat1D | float],
    y0: onp.ToFloat1D | float,
    t: onp.ToFloat1D,
    args: tuple[*_Ts],
    Dfun: _FuncTY[*_Ts, onp.ToFloat2D] | None = None,
    col_deriv: onp.ToBool = 0,
    *,
    full_output: Literal[True, 1],
    ml: int | None = None,
    mu: int | None = None,
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
    printmessg: onp.ToBool = 0,
    tfirst: Literal[True, 1],
) -> tuple[onp.Array2D[np.float64], _InfoDict]: ...
