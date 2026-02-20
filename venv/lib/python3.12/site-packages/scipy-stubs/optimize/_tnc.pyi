from collections.abc import Callable, Sequence
from typing import Final, Literal, TypeAlias

import numpy as np
import optype.numpy as onp

__all__ = ["fmin_tnc"]

_ReturnCode: TypeAlias = Literal[-1, 0, 1, 2, 3, 4, 5, 6, 7]

MSG_NONE: Final = 0
MSG_ITER: Final = 1
MSG_INFO: Final = 2
MSG_VERS: Final = 4
MSG_EXIT: Final = 8
MSG_ALL: Final = 15
MSGS: Final[dict[Literal[0, 1, 2, 4, 8, 15], str]]

INFEASIBLE: Final = -1
LOCALMINIMUM: Final = 0
FCONVERGED: Final = 1
XCONVERGED: Final = 2
MAXFUN: Final = 3
LSFAIL: Final = 4
CONSTANT: Final = 5
NOPROGRESS: Final = 6
USERABORT: Final = 7
RCSTRINGS: Final[dict[_ReturnCode, str]]

def fmin_tnc(
    func: Callable[..., onp.ToFloat] | Callable[..., tuple[onp.ToFloat, onp.ToFloat]],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: Callable[..., onp.ToFloat] | None = None,
    args: tuple[object, ...] = (),
    approx_grad: int = 0,
    bounds: Sequence[tuple[float | None, float | None]] | None = None,
    epsilon: float = 1e-08,
    scale: onp.ToFloat | onp.ToFloat1D | None = None,
    offset: onp.ToFloat | onp.ToFloat1D | None = None,
    messages: int = 15,
    maxCGit: int = -1,
    maxfun: int | None = None,
    eta: float = -1,
    stepmx: float = 0,
    accuracy: float = 0,
    fmin: float = 0,
    ftol: float = -1,
    xtol: float = -1,
    pgtol: float = -1,
    rescale: float = -1,
    disp: bool | None = None,
    callback: Callable[[onp.Array1D[np.float64]], None] | None = None,
) -> tuple[onp.Array1D[np.float64], int, _ReturnCode]: ...
