from typing import Literal

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["lsqr"]

###

type _Real = np.bool | npc.integer | npc.floating
type _ToRealMatrix = onp.CanArrayND[_Real] | _spbase[_Real] | LinearOperator[_Real]

type _IStop = Literal[0, 1, 2, 3, 4, 5, 6, 7]

###

def lsqr(
    A: _ToRealMatrix,
    b: onp.ToFloat1D,
    damp: float = 0.0,
    atol: float = 1e-6,
    btol: float = 1e-6,
    conlim: float = 1e8,
    iter_lim: int | None = None,
    show: bool = False,
    calc_var: bool = False,
    x0: onp.ToFloat1D | None = None,
) -> tuple[
    onp.Array1D[np.float64],  # x
    _IStop,  # istop
    int,  # itn
    float,  # r1norm
    float,  # r2norm
    float,  # anorm
    float,  # acond
    float,  # arnorm
    float,  # xnorm
    onp.Array1D[np.float64],  # var
]: ...
