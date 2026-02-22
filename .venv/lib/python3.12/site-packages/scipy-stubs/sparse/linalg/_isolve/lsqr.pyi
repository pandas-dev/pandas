from typing import Literal, TypeAlias

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["lsqr"]

_Real: TypeAlias = np.bool_ | npc.integer | npc.floating
_ToRealMatrix: TypeAlias = onp.CanArrayND[_Real] | _spbase[_Real] | LinearOperator[_Real]

_IStop: TypeAlias = Literal[0, 1, 2, 3, 4, 5, 6, 7]

###

def lsqr(
    A: _ToRealMatrix,
    b: onp.ToFloat1D,
    damp: float | npc.floating = 0.0,
    atol: float | npc.floating = 1e-6,
    btol: float | npc.floating = 1e-6,
    conlim: onp.ToFloat = 1e8,
    iter_lim: int | None = None,
    show: onp.ToBool = False,
    calc_var: onp.ToBool = False,
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
