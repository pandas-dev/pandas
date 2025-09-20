from collections.abc import Callable
from typing import Final, Generic, Never, TypeAlias
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from .base import DenseOutput, OdeSolver
from scipy.sparse import sparray, spmatrix

###

_SCT_co = TypeVar("_SCT_co", covariant=True, bound=npc.inexact, default=np.float64 | np.complex128)

_LU: TypeAlias = tuple[onp.ArrayND[npc.inexact], onp.ArrayND[npc.integer]]
_FuncLU: TypeAlias = Callable[[onp.ArrayND[np.float64]], _LU] | Callable[[onp.ArrayND[np.complex128]], _LU]
_FuncSolveLU: TypeAlias = Callable[[_LU, onp.ArrayND], onp.ArrayND[npc.inexact]]

_ToJac: TypeAlias = onp.ToComplex2D | spmatrix | sparray

###

MAX_ORDER: Final = 5
NEWTON_MAXITER: Final = 4
MIN_FACTOR: Final = 0.2
MAX_FACTOR: Final = 10

class BDF(OdeSolver, Generic[_SCT_co]):
    max_step: float
    h_abs: float
    h_abs_old: float | None
    error_norm_old: None
    newton_tol: float
    jac_factor: onp.ArrayND[np.float64] | None  # 1d

    LU: _LU
    lu: _FuncLU
    solve_lu: _FuncSolveLU

    I: onp.ArrayND[_SCT_co]
    error_const: onp.ArrayND[np.float64]
    gamma: onp.ArrayND[np.float64]
    alpha: onp.ArrayND[np.float64]
    D: onp.ArrayND[np.float64]
    order: int
    n_equal_steps: int

    def __init__(
        self,
        /,
        fun: Callable[[float, onp.Array1D[_SCT_co]], onp.ToComplex1D],
        t0: onp.ToFloat,
        y0: onp.Array1D[_SCT_co] | onp.ToComplexND,
        t_bound: onp.ToFloat,
        max_step: onp.ToFloat = ...,
        rtol: onp.ToFloat = 0.001,
        atol: onp.ToFloat = 1e-06,
        jac: _ToJac | Callable[[float, onp.ArrayND[_SCT_co]], _ToJac] | None = None,
        jac_sparsity: _ToJac | None = None,
        vectorized: bool = False,
        first_step: onp.ToFloat | None = None,
        **extraneous: Never,
    ) -> None: ...

class BdfDenseOutput(DenseOutput):
    order: int
    t_shift: onp.ArrayND[np.float64]
    denom: onp.ArrayND[np.float64]
    D: onp.ArrayND[np.float64]
    def __init__(self, /, t_old: float, t: float, h: float, order: int, D: onp.ArrayND[np.float64]) -> None: ...

def compute_R(order: int, factor: float) -> onp.ArrayND[np.float64]: ...
def change_D(D: onp.ArrayND[np.float64], order: int, factor: float) -> None: ...
def solve_bdf_system(
    fun: Callable[[float, onp.ArrayND[_SCT_co]], onp.ToComplex1D],
    t_new: onp.ToFloat,
    y_predict: onp.ArrayND[_SCT_co],
    c: float,
    psi: onp.ArrayND[np.float64],
    LU: _FuncLU,
    solve_lu: _FuncSolveLU,
    scale: onp.ArrayND[np.float64],
    tol: float,
) -> tuple[bool, int, onp.ArrayND[_SCT_co], onp.ArrayND[_SCT_co]]: ...
