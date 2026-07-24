from collections.abc import Callable
from typing import Any, Final, Generic, Never, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from .base import DenseOutput, OdeSolver
from scipy.sparse import csc_matrix, sparray, spmatrix

###

_InexactT = TypeVar("_InexactT", bound=np.float64 | np.complex128, default=np.float64 | Any)

type _LU = tuple[onp.ArrayND[npc.inexact], onp.ArrayND[npc.integer]]
type _FuncLU = Callable[[onp.ArrayND[np.float64]], _LU] | Callable[[onp.ArrayND[np.complex128]], _LU]
type _FuncSolveLU = Callable[[_LU, onp.ArrayND], onp.ArrayND[npc.inexact]]

type _Sparse2D[NumberT: npc.number] = spmatrix[NumberT] | sparray[NumberT, tuple[int, int]]
type _ArrayOrCSC[NumberT: npc.number] = onp.Array2D[NumberT] | csc_matrix[NumberT]

type _ToJacReal = onp.ToFloat2D | _Sparse2D[npc.floating | npc.integer]
type _ToJacComplex = onp.ToComplex2D | _Sparse2D[npc.number]

###

MAX_ORDER: Final = 5
NEWTON_MAXITER: Final = 4
MIN_FACTOR: Final = 0.2
MAX_FACTOR: Final = 10

class BDF(OdeSolver[_InexactT], Generic[_InexactT]):
    max_step: float
    rtol: float
    atol: float
    h_abs: float
    h_abs_old: float | None
    error_norm_old: None
    newton_tol: float

    jac_factor: onp.Array1D[np.float64] | None
    jac: Callable[[float, onp.ArrayND[_InexactT]], _ArrayOrCSC[_InexactT]] | None

    J: _ArrayOrCSC[_InexactT]
    I: _ArrayOrCSC[_InexactT]
    D: onp.Array2D[_InexactT]

    LU: _LU | None
    lu: _FuncLU
    solve_lu: _FuncSolveLU

    gamma: onp.Array1D[np.float64]
    alpha: onp.Array1D[np.float64]
    error_const: onp.Array1D[np.float64]

    order: int | np.intp
    n_equal_steps: int

    @overload
    def __init__(
        self: BDF[np.float64],
        /,
        fun: Callable[[float, onp.ArrayND[np.float64]], onp.ToFloatND],
        t0: float,
        y0: onp.ToFloatND,
        t_bound: float,
        max_step: float = ...,  # = np.inf
        rtol: float = 1e-3,
        atol: float = 1e-6,
        jac: _ToJacReal | Callable[[float, onp.ArrayND[np.float64]], _ToJacReal] | None = None,
        jac_sparsity: _ToJacReal | None = None,
        vectorized: bool = False,
        first_step: float | None = None,
        **extraneous: Never,
    ) -> None: ...
    @overload
    def __init__(
        self: BDF[np.complex128],
        /,
        fun: Callable[[float, onp.ArrayND[np.complex128]], onp.ToComplexND],
        t0: float,
        y0: onp.ToJustComplexND,
        t_bound: float,
        max_step: float = ...,  # = np.inf
        rtol: float = 1e-3,
        atol: float = 1e-6,
        jac: _ToJacComplex | Callable[[float, onp.ArrayND[np.complex128]], _ToJacComplex] | None = None,
        jac_sparsity: _ToJacComplex | None = None,
        vectorized: bool = False,
        first_step: float | None = None,
        **extraneous: Never,
    ) -> None: ...

class BdfDenseOutput(DenseOutput[np.float64]):
    order: int
    t_shift: onp.ArrayND[np.float64]
    denom: onp.ArrayND[np.float64]
    D: onp.ArrayND[np.float64]
    def __init__(self, /, t_old: float, t: float, h: float, order: int, D: onp.ArrayND[np.float64]) -> None: ...

def compute_R(order: int, factor: float) -> onp.ArrayND[np.float64]: ...
def change_D(D: onp.ArrayND[np.float64], order: int, factor: float) -> None: ...
def solve_bdf_system(
    fun: Callable[[float, onp.ArrayND[_InexactT]], onp.ToComplex1D],
    t_new: onp.ToFloat,
    y_predict: onp.ArrayND[_InexactT],
    c: float,
    psi: onp.ArrayND[np.float64],
    LU: _FuncLU,
    solve_lu: _FuncSolveLU,
    scale: onp.ArrayND[np.float64],
    tol: float,
) -> tuple[bool, int, onp.ArrayND[_InexactT], onp.ArrayND[_InexactT]]: ...
