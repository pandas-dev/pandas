from collections.abc import Callable
from typing import Final, Never, TypeAlias

import numpy as np
import optype.numpy as onp

from .base import DenseOutput, OdeSolver
from scipy.sparse import sparray, spmatrix

_LU: TypeAlias = tuple[onp.ArrayND[np.float64], onp.ArrayND[np.intp]]
_FuncSolveLU: TypeAlias = Callable[[_LU, onp.ArrayND[np.float64]], onp.ArrayND[np.float64]]
_ToJac: TypeAlias = onp.ToArray2D | spmatrix | sparray

###

S6: Final[float] = ...

C: Final[onp.ArrayND[np.float64]] = ...
E: Final[onp.ArrayND[np.float64]] = ...

MU_REAL: Final[float] = ...
MU_COMPLEX: Final[complex] = ...

T: Final[onp.ArrayND[np.float64]] = ...
TI: Final[onp.ArrayND[np.float64]] = ...
TI_REAL: Final[onp.ArrayND[np.float64]] = ...
TI_COMPLEX: Final[onp.ArrayND[np.complex128]] = ...

P: Final[onp.ArrayND[np.float64]] = ...

NEWTON_MAXITER: Final = 6
MIN_FACTOR: Final = 0.2
MAX_FACTOR: Final = 10

class Radau(OdeSolver[np.float64]):
    max_step: float
    h_abs: float
    h_abs_old: float | None
    error_norm_old: float | None
    newton_tol: float
    jac_factor: onp.ArrayND[np.float64] | None  # 1d
    current_jac: bool

    LU_real: _LU
    LU_complex: _LU
    lu: Callable[[onp.ArrayND[np.float64]], _LU]
    solve_lu: _FuncSolveLU

    y_old: onp.ArrayND[np.float64] | None
    f: onp.ArrayND[np.float64]
    I: onp.ArrayND[np.float64]
    Z: onp.ArrayND[np.float64] | None
    sol: RadauDenseOutput | None

    def __init__(
        self,
        /,
        fun: Callable[[float, onp.Array1D[np.float64]], onp.ToArray1D],
        t0: onp.ToFloat,
        y0: onp.ToArray1D,
        t_bound: onp.ToFloat,
        max_step: onp.ToFloat = ...,
        rtol: onp.ToFloat = 0.001,
        atol: onp.ToFloat = 1e-06,
        jac: _ToJac | Callable[[float, onp.Array1D[np.float64]], _ToJac] | None = None,
        jac_sparsity: _ToJac | None = None,
        vectorized: bool = False,
        first_step: onp.ToFloat | None = None,
        **extraneous: Never,
    ) -> None: ...

class RadauDenseOutput(DenseOutput[np.float64]):
    order: int
    h: float
    Q: onp.ArrayND[np.float64]
    y_old: onp.ArrayND[np.float64]

    def __init__(self, /, t_old: float, t: float, y_old: onp.ArrayND[np.float64], Q: onp.ArrayND[np.float64]) -> None: ...

def solve_collocation_system(
    fun: Callable[[float, onp.Array1D[np.float64]], onp.ToFloat1D],
    t: onp.ToFloat,
    y: onp.ArrayND[np.float64],
    h: onp.ToFloat,
    Z0: onp.ArrayND[np.float64],
    scale: onp.ArrayND[np.float64],
    tol: onp.ToFloat,
    LU_real: _LU,
    LU_complex: _LU,
    solve_lu: _FuncSolveLU,
) -> tuple[bool, int, onp.Array2D[np.float64], float | None]: ...
def predict_factor(
    h_abs: onp.ToFloat, h_abs_old: onp.ToFloat, error_norm: onp.ToFloat, error_norm_old: onp.ToFloat
) -> onp.ToFloat: ...
