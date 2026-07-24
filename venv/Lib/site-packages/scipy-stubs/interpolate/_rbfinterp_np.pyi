from _typeshed import Unused

import numpy as np
import optype.numpy as onp

from scipy.linalg.blas import _FortranFunction  # type-check-only

###

type _Int2D = onp.Array2D[np.int64]
type _Float1D = onp.Array1D[np.float64]
type _Float2D = onp.Array2D[np.float64]

###

dsysv: _FortranFunction = ...

def _build_evaluation_coefficients(
    x: _Float2D, y: _Float2D, kernel: str, epsilon: float, powers: _Int2D, shift: _Float1D, scale: _Float1D, xp: Unused
) -> _Float2D: ...

#
def polynomial_matrix(x: _Float2D, powers: _Int2D, xp: Unused) -> _Float2D: ...

#
def _monomial_powers(ndim: int, degree: int, xp: Unused) -> _Int2D: ...

#
def _build_system(
    y: _Float2D, d: _Float2D, smoothing: _Float1D, kernel: str, epsilon: float, powers: _Int2D, xp: Unused
) -> tuple[_Float2D, _Float2D, _Float1D, _Float1D]: ...

#
def _build_and_solve_system(
    y: _Float2D, d: _Float2D, smoothing: _Float1D, kernel: str, epsilon: float, powers: _Int2D, xp: Unused
) -> tuple[_Float2D, _Float1D, _Float1D]: ...

#
def compute_interpolation(
    x: _Float2D,
    y: _Float2D,
    kernel: str,
    epsilon: float,
    powers: _Int2D,
    shift: _Float1D,
    scale: _Float1D,
    coeffs: _Float2D,
    xp: Unused,
) -> _Float2D: ...
