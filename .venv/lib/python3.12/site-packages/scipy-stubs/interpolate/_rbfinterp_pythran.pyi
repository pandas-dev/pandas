from typing import TypeAlias

import numpy as np
import optype.numpy as onp

_Int2D: TypeAlias = onp.Array2D[np.int64]
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]

###

__pythran__: tuple[str, str]

def _pythran_polynomial_matrix(x: _Float2D, powers: _Int2D, /) -> _Float2D: ...
def _pythran_build_system(
    y: _Float2D, d: _Float2D, smoothing: _Float1D, kernel: str, epsilon: float, powers: _Int2D, /
) -> tuple[_Float2D, _Float2D, _Float1D, _Float1D]: ...
def _build_evaluation_coefficients(
    x: _Float2D, y: _Float2D, kernel: str, epsilon: float, powers: _Int2D, shift: _Float1D, scale: _Float1D
) -> _Float2D: ...
