from typing import overload

import numpy as np
import optype.numpy as onp

from ._interface import LinearOperator
from scipy.sparse._base import _spbase

__all__ = ["onenormest"]

###

type _Float1D = onp.Array1D[np.float64]
type _ToMatrix = onp.ToComplex2D | LinearOperator | _spbase

#

@overload  # compute_v: falsy, compute_w: falsy
def onenormest(
    A: _ToMatrix, t: int = 2, itmax: int = 5, compute_v: onp.ToFalse = False, compute_w: onp.ToFalse = False
) -> np.float64: ...
@overload  # compute_v: falsy, compute_w: truthy  (positional)
def onenormest(
    A: _ToMatrix, t: int, itmax: int, compute_v: onp.ToFalse, compute_w: onp.ToTrue
) -> tuple[np.float64, _Float1D]: ...
@overload  # compute_v: falsy, compute_w: truthy  (keyword)
def onenormest(
    A: _ToMatrix, t: int = 2, itmax: int = 5, compute_v: onp.ToFalse = False, *, compute_w: onp.ToTrue
) -> tuple[np.float64, _Float1D]: ...
@overload  # compute_v: truthy  (positional), compute_w: falsy
def onenormest(
    A: _ToMatrix, t: int, itmax: int, compute_v: onp.ToTrue, compute_w: onp.ToFalse = False
) -> tuple[np.float64, _Float1D]: ...
@overload  # compute_v: truthy  (keyword), compute_w: falsy
def onenormest(
    A: _ToMatrix, t: int = 2, itmax: int = 5, *, compute_v: onp.ToTrue, compute_w: onp.ToFalse = False
) -> tuple[np.float64, _Float1D]: ...
@overload  # compute_v: truthy  (positional), compute_w: truthy
def onenormest(
    A: _ToMatrix, t: int, itmax: int, compute_v: onp.ToTrue, compute_w: onp.ToTrue
) -> tuple[np.float64, _Float1D, _Float1D]: ...
@overload  # compute_v: truthy  (keyword), compute_w: truthy
def onenormest(
    A: _ToMatrix, t: int = 2, itmax: int = 5, *, compute_v: onp.ToTrue, compute_w: onp.ToTrue
) -> tuple[np.float64, _Float1D, _Float1D]: ...
