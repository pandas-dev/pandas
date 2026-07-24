from collections.abc import Sequence
from typing import Final

import numpy as np
import optype.numpy as onp

from ._ndbspline import NdBSpline
from scipy.sparse import csr_array

type _Float1D = onp.Array1D[np.float64]
type _Float2D = onp.Array2D[np.float64]

type _ToBBox = Sequence[onp.ToFloat | None]

###

TOL: Final = 0.001

class PackedMatrix:
    a: _Float2D
    offset: onp.Array1D[np.int64]
    nc: int

    @property
    def shape(self, /) -> tuple[int, int]: ...

    #
    def __init__(self, /, a: _Float2D, offset: onp.Array1D[np.int64], nc: int) -> None: ...
    def todense(self, /) -> _Float2D: ...
    def tocsr(self, /, k: int, m: int, len_t: int) -> csr_array[np.float64, tuple[int, int]]: ...

class F:
    Ax: PackedMatrix
    Dx: PackedMatrix
    Ay: PackedMatrix
    Dy: PackedMatrix
    Q: _Float2D
    kx: int
    tx: _Float1D
    x_x: _Float1D
    ky: int
    ty: _Float1D
    x_y: _Float1D
    z: _Float2D

    C: _Float2D  # set in __call__
    fp: float  # set in __call__

    def __init__(
        self,
        /,
        Ax: PackedMatrix,
        Dx: PackedMatrix,
        Ay: PackedMatrix,
        Dy: PackedMatrix,
        Q: _Float2D,
        kx: int,
        tx: _Float1D,
        x_x: _Float1D,
        ky: int,
        ty: _Float1D,
        x_y: _Float1D,
        z: _Float2D,
    ) -> None: ...
    def __call__(self, /, p: float) -> float: ...

#
def _ndbspline_call_like_bivariate[ScalarT: np.float64 | np.complex128](
    ndbs: NdBSpline[ScalarT], x: onp.ToFloat1D, y: onp.ToFloat1D, dx: int = 0, dy: int = 0, grid: bool = True
) -> onp.ArrayND[ScalarT]: ...

#
def return_NdBSpline(fp: float, tck: tuple[_Float1D, _Float1D, _Float2D], degrees: tuple[int, int]) -> NdBSpline[np.float64]: ...

#
def _stack_augmented_fitpack(
    A: PackedMatrix, D: PackedMatrix, nc: int, k: int, p: float
) -> tuple[_Float2D, onp.Array1D[np.int64], int]: ...

#
def _solve_2d_fitpack(
    Ax: PackedMatrix,
    Ay: PackedMatrix,
    Q: _Float2D,
    p: float,
    kx: int,
    tx: _Float1D,
    x_x: _Float1D,
    ky: int,
    ty: _Float1D,
    x_y: _Float1D,
    z: _Float2D,
    Dx: PackedMatrix | None = None,
    Dy: PackedMatrix | None = None,
) -> tuple[_Float2D, float, _Float2D]: ...

#
def _p_search_hit_s(
    Ax: PackedMatrix,
    Dx: PackedMatrix,
    Ay: PackedMatrix,
    Dy: PackedMatrix,
    Q: _Float2D,
    kx: int,
    tx: _Float1D,
    x_x: _Float1D,
    ky: int,
    ty: _Float1D,
    x_y: _Float1D,
    z: _Float2D,
    s: float,
    fp0: float | None,
    *,
    p_init: float = 1.0,
    tol_rel: float = 1e-3,
    maxit: int = 40,
) -> tuple[float, _Float2D, float]: ...

#
def _apply_bbox_grid(
    x: _Float1D, y: _Float1D, Z: _Float2D, bbox: _ToBBox
) -> tuple[_Float1D, _Float1D, _Float2D, slice | onp.Array1D[np.intp], slice | onp.Array1D[np.intp]]: ...

#
def _build_design_matrices(
    x: _Float1D, y: _Float1D, z: _Float2D, tx: _Float1D, ty: _Float1D, kx: int, ky: int
) -> tuple[PackedMatrix, PackedMatrix, _Float2D]: ...

#
def _initialise_knots(m: int, xb: float, xe: float, k: int, nest: int | None = None) -> tuple[_Float1D, int, int, int]: ...

#
def _add_knots(
    x: _Float1D,
    k: int,
    s: float,
    t: _Float1D,
    nmin: int,
    nmax: int,
    nest: int,
    fp: float,
    fpold: float | None,
    residuals: _Float1D,
    nplus: int | None,
) -> tuple[_Float1D, int]: ...

#
def _regrid_fitpack(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    Z: onp.ToFloat2D,
    *,
    kx: int = 3,
    ky: int = 3,
    s: float = 0.0,
    maxit: int = 50,
    nestx: int | None = None,
    nesty: int | None = None,
    bbox: _ToBBox | None = None,
) -> NdBSpline[np.float64]: ...

#
def _regrid(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    z: onp.ToFloat2D,
    *,
    bbox: _ToBBox | None = None,
    kx: int = 3,
    ky: int = 3,
    s: float = 0.0,
    maxit: int = 50,
) -> NdBSpline[np.float64]: ...
