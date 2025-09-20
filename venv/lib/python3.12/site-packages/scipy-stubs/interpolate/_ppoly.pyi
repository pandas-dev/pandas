# defined in scipy/interpolate/_ppoly.pyx

import numpy as np
import optype.numpy as onp

def _croots_poly1(c: onp.Array1D[np.float64], w: onp.Array1D[np.complex128], y: float = 0.0) -> None: ...  # undocumented
def evaluate(
    c: onp.Array3D[np.float64 | np.complex128],
    x: onp.Array1D[np.float64],
    xp: onp.Array1D[np.float64],
    dx: int,
    extrapolate: bool,
    out: onp.Array2D[np.float64 | np.complex128],
) -> None: ...  # undocumented
def evaluate_bernstein(
    c: onp.Array3D[np.float64 | np.complex128],
    x: onp.Array1D[np.float64],
    xp: onp.Array1D[np.float64],
    nu: int,
    extrapolate: bool,
    out: onp.Array2D[np.float64 | np.complex128],
) -> None: ...  # undocumented
def evaluate_nd(
    c: onp.Array3D[np.float64 | np.complex128],
    xs: tuple[onp.Array1D[np.float64], ...],
    ks: onp.Array1D[np.int32 | np.int64],
    xp: onp.Array2D[np.float64],
    dx: onp.Array2D[np.int32 | np.int64],
    extrapolate: bool,
    out: onp.Array2D[np.float64 | np.complex128],
) -> None: ...  # undocumented
def fix_continuity(c: onp.Array3D[np.float64 | np.complex128], x: onp.Array1D[np.float64], order: int) -> None: ...
def integrate(
    c: onp.Array3D[np.float64 | np.complex128],
    x: onp.Array1D[np.float64],
    a: float,
    b: float,
    extrapolate: bool,
    out: onp.Array1D[np.float64 | np.complex128],
) -> None: ...  # undocumented
def real_roots(
    c: onp.Array3D[np.float64], x: onp.Array1D[np.float64], y: float, report_discont: bool, extrapolate: bool
) -> list[onp.Array1D[np.float64]]: ...  # undocumented
