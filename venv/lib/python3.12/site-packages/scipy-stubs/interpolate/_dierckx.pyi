# defined in scipy/interpolate/src/_dierckxmodule.cc

import numpy as np
import optype.numpy as onp

def _coloc(
    x: onp.Array1D[np.float64], t: onp.Array1D[np.float64], k: int, ab: onp.Array2D[np.float64], offset: int = 0, /
) -> None: ...  # undocumented
def _coloc_nd(
    xvals: onp.Array2D[np.float64],
    _t: onp.Array2D[np.float64],
    len_t: onp.Array1D[np.int64],
    k: onp.Array1D[np.int64],
    _indices_k1d: onp.Array2D[np.int64],
    _cstrides: onp.Array1D[np.int64],
    /,
) -> None: ...  # undocumented
def _norm_eq_lsq(
    x: onp.Array1D[np.float64],
    t: onp.Array1D[np.float64],
    k: int,
    y: onp.Array2D[np.float64],
    w: onp.Array1D[np.float64],
    ab: onp.Array2D[np.float64],
    rhs: onp.Array2D[np.float64],
    /,
) -> None: ...  # undocumented

#
def data_matrix(
    x: onp.Array1D[np.float64], t: onp.Array1D[np.float64], k: int, w: onp.Array1D[np.float64], extrapolate: bool = False, /
) -> onp.Array2D[np.float64]: ...  # undocumented
def evaluate_all_bspl(
    t: onp.Array1D[np.float64], k: int, xval: float, m: int, mu: int = 0, /
) -> onp.Array1D[np.float64]: ...  # undocumented
def evaluate_ndbspline(
    xi: onp.Array2D[np.float64],
    t: onp.Array2D[np.float64],
    len_t: onp.Array1D[np.int64],
    k: onp.Array1D[np.int64],
    nu: onp.Array1D[np.int64],
    extrapolate: bool,
    c1r: onp.Array1D[np.float64],
    num_c_tr: int,
    strides_c1r: onp.Array1D[np.int64],
    indices_k1d: onp.Array2D[np.int64],
    /,
) -> onp.Array2D[np.float64]: ...  # undocumented
def evaluate_spline(
    t: onp.Array1D[np.float64], c: onp.Array2D[np.float64], k: int, xp: onp.Array1D[np.float64], nu: int, extrapolate: bool
) -> onp.Array2D[np.float64]: ...  # undocumented
def find_interval(t: onp.Array1D[np.float64], k: int, xval: float, prev_l: int, extrapolate: bool) -> int: ...  # undocumented
def fpback(R: onp.Array2D[np.float64], nc: int, y: onp.Array2D[np.float64]) -> onp.Array2D[np.float64]: ...  # undocumented
def fpknot(
    x: onp.Array1D[np.float64], t: onp.Array1D[np.float64], k: int, residuals: onp.Array2D[np.float64]
) -> float: ...  # undocumented
def qr_reduce(
    a: onp.Array2D[np.float64], offset: int, nc: int, y: onp.Array2D[np.float64], startrow: int = 1
) -> None: ...  # undocumented
