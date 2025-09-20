# defined in scipy/interpolate/src/_fitpackmodule.c

from _typeshed import Incomplete

import numpy as np
import optype.numpy as onp

def _insert(
    iopt: int, t: onp.Array1D[np.float64], c: onp.Array2D[np.float64], k: int, x: onp.Array1D[np.float64], m: int, /
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.float64], int]: ...
def _parcur(
    x: onp.Array1D[np.float64],
    w: onp.Array1D[np.float64],
    u: onp.Array1D[np.float64],
    ub: onp.Array1D[np.float64],
    ue: onp.Array1D[np.float64],
    k: int,
    iopt: int,
    ipar: int,
    s: onp.Array1D[np.float64],
    t: onp.Array1D[np.float64],
    nest: int,
    wrk: onp.Array1D[np.float64],
    iwrk: int,
    per: float,
    /,
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.float64], dict[str, Incomplete]]: ...
def _surfit(
    x: onp.Array1D[np.float64],
    y: onp.Array1D[np.float64],
    z: onp.Array1D[np.float64],
    w: onp.Array1D[np.float64],
    xb: float,
    xe: float,
    yb: float,
    ye: float,
    kx: int,
    ky: int,
    iopt: int,
    s: onp.Array1D[np.float64],
    eps: float,
    tx: onp.Array1D[np.float64],
    ty: onp.Array1D[np.float64],
    nxest: int,
    nyest: int,
    wrk: onp.Array1D[np.float64],
    lwrk1: int,
    lwrk2: int,
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.float64], onp.Array1D[np.float64], dict[str, Incomplete]]: ...
