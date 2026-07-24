# defined in scipy/interpolate/src/_fitpackmodule.c

from typing import Any, Literal as L, final

import numpy as np
import optype.numpy as onp

type _Float1D = onp.Array1D[np.float64]

###

@final
class error(Exception): ...

#
def bispeu(
    tx: onp.ToFloat1D,
    ty: onp.ToFloat1D,
    c: onp.ToFloat1D,  # len(c) == (len(x) - kx - 1) * (len(y) - ky - 1)
    kx: int,  # 0 <= kx < len(x)
    ky: int,  # 0 <= ky < len(y)
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    /,
) -> tuple[_Float1D, int]: ...

#
def bispev(
    tx: onp.ToFloat1D,
    ty: onp.ToFloat1D,
    c: onp.ToFloat1D,  # len(c) == (len(x) - kx - 1) * (len(y) - ky - 1)
    kx: int,  # 0 <= kx < len(x)
    ky: int,  # 0 <= ky < len(y)
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    /,
) -> tuple[_Float1D, int]: ...

#
def splev(
    t: onp.ToFloat1D,
    c: onp.ToFloat1D,  # len(c) >= len(t) - k - 1
    k: int,  # 0 <= k < len(c)
    x: onp.ToFloat1D,
    e: L[0, 1, 2, 3] = 0,
    /,
) -> tuple[_Float1D, int]: ...

#
def splder(
    t: onp.ToFloat1D,
    c: onp.ToFloat1D,  # len(c) >= len(t) - k - 1
    k: int,  # 0 <= k < len(c)
    x: onp.ToFloat1D,
    nu: int = 1,  # 0 <= nu <= k < len(c)
    e: L[0, 1, 2, 3] = 0,
    /,
) -> tuple[_Float1D, int]: ...

#
def splint(
    t: onp.ToFloat1D,
    c: onp.ToFloat1D,  # len(c) >= len(t) - k - 1
    k: int,  # 0 <= k < len(c)
    a: onp.ToFloat | onp.ToFloatND,
    b: onp.ToFloat | onp.ToFloatND,
    /,
) -> tuple[float, _Float1D]: ...

#
def sproot(
    t: onp.ToFloat1D,  # len(t) >= 8
    c: onp.ToFloat1D,  # len(c) >= len(t) - 4
    mest: int = -1,  # mest = 3 * (len(t) - 7)
    /,
) -> tuple[_Float1D, int, int]: ...

#
def spalde(
    t: onp.ToFloat1D,  # len(t) >= 8
    c: onp.ToFloat1D,  # len(c) >= len(t) - 4
    k1: int,
    x: onp.ToFloat | onp.ToFloatND,
    /,
) -> tuple[_Float1D, int]: ...

#
def fpchec(x: onp.ToFloat1D, t: onp.ToFloat1D, k: int, /) -> int: ...

#
def parder(
    tx: onp.ToFloat1D,
    ty: onp.ToFloat1D,
    c: onp.ToFloat1D,
    kx: int,
    ky: int,
    nux: int,
    nuy: int,
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    /,
) -> tuple[onp.Array2D[np.float64], int]: ...

#
def pardtc(
    tx: onp.ToFloat1D,
    ty: onp.ToFloat1D,
    c: onp.ToFloat1D,
    kx: L[1, 2, 3, 4, 5],
    ky: L[1, 2, 3, 4, 5],
    nux: L[0, 1, 2, 3, 4],  # 0 <= nux < kx
    nuy: L[0, 1, 2, 3, 4],  # 0 <= nux < ky
    /,
) -> tuple[_Float1D, int]: ...

#
def pardeu(
    tx: onp.ToFloat1D,
    ty: onp.ToFloat1D,
    c: onp.ToFloat1D,
    kx: int,
    ky: int,
    nux: int,
    nuy: int,
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    /,
) -> tuple[_Float1D, int]: ...

#
def dblint(
    tx: onp.ToFloat1D,
    ty: onp.ToFloat1D,
    c: onp.ToFloat1D,
    # requires c.shape[0] == (nx - kx - 1) * (ny - ky - 1)
    kx: int,
    ky: int,
    xb: float,
    xe: float,
    yb: float,
    ye: float,
    /,
) -> float: ...

#
def curfit(
    iopt: int,
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    w: onp.ToFloat1D,
    xb: float,
    xe: float,
    k: int,
    s: float,
    t: onp.ToFloat1D,
    wrk: onp.ToFloat1D,
    iwrk: onp.ToInt1D,
    /,
) -> tuple[int, _Float1D, _Float1D, float, int]: ...

#
def percur(
    iopt: int,
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    w: onp.ToFloat1D,
    k: int,
    s: float,
    # nest: int,
    t: onp.ToFloat1D,
    wrk: onp.ToFloat1D,
    iwrk: onp.ToInt1D,
    /,
) -> tuple[int, _Float1D, _Float1D, float, int]: ...

#
def surfit_lsq(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    z: onp.ToFloatND,
    w: onp.ToFloatND,
    xb: float,
    xe: float,
    yb: float,
    ye: float,
    kx: int,
    ky: int,
    s: float,
    nxest: int,
    nyest: int,
    nmax: int,
    eps: float,
    tx: onp.Array1D[np.float64],
    ty: onp.Array1D[np.float64],
    /,
) -> tuple[_Float1D, float, int]: ...

#
def surfit_smth(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    z: onp.ToFloatND,
    w: onp.ToFloatND,
    xb: float,
    xe: float,
    yb: float,
    ye: float,
    kx: int,
    ky: int,
    s: float,
    nxest: int,
    nyest: int,
    nmax: int,
    eps: float,
    /,
) -> tuple[int, _Float1D, int, _Float1D, _Float1D, float, int]: ...

#
def surfit(
    iopt: int,
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    z: onp.ToFloatND,
    w: onp.ToFloatND,
    xb: float,
    xe: float,
    yb: float,
    ye: float,
    kx: int,
    ky: int,
    s: float,
    nxest: int,
    nyest: int,
    nmax: int,
    eps: float,
    wrk: onp.ToFloat1D,
    /,
) -> tuple[int, _Float1D, int, _Float1D, _Float1D, float, int, _Float1D]: ...

#
def regrid(
    iopt: int,
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    z: onp.ToFloat1D,
    xb: float,
    xe: float,
    yb: float,
    ye: float,
    kx: int,
    ky: int,
    s: float,
    nxest: int,
    nyest: int,
    maxit: int,
    wrk: onp.ToFloat1D,
    iwrk: onp.ToInt1D,
    /,
) -> tuple[int, _Float1D, int, _Float1D, _Float1D, float, int]: ...

#
def sphere(
    iopt: int,
    teta: onp.ToFloat1D,
    phi: onp.ToFloat1D,
    r: onp.ToFloat1D,
    w: onp.ToFloat1D,
    s: float,
    ntest: int,
    npest: int,
    tt: _Float1D | None,
    tp: _Float1D | None,
    eps: float,
    wrk1: onp.ToFloat1D,
    wrk2: onp.ToFloat1D,
    iwrk: onp.ToInt1D,
    /,
) -> tuple[int, _Float1D, int, _Float1D, _Float1D, float, int]: ...

#
def spgrid(
    iopt: onp.ToInt1D,  # (3,)
    ider: onp.ToInt1D,  # (4,)
    u: onp.ToFloat1D,
    v: onp.ToFloat1D,
    r: onp.ToFloat1D,
    r0: float,
    r1: float,
    s: float,
    nuest: int,
    nvest: int,
    /,
) -> tuple[int, _Float1D, int, _Float1D, _Float1D, float, int]: ...

#
def parcur(
    x: onp.ToFloat1D,
    w: onp.ToFloat1D,
    u: onp.ToFloat1D,
    ub: float,
    ue: float,
    k: int,
    task: int,
    ipar: int,
    s: float,
    t: onp.ToFloat1D,
    nest: int,
    wrk: onp.ToFloat1D,
    iwrk: onp.ToInt1D,
    per: int,
    /,
) -> tuple[_Float1D, _Float1D, dict[str, Any]]: ...

#
def insert(iopt: int, t: onp.ToFloat1D, c: onp.ToFloat1D, k: int, x: float, nest: int, /) -> tuple[_Float1D, _Float1D, int]: ...
