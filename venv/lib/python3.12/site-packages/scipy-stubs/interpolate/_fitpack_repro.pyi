from collections.abc import Callable, Generator
from typing import Final, Literal, TypeAlias, type_check_only

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._bsplines import BSpline

_Float: TypeAlias = float | np.float32 | np.float64 | np.longdouble
_Float64: TypeAlias = float | np.float64

@type_check_only
class _RootRatiBunch(Bunch):
    root: _Float
    converged: bool
    iterations: int
    ier: Literal[0, 2, 3]

###

# private (yet exported) api

TOL: Final = 0.001
MAXIT: Final = 20

class F:  # undocumented
    x: onp.Array1D[np.float64]
    y: onp.Array1D[np.float64]
    t: onp.Array1D[np.float64]
    k: int
    w: onp.Array1D[np.float64] | None
    s: _Float

    YY: onp.Array2D[np.float64]
    AA: onp.Array2D[np.float64]
    offset: onp.Array1D[np.int64]

    nc: int
    b: onp.Array2D[np.float64]

    spl: BSpline  # set in __call__

    def __init__(
        self,
        /,
        x: onp.Array1D[np.float64],
        y: onp.Array1D[np.float64],
        t: onp.Array1D[np.float64],
        k: int,
        s: _Float64,
        w: onp.Array1D[np.float64] | None = None,
        *,
        R: onp.Array2D[np.float64] | None = None,
        Y: onp.Array2D[np.float64] | None = None,
    ) -> None: ...
    def __call__(self, /, p: onp.ToFloat) -> _Float: ...

class Bunch: ...  # undocumented

#
def add_knot(
    x: onp.Array1D[np.float64], t: onp.Array1D[np.float64], k: int, residuals: onp.Array1D[np.float64]
) -> onp.Array1D[np.float64]: ...  # undocumented

#
def prodd(t: onp.Array1D[np.float64], i: int, j: int, k: int) -> _Float64: ...  # undocumented

#
def disc(t: onp.ArrayND[npc.floating], k: int) -> tuple[onp.Array2D[np.float64], onp.Array1D[np.int64], int]: ...  # undocumented

#
def fprati(
    p1: onp.ToFloat, f1: onp.ToFloat, p2: onp.ToFloat, f2: onp.ToFloat, p3: onp.ToFloat, f3: onp.ToFloat
) -> _Float: ...  # undocumented

#
def root_rati(
    f: Callable[[float], onp.ToFloat],
    p0: onp.ToFloat,
    bracket: tuple[tuple[onp.ToFloat, onp.ToFloat], tuple[onp.ToFloat, onp.ToFloat]],
    acc: onp.ToFloat,
) -> _RootRatiBunch: ...  # undocumented

# public api

def generate_knots(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    *,
    w: onp.ToFloat1D | None = None,
    xb: onp.ToFloat | None = None,
    xe: onp.ToFloat | None = None,
    k: op.JustInt = 3,
    s: onp.ToFloat = 0,
    nest: op.JustInt | None = None,
) -> Generator[onp.Array1D[np.float64]]: ...

#
def make_splrep(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    *,
    w: onp.ToFloat1D | None = None,
    xb: onp.ToFloat | None = None,
    xe: onp.ToFloat | None = None,
    k: op.JustInt = 3,
    s: onp.ToFloat = 0,
    t: onp.ToFloat1D | None = None,
    nest: op.JustInt | None = None,
) -> BSpline: ...

#
def make_splprep(
    x: onp.ToFloat2D,
    *,
    w: onp.ToFloat1D | None = None,
    u: onp.ToFloat1D | None = None,
    ub: onp.ToFloat | None = None,
    ue: onp.ToFloat | None = None,
    k: op.JustInt = 3,
    s: onp.ToFloat = 0,
    t: onp.ToFloat1D | None = None,
    nest: op.JustInt | None = None,
) -> tuple[BSpline, onp.Array1D[npc.floating]]: ...
