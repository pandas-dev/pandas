# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from ._bsplines import BSpline as _BSpline

__all__ = [
    "BSpline",
    "bisplev",
    "bisplrep",
    "insert",
    "spalde",
    "splantider",
    "splder",
    "splev",
    "splint",
    "splprep",
    "splrep",
    "sproot",
]

# _bsplines
@deprecated("will be removed in SciPy v2.0.0")
class BSpline(_BSpline): ...

# _fitpack_impl
@deprecated("will be removed in SciPy v2.0.0")
def splprep(
    x: object,
    w: object = None,
    u: object = None,
    ub: object = None,
    ue: object = None,
    k: object = 3,
    task: object = 0,
    s: object = None,
    t: object = None,
    full_output: object = 0,
    nest: object = None,
    per: object = 0,
    quiet: object = 1,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def splrep(
    x: object,
    y: object,
    w: object = None,
    xb: object = None,
    xe: object = None,
    k: object = 3,
    task: object = 0,
    s: object = None,
    t: object = None,
    full_output: object = 0,
    per: object = 0,
    quiet: object = 1,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def splev(x: object, tck: object, der: object = 0, ext: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def splint(a: object, b: object, tck: object, full_output: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sproot(tck: object, mest: object = 10) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spalde(x: object, tck: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def insert(x: object, tck: object, m: object = 1, per: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def splder(tck: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def splantider(tck: object, n: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bisplrep(
    x: object,
    y: object,
    z: object,
    w: object = None,
    xb: object = None,
    xe: object = None,
    yb: object = None,
    ye: object = None,
    kx: object = 3,
    ky: object = 3,
    task: object = 0,
    s: object = None,
    eps: object = 1e-16,
    tx: object = None,
    ty: object = None,
    full_output: object = 0,
    nxest: object = None,
    nyest: object = None,
    quiet: object = 1,
) -> object: ...
def bisplev(x: object, y: object, tck: object, dx: object = 0, dy: object = 0) -> object: ...
