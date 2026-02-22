# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["line_search"]

@deprecated("will be removed in SciPy v2.0.0")
def line_search(
    f: object,
    fprime: object,
    xk: object,
    pk: object,
    gfk: object = None,
    old_fval: object = None,
    old_old_fval: object = None,
    args: object = (),
    c1: object = 1e-4,
    c2: object = 0.9,
    amax: object = 50,
    amin: object = 1e-8,
    xtol: object = 1e-14,
) -> object: ...
