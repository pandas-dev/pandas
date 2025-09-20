# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["line_search"]

@deprecated("will be removed in SciPy v2.0.0")
def line_search(
    f: object,
    fprime: object,
    xk: object,
    pk: object,
    gfk: object = ...,
    old_fval: object = ...,
    old_old_fval: object = ...,
    args: object = (),
    c1: object = ...,
    c2: object = ...,
    amax: object = ...,
    amin: object = ...,
    xtol: object = ...,
) -> object: ...
