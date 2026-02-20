# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["chirp", "gausspulse", "sawtooth", "square", "sweep_poly", "unit_impulse"]

@deprecated("will be removed in SciPy v2.0.0")
def sawtooth(t: object, width: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def square(t: object, duty: object = 0.5) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gausspulse(
    t: object,
    fc: object = 1_000,
    bw: object = 0.5,
    bwr: object = -6,
    tpr: object = -60,
    retquad: object = False,
    retenv: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chirp(
    t: object,
    f0: object,
    t1: object,
    f1: object,
    method: object = "linear",
    phi: object = 0,
    vertex_zero: object = True,
    *,
    complex: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sweep_poly(t: object, poly: object, phi: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def unit_impulse(shape: object, idx: object = None, dtype: object = ...) -> object: ...
