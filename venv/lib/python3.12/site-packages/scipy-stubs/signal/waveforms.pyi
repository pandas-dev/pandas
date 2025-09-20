# This module is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

__all__ = ["chirp", "gausspulse", "sawtooth", "square", "sweep_poly", "unit_impulse"]

@deprecated("will be removed in SciPy v2.0.0")
def sawtooth(t: object, width: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def square(t: object, duty: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gausspulse(
    t: object,
    fc: object = ...,
    bw: object = ...,
    bwr: object = ...,
    tpr: object = ...,
    retquad: object = ...,
    retenv: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chirp(
    t: object,
    f0: object,
    t1: object,
    f1: object,
    method: object = ...,
    phi: object = ...,
    vertex_zero: object = ...,
    *,
    complex: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sweep_poly(t: object, poly: object, phi: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def unit_impulse(shape: object, idx: object = ..., dtype: object = ...) -> object: ...
