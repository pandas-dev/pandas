# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

__all__ = ["OptimizeResult", "fmin_tnc", "zeros"]

@deprecated("will be removed in SciPy v2.0.0")
def zeros(shape: object, dtype: object = ..., order: object = ..., *, device: object = ..., like: object = ...) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(Any): ...

@deprecated("will be removed in SciPy v2.0.0")
def fmin_tnc(
    func: object,
    x0: object,
    fprime: object = ...,
    args: object = ...,
    approx_grad: object = ...,
    bounds: object = ...,
    epsilon: object = ...,
    scale: object = ...,
    offset: object = ...,
    messages: object = ...,
    maxCGit: object = ...,
    maxfun: object = ...,
    eta: object = ...,
    stepmx: object = ...,
    accuracy: object = ...,
    fmin: object = ...,
    ftol: object = ...,
    xtol: object = ...,
    pgtol: object = ...,
    rescale: object = ...,
    disp: object = ...,
    callback: object = ...,
) -> object: ...
