# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["OptimizeResult", "fmin_tnc", "zeros"]

@deprecated("will be removed in SciPy v2.0.0")
def zeros(shape: object, dtype: object = ..., order: object = ..., *, device: object = ..., like: object = ...) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(_OptimizeResult): ...

@deprecated("will be removed in SciPy v2.0.0")
def fmin_tnc(
    func: object,
    x0: object,
    fprime: object = None,
    args: object = (),
    approx_grad: object = 0,
    bounds: object = None,
    epsilon: object = 1e-8,
    scale: object = None,
    offset: object = None,
    messages: object = 15,
    maxCGit: object = -1,
    maxfun: object = None,
    eta: object = -1,
    stepmx: object = 0,
    accuracy: object = 0,
    fmin: object = 0,
    ftol: object = -1,
    xtol: object = -1,
    pgtol: object = -1,
    rescale: object = -1,
    disp: object = None,
    callback: object = None,
) -> object: ...
