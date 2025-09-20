# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

__all__ = ["LbfgsInvHessProduct", "OptimizeResult", "fmin_l_bfgs_b", "zeros"]

@deprecated("will be removed in SciPy v2.0.0")
def zeros(shape: object, dtype: object = ..., order: object = ..., *, device: object = ..., like: object = ...) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(Any): ...

@deprecated("will be removed in SciPy v2.0.0")
def fmin_l_bfgs_b(
    func: object,
    x0: object,
    fprime: object = ...,
    args: object = ...,
    approx_grad: object = ...,
    bounds: object = ...,
    m: object = ...,
    factr: object = ...,
    pgtol: object = ...,
    epsilon: object = ...,
    iprint: object = ...,
    maxfun: object = ...,
    maxiter: object = ...,
    disp: object = ...,
    callback: object = ...,
    maxls: object = ...,
) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class LbfgsInvHessProduct:
    def __init__(self, /, sk: object, yk: object) -> None: ...
    def todense(self, /) -> object: ...
