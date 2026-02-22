# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["LbfgsInvHessProduct", "OptimizeResult", "fmin_l_bfgs_b", "zeros"]

@deprecated("will be removed in SciPy v2.0.0")
def zeros(shape: object, dtype: object = ..., order: object = ..., *, device: object = ..., like: object = ...) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(_OptimizeResult): ...

@deprecated("will be removed in SciPy v2.0.0")
def fmin_l_bfgs_b(
    func: object,
    x0: object,
    fprime: object = None,
    args: object = (),
    approx_grad: object = 0,
    bounds: object = None,
    m: object = 10,
    factr: object = 1e7,
    pgtol: object = 1e-5,
    epsilon: object = 1e-8,
    iprint: object = ...,
    maxfun: object = 15_000,
    maxiter: object = 15_000,
    disp: object = ...,
    callback: object = None,
    maxls: object = 20,
) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class LbfgsInvHessProduct:
    def __init__(self, /, sk: object, yk: object) -> None: ...
    def todense(self, /) -> object: ...
