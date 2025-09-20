# This file is not meant for public use and will be removed in SciPy v2.0.0.

from collections.abc import Callable
from typing import Any
from typing_extensions import deprecated

__all__ = ["OptimizeResult", "fmin_slsqp", "slsqp"]

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(Any): ...

@deprecated("will be removed in SciPy v2.0.0")
def fmin_slsqp(
    func: object,
    x0: object,
    eqcons: object = ...,
    f_eqcons: object = ...,
    ieqcons: object = ...,
    f_ieqcons: object = ...,
    bounds: object = ...,
    fprime: object = ...,
    fprime_eqcons: object = ...,
    fprime_ieqcons: object = ...,
    args: object = ...,
    iter: object = ...,
    acc: object = ...,
    iprint: object = ...,
    disp: object = ...,
    full_output: object = ...,
    epsilon: object = ...,
    callback: object = ...,
) -> Any: ...

slsqp: Callable[..., Any] = ...  # deprecated
