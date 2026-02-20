# This file is not meant for public use and will be removed in SciPy v2.0.0.

from collections.abc import Callable
from typing import Any
from typing_extensions import deprecated

from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["OptimizeResult", "fmin_slsqp", "slsqp"]

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(_OptimizeResult): ...

@deprecated("will be removed in SciPy v2.0.0")
def fmin_slsqp(
    func: object,
    x0: object,
    eqcons: object = (),
    f_eqcons: object = None,
    ieqcons: object = (),
    f_ieqcons: object = None,
    bounds: object = (),
    fprime: object = None,
    fprime_eqcons: object = None,
    fprime_ieqcons: object = None,
    args: object = (),
    iter: object = 100,
    acc: object = 1e-6,
    iprint: object = 1,
    disp: object = None,
    full_output: object = 0,
    epsilon: object = ...,
    callback: object = None,
) -> Any: ...

slsqp: Callable[..., Any] = ...  # deprecated
