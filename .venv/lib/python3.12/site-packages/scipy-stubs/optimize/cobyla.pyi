# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["OptimizeResult", "fmin_cobyla"]

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(_OptimizeResult): ...

@deprecated("will be removed in SciPy v2.0.0")
def fmin_cobyla(
    func: object,
    x0: object,
    cons: object,
    args: object = (),
    consargs: object = None,
    rhobeg: object = 1.0,
    rhoend: object = 0.0001,
    maxfun: object = 1000,
    disp: object = None,
    catol: object = 0.0002,
    *,
    callback: object = None,
) -> object: ...
