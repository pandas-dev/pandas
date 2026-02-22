# This file is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = ["ODEintWarning", "odeint"]

@deprecated("will be removed in SciPy 2.0.0.")
class ODEintWarning(Warning): ...

@deprecated("will be removed in SciPy v2.0.0")
def odeint(
    func: object,
    y0: object,
    t: object,
    args: object = (),
    Dfun: object = None,
    col_deriv: object = 0,
    full_output: object = 0,
    ml: object = None,
    mu: object = None,
    rtol: object = None,
    atol: object = None,
    tcrit: object = None,
    h0: object = 0.0,
    hmax: object = 0.0,
    hmin: object = 0.0,
    ixpr: object = 0,
    mxstep: object = 0,
    mxhnil: object = 0,
    mxordn: object = 12,
    mxords: object = 5,
    printmessg: object = 0,
    tfirst: object = False,
) -> object: ...
