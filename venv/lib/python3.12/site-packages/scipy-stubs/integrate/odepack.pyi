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
    args: object = ...,
    Dfun: object = ...,
    col_deriv: object = ...,
    full_output: object = ...,
    ml: object = ...,
    mu: object = ...,
    rtol: object = ...,
    atol: object = ...,
    tcrit: object = ...,
    h0: object = ...,
    hmax: object = ...,
    hmin: object = ...,
    ixpr: object = ...,
    mxstep: object = ...,
    mxhnil: object = ...,
    mxordn: object = ...,
    mxords: object = ...,
    printmessg: object = ...,
    tfirst: object = ...,
) -> object: ...
