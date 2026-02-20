# This file is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = ["IntegrationWarning", "dblquad", "nquad", "quad", "tplquad"]

@deprecated("will be removed in SciPy v2.0.0")
class IntegrationWarning(UserWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def dblquad(
    func: object,
    a: object,
    b: object,
    gfun: object,
    hfun: object,
    args: object = (),
    epsabs: object = 1.49e-8,
    epsrel: object = 1.49e-8,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def nquad(func: object, ranges: object, args: object = None, opts: object = None, full_output: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def quad(
    func: object,
    a: object,
    b: object,
    args: object = (),
    full_output: object = 0,
    epsabs: object = 1.49e-8,
    epsrel: object = 1.49e-8,
    limit: object = 50,
    points: object = None,
    weight: object = None,
    wvar: object = None,
    wopts: object = None,
    maxp1: object = 50,
    limlst: object = 50,
    complex_func: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tplquad(
    func: object,
    a: object,
    b: object,
    gfun: object,
    hfun: object,
    qfun: object,
    rfun: object,
    args: object = (),
    epsabs: object = 1.49e-8,
    epsrel: object = 1.49e-8,
) -> object: ...
