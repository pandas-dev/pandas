# This file is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = ["IntegrationWarning", "dblquad", "nquad", "quad", "tplquad"]

@deprecated("will be removed in SciPy v2.0.0")
class IntegrationWarning(UserWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def dblquad(
    func: object, a: object, b: object, gfun: object, hfun: object, args: object = ..., epsabs: object = ..., epsrel: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def nquad(func: object, ranges: object, args: object = ..., opts: object = ..., full_output: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def quad(
    func: object,
    a: object,
    b: object,
    args: object = ...,
    full_output: object = ...,
    epsabs: object = ...,
    epsrel: object = ...,
    limit: object = ...,
    points: object = ...,
    weight: object = ...,
    wvar: object = ...,
    wopts: object = ...,
    maxp1: object = ...,
    limlst: object = ...,
    complex_func: object = ...,
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
    args: object = ...,
    epsabs: object = ...,
    epsrel: object = ...,
) -> object: ...
