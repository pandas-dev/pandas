# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

__all__ = [
    "OptimizeResult",
    "OptimizeWarning",
    "approx_fprime",
    "bracket",
    "brent",
    "brute",
    "check_grad",
    "fmin",
    "fmin_bfgs",
    "fmin_cg",
    "fmin_ncg",
    "fmin_powell",
    "fminbound",
    "golden",
    "line_search",
    "rosen",
    "rosen_der",
    "rosen_hess",
    "rosen_hess_prod",
    "show_options",
    "zeros",
]

@deprecated("will be removed in SciPy v2.0.0")
def zeros(shape: object, dtype: object = ..., order: object = ..., *, device: object = ..., like: object = ...) -> object: ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(Any): ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeWarning(UserWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def approx_fprime(xk: object, f: object, epsilon: object = ..., *args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bracket(
    func: object, xa: object = ..., xb: object = ..., args: object = ..., grow_limit: object = ..., maxiter: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brent(
    func: object, args: object = ..., brack: object = ..., tol: object = ..., full_output: object = ..., maxiter: object = ...
) -> float | object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brute(
    func: object,
    ranges: object,
    args: object = ...,
    Ns: object = ...,
    full_output: object = ...,
    finish: object = ...,
    disp: object = ...,
    workers: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def check_grad(
    func: object,
    grad: object,
    x0: object,
    *args: object,
    epsilon: object = ...,
    direction: object = ...,
    rng: object = None,
    seed: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin(
    func: object,
    x0: object,
    args: object = ...,
    xtol: object = ...,
    ftol: object = ...,
    maxiter: object = ...,
    maxfun: object = ...,
    full_output: object = ...,
    disp: object = ...,
    retall: object = ...,
    callback: object = ...,
    initial_simplex: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_bfgs(
    f: object,
    x0: object,
    fprime: object = ...,
    args: object = ...,
    gtol: object = ...,
    norm: object = ...,
    epsilon: object = ...,
    maxiter: object = ...,
    full_output: object = ...,
    disp: object = ...,
    retall: object = ...,
    callback: object = ...,
    xrtol: object = ...,
    c1: object = ...,
    c2: object = ...,
    hess_inv0: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_cg(
    f: object,
    x0: object,
    fprime: object = ...,
    args: object = ...,
    gtol: object = ...,
    norm: object = ...,
    epsilon: object = ...,
    maxiter: object = ...,
    full_output: object = ...,
    disp: object = ...,
    retall: object = ...,
    callback: object = ...,
    c1: object = ...,
    c2: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_ncg(
    f: object,
    x0: object,
    fprime: object,
    fhess_p: object = ...,
    fhess: object = ...,
    args: object = ...,
    avextol: object = ...,
    epsilon: object = ...,
    maxiter: object = ...,
    full_output: object = ...,
    disp: object = ...,
    retall: object = ...,
    callback: object = ...,
    c1: object = ...,
    c2: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_powell(
    func: object,
    x0: object,
    args: object = ...,
    xtol: object = ...,
    ftol: object = ...,
    maxiter: object = ...,
    maxfun: object = ...,
    full_output: object = ...,
    disp: object = ...,
    retall: object = ...,
    callback: object = ...,
    direc: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fminbound(
    func: object,
    x1: object,
    x2: object,
    args: object = ...,
    xtol: object = ...,
    maxfun: object = ...,
    full_output: object = ...,
    disp: object = ...,
) -> float | object: ...
@deprecated("will be removed in SciPy v2.0.0")
def golden(
    func: object, args: object = ..., brack: object = ..., tol: object = ..., full_output: object = ..., maxiter: object = ...
) -> float | object: ...
@deprecated("will be removed in SciPy v2.0.0")
def line_search(
    f: object,
    myfprime: object,
    xk: object,
    pk: object,
    gfk: object = ...,
    old_fval: object = ...,
    old_old_fval: object = ...,
    args: object = (),
    c1: object = ...,
    c2: object = ...,
    amax: object = ...,
    extra_condition: object = ...,
    maxiter: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rosen(x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rosen_der(x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rosen_hess(x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rosen_hess_prod(x: object, p: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def show_options(solver: object = ..., method: object = ..., disp: object = ...) -> str: ...
