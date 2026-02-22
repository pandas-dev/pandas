# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from ._optimize import OptimizeResult as _OptimizeResult

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
class OptimizeResult(_OptimizeResult): ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeWarning(UserWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def approx_fprime(xk: object, f: object, epsilon: object = ..., *args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bracket(
    func: object, xa: object = 0.0, xb: object = 1.0, args: object = (), grow_limit: object = 110.0, maxiter: object = 1000
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brent(
    func: object, args: object = (), brack: object = None, tol: object = 1.48e-08, full_output: object = 0, maxiter: object = 500
) -> float | object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brute(
    func: object,
    ranges: object,
    args: object = (),
    Ns: object = 20,
    full_output: object = 0,
    finish: object = ...,
    disp: object = False,
    workers: object = 1,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def check_grad(
    func: object,
    grad: object,
    x0: object,
    *args: object,
    epsilon: object = ...,
    direction: object = "all",
    rng: object = None,
    seed: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin(
    func: object,
    x0: object,
    args: object = (),
    xtol: object = 0.0001,
    ftol: object = 0.0001,
    maxiter: object = None,
    maxfun: object = None,
    full_output: object = 0,
    disp: object = 1,
    retall: object = 0,
    callback: object = None,
    initial_simplex: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_bfgs(
    f: object,
    x0: object,
    fprime: object = None,
    args: object = (),
    gtol: object = 1e-05,
    norm: object = ...,
    epsilon: object = ...,
    maxiter: object = None,
    full_output: object = 0,
    disp: object = 1,
    retall: object = 0,
    callback: object = None,
    xrtol: object = 0,
    c1: object = 0.0001,
    c2: object = 0.9,
    hess_inv0: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_cg(
    f: object,
    x0: object,
    fprime: object = None,
    args: object = (),
    gtol: object = 1e-05,
    norm: object = ...,
    epsilon: object = ...,
    maxiter: object = None,
    full_output: object = 0,
    disp: object = 1,
    retall: object = 0,
    callback: object = None,
    c1: object = 0.0001,
    c2: object = 0.4,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_ncg(
    f: object,
    x0: object,
    fprime: object,
    fhess_p: object = None,
    fhess: object = None,
    args: object = (),
    avextol: object = 1e-05,
    epsilon: object = ...,
    maxiter: object = None,
    full_output: object = 0,
    disp: object = 1,
    retall: object = 0,
    callback: object = None,
    c1: object = 0.0001,
    c2: object = 0.9,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fmin_powell(
    func: object,
    x0: object,
    args: object = (),
    xtol: object = 0.0001,
    ftol: object = 0.0001,
    maxiter: object = None,
    maxfun: object = None,
    full_output: object = 0,
    disp: object = 1,
    retall: object = 0,
    callback: object = None,
    direc: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fminbound(
    func: object,
    x1: object,
    x2: object,
    args: object = (),
    xtol: object = 1e-05,
    maxfun: object = 500,
    full_output: object = 0,
    disp: object = 1,
) -> float | object: ...
@deprecated("will be removed in SciPy v2.0.0")
def golden(
    func: object, args: object = (), brack: object = None, tol: object = ..., full_output: object = 0, maxiter: object = 5000
) -> float | object: ...
@deprecated("will be removed in SciPy v2.0.0")
def line_search(
    f: object,
    myfprime: object,
    xk: object,
    pk: object,
    gfk: object = None,
    old_fval: object = None,
    old_old_fval: object = None,
    args: object = (),
    c1: object = 0.0001,
    c2: object = 0.9,
    amax: object = None,
    extra_condition: object = None,
    maxiter: object = 10,
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
def show_options(solver: object = None, method: object = None, disp: object = True) -> str: ...
