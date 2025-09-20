# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing import Any
from typing_extensions import deprecated

__all__ = ["OptimizeResult", "OptimizeWarning", "curve_fit", "fixed_point", "fsolve", "least_squares", "leastsq", "zeros"]

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(Any): ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeWarning(UserWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def zeros(shape: object, dtype: object = ..., order: object = ..., *, device: object = ..., like: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def least_squares(
    fun: object,
    x0: object,
    jac: object = ...,
    bounds: object = ...,
    method: object = ...,
    ftol: object = ...,
    xtol: object = ...,
    gtol: object = ...,
    x_scale: object = ...,
    loss: object = ...,
    f_scale: object = ...,
    diff_step: object = ...,
    tr_solver: object = ...,
    tr_options: object = ...,
    jac_sparsity: object = ...,
    max_nfev: object = ...,
    verbose: object = ...,
    args: tuple[object, ...] = (),
    kwargs: dict[str, object] | None = None,
    callback: object | None = None,
    workers: object | None = None,
) -> OptimizeResult: ...  # pyright: ignore[reportDeprecated]
@deprecated("will be removed in SciPy v2.0.0")
def fsolve(
    func: object,
    x0: object,
    args: object = ...,
    fprime: object = ...,
    full_output: object = ...,
    col_deriv: object = ...,
    xtol: object = ...,
    maxfev: object = ...,
    band: object = ...,
    epsfcn: object = ...,
    factor: object = ...,
    diag: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def leastsq(
    func: object,
    x0: object,
    args: object = ...,
    Dfun: object = ...,
    full_output: object = ...,
    col_deriv: object = ...,
    ftol: object = ...,
    xtol: object = ...,
    gtol: object = ...,
    maxfev: object = ...,
    epsfcn: object = ...,
    factor: object = ...,
    diag: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def curve_fit(
    f: object,
    xdata: object,
    ydata: object,
    p0: object = ...,
    sigma: object = ...,
    absolute_sigma: object = ...,
    check_finite: object = ...,
    bounds: object = ...,
    method: object = ...,
    jac: object = ...,
    *,
    full_output: object = ...,
    nan_policy: object = ...,
    **kwargs: object,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fixed_point(
    func: object, x0: object, args: object = ..., xtol: object = ..., maxiter: object = ..., method: object = ...
) -> object: ...
