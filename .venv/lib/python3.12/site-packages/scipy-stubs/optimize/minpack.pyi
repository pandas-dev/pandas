# This file is not meant for public use and will be removed in SciPy v2.0.0.

from typing_extensions import deprecated

from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["OptimizeResult", "OptimizeWarning", "curve_fit", "fixed_point", "fsolve", "least_squares", "leastsq", "zeros"]

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeResult(_OptimizeResult): ...

@deprecated("will be removed in SciPy v2.0.0")
class OptimizeWarning(UserWarning): ...

@deprecated("will be removed in SciPy v2.0.0")
def zeros(shape: object, dtype: object = ..., order: object = ..., *, device: object = ..., like: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def least_squares(
    fun: object,
    x0: object,
    jac: object = "2-point",
    bounds: object = ...,
    method: object = "trf",
    ftol: object = 1e-8,
    xtol: object = 1e-8,
    gtol: object = 1e-8,
    x_scale: object = None,
    loss: object = "linear",
    f_scale: object = 1.0,
    diff_step: object = None,
    tr_solver: object = None,
    tr_options: object = None,
    jac_sparsity: object = None,
    max_nfev: object = None,
    verbose: object = 0,
    args: tuple[object, ...] = (),
    kwargs: dict[str, object] | None = None,
    callback: object | None = None,
    workers: object | None = None,
) -> OptimizeResult: ...  # pyright: ignore[reportDeprecated]  # ty: ignore[deprecated]
@deprecated("will be removed in SciPy v2.0.0")
def fsolve(
    func: object,
    x0: object,
    args: object = (),
    fprime: object = None,
    full_output: object = 0,
    col_deriv: object = 0,
    xtol: object = ...,
    maxfev: object = 0,
    band: object = None,
    epsfcn: object = None,
    factor: object = 100,
    diag: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def leastsq(
    func: object,
    x0: object,
    args: object = (),
    Dfun: object = None,
    full_output: object = False,
    col_deriv: object = False,
    ftol: object = ...,
    xtol: object = ...,
    gtol: object = 0.0,
    maxfev: object = 0,
    epsfcn: object = None,
    factor: object = 100,
    diag: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def curve_fit(
    f: object,
    xdata: object,
    ydata: object,
    p0: object = None,
    sigma: object = None,
    absolute_sigma: object = False,
    check_finite: object = None,
    bounds: object = ...,
    method: object = None,
    jac: object = None,
    *,
    full_output: object = False,
    nan_policy: object = None,
    **kwargs: object,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fixed_point(
    func: object, x0: object, args: object = (), xtol: object = 1e-08, maxiter: object = 500, method: object = "del2"
) -> object: ...
