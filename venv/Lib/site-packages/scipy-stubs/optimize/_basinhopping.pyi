from collections.abc import Callable
from typing import Concatenate, Literal, Protocol, TypeVar, overload, type_check_only

import numpy as np
import optype.numpy as onp

from ._minimize import OptimizeResult as _MinimizeResult
from ._optimize import OptimizeResult as _OptimizeResult
from ._typing import MinimizerKwargs, MinimizerKwargsJac

__all__ = ["basinhopping"]

###

type _Float = float | np.float64
type _Float1D = onp.Array1D[np.float64]

type _CallbackFun[FloatT: onp.ToFloat | onp.ToFloatND] = Callable[[_Float1D, FloatT, bool], bool | None]

_FT_contra = TypeVar("_FT_contra", bound=onp.ToFloat | onp.ToFloatND, contravariant=True)

@type_check_only
class _AcceptTestFun(Protocol[_FT_contra]):
    def __call__(
        self, /, *, f_new: _FT_contra, x_new: onp.ToFloat1D, f_old: _FT_contra, x_old: onp.ToFloat1D
    ) -> bool | Literal["force accept"]: ...

@type_check_only
class OptimizeResult(_OptimizeResult):
    minimization_failures: int
    lowest_optimization_result: _MinimizeResult

    x: _Float1D
    fun: _Float

    message: list[str]
    nit: int
    success: bool

    nfev: int  # might not exist
    njev: int  # might not exist
    nhev: int  # might not exist

###

@overload
def basinhopping(
    func: Callable[Concatenate[_Float1D, ...], onp.ToFloat],
    x0: onp.ToFloat1D,
    niter: int = 100,
    T: onp.ToFloat = 1.0,
    stepsize: onp.ToFloat = 0.5,
    minimizer_kwargs: MinimizerKwargs | None = None,
    take_step: Callable[[_Float1D], onp.ToFloat] | None = None,
    accept_test: _AcceptTestFun[onp.ToFloat] | None = None,
    callback: _CallbackFun[float] | _CallbackFun[np.float64] | None = None,
    interval: int = 50,
    disp: bool = False,
    niter_success: int | None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
    target_accept_rate: onp.ToFloat = 0.5,
    stepwise_factor: onp.ToFloat = 0.9,
) -> OptimizeResult: ...
@overload  # minimizer_kwargs={"jac": True, ...}
def basinhopping(
    func: Callable[Concatenate[_Float1D, ...], tuple[onp.ToFloat, onp.ToFloat1D]],
    x0: onp.ToFloat1D,
    niter: int = 100,
    T: onp.ToFloat = 1.0,
    stepsize: onp.ToFloat = 0.5,
    *,
    minimizer_kwargs: MinimizerKwargsJac,
    take_step: Callable[[_Float1D], onp.ToFloat] | None = None,
    accept_test: _AcceptTestFun[onp.ToFloat] | None = None,
    callback: _CallbackFun[float] | _CallbackFun[np.float64] | None = None,
    interval: int = 50,
    disp: bool = False,
    niter_success: int | None = None,
    rng: onp.random.ToRNG | None = None,
    seed: onp.random.ToRNG | None = None,
    target_accept_rate: onp.ToFloat = 0.5,
    stepwise_factor: onp.ToFloat = 0.9,
) -> OptimizeResult: ...
