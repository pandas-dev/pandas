from collections.abc import Callable
from typing import Concatenate, Literal, overload

import numpy as np
import optype.numpy as onp

from ._constraints import Bounds
from ._typing import MinimizerKwargs, MinimizerKwargsJac
from scipy.optimize import OptimizeResult as _OptimizeResult

__all__ = ["dual_annealing"]

###

type _ToBounds = Bounds | tuple[onp.ToFloat1D, onp.ToFloat1D] | onp.ToFloat2D
type _CallbackFun = Callable[[onp.Array1D[np.float64], float, Literal[0, 1, 2]], bool | None]

###

class OptimizeResult(_OptimizeResult):
    success: bool
    status: int
    x: onp.Array1D[np.float64]
    fun: float
    nit: int
    nfev: int
    njev: int
    nhev: int
    message: list[str]

@overload
def dual_annealing(
    func: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat],
    bounds: _ToBounds,
    args: tuple[object, ...] = (),
    maxiter: int = 1_000,
    minimizer_kwargs: MinimizerKwargs | None = None,
    initial_temp: onp.ToFloat = 5_230.0,
    restart_temp_ratio: onp.ToFloat = 2e-5,
    visit: onp.ToFloat = 2.62,
    accept: onp.ToFloat = -5.0,
    maxfun: onp.ToFloat = 10_000_000.0,
    rng: onp.random.ToRNG | None = None,
    no_local_search: bool = False,
    callback: _CallbackFun | None = None,
    x0: onp.ToFloat1D | None = None,
    *,
    seed: onp.random.ToRNG | None = None,
) -> _OptimizeResult: ...
@overload  # minimizer_kwargs={"jac": True, ...}
def dual_annealing(
    func: Callable[Concatenate[onp.Array1D[np.float64], ...], tuple[onp.ToFloat, onp.ToFloat1D]],
    bounds: _ToBounds,
    args: tuple[object, ...] = (),
    maxiter: int = 1_000,
    *,
    minimizer_kwargs: MinimizerKwargsJac,
    initial_temp: onp.ToFloat = 5_230.0,
    restart_temp_ratio: onp.ToFloat = 2e-5,
    visit: onp.ToFloat = 2.62,
    accept: onp.ToFloat = -5.0,
    maxfun: onp.ToFloat = 10_000_000.0,
    rng: onp.random.ToRNG | None = None,
    no_local_search: bool = False,
    callback: _CallbackFun | None = None,
    x0: onp.ToFloat1D | None = None,
    seed: onp.random.ToRNG | None = None,
) -> _OptimizeResult: ...
