from collections.abc import Callable, Iterable
from typing import Concatenate, Literal, TypeAlias, TypeVar

import numpy as np
import optype.numpy as onp

from ._constraints import Bounds, LinearConstraint, NonlinearConstraint
from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["differential_evolution"]

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]

_StrategyName: TypeAlias = Literal[
    "best1bin",
    "best1exp",
    "best2exp",
    "best2bin",
    "rand1bin",
    "rand1exp",
    "rand2bin",
    "rand2exp",
    "randtobest1bin",
    "randtobest1exp",
    "currenttobest1bin",
    "currenttobest1exp",
]

_S = TypeVar("_S")
_T = TypeVar("_T")

###

class OptimizeResult(_OptimizeResult):
    x: _Float1D
    fun: float | np.float64
    population: _Float2D
    population_energies: _Float1D
    jac: _Float2D  # only if `polish=True`

    success: bool
    message: str
    nit: int
    nfev: int

def differential_evolution(
    func: Callable[Concatenate[_Float1D, ...], onp.ToFloat],
    bounds: tuple[onp.ToFloat | onp.ToFloat1D, onp.ToFloat | onp.ToFloat1D] | Bounds,
    args: tuple[object, ...] = (),
    strategy: _StrategyName | Callable[[int, _Float2D, np.random.Generator], onp.ToFloat1D] = "best1bin",
    maxiter: onp.ToJustInt = 1000,
    popsize: onp.ToJustInt = 15,
    tol: onp.ToFloat = 0.01,
    mutation: onp.ToFloat | tuple[onp.ToFloat, onp.ToFloat] = (0.5, 1),
    recombination: onp.ToFloat = 0.7,
    rng: onp.random.ToRNG | None = None,
    callback: Callable[[OptimizeResult], None] | Callable[[_Float1D, onp.ToFloat], None] | None = None,
    disp: onp.ToBool = False,
    polish: onp.ToBool = True,
    init: onp.ToFloat2D | Literal["sobol", "halton", "random", "latinhypercube"] = "latinhypercube",
    atol: onp.ToFloat = 0,
    updating: Literal["immediate", "deferred"] = "immediate",
    workers: onp.ToJustInt | Callable[[Callable[[_S], _T], Iterable[_S]], Iterable[_T]] = 1,
    constraints: NonlinearConstraint | LinearConstraint | Bounds | tuple[()] = (),
    x0: onp.ToArray1D | None = None,
    *,
    integrality: onp.ToBool1D | None = None,
    vectorized: onp.ToBool = False,
    seed: onp.random.ToRNG | None = None,
) -> OptimizeResult: ...
