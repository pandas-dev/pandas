from collections.abc import Callable, Iterable
from typing import Concatenate, Literal, Protocol, type_check_only

import numpy as np
import optype.numpy as onp

from ._constraints import Bounds, LinearConstraint, NonlinearConstraint
from ._optimize import OptimizeResult as _OptimizeResult

__all__ = ["differential_evolution"]

###

type _StrategyName = Literal[
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

@type_check_only
class _DoesMap(Protocol):
    def __call__[VT, RT](self, func: Callable[[VT], RT], iterable: Iterable[VT], /) -> Iterable[RT]: ...

###

class OptimizeResult(_OptimizeResult):
    x: onp.Array1D[np.float64]
    fun: float | np.float64
    population: onp.Array2D[np.float64]
    population_energies: onp.Array1D[np.float64]
    jac: onp.Array2D[np.float64]  # only if `polish=True`

    success: bool
    message: str
    nit: int
    nfev: int

def differential_evolution(
    func: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat],
    bounds: tuple[onp.ToFloat | onp.ToFloat1D, onp.ToFloat | onp.ToFloat1D] | Bounds,
    args: tuple[object, ...] = (),
    strategy: _StrategyName | Callable[[int, onp.Array2D[np.float64], np.random.Generator], onp.ToFloat1D] = "best1bin",
    maxiter: onp.ToJustInt = 1000,
    popsize: onp.ToJustInt = 15,
    tol: onp.ToFloat = 0.01,
    mutation: onp.ToFloat | tuple[onp.ToFloat, onp.ToFloat] = (0.5, 1),
    recombination: onp.ToFloat = 0.7,
    rng: onp.random.ToRNG | None = None,
    callback: Callable[[OptimizeResult], None] | Callable[[onp.Array1D[np.float64], onp.ToFloat], None] | None = None,
    disp: bool = False,
    polish: bool = True,
    init: onp.ToFloat2D | Literal["sobol", "halton", "random", "latinhypercube"] = "latinhypercube",
    atol: onp.ToFloat = 0,
    updating: Literal["immediate", "deferred"] = "immediate",
    workers: int | _DoesMap = 1,
    constraints: NonlinearConstraint | LinearConstraint | Bounds | tuple[()] = (),
    x0: onp.ToArray1D | None = None,
    *,
    integrality: onp.ToBool1D | None = None,
    vectorized: bool = False,
    seed: onp.random.ToRNG | None = None,
) -> OptimizeResult: ...
