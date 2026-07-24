from collections.abc import Callable, Iterable, Sequence
from typing import Concatenate, Literal, Protocol, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp

from ._constraints import Bounds
from ._optimize import OptimizeResult as _OptimizeResult
from ._typing import Constraints, MinimizerKwargs, MinimizerKwargsJac

__all__ = ["shgo"]

###

type _Fun1D[_RT] = Callable[Concatenate[onp.Array1D[np.float64], ...], _RT]

type _ToBounds = tuple[onp.ToFloat | onp.ToFloat1D, onp.ToFloat | onp.ToFloat1D] | Bounds
type _SamplingMethod = Callable[[int, int], onp.ToFloat2D] | Literal["simplicial", "halton", "sobol"]

@type_check_only
class _SHGOOptions(TypedDict, total=False):
    f_min: float
    f_tol: float
    maxiter: int
    maxfev: int
    maxev: int
    mmaxtime: float
    minhgrd: int
    symmetry: Sequence[int] | bool
    jac: _Fun1D[onp.ToFloat1D] | bool  # gradient
    hess: _Fun1D[onp.ToFloat2D]
    hessp: Callable[Concatenate[onp.Array1D[np.float64], onp.Array1D[np.float64], ...], onp.ToFloat1D]
    minimize_every_iter: bool
    local_iter: int
    infty_constraints: bool
    disp: bool

@type_check_only
class _DoesMap(Protocol):
    def __call__[VT, RT](self, func: Callable[[VT], RT], iterable: Iterable[VT], /) -> Iterable[RT]: ...

###

class OptimizeResult(_OptimizeResult):
    x: onp.Array1D[np.float64]
    xl: onp.Array2D[np.float64]
    fun: np.float64
    funl: onp.Array1D[np.float64]
    success: bool
    message: str
    nfev: int
    nlfev: int
    nljev: int  # undocumented
    nlhev: int  # undocumented
    nit: int

@overload
def shgo(
    func: _Fun1D[onp.ToFloat],
    bounds: _ToBounds,
    args: tuple[object, ...] = (),
    constraints: Constraints | None = None,
    n: int = 100,
    iters: int = 1,
    callback: Callable[[onp.Array1D[np.float64]], None] | None = None,
    minimizer_kwargs: MinimizerKwargs | None = None,
    options: _SHGOOptions | None = None,
    sampling_method: _SamplingMethod = "simplicial",
    *,
    workers: int | _DoesMap = 1,
) -> OptimizeResult: ...
@overload  # minimizer_kwargs={"jac": True, ...}
def shgo(
    func: _Fun1D[tuple[onp.ToFloat, onp.ToFloat1D]],
    bounds: _ToBounds,
    args: tuple[object, ...] = (),
    constraints: Constraints | None = None,
    n: int = 100,
    iters: int = 1,
    callback: Callable[[onp.Array1D[np.float64]], None] | None = None,
    *,
    minimizer_kwargs: MinimizerKwargsJac,
    options: _SHGOOptions | None = None,
    sampling_method: _SamplingMethod = "simplicial",
    workers: int | _DoesMap = 1,
) -> OptimizeResult: ...
