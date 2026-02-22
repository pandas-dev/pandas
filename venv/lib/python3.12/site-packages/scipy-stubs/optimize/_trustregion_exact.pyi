from collections.abc import Callable
from typing import ClassVar, Concatenate, Final, TypedDict
from typing_extensions import Unpack, override

import numpy as np
import optype.numpy as onp

from ._minimize import OptimizeResult
from ._trustregion import BaseQuadraticSubproblem

__all__ = ["IterativeSubproblem", "_minimize_trustregion_exact", "estimate_smallest_singular_value", "singular_leading_submatrix"]

class _TrustRegionOptions(TypedDict, total=False):
    initial_trust_radius: onp.ToFloat
    max_trust_radius: onp.ToFloat
    eta: onp.ToFloat
    gtol: onp.ToFloat

###

def _minimize_trustregion_exact(
    fun: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat],
    x0: onp.ToFloat1D,
    args: tuple[object, ...] = (),
    jac: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat1D] | None = None,
    hess: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat2D] | None = None,
    **trust_region_options: Unpack[_TrustRegionOptions],
) -> OptimizeResult: ...
def estimate_smallest_singular_value(U: onp.ToFloat2D) -> tuple[float | np.float64, onp.Array1D[np.float64]]: ...
def gershgorin_bounds(H: onp.ToFloat2D) -> tuple[float | np.float64, float | np.float64]: ...
def singular_leading_submatrix(
    A: onp.ToFloat2D, U: onp.ToFloat2D, k: onp.ToJustInt
) -> tuple[float | np.float64, onp.Array1D[np.float64]]: ...

class IterativeSubproblem(BaseQuadraticSubproblem):
    UPDATE_COEFF: ClassVar[float] = 0.01
    EPS: ClassVar[float | np.float64] = ...
    MAXITER_DEFAULT: ClassVar[float] = 25

    CLOSE_TO_ZERO: Final[float | np.float64]
    dimension: Final[int]
    previous_tr_radius: int
    niter: int
    k_easy: float | np.float64
    k_hard: float | np.float64
    hess_inf: float | np.float64
    hess_fro: float | np.float64
    lambda_lb: float | np.float64 | None  # intially `None`
    lambda_current: float | np.float64  # set in `solve()`

    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        fun: Callable[[onp.Array1D[np.float64]], onp.ToFloat],
        jac: Callable[[onp.Array1D[np.float64]], onp.ToFloat1D],
        hess: Callable[[onp.Array1D[np.float64]], onp.ToFloat2D],
        hessp: Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat1D] | None = None,
        k_easy: onp.ToFloat = 0.1,
        k_hard: onp.ToFloat = 0.2,
        maxiter: float | None = None,
    ) -> None: ...
    @override
    # pyrefly: ignore[bad-param-name-override]
    def solve(self, /, tr_radius: onp.ToFloat) -> tuple[onp.ArrayND[np.float64], bool]: ...  # pyright: ignore[reportIncompatibleMethodOverride] # ty: ignore[invalid-method-override]
