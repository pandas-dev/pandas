from collections.abc import Callable
from typing import Concatenate, TypedDict, Unpack, type_check_only

import numpy as np
import optype.numpy as onp

from ._minimize import OptimizeResult

__all__ = ["_minimize_trust_krylov"]

@type_check_only
class _TrustRegionOptions(TypedDict, total=False):
    initial_trust_radius: onp.ToFloat
    max_trust_radius: onp.ToFloat
    eta: onp.ToFloat
    gtol: onp.ToFloat
    maxiter: onp.ToJustInt
    disp: bool

###

def _minimize_trust_krylov(
    fun: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat],
    x0: onp.ToFloat1D,
    args: tuple[object, ...] = (),
    jac: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat1D] | None = None,
    hess: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat2D] | None = None,
    hessp: Callable[Concatenate[onp.Array1D[np.float64], onp.Array1D[np.float64], ...], onp.ToFloat1D] | None = None,
    inexact: bool = True,
    **trust_region_options: Unpack[_TrustRegionOptions],
) -> OptimizeResult: ...
