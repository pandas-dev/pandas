from collections.abc import Callable

import numpy as np
import optype.numpy as onp

__all__: list[str] = []

class BaseQuadraticSubproblem:
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        fun: Callable[[onp.Array1D[np.float64]], onp.ToFloat],
        jac: Callable[[onp.Array1D[np.float64]], onp.ToFloat1D],
        hess: Callable[[onp.Array1D[np.float64]], onp.ToFloat2D] | None = None,
        hessp: Callable[[onp.Array1D[np.float64], onp.Array1D[np.float64]], onp.ToFloat1D] | None = None,
    ) -> None: ...
    def __call__(self, /, p: onp.ToFloat1D) -> float | np.float64: ...

    #
    @property
    def fun(self, /) -> float | np.float64: ...
    @property
    def jac_mag(self, /) -> float | np.float64: ...
    @property
    def jac(self, /) -> onp.Array1D[np.float64]: ...
    @property
    def hess(self, /) -> onp.Array2D[np.float64]: ...

    #
    def hessp(self, /, p: onp.ToFloat1D) -> onp.Array1D[np.float64]: ...
    def get_boundaries_intersections(
        self, /, z: onp.ToArray1D, d: onp.ToArray1D, trust_radius: onp.ToFloat
    ) -> list[float | np.float64]: ...  # list of size 2

    #
    def solve(self, /, trust_radius: onp.ToFloat) -> tuple[onp.Array1D[np.float64], bool]: ...
