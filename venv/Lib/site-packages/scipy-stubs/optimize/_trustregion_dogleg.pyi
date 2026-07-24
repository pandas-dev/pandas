from typing import override

import numpy as np
import optype.numpy as onp

from ._trustregion import BaseQuadraticSubproblem

__all__: list[str] = []

class DoglegSubproblem(BaseQuadraticSubproblem):
    def cauchy_point(self, /) -> onp.Array1D[np.float64]: ...
    def newton_point(self, /) -> onp.Array1D[np.float64]: ...

    #
    @override
    def solve(self, /, trust_radius: onp.ToFloat) -> tuple[onp.Array1D[np.float64], bool]: ...
